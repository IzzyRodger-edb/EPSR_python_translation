#!/usr/bin/env python3
"""Run one test case with the Fortran reference programs and with the Python translation,
and compare every file after every program call.

usage: tools/run_ex_py.py CASE_DIR [--fortran-build DIR] [--work DIR] [--mode resume|chain]
                                    [--passes N] [--steps 1,3] [--max-niter N] [--timeout S]

CASE_DIR holds the inputs plus either ``commands.txt`` (one program call per line,
optionally ``< stdin_file``; ``repeat N`` chains the list N times), as in test_problems/,
or ``steps.txt`` (``run <program> args [< stdin]``, ``set <file> <keyword> <value>``,
``delete <files>``), as in validation/.

Modes (PLAN.md section 6):

* ``resume`` (default): before each program step the Python folder is reset to the
  Fortran's folder *before that step*, so every step is tested on its own;
* ``chain``: the Python runs every step from the raw inputs, like the Fortran.

``--steps 1,3`` runs the Python only for those steps (the Fortran runs all of them).

``--max-niter N`` (resume mode): an ``epsr`` step whose input file asks for more than N
iterations is compared on a side copy of the folder with ``niter`` set to N, run by both
the Fortran and the Python; the main Fortran chain still runs the full step.

Prints one line per step with IDENTICAL/DIFFERS and the Fortran and Python wall times, and
writes ``<work>/summary.json``.  Exit status 0 if every compared step was identical.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SKIP = {"commands.txt", "README.md", "fortran_reference.tar.gz", "steps.txt"}


def copy_inputs(src, dst):
    os.makedirs(dst, exist_ok=True)
    for f in os.listdir(src):
        p = os.path.join(src, f)
        if os.path.isfile(p) and f not in SKIP and not f.endswith((".pdf", ".png", ".py")):
            shutil.copy2(p, os.path.join(dst, f))


def snapshot(src, dst):
    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


def read_steps(case):
    """[(kind, words, stdin)]: kind is run (words = program + args), set or delete."""
    steps = []
    if os.path.exists(os.path.join(case, "commands.txt")):
        lines, repeat = [], 1
        with open(os.path.join(case, "commands.txt")) as fh:
            for ln in fh:
                ln = ln.split("#", 1)[0].strip()
                if not ln:
                    continue
                if ln.startswith("repeat"):
                    repeat = int(ln.split()[1])
                    continue
                lines.append("run " + ln)
        lines = lines * repeat
    else:
        with open(os.path.join(case, "steps.txt")) as fh:
            lines = [ln.split("#", 1)[0].strip() for ln in fh]
    for ln in lines:
        if not ln:
            continue
        kind, _, rest = ln.partition(" ")
        stdin = None
        if kind == "run" and "<" in rest:
            rest, stdin = rest.split("<", 1)
            rest, stdin = rest.strip(), stdin.strip()
        steps.append((kind, rest.split(), stdin))
    return steps


def get_niter(folder, inp):
    try:
        with open(os.path.join(folder, inp), encoding="latin-1") as fh:
            for ln in fh:
                if ln[:12].strip() == "niter":
                    return int(ln[12:].split()[0])
    except (OSError, ValueError, IndexError):
        pass
    return None


def set_keyword(folder, fname, kw, *value):
    subprocess.run([sys.executable, os.path.join(ROOT, "validation", "set_keyword.py"),
                    os.path.join(folder, fname), kw] + list(value), check=True)


def edit(kind, words, folder):
    if kind == "set":
        set_keyword(folder, words[0], words[1], *words[2:])
    elif kind == "delete":
        for f in words:
            try:
                os.remove(os.path.join(folder, f))
            except FileNotFoundError:
                pass


def run(cmd, cwd, stdin_path, log, timeout, env):
    t = time.perf_counter()
    with open(stdin_path or os.devnull, "rb") as fi, open(log, "wb") as fo:
        try:
            p = subprocess.run(cmd, cwd=cwd, stdin=fi, stdout=fo, stderr=subprocess.STDOUT, timeout=timeout, env=env)
            rc = p.returncode
        except subprocess.TimeoutExpired:
            rc = "timeout"
    return rc, time.perf_counter() - t


def compare(ref, cand):
    p = subprocess.run([sys.executable, os.path.join(HERE, "compare_outputs.py"), ref, cand, "--quiet"],
                       capture_output=True, text=True)
    # both folders are complete copies, so a file only in one of them is a difference too
    lines = [ln for ln in p.stdout.splitlines() if ln.strip() and ln.strip() not in ("SOME FILES DIFFER", "ALL IDENTICAL")]
    return p.returncode == 0, lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--fortran-build", default=os.path.join(ROOT, "fortran_build"))
    ap.add_argument("--work")
    ap.add_argument("--mode", choices=["resume", "chain"], default="resume")
    ap.add_argument("--passes", type=int, help="stop after N program steps")
    ap.add_argument("--steps", help="comma-separated program-step numbers to run in Python")
    ap.add_argument("--max-niter", type=int)
    ap.add_argument("--timeout", type=float, default=36000)
    a = ap.parse_args()
    case = os.path.abspath(a.case)
    fb = os.path.abspath(a.fortran_build)
    work = os.path.abspath(a.work or tempfile.mkdtemp(prefix="epsrpy-"))
    if os.path.exists(work):
        shutil.rmtree(work)
    fdir, pdir, logs, pristine = (os.path.join(work, d) for d in ("fortran", "python", "logs", "pristine"))
    for d in (fdir, pdir, pristine):
        copy_inputs(case, d)
    os.makedirs(logs)
    steps = read_steps(case)
    only = {int(s) for s in a.steps.split(",")} if a.steps else None
    fenv = dict(os.environ, OMP_NUM_THREADS="1")
    penv = dict(os.environ)
    penv["PYTHONPATH"] = os.path.join(ROOT, "src") + os.pathsep + penv.get("PYTHONPATH", "")
    try:
        import resource

        resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
    except (ValueError, OSError):
        pass
    results = []
    all_ok = True
    n = 0
    for kind, words, stdin in steps:
        if kind != "run":
            edit(kind, words, fdir)
            if a.mode == "chain":
                edit(kind, words, pdir)
            continue
        n += 1
        if a.passes and n > a.passes:
            break
        prog, args = words[0], words[1:]
        stdin_path = os.path.join(pristine, stdin) if stdin else None
        before = os.path.join(work, "before")
        snapshot(fdir, before)
        frc, ft = run([os.path.join(fb, prog), "./", prog] + args, fdir, stdin_path, os.path.join(logs, f"f{n}.log"),
                      a.timeout, fenv)
        res = {"step": n, "program": prog, "args": args, "fortran_rc": frc, "fortran_s": round(ft, 3)}
        if only is None or n in only:
            ref = fdir
            short = None
            if a.mode == "resume":
                shutil.rmtree(pdir)
                shutil.copytree(before, pdir)
                if a.max_niter and prog == "epsr" and args:
                    ni = get_niter(pdir, args[0])
                    if ni is not None and ni > a.max_niter:
                        short = a.max_niter
                        set_keyword(pdir, args[0], "niter", str(short))
                        ref = os.path.join(work, "short_ref")
                        snapshot(pdir, ref)
                        _, ft = run([os.path.join(fb, prog), "./", prog] + args, ref, stdin_path,
                                    os.path.join(logs, f"f{n}_short.log"), a.timeout, fenv)
                        res.update(niter_full=ni, niter_compared=short, fortran_short_s=round(ft, 3))
            prc, pt = run([sys.executable, "-m", "epsr26.run", prog, "./", prog] + args, pdir, stdin_path,
                          os.path.join(logs, f"p{n}.log"), a.timeout, penv)
            ok, diffs = compare(ref, pdir)
            res.update(python_rc=prc, python_s=round(pt, 3), identical=ok, diffs=diffs[:40],
                       slowdown=round(pt / ft, 1) if ft > 0 else None)
            all_ok &= ok
            status = "IDENTICAL" if ok else "DIFFERS"
            extra = f" [niter {short} of {res['niter_full']}]" if short else ""
            print(f"step {n:3d} {prog:14s} {status:9s} fortran {ft:8.2f} s  python {pt:9.2f} s"
                  f"  (x{res['slowdown']}){extra}" + ("" if prc == 0 else f"  python exit {prc}"), flush=True)
            if not ok:
                for d in diffs[:12]:
                    print("      " + d)
        results.append(res)
    with open(os.path.join(work, "summary.json"), "w") as fh:
        json.dump({"case": case, "mode": a.mode, "steps": results}, fh, indent=1)
    print(f"work folder: {work}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
