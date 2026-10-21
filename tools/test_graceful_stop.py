#!/usr/bin/env python3
"""End-to-end check of Ctrl+C during `epsr` in the Python shell (PLAN.md 4.3).

usage: tools/test_graceful_stop.py WORK_DIR

A copy of EPSR-test-ex1 is set to niter 3 and `epsr` is started from the shell.  Ctrl+C is
pressed during iteration 1.  The shell must finish iteration 1, write every result file,
return to the prompt and exit cleanly with `bye`; the folder must then be byte-identical to
a Fortran run of the same inputs with niter 1, and no child process may remain.
"""

import os
import shutil
import subprocess
import sys
import time

import pexpect

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CASE = os.path.join(ROOT, "test_problems", "200DMSO400Water", "EPSR-test-ex1")
INP = "200DMSO400Water.EPSR.inp"


def copy_case(dst):
    os.makedirs(dst)
    for f in os.listdir(CASE):
        p = os.path.join(CASE, f)
        if os.path.isfile(p) and f not in ("commands.txt", "README.md", "fortran_reference.tar.gz") and not f.endswith(".pdf"):
            shutil.copy2(p, dst)


def main():
    work = os.path.abspath(sys.argv[1])
    if os.path.exists(work):
        shutil.rmtree(work)
    py, fo = os.path.join(work, "python"), os.path.join(work, "fortran")
    copy_case(py)
    copy_case(fo)
    setkw = [sys.executable, os.path.join(ROOT, "validation", "set_keyword.py")]
    subprocess.run(setkw + [os.path.join(py, INP), "niter", "3"], check=True)
    subprocess.run(setkw + [os.path.join(fo, INP), "niter", "1"], check=True)
    # Fortran reference: niter 1
    subprocess.run(f"ulimit -s unlimited; OMP_NUM_THREADS=1 {ROOT}/fortran_build/epsr ./ epsr {INP} > run.log 2>&1",
                   shell=True, cwd=fo, check=True)
    # Python shell: niter 3, Ctrl+C during the first iteration
    env = dict(os.environ, PYTHONUNBUFFERED="1", TERM="dumb")
    env["PYTHONPATH"] = os.path.join(ROOT, "src") + os.pathsep + env.get("PYTHONPATH", "")
    log = open(os.path.join(work, "shell.log"), "w")
    c = pexpect.spawn(sys.executable, ["-m", "epsr26"], cwd=py, env=env, encoding="utf-8", timeout=7200, logfile=log)
    c.expect("EPSRshell> ")
    t0 = time.time()
    c.sendline(f"epsr {INP}")
    c.expect("Start molfit iterations")  # well inside iteration 1
    c.sendintr()
    c.expect("EPSR will stop after the current iteration")
    c.expect("this is the last iteration")
    c.expect("EPSRshell> ")
    t1 = time.time()
    c.sendline("bye")
    c.expect(pexpect.EOF)
    c.close()
    log.close()
    print(f"python run with Ctrl+C: {t1 - t0:.0f} s; shell exit status {c.exitstatus}")
    # compare (the Python folder also has the shell's runflag.txt / FChistory.txt / shell.log)
    for f in ("runflag.txt", "FChistory.txt"):
        try:
            os.remove(os.path.join(py, f))
        except FileNotFoundError:
            pass
    # The saved .EPSR.inp keeps the niter the user asked for (3); that is intended. Compare
    # everything else, and check that niter is the only difference in that file.
    inp_py = open(os.path.join(py, INP), encoding="latin-1").read().split("\n")
    inp_fo = open(os.path.join(fo, INP), encoding="latin-1").read().split("\n")
    diff = [(a, b) for a, b in zip(inp_fo, inp_py) if a != b]
    niter_only = len(inp_py) == len(inp_fo) and len(diff) == 1 and diff[0][0].startswith("niter") \
        and diff[0][1].replace("niter       3", "niter       1") == diff[0][0]
    print("saved .EPSR.inp differs only by the requested niter:", niter_only)
    shutil.copy2(os.path.join(fo, INP), os.path.join(work, "fortran_inp"))
    shutil.copy2(os.path.join(py, INP), os.path.join(work, "python_inp"))
    shutil.copy2(os.path.join(fo, INP), os.path.join(py, INP))
    p = subprocess.run([sys.executable, os.path.join(HERE, "compare_outputs.py"), fo, py, "--quiet"],
                       capture_output=True, text=True)
    out = [ln for ln in p.stdout.splitlines() if "run.log" not in ln]
    print("\n".join(out))
    leftovers = subprocess.run(["pgrep", "-f", f"epsr26.*{work}"], capture_output=True, text=True).stdout.split()
    print("residual processes:", leftovers or "none")
    ok = p.returncode == 0 and niter_only and c.exitstatus == 0 and not leftovers
    print("GRACEFUL STOP:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
