#!/usr/bin/env python3
"""For a test case where the Python deliberately differs from the Fortran (a bug fixed as in
docs/FORTRAN_BUGS.md section 1), check every epsr pass against the Rust port instead.

usage: tools/check_fixed_case.py CASE_DIR WORK_DIR

The Fortran chain (fortran_build/epsr) provides the folder before each pass, exactly as in
tools/run_ex_py.py resume mode; from each such folder the Rust (with the same fixes) and
the Python run one pass and all their files are compared.  Exit 0 if all passes agree.
"""

import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import run_ex_py as rx  # noqa: E402


def main():
    case, work = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    if os.path.exists(work):
        shutil.rmtree(work)
    chain = os.path.join(work, "fortran_chain")
    rx.copy_inputs(case, chain)
    steps = [s for s in rx.read_steps(case) if s[0] == "run"]
    ok_all = True
    for n, (_, words, _) in enumerate(steps, 1):
        prog, args = words[0], words[1:]
        assert prog == "epsr", "only epsr steps are supported"
        before = os.path.join(work, f"before{n}")
        shutil.copytree(chain, before)
        res = subprocess.run([sys.executable, os.path.join(HERE, "compare_with_rust.py"), before, args[0],
                              os.path.join(work, f"pass{n}")], capture_output=True, text=True)
        verdict = "IDENTICAL to Rust" if res.returncode == 0 else "DIFFERS from Rust"
        ok_all &= res.returncode == 0
        print(f"pass {n}: Python {verdict}", flush=True)
        if res.returncode:
            print(res.stdout[-3000:])
        subprocess.run(f"ulimit -s unlimited; OMP_NUM_THREADS=1 {ROOT}/fortran_build/epsr ./ epsr {' '.join(args)}"
                       " > /dev/null 2>&1", shell=True, cwd=chain)
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
