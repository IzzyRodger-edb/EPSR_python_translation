#!/usr/bin/env python3
"""Compare the Python with the Rust port (EPSRpy epsr-rs) on one epsr step of a test case.

usage: tools/compare_with_rust.py START_DIR INPUT_FILE WORK_DIR [--rust PATH]

For the steps where the Python deliberately differs from the Fortran (the bug fixes of
docs/FORTRAN_BUGS.md section 1, which the Rust port has too), the Rust is the reference:
both programs run one `epsr` step from a copy of START_DIR (e.g. a run_ex_py.py
`pristine/` or `before/` folder) and every file is compared.
"""

import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def main():
    start, inp, work = sys.argv[1:4]
    rust = os.path.expanduser("~/EPSRpy/rust/target/release/epsr-rs")
    if "--rust" in sys.argv:
        rust = sys.argv[sys.argv.index("--rust") + 1]
    if os.path.exists(work):
        shutil.rmtree(work)
    r, p = os.path.join(work, "rust"), os.path.join(work, "python")
    shutil.copytree(start, r)
    shutil.copytree(start, p)
    with open(os.path.join(work, "rust.log"), "w") as fo:
        subprocess.run([rust, "./", "epsr", inp, "--quiet"], cwd=r, stdout=fo, stderr=subprocess.STDOUT)
    env = dict(os.environ)
    env["PYTHONPATH"] = os.path.join(ROOT, "src") + os.pathsep + env.get("PYTHONPATH", "")
    with open(os.path.join(work, "python.log"), "w") as fo:
        subprocess.run([sys.executable, "-m", "epsr26.run", "epsr", "./", "epsr", inp], cwd=p, stdout=fo,
                       stderr=subprocess.STDOUT, env=env)
    res = subprocess.run([sys.executable, os.path.join(HERE, "compare_outputs.py"), r, p, "--quiet"])
    return res.returncode


if __name__ == "__main__":
    sys.exit(main())
