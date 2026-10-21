"""ran1/initran1/saveran1 (library/auxiliary_routines.py) against gfortran.

The reference (data/ran1_reference.txt) was produced by data/ran1_driver.f90 linked with
fortran_src/auxiliary_routines.f; see the driver header for the build command.  Each
block is 1000 deviates (bit patterns of the real*4 values) from one seed, then the
state line written by saveran1; the last block continues from a state restored by
initran1.
"""

import pathlib

import numpy as np

from epsr26 import common
from epsr26.common import ran1com
from epsr26.fortran import units
from epsr26.library import auxiliary_routines

DATA = pathlib.Path(__file__).parent / "data"


def _blocks():
    lines = (DATA / "ran1_reference.txt").read_text().splitlines()
    blocks = []
    i = 0
    while i < len(lines):
        head = lines[i]
        vals = [int(v) for v in lines[i + 1 : i + 1001]]
        state = lines[i + 1001]
        blocks.append((head, vals, state))
        i += 1002
    return blocks


def _bits(x):
    return int(np.float32(x).view(np.int32))


def _run(n):
    return [_bits(auxiliary_routines.ran1()) for _ in range(n)]


def _state_line(tmp_path):
    f = tmp_path / "state.txt"
    units.open(11, str(f), status="unknown")
    auxiliary_routines.saveran1(11)
    units.close(11)
    return f.read_text().rstrip("\n")


def test_ran1_matches_gfortran(tmp_path):
    common.reset()
    units.reset()
    blocks = _blocks()
    assert len(blocks) == 7
    for head, vals, state in blocks[:-1]:
        seed = int(head.split()[1])
        ran1com.mcall = 0
        ran1com.idum = seed
        ran1com.iy = 0
        ran1com.iv[:] = 0
        assert _run(1000) == vals, head
        assert _state_line(tmp_path) == state, head
    # continue from a state saved and restored with saveran1/initran1
    f = tmp_path / "saved.txt"
    units.open(10, str(f), status="unknown")
    auxiliary_routines.saveran1(10)
    units.close(10)
    ran1com.mcall = 0
    ran1com.idum = 0
    ran1com.iy = 0
    ran1com.iv[:] = 0
    units.open(10, str(f), status="old")
    auxiliary_routines.initran1(10)
    units.close(10)
    head, vals, state = blocks[-1]
    assert head == "continued"
    assert _run(1000) == vals
    assert _state_line(tmp_path) == state


def test_initran1_error_resets(tmp_path):
    common.reset()
    units.reset()
    f = tmp_path / "empty.txt"
    f.write_text("")
    ran1com.mcall = 5
    ran1com.idum = 7
    ran1com.iy = 9
    units.open(10, str(f), status="old")
    auxiliary_routines.initran1(10)
    units.close(10)
    assert (ran1com.mcall, ran1com.idum, ran1com.iy) == (0, 0, 0)
