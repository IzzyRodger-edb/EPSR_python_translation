"""Generated from fortran_src/weights.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``weights.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401
from .epsr_dimension import mfile, mpairs  # noqa: F401

mqinv = 15

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of weights.inc."""
    global cwt, ciwt, cwtpot, cwtt, carray, carrayold, cwttq, cwtav
    global nxfile, xfile, nqinv, qinv
    # common /weights_d/
    cwt = dble_array(mfile, mpairs)
    ciwt = dble_array(mfile, mpairs)
    cwtpot = dble_array(mpairs, mpairs)
    cwtt = dble_array(mpairs, 850)
    carray = dble_array(850)
    carrayold = dble_array(850, 850)
    cwttq = dble_array(mqinv, mpairs, 850)
    # common /weights_r/
    cwtav = f32(0.0)
    # common /weights_i/
    nxfile = 0
    xfile = int_array(mfile)
    nqinv = 0
    qinv = int_array(mqinv)


reset()
