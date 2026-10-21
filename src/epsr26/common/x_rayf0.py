"""Generated from fortran_src/X-rayf0.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``x_rayf0.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401
from .epsr_dimension import mcomp, mpairs, mpts  # noqa: F401

mxfile = 2

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of X-rayf0.inc."""
    global f0x, fm0x, weightx, backx, datapoint, datanorm, normxx, iafformaffx
    # common /xrayf0_r/
    f0x = dble_array(mpts, mcomp, mxfile)
    fm0x = dble_array(mpts, mcomp, mxfile)
    weightx = dble_array(mpts, mpairs, mxfile)
    backx = dble_array(mpts, mxfile)
    # common /xrayf0_i/
    datapoint = int_array(mxfile)
    datanorm = int_array(mxfile)
    normxx = int_array(mxfile)
    iafformaffx = int_array(mxfile)


reset()
