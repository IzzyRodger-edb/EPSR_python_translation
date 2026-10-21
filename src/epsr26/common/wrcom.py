"""Generated from fortran_src/wrcom.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``wrcom.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401

nwrmax = 100
nform = 201  # (nwrmax*2+1)

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of wrcom.inc."""
    global xo, yo, eo, ext, titleo
    # common /wrcom/
    xo = real_array(10000)
    yo = real_array(10000, nwrmax)
    eo = real_array(10000, nwrmax)
    # common /wrcomc/
    ext = ''
    titleo = char_array(nwrmax)


reset()
