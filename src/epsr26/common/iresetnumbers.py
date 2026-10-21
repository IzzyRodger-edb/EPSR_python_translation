"""Generated from fortran_src/iresetnumbers.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``iresetnumbers.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of iresetnumbers.inc."""
    global iresetref, iinitref, iqminref, iroverlap, irminfac, incoeffp, iereqfirst, iereqlast
    # common /iresetnumbers/
    iresetref = 0
    iinitref = 0
    iqminref = 0
    iroverlap = 0
    irminfac = 0
    incoeffp = 0
    iereqfirst = 0
    iereqlast = 0


reset()
