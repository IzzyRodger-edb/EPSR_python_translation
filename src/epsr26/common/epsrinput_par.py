"""Generated from fortran_src/epsrinput_par.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``epsrinput_par.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of epsrinput_par.inc."""
    global nisum, nesum, nesumt, nsum, nsums, roverlap, rminfac, rminpt
    global rmaxpt, widthptc, rbroad, expecf, expecfpwr, rcharge, ensum, ensum2
    global prsum, prsum2, lorchqwidth, lorchqwidthq, rminbragg, rmaxgr
    # common /epsrinput_par/
    nisum = 0
    nesum = 0
    nesumt = 0
    nsum = 0
    nsums = 0
    roverlap = f32(0.0)
    rminfac = f32(0.0)
    rminpt = f32(0.0)
    rmaxpt = f32(0.0)
    widthptc = f32(0.0)
    rbroad = f32(0.0)
    expecf = f32(0.0)
    expecfpwr = f32(0.0)
    rcharge = f32(0.0)
    ensum = f32(0.0)
    ensum2 = f32(0.0)
    prsum = f32(0.0)
    prsum2 = f32(0.0)
    lorchqwidth = f32(0.0)
    lorchqwidthq = f32(0.0)
    rminbragg = f32(0.0)
    rmaxgr = f32(0.0)


reset()
