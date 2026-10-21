"""Generated from fortran_src/lattice.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``lattice.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of lattice.inc."""
    global alat, blat, clat, phib, thetac, phic, phibr, phicr
    global thetacr, alatx, blatx, clatx, blaty, claty, clatz, bxax
    global cxax, cyby, alatc, blatc, clatc, twopi, cubic, absolute
    # common /lattice_r/
    alat = f32(0.0)
    blat = f32(0.0)
    clat = f32(0.0)
    phib = f32(0.0)
    thetac = f32(0.0)
    phic = f32(0.0)
    phibr = f32(0.0)
    phicr = f32(0.0)
    thetacr = f32(0.0)
    alatx = f32(0.0)
    blatx = f32(0.0)
    clatx = f32(0.0)
    blaty = f32(0.0)
    claty = f32(0.0)
    clatz = f32(0.0)
    bxax = f32(0.0)
    cxax = f32(0.0)
    cyby = f32(0.0)
    alatc = real_array(3)
    blatc = real_array(3)
    clatc = real_array(3)
    twopi = f32(0.0)
    # common /lattice_l/
    cubic = False
    absolute = False


reset()
