"""Generated from fortran_src/sharm_sdf.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``sharm_sdf.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401
from .sharm import mcentre  # noqa: F401

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of sharm_sdf.inc."""
    global naxis1, naxis2, pi, piconv, a1, b1, c1, pdif
    global perp1, alab, blab, clab, at, bt, ct, pt
    global phim1, thetam1, chim1, phim2, thetam2, chim2, thetallab, phillab
    global phispread, thetaspread, phinew, thetanew, chinew, atomno
    # common /sharm_sdf_c/
    naxis1 = ''
    naxis2 = ''
    # common /sharm_sdf_r/
    pi = f32(0.0)
    piconv = f32(0.0)
    a1 = real_array(3)
    b1 = real_array(3)
    c1 = real_array(3)
    pdif = real_array(3)
    perp1 = real_array(3)
    alab = real_array(3)
    blab = real_array(3)
    clab = real_array(3)
    at = real_array(3)
    bt = real_array(3)
    ct = real_array(3)
    pt = real_array(3)
    phim1 = f32(0.0)
    thetam1 = f32(0.0)
    chim1 = f32(0.0)
    phim2 = f32(0.0)
    thetam2 = f32(0.0)
    chim2 = f32(0.0)
    thetallab = f32(0.0)
    phillab = f32(0.0)
    phispread = f32(0.0)
    thetaspread = f32(0.0)
    phinew = f32(0.0)
    thetanew = f32(0.0)
    chinew = f32(0.0)
    # common /sharm_sdf_i/
    atomno = int_array(mcentre)


reset()
