"""Generated from fortran_src/coord.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``coord.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401

mdist = 100
mz = 1000

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of coord.inc."""
    global titler, titledat, ext, atom1, atom2, atomc, atomz, atomr
    global damin, damax, rminc, rmaxc, coordn, zstep, zhalf, z
    global gofz, gzt, gzt2, fract2, nz, nsumt, ltyp1, ltyp2
    global ltypc, ndist, ndist2, dorhoz, dorhor
    # common /coord_c/
    titler = ''
    titledat = char_array(mdist)
    ext = ''
    atom1 = char_array(mdist)
    atom2 = char_array(mdist)
    atomc = ''
    atomz = ''
    atomr = ''
    # common /coord_r/
    damin = real_array(mdist)
    damax = real_array(mdist)
    rminc = real_array(mdist)
    rmaxc = real_array(mdist)
    coordn = real_array(mdist, 2)
    zstep = f32(0.0)
    zhalf = f32(0.0)
    z = real_array(mz)
    gofz = real_array(mz, mdist)
    gzt = real_array(mz, mdist)
    gzt2 = real_array(mz, mdist)
    fract2 = f32(0.0)
    # common /coord_i/
    nz = 0
    nsumt = 0
    ltyp1 = int_array(mdist)
    ltyp2 = int_array(mdist)
    ltypc = 0
    ndist = 0
    ndist2 = 0
    # common /coord_l/
    dorhoz = False
    dorhor = False


reset()
