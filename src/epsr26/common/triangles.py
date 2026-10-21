"""Generated from fortran_src/triangles.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``triangles.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401
from .epsr_dimension import mfile  # noqa: F401

mdist = 30  # mfile
mz = 200

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of triangles.inc."""
    global titler, titledat, ext, atom1, atom2, atom3, atomc, atomz
    global atomr, damin, damax, rminc, rmaxc, zstep, zhalf, z
    global gofz, gzt, gzt2, nz, nsumt, ltyp, ltyp1, ltyp2
    global ltyp3, ltypc, ndist, dorhoz, dorhor
    # common /triangles_c/
    titler = ''
    titledat = char_array(mdist)
    ext = ''
    atom1 = char_array(mdist)
    atom2 = char_array(mdist)
    atom3 = char_array(mdist)
    atomc = ''
    atomz = ''
    atomr = ''
    # common /triangles_r/
    damin = real_array(mdist, 2)
    damax = real_array(mdist, 2)
    rminc = real_array(mdist)
    rmaxc = real_array(mdist)
    zstep = f32(0.0)
    zhalf = f32(0.0)
    z = real_array(mz)
    gofz = real_array(mz, mdist)
    gzt = real_array(mz, mdist)
    gzt2 = real_array(mz, mdist)
    # common /triangles_i/
    nz = 0
    nsumt = 0
    ltyp = int_array(mdist)
    ltyp1 = int_array(mdist)
    ltyp2 = int_array(mdist)
    ltyp3 = int_array(mdist)
    ltypc = 0
    ndist = 0
    # common /triangles_l/
    dorhoz = False
    dorhor = False


reset()
