"""Generated from fortran_src/sharm.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``sharm.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401

mcentre = 20
mcentre2 = 190  # mcentre*(mcentre-1)/2

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of sharm.inc."""
    global titler, ext, nzt1, nzt2, nyt1, nyt2, vol, zstep
    global zhalf, zmax, zmaxd, distfirst, fractfirst, distsecond, fractsecond, distcentre
    global fractcentre, xtemp, ytemp, ztemp, nz, nsumt, ncoeffs, ncent1
    global ncent2, ncz1, ncy1, ncz2, ncy2, lztyp1, lztyp2, lytyp1
    global lytyp2, lctyp1, lctyp2, lctyp, natofirst, npairfirst, typefirst, natosecond
    global npairsecond, typesecond, natomol, npairmol
    # common /sharm_c/
    titler = ''
    ext = ''
    nzt1 = ''
    nzt2 = ''
    nyt1 = ''
    nyt2 = ''
    # common /sharm_r/
    vol = f32(0.0)
    zstep = f32(0.0)
    zhalf = f32(0.0)
    zmax = f32(0.0)
    zmaxd = f32(0.0)
    distfirst = real_array(mcentre2)
    fractfirst = f32(0.0)
    distsecond = real_array(mcentre2)
    fractsecond = f32(0.0)
    distcentre = real_array(mcentre2)
    fractcentre = f32(0.0)
    xtemp = real_array(mcentre)
    ytemp = real_array(mcentre)
    ztemp = real_array(mcentre)
    # common /sharm_i/
    nz = 0
    nsumt = 0
    ncoeffs = 0
    ncent1 = 0
    ncent2 = 0
    ncz1 = 0
    ncy1 = 0
    ncz2 = 0
    ncy2 = 0
    lztyp1 = int_array(mcentre)
    lztyp2 = int_array(mcentre)
    lytyp1 = int_array(mcentre)
    lytyp2 = int_array(mcentre)
    lctyp1 = int_array(mcentre)
    lctyp2 = int_array(mcentre)
    lctyp = int_array(mcentre)
    natofirst = 0
    npairfirst = 0
    typefirst = int_array(mcentre)
    natosecond = 0
    npairsecond = 0
    typesecond = int_array(mcentre)
    natomol = 0
    npairmol = 0


reset()
