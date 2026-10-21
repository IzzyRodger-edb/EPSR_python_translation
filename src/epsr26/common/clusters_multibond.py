"""Generated from fortran_src/clusters_multibond.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``clusters_multibond.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401
from .epsr_dimension import mfile  # noqa: F401

mdist = 30  # mfile
mz = 6000
mpix = 201
maxncbonds = 10

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of clusters_multibond.inc."""
    global titler, titledat, ext, atom1, atom2, atom1s, atom2s, zstep
    global zhalf, z, gofz, gzt, gzt2, fract2, voidmax, rmins
    global rmaxs, damin, damax, nz, nsumt, ltyp, ltyp1, ndist
    global npix, numbdefs
    # common /clusters_c/
    titler = ''
    titledat = char_array(mdist)
    ext = ''
    atom1 = char_array(mdist, maxncbonds)
    atom2 = char_array(mdist, maxncbonds)
    atom1s = char_array(mdist, maxncbonds)
    atom2s = char_array(mdist, maxncbonds)
    # common /clusters_r/
    zstep = f32(0.0)
    zhalf = f32(0.0)
    z = real_array(mz)
    gofz = real_array(mz, mdist)
    gzt = real_array(mz, mdist)
    gzt2 = real_array(mz, mdist)
    fract2 = real_array(mdist)
    voidmax = f32(0.0)
    rmins = real_array(mdist, maxncbonds)
    rmaxs = real_array(mdist, maxncbonds)
    damin = real_array(mdist)
    damax = real_array(mdist)
    # common /clusters_i/
    nz = 0
    nsumt = 0
    ltyp = int_array(mdist, maxncbonds)
    ltyp1 = int_array(mdist, maxncbonds)
    ndist = 0
    npix = 0
    numbdefs = int_array(mdist)


reset()
