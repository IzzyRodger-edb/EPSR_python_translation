"""Generated from fortran_src/chains.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``chains.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401
from .epsr_dimension import mmol  # noqa: F401

mdist = 20
mz = 200
mchain = 200
mbond = 30
mring = 25

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of chains.inc."""
    global titler, titledat, ext, atom1, atom2, atom3, dff1, dff2
    global dfs1, dfs2, dft1, dft2, chainlen, zstep, zhalf, z
    global gofz, gzt, gzt2, ang213min, ang213max, nz, nsumt, lftype
    global lstype, lttype, ndist, ndist3, neighb, nneighb, icount, nrlab
    global rneighb, spvalue, spvaluet, spseq, sptest, spgrad, spmin, spmax
    global rptest, rneighb1, rneighb2, nwordsave, nodd, sptrue1, sptrue2, carryon
    # common /chains_c/
    titler = ''
    titledat = char_array(mdist)
    ext = ''
    atom1 = char_array(mdist, mdist)
    atom2 = char_array(mdist, mdist)
    atom3 = char_array(mdist, mdist)
    # common /chains_r/
    dff1 = real_array(mdist, mdist)
    dff2 = real_array(mdist, mdist)
    dfs1 = real_array(mdist, mdist)
    dfs2 = real_array(mdist, mdist)
    dft1 = real_array(mdist, mdist)
    dft2 = real_array(mdist, mdist)
    chainlen = real_array(mdist, 2)
    zstep = f32(0.0)
    zhalf = f32(0.0)
    z = real_array(mz)
    gofz = real_array(mz, mdist)
    gzt = real_array(mz, mdist)
    gzt2 = real_array(mz, mdist)
    ang213min = real_array(mdist, mdist)
    ang213max = real_array(mdist, mdist)
    # common /chains_i/
    nz = 0
    nsumt = 0
    lftype = int_array(mdist)
    lstype = int_array(mdist)
    lttype = int_array(mdist)
    ndist = 0
    ndist3 = 0
    neighb = int_array(mmol, mbond)
    nneighb = int_array(mmol)
    icount = 0
    nrlab = int_array((0, mz))
    rneighb = int_array(mz)
    spvalue = int_array(mmol, mring)
    spvaluet = int_array(mmol)
    spseq = int_array((0, mchain))
    sptest = 0
    spgrad = 0
    spmin = 0
    spmax = 0
    rptest = 0
    rneighb1 = int_array(mmol)
    rneighb2 = int_array(mmol)
    nwordsave = int_array(mdist)
    # common /chains_l/
    nodd = False
    sptrue1 = False
    sptrue2 = False
    carryon = logical_array(mmol)


reset()
