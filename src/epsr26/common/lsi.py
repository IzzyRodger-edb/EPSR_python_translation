"""Generated from fortran_src/lsi.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``lsi.<name>``.

Declared here but in no COMMON block (locals of each including routine, see
``LOCALS``): bneighb.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401
from .epsr_dimension import mmol  # noqa: F401

mdist = 30
mz = 200
mchain = 200
mbond = 20
mring = 25

LOCALS = {'bneighb': 'logical(mmol, mbond)'}


def reset():
    """(Re)create the COMMON storage of lsi.inc."""
    global titler, titledat, ext, atom1, atom2, atom3, dff1, dff2
    global dfs1, dfs2, dft1, dft2, chainlen, zstep, zhalf, z
    global gofz, gzt, gzt2, ang213min, ang213max, lsimin, lsimax, rtrimax
    global dneighb, nz, nsumt, lftype, lstype, lttype, ndist, ndist3
    global neighb, nneighb, ia1min, ia2min, iwmin, nwordsave
    # common /lsi_c/
    titler = ''
    titledat = char_array(mdist)
    ext = ''
    atom1 = char_array(mdist, mdist)
    atom2 = char_array(mdist, mdist)
    atom3 = char_array(mdist, mdist)
    # common /lsi_r/
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
    lsimin = f32(0.0)
    lsimax = f32(0.0)
    rtrimax = f32(0.0)
    dneighb = real_array(mmol, mbond)
    # common /lsi_i/
    nz = 0
    nsumt = 0
    lftype = int_array(mdist)
    lstype = int_array(mdist)
    lttype = int_array(mdist)
    ndist = 0
    ndist3 = 0
    neighb = int_array(mmol, mbond)
    nneighb = int_array(mmol)
    ia1min = int_array(mmol, mbond)
    ia2min = int_array(mmol, mbond)
    iwmin = int_array(mmol, mbond)
    nwordsave = int_array(mdist)


reset()
