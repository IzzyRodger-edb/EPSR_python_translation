"""Generated from fortran_src/torangles.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``torangles.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401
from .epsr_dimension import mfile  # noqa: F401

mdist = 30  # mfile
mz = 200

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of torangles.inc."""
    global titler, titledat, ext, atom1, atom2, atom3, damin, damax
    global zstep, zhalf, z, gofz, gzt, gzt2, nz, nsumt
    global ltyp, ltyp1, ltyp2, ltyp3, ndist, nfirst
    # common /torangles_c/
    titler = ''
    titledat = char_array(mdist)
    ext = ''
    atom1 = char_array(mdist)
    atom2 = char_array(mdist)
    atom3 = char_array(mdist)
    # common /torangles_r/
    damin = real_array(mdist, 2)
    damax = real_array(mdist, 2)
    zstep = f32(0.0)
    zhalf = f32(0.0)
    z = real_array(mz)
    gofz = real_array(mz, mdist)
    gzt = real_array(mz, mdist)
    gzt2 = real_array(mz, mdist)
    # common /torangles_i/
    nz = 0
    nsumt = 0
    ltyp = int_array(mdist)
    ltyp1 = int_array(mdist)
    ltyp2 = int_array(mdist, 2)
    ltyp3 = int_array(mdist, 2)
    ndist = 0
    nfirst = int_array(mdist)


reset()
