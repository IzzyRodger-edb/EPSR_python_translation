"""Generated from fortran_src/partials.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``partials.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401
from .epsr_dimension import mpairs  # noqa: F401

mdist = 820  # mpairs
mz = 1000

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of partials.inc."""
    global titler, titledat, ext, atom1, atom2, atomc, atomz, atomr
    global zstep, zhalf, z, gofz, gzt, gzt2, rnormalise, nz
    global nsumt, ltyp1, ltyp2, ltypc, ndist, dorhoz, dorhor
    # common /partials_c/
    titler = ''
    titledat = char_array(mdist)
    ext = ''
    atom1 = char_array(mdist)
    atom2 = char_array(mdist)
    atomc = ''
    atomz = ''
    atomr = ''
    # common /partials_r/
    zstep = f32(0.0)
    zhalf = f32(0.0)
    z = real_array(mz)
    gofz = real_array(mz, mdist)
    gzt = real_array(mz, mdist)
    gzt2 = real_array(mz, mdist)
    rnormalise = f32(0.0)
    # common /partials_i/
    nz = 0
    nsumt = 0
    ltyp1 = int_array(mdist)
    ltyp2 = int_array(mdist)
    ltypc = 0
    ndist = 0
    # common /partials_l/
    dorhoz = False
    dorhor = False


reset()
