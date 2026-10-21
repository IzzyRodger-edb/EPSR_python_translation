"""Generated from fortran_src/ran1com.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``ran1com.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401

ia = 16807
im = 2147483647
am = f32(4.656612873077393e-10)  # 1./im
iq = 127773
ir = 2836
ntab = 32
ndiv = 67108864  # 1+(im-1)/ntab
eps = f32(1.199999957179898e-07)  # 1.2e-7
rnmx = f32(0.9999998807907104)  # 1.-eps

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of ran1com.inc."""
    global mcall, idum, iy, iv
    # common /ran1com/
    mcall = 0
    idum = 0
    iy = 0
    iv = int_array(ntab)


reset()
