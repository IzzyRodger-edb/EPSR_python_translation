"""Generated from fortran_src/version.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``version.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401

date = '26: 2017-03-10'
ilendate = 14

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of version.inc."""
    global version, ilenversion
    # common /version_c/
    version = ''
    # common /version_i/
    ilenversion = 0


reset()
