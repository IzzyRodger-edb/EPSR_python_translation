"""Generated from fortran_src/directory_name.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``directory_name.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of directory_name.inc."""
    global directory_name, home_directory, nchardir, ncharhomedir
    # common /dir_c/
    directory_name = ''
    home_directory = ''
    # common /dir_i/
    nchardir = 0
    ncharhomedir = 0


reset()
