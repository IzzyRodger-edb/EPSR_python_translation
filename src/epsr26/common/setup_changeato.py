"""Generated from fortran_src/setup_changeato.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``setup_changeato.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of setup_changeato.inc."""
    global bond, label, addbond, removebond, density, temperature
    # common /setup_changeato_c/
    bond = ''
    label = ''
    addbond = ''
    removebond = ''
    # common /setup_changeato_r/
    density = f32(0.0)
    temperature = f32(0.0)


reset()
