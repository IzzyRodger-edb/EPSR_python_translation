"""Generated from fortran_src/epsr_dimension.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``epsr_dimension.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401

mqd = 5000
mpz = 3000
mpts = 2001
mfile = 30
mcomp = 40
mpairs = 820  # ((mcomp*(mcomp+1))/2)
matom = 150000
near = 200
mmol = 50000

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of epsr_dimension.inc."""
    pass


reset()
