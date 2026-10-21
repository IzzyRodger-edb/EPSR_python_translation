"""Fortran module data of ``module pfroutines`` (fortran_src/pfroutines_mod.f90).

Hand-written (TRANSLATION_GUIDE 5.2).  The routines (``getpfprobability``, ``lngamma``,
``getrratio``) are in ``library/pfroutines_mod.py``.

Both arrays are allocated ``(0:nr)``, i.e. with ``real_array((0, nr))``: index 0 is
element 0, and the Fortran ``size(a)`` is ``a.shape[0]``.
"""

import sys

from . import register


def reset():
    """(Re)create the module variables of pfroutines."""
    global rratio, intenprobfunction
    rratio = None  # real, dimension(:), allocatable
    intenprobfunction = None  # real, dimension(:), allocatable


register(sys.modules[__name__])
reset()
