"""Translation of fortran_src/molecom_mod.f90: routines of ``module molecom_mod``.

The module data is in ``common/molecom_mod.py`` (imported here as ``molecom_mod``).
A caller that needs both must alias one of them, e.g.
``from . import molecom_mod as molecom_mod_routines``.
"""

from ..common import molecom_mod
from . import reallocation_routines
from ..trace import routine


@routine
def testsize_rotgrps(newsizein=None):
    """molecom_mod.f90:testsize_rotgrps -- grow the rotational-group arrays to newsize.

    newsizein is ``optional``; when absent the size is ngrrot.
    """
    # use reallocation_routines
    if newsizein is not None:
        newsize = newsizein
    else:
        newsize = molecom_mod.ngrrot
    while newsize > molecom_mod.ngrrotmax:
        molecom_mod.ngrrotmax = newsize
        # len(rotgrp) = 3 (character(len=3))
        molecom_mod.rotgrp = reallocation_routines.reallocate1d_c(molecom_mod.rotgrp, 3, molecom_mod.ngrrotmax)
        molecom_mod.labaxi = reallocation_routines.reallocate2d_i(molecom_mod.labaxi, molecom_mod.ngrrotmax, 2)
        molecom_mod.nrot = reallocation_routines.reallocate1d_i(molecom_mod.nrot, molecom_mod.ngrrotmax)
        molecom_mod.rotf = reallocation_routines.reallocate1d_i(molecom_mod.rotf, molecom_mod.ngrrotmax)
        molecom_mod.rotl = reallocation_routines.reallocate1d_i(molecom_mod.rotl, molecom_mod.ngrrotmax)
        molecom_mod.mrot = reallocation_routines.reallocate1d_i(molecom_mod.mrot, molecom_mod.ngrrotmax)
    return


@routine
def testsize_labrot(newsizein=None):
    """molecom_mod.f90:testsize_labrot -- grow labrot to newsize (default nrottot)."""
    # use reallocation_routines
    if newsizein is not None:
        newsize = newsizein
    else:
        newsize = molecom_mod.nrottot
    # We only reallocate labrot if the new size is bigger than the old size
    while newsize > molecom_mod.nrottotmax:
        molecom_mod.nrottotmax = newsize
        molecom_mod.labrot = reallocation_routines.reallocate1d_i(molecom_mod.labrot, molecom_mod.nrottotmax)
    return
