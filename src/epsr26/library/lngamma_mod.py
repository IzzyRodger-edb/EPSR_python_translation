"""Translation of fortran_src/lngamma_mod.f90: routines of ``module lngamma_mod``.

The module data (``mgamma``, ``gamln``) is in ``common/lngamma_mod.py``.  Callers that
also need the data module must import one of the two under another name, e.g.
``from ..common import lngamma_mod as lngamma_mod_data``.
"""

import math

from ..common import lngamma_mod
from ..fortran import dble_array
from ..trace import routine


@routine
def fctrln(nterms):
    """lngamma_mod.f90:fctrln -- tabulate gamln(i) = log((i-1)!) for i = 1..nterms+1."""
    # integer nterms,ncalc,i
    # real(kind=8) g,a
    if lngamma_mod.gamln is not None:
        lngamma_mod.gamln = None
    lngamma_mod.mgamma = nterms + 1
    lngamma_mod.gamln = dble_array(lngamma_mod.mgamma)
    ncalc = lngamma_mod.mgamma
    g = 0.0
    a = 0.0
    gamln = lngamma_mod.gamln
    for i in range(1, ncalc + 1):
        gamln[i] = g
        a = a + 1.0
        g = g + math.log(a)
    return


@routine
def lngamma(ic):
    """lngamma_mod.f90:lngamma -- log(ic!) (real*8) from the gamln table.

    Returns lngamma.
    """
    if ic > lngamma_mod.mgamma or lngamma_mod.gamln is None:
        fctrln(ic)
    lngamma_ = float(lngamma_mod.gamln[ic + 1])
    return lngamma_
