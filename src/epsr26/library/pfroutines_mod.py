"""Translation of fortran_src/pfroutines_mod.f90: routines of ``module pfroutines``.

The module data (``rratio``, ``intenprobfunction``, both ``(0:nr)``) is in
``common/pfroutines.py``.

The file has no ``implicit none``: ``r``, ``arg`` (getpfprobability) and ``diff``
(getrratio) are implicit real*4.
"""

import ctypes
import ctypes.util
import math

from ..common import pfroutines
from ..fortran import expf, f32, ifix, logf, lw, r4, real_array, sqrtf, units
from ..trace import routine

# gfortran's GAMMA of a real*8 is libm tgamma; Python's math.gamma is CPython's own
# Lanczos implementation and need not agree to the last bit, so call libm directly.
_libm = ctypes.CDLL(ctypes.util.find_library("m"))
_libm.tgamma.restype = ctypes.c_double
_libm.tgamma.argtypes = [ctypes.c_double]


def _tgamma(x):
    """``gamma(x)`` for real*8 x (glibc tgamma)."""
    return float(_libm.tgamma(float(x)))


@routine
def getpfprobability(rfract):
    """pfroutines_mod.f90:getpfprobability -- cumulative Poisson-like distribution of r ratios."""
    # real, intent(in) :: rfract
    nr = 50  # Must be greater than 1
    if pfroutines.intenprobfunction is None or pfroutines.intenprobfunction.shape[0] < nr + 1:
        # size(intenprobfunction,1) of an (0:nr) array is shape[0]
        if pfroutines.intenprobfunction is not None:
            pfroutines.rratio = None
            pfroutines.intenprobfunction = None
        pfroutines.rratio = real_array((0, nr))
        pfroutines.intenprobfunction = real_array((0, nr))
    rratio = pfroutines.rratio
    intenprobfunction = pfroutines.intenprobfunction
    rsigma = max((r4(1.0) / (rfract * rfract) - r4(1.0)), f32(0.0))  # Keep the widths sensible: rfract cannot be > 1.
    gammacons = lngamma(rsigma)  # log of gamma(rsigma+1)
    # Set up the probability function - assumed Poisson in form for r > 0.
    fullwidth = rsigma + r4(1.0) + r4(10.0) * sqrtf(rsigma + r4(1.0))  # Full width of Poisson function in normalised units = position of max + 10 stds.
    rstep = fullwidth / f32(nr)
    units.write(6, lw("rsigma, gammacons, fullwidth, rstep: ", rsigma, gammacons, fullwidth, rstep))
    rratio[0] = 0.0
    intenprobfunction[0] = 0.0
    for ir in range(1, nr + 1):
        r = f32(ir) * rstep
        arg = rsigma * logf(r) - r - gammacons
        if abs(arg) < 30.0:
            probfunction = expf(arg)
        else:
            probfunction = f32(0.0)
        rratio[ir] = r / (rsigma + r4(1.0))
        intenprobfunction[ir] = intenprobfunction[ir - 1] + probfunction * rstep  # /(rratio(ir)**3)
    # Final renormalisation to correct round-off errors
    for ir in range(1, nr + 1):
        intenprobfunction[ir] = intenprobfunction[ir] / intenprobfunction[nr]
    return


@routine
def lngamma(argin):
    """pfroutines_mod.f90:lngamma -- log of gamma(argin+1.0) (real*4; 0 if argin < 0).

    Returns lngamma.
    """
    # real argin,lngamma
    # real*8 arg,gammaval,sumlngamma
    lngamma_ = f32(0.0)
    if argin < 0.0:
        return lngamma_  # return 0 if illegal argument
    narg = ifix(argin)  # Integer part of arg
    arg = float(argin - f32(narg) + r4(1.0))
    gammaval = _tgamma(arg)
    sumlngamma = math.log(gammaval)
    #      write(6,*) 'lngamma> ',argin,narg,arg,gammaval,sumlngamma
    if narg > 0:
        for i in range(1, narg + 1):
            sumlngamma = sumlngamma + math.log(arg)
            arg = arg + 1.0
    #      write(6,*) 'lngamma> ',argin,narg,arg,gammaval,sumlngamma
    lngamma_ = f32(sumlngamma)
    return lngamma_


@routine
def getrratio():
    """pfroutines_mod.f90:getrratio -- sample an r ratio from intenprobfunction with ran1.

    Returns getrratio.
    """
    from . import auxiliary_routines  # imported here: auxiliary_routines.f belongs to another package

    rratio = pfroutines.rratio
    intenprobfunction = pfroutines.intenprobfunction
    ranrratio = auxiliary_routines.ran1()
    i = 0
    found = False
    # size(rratio) of an (0:nr) array is shape[0]
    while i < rratio.shape[0] - 1 and not found:
        i = i + 1
        found = ranrratio < intenprobfunction[i]
    diff = intenprobfunction[i] - intenprobfunction[i - 1]
    grad = f32(0.0)
    if diff > 0.0:
        grad = (rratio[i] - rratio[i - 1]) / diff
    getrratio_ = rratio[i - 1] + (ranrratio - intenprobfunction[i - 1]) * grad
    return getrratio_
