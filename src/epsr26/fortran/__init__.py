"""Fortran compatibility layer: no EPSR physics, only the parts of the Fortran language
and of the gfortran runtime that EPSR's results depend on.

==============  ==========================================================
module          purpose
==============  ==========================================================
``fmath``       real*4 type, libm intrinsics, powers, nint, integer mod/div
``fformat``     FORMAT interpreter (``fmt``) and list-directed output (``lw``)
``listio``      list-directed and formatted input
``units``       I/O units: open/close/read/write/backspace/rewind
``text``        fixed-length character helpers
``arrays``      1-based column-major arrays
``runtime``     getarg/iargc, clocks, date
``stop``        FortranStop, FortranEOF, FortranIOError
==============  ==========================================================

Translated modules usually start with::

    from epsr26.fortran import *          # f32, r4, fmt, lw, Ch, units, fstr, ...
"""

from . import runtime, units
from .arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array, resize
from .fformat import I8, Ch, fmt, lw
from .fmath import (
    PI4ATAN,
    acosf,
    aint,
    amod,
    anint,
    asinf,
    atan2d_intrinsic,
    atan2f,
    atanf,
    cosd_intrinsic,
    cosf,
    coshf,
    dsign,
    dsqrt,
    erfcf,
    erff,
    expf,
    f32,
    i32,
    idiv,
    ifix,
    imod,
    isign,
    log10f,
    logf,
    nint,
    pow2,
    pow3,
    pow4,
    powf,
    powi,
    powi_d,
    r4,
    sign,
    sind_intrinsic,
    sinf,
    sinhf,
    sqrtf,
    tanf,
    tanhf,
    wrap32,
)
from .listio import parse_dble, parse_int, parse_real, read_internal
from .stop import FortranEOF, FortranIOError, FortranStop, iostat_of
from .text import adjustl, fstr, findex, findex_padded, len_trim, pad, setsub, substr

__all__ = [n for n in dir() if not n.startswith("_")]
