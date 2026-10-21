"""Fortran numeric types and intrinsics with gfortran's exact results.

Types
-----
* ``real`` (real*4)          -> ``numpy.float32`` scalars (``f32``)
* ``real*8``/``double``      -> Python ``float``
* ``integer``                -> Python ``int`` (``numpy.int32`` where it comes out of an array)

Traps (read these before translating):

* ``numpy.int32 * numpy.float32`` gives **float64**.  Fortran converts the integer to
  ``real`` first, so write ``f32(i) * x``.  A plain Python ``int`` is "weak" under NEP 50
  and is converted to float32, but an element of an integer array is not.
* ``float32 op python_float`` gives float32 (the Python float is weak).  That matches a
  Fortran single-precision literal (``x + 0.5``), but NOT a ``real*8`` variable:
  ``real*8 d; d + x`` must be written ``d + float(x)``.
* A constant expression of literals is folded by gfortran in single precision:
  ``1.0/3.0`` -> ``r4(1.0) / r4(3.0)`` (or ``f32(1.0)/f32(3.0)``), never ``f32(1.0/3.0)``.
* ``x**2`` -> ``x*x``; ``x**3`` -> ``x*x*x``; ``x**4`` -> ``pow4(x)``; ``x**n`` with an
  integer *variable* -> ``powi(x, n)`` (libgcc ``__powisf2``); ``x**y`` real -> ``powf``.

Single-precision transcendentals call glibc ``libm`` through ctypes, because NumPy's
float32 ufuncs do not use glibc.  ``sqrt`` is correctly rounded, so ``f32(math.sqrt(x))``
equals ``sqrtf(x)`` exactly and is used for speed.  Double precision uses ``math.*``,
which calls glibc already.
"""

import ctypes
import ctypes.util
import math

import numpy as np

np.seterr(all="ignore")  # Fortran does not trap on overflow / division by zero

f32 = np.float32
i32 = np.int32


def r4(x):
    """A single-precision literal or constant: ``r4(0.1)`` is Fortran ``0.1``."""
    return np.float32(x)


# --------------------------------------------------------------------------- libm
_libm = ctypes.CDLL(ctypes.util.find_library("m"))


def _f1(name):
    fn = getattr(_libm, name)
    fn.restype = ctypes.c_float
    fn.argtypes = [ctypes.c_float]
    return fn


def _f2(name):
    fn = getattr(_libm, name)
    fn.restype = ctypes.c_float
    fn.argtypes = [ctypes.c_float, ctypes.c_float]
    return fn


_expf = _f1("expf")
_logf = _f1("logf")
_log10f = _f1("log10f")
_sinf = _f1("sinf")
_cosf = _f1("cosf")
_tanf = _f1("tanf")
_asinf = _f1("asinf")
_acosf = _f1("acosf")
_atanf = _f1("atanf")
_sinhf = _f1("sinhf")
_coshf = _f1("coshf")
_tanhf = _f1("tanhf")
_erff = _f1("erff")
_erfcf = _f1("erfcf")
_lgammaf = _f1("lgammaf")
_atan2f = _f2("atan2f")
_powf = _f2("powf")
_fmodf = _f2("fmodf")


def expf(x):
    return f32(_expf(x))


def logf(x):
    return f32(_logf(x))


def log10f(x):
    return f32(_log10f(x))


def sinf(x):
    return f32(_sinf(x))


def cosf(x):
    return f32(_cosf(x))


def tanf(x):
    return f32(_tanf(x))


def asinf(x):
    return f32(_asinf(x))


def acosf(x):
    return f32(_acosf(x))


def atanf(x):
    return f32(_atanf(x))


def atan2f(y, x):
    return f32(_atan2f(y, x))


def sinhf(x):
    return f32(_sinhf(x))


def coshf(x):
    return f32(_coshf(x))


def tanhf(x):
    return f32(_tanhf(x))


def erff(x):
    return f32(_erff(x))


def erfcf(x):
    return f32(_erfcf(x))


def lgammaf(x):
    return f32(_lgammaf(x))


def powf(x, y):
    """``x**y`` with a real exponent (both real*4)."""
    return f32(_powf(x, y))


def sqrtf(x):
    """``sqrt`` of a real*4 (correctly rounded, equal to ``sqrtss``)."""
    x = float(x)
    if x < 0.0:
        return f32("nan")
    return f32(math.sqrt(x))


def dsqrt(x):
    x = float(x)
    if x < 0.0:
        return math.nan
    return math.sqrt(x)


def amod(a, b):
    """Fortran ``mod(a, b)`` for reals (exact, = fmodf)."""
    return f32(_fmodf(a, b))


# --------------------------------------------------------------------------- powers
def pow2(x):
    return x * x


def pow3(x):
    return x * x * x


def pow4(x):
    x2 = x * x
    return x2 * x2


def powi(x, m):
    """``x**m`` for real*4 x and integer variable m: port of libgcc ``__powisf2``."""
    m = int(m)
    n = -m if m < 0 else m
    y = x if n % 2 == 1 else f32(1.0)
    while True:
        n >>= 1
        if n == 0:
            break
        x = x * x
        if n % 2 == 1:
            y = y * x
    if m < 0:
        return f32(1.0) / y
    return y


def powi_d(x, m):
    """``x**m`` for real*8 x and integer variable m: port of libgcc ``__powidf2``."""
    m = int(m)
    n = -m if m < 0 else m
    x = float(x)
    y = x if n % 2 == 1 else 1.0
    while True:
        n >>= 1
        if n == 0:
            break
        x = x * x
        if n % 2 == 1:
            y = y * x
    if m < 0:
        return 1.0 / y
    return y


# --------------------------------------------------------------------------- integer intrinsics
def nint(x):
    """Fortran ``NINT``: round half away from zero.

    gfortran calls ``lroundf``/``lround`` (64-bit) and keeps the low 32 bits, so a value
    outside the int64 range (or NaN) gives ``LONG_MIN`` -> 0, and one outside int32 wraps.
    """
    x = float(x)
    if not (-9.2233720368547758e18 <= x < 9.2233720368547758e18):
        return 0
    t = math.trunc(x)
    frac = x - t  # exact
    if frac >= 0.5:
        t += 1
    elif frac <= -0.5:
        t -= 1
    if -2147483648 <= t <= 2147483647:
        return t
    return wrap32(t)


def anint(x):
    """Fortran ``ANINT`` of a real*4 (``roundf``)."""
    return f32(nint(x)) if math.isfinite(float(x)) else x


def aint(x):
    """Fortran ``AINT`` of a real*4 (truncation, result real)."""
    return f32(math.trunc(float(x))) if math.isfinite(float(x)) else x


def ifix(x):
    """Fortran ``INT``/``IFIX`` of a real: truncation toward zero.

    Out of the int32 range (or NaN) the x86 conversion gives -2147483648.
    """
    x = float(x)
    if not (-2147483649.0 < x < 2147483648.0):
        return -2147483648
    return math.trunc(x)


def idiv(a, b):
    """Fortran integer division ``a/b`` (truncates toward zero)."""
    a = int(a)
    b = int(b)
    q = abs(a) // abs(b)
    return q if (a >= 0) == (b >= 0) else -q


def imod(a, b):
    """Fortran ``mod(a, b)`` for integers (result has the sign of a)."""
    a = int(a)
    b = int(b)
    return a - idiv(a, b) * b


def wrap32(i):
    """32-bit two's-complement wrap-around of an integer result."""
    i = int(i) & 0xFFFFFFFF
    return i - 0x100000000 if i >= 0x80000000 else i


def sign(a, b):
    """Fortran ``SIGN(a, b)`` for reals: ``|a|`` with the sign of b."""
    return f32(math.copysign(abs(float(a)), float(b)))


def isign(a, b):
    """Fortran ``SIGN(a, b)`` for integers."""
    a = abs(int(a))
    return a if int(b) >= 0 else -a


def dsign(a, b):
    return math.copysign(abs(float(a)), float(b))


# --------------------------------------------------------------------------- gfortran SIND/COSD
# A Fortran unit that calls sind/cosd without declaring them EXTERNAL gets libgfortran's
# intrinsics, which differ from EPSR's own sind/cosd (auxiliary_routines.f).  Port of
# libgfortran (via the Rust fortran/fmath.rs).
_PIO180H = f32(0.0174560546875)
_PIO180L = f32(-2.7621674689726206e-06)


def _fmaf(a, b, c):
    # fused multiply-add of float32 values, exact via double (24+24 bit product fits)
    # followed by one rounding: a*b exact in double, then + c rounded in double could
    # double-round; use the exact fraction arithmetic instead.
    from fractions import Fraction

    return f32(float(Fraction(float(a)) * Fraction(float(b)) + Fraction(float(c))))


def _d2r(t):
    return _fmaf(t, _PIO180H, t * _PIO180L)


def _integral(z):
    n = math.trunc(float(z))
    if f32(n) - z == f32(0.0):
        return n
    return None


def sind_intrinsic(x):
    x = f32(x)
    if not math.isfinite(float(x)):
        return x - x
    ax = abs(x)
    if ax < f32(0.03125):
        return _fmaf(x, _PIO180H, x * _PIO180L)
    s = f32(math.copysign(1.0, float(x)))
    z = amod(ax, f32(360.0))
    n = _integral(z)
    if n is not None and n % 30 == 0:
        if n % 180 == 0:
            if n == 180:
                s = -s
            return s * f32(0.0)
        if n % 90 == 0:
            return s if n == 90 else -s
        if n % 60 == 0:
            return s * f32(-0.8660253882408142) if n > 179 else s * f32(0.8660253882408142)
        return s * f32(-0.5) if n > 179 else s * f32(0.5)
    if z <= 180.0:
        if z <= 90.0:
            if z <= 45.0:
                return s * sinf(_d2r(z))
            return s * cosf(_d2r(f32(90.0) - z))
        if z <= 135.0:
            return s * cosf(_d2r(z - f32(90.0)))
        return s * sinf(_d2r(f32(180.0) - z))
    s = -s
    if z <= 270.0:
        if z <= 225.0:
            return s * sinf(_d2r(z - f32(180.0)))
        return s * cosf(_d2r(f32(270.0) - z))
    if z <= 315.0:
        return s * cosf(_d2r(z - f32(270.0)))
    return s * sinf(_d2r(f32(360.0) - z))


def cosd_intrinsic(x):
    x = f32(x)
    if not math.isfinite(float(x)):
        return x - x
    ax = abs(x)
    if ax <= f32(0.0078125):
        return f32(1.0) if x == 0.0 else f32(1.0) - f32(7.888609052210118e-31)
    z = amod(ax, f32(360.0))
    n = _integral(z)
    if n is not None and n % 30 == 0:
        if n % 180 == 0:
            return f32(-1.0) if n == 180 else f32(1.0)
        if n % 90 == 0:
            return f32(0.0)
        if n % 60 == 0:
            return f32(0.5) if n in (60, 300) else f32(-0.5)
        return f32(0.8660253882408142) if n in (30, 330) else f32(-0.8660253882408142)
    if z <= 180.0:
        if z <= 90.0:
            if z <= 45.0:
                return cosf(_d2r(z))
            return sinf(_d2r(f32(90.0) - z))
        if z <= 135.0:
            return -sinf(_d2r(z - f32(90.0)))
        return -cosf(_d2r(f32(180.0) - z))
    if z <= 270.0:
        if z <= 225.0:
            return -cosf(_d2r(z - f32(180.0)))
        return -sinf(_d2r(f32(270.0) - z))
    if z <= 315.0:
        return sinf(_d2r(z - f32(270.0)))
    return cosf(_d2r(f32(360.0) - z))


def atan2d_intrinsic(y, x):
    """gfortran intrinsic ``ATAN2D`` for real(4): ``atan2f(y, x) * (180.0/pi)`` with the
    factor folded in single precision (57.2957764), not a correctly rounded result
    (checked against gfortran 15 on 2000 points, package F notes)."""
    return atan2f(y, x) * (f32(180.0) / f32(math.pi))


#: The value gfortran folds ``4.0*atan(1.0)`` to (single precision).
PI4ATAN = f32(math.pi)
