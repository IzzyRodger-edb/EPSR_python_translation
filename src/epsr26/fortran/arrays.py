"""Fortran arrays as NumPy arrays with one extra element per dimension.

``a = real_array(n, m)`` is Fortran ``real a(n,m)``: index 0 is unused so that
``a[i, j]`` is ``a(i,j)``.  Arrays are column-major (``order='F'``) so that a column
``a[:, j]`` is a contiguous view, which is what Fortran sequence association passes when
it hands ``a(1,j)`` to a routine expecting a 1-D array.

A dimension may also be given as a ``(lower, upper)`` pair when the lower bound is 0,
e.g. ``real_array((0, n))`` for ``a(0:n)``.  Negative lower bounds are not supported by
these helpers; translate them with an explicit offset.

Do not use NumPy reductions or vectorised arithmetic on these arrays in translated code
(see ``fmath.py``).  Whole-array assignment of a constant (``a[:] = 0``) is fine.
"""

import numpy as np


def _shape(dims):
    shape = []
    for d in dims:
        if isinstance(d, tuple):
            lo, hi = d
            if lo != 0:
                raise ValueError("only lower bound 0 is supported")
            shape.append(int(hi) + 1)
        else:
            shape.append(int(d) + 1)
    return tuple(shape)


def real_array(*dims):
    """``real a(dims)`` (real*4), zero-initialised."""
    return np.zeros(_shape(dims), dtype=np.float32, order="F")


def dble_array(*dims):
    """``real*8 a(dims)``, zero-initialised."""
    return np.zeros(_shape(dims), dtype=np.float64, order="F")


def int_array(*dims):
    """``integer a(dims)`` (32-bit), zero-initialised."""
    return np.zeros(_shape(dims), dtype=np.int32, order="F")


def int8_array(*dims):
    """``integer*8 a(dims)``."""
    return np.zeros(_shape(dims), dtype=np.int64, order="F")


def logical_array(*dims):
    """``logical a(dims)``, all ``.false.``."""
    return np.zeros(_shape(dims), dtype=bool, order="F")


def char_array(*dims):
    """``character*n a(dims)``: an object array of unpadded strings, all blank ('')."""
    a = np.empty(_shape(dims), dtype=object, order="F")
    a.fill("")
    return a


def resize(a, *dims):
    """Reallocate ``a`` to the new dims keeping the overlapping elements
    (the ``reallocation_routines.f90`` pattern: allocate temp, copy, move_alloc)."""
    shape = _shape(dims)
    b = np.zeros(shape, dtype=a.dtype, order="F")
    if a.dtype == object:
        b.fill("")
    sl = tuple(slice(0, min(x, y)) for x, y in zip(a.shape, shape))
    b[sl] = a[sl]
    return b
