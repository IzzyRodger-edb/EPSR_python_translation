"""Fortran fixed-length ``character`` semantics.

Convention: a ``character*n`` variable is stored as a Python ``str`` **without trailing
blanks** (Fortran comparisons ignore trailing blanks, so ``==`` then behaves like
``.eq.``).  The declared length matters only in a few places, which use these helpers:

==============================  =============================================
Fortran                         Python
==============================  =============================================
``s = expr`` (character*n)      ``s = fstr(expr, n)``
``s(i:j)``                      ``substr(s, i, j)``   (blank-padded)
``s(i:j) = t``                  ``s = setsub(s, i, j, t, n)``
``index(s, t)``                 ``findex(s, t)``      (on the padded value if t has blanks)
``len_trim(s)``                 ``len(s)``
``len(s)``                      the declared length n
``s // t``                      ``s + t`` (or ``substr`` if s's padding matters)
``write(u,'(a)') s``            ``fmt('(a)', Ch(s, n))``  (padding reaches the output)
``ichar(c)``/``char(i)``        ``ord(c)``/``chr(i)``
==============================  =============================================
"""


def fstr(s, n):
    """Value of a ``character*n`` variable after ``var = s``: truncated to n, unpadded."""
    if len(s) > n:
        s = s[:n]
    return s.rstrip(" ")


def pad(s, n):
    """The full blank-padded ``character*n`` value."""
    if len(s) >= n:
        return s[:n]
    return s + " " * (n - len(s))


def substr(s, i, j):
    """Fortran substring ``s(i:j)`` (1-based, inclusive) of a blank-padded variable."""
    i = int(i)
    j = int(j)
    if j < i:
        return ""
    if i < 1:
        i = 1
    t = s[i - 1 : j]
    want = j - i + 1
    if len(t) < want:
        t = t + " " * (want - len(t))
    return t


def setsub(s, i, j, t, n):
    """``s(i:j) = t`` for a ``character*n`` variable s; returns the new (unpadded) value."""
    full = pad(s, n)
    piece = pad(t, j - i + 1)
    return fstr(full[: i - 1] + piece + full[j:], n)


def findex(s, sub, back=False):
    """Fortran ``index(s, sub [, back])``: 1-based position, 0 if absent."""
    p = s.rfind(sub) if back else s.find(sub)
    return p + 1


def findex_padded(s, sub, n, back=False):
    """``index`` on the blank-padded ``character*n`` value of s (e.g. ``index(s,' ')``)."""
    return findex(pad(s, n), sub, back)


def len_trim(s):
    return len(s.rstrip(" "))


def adjustl(s, n):
    return fstr(s.lstrip(" "), n)


def lowercase(s):
    return s.lower()
