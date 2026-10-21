"""Fortran formatted and list-directed *output*, byte-compatible with gfortran.

Port of the Rust ``fortran/format.rs``.

* ``fmt(format, *values)`` interprets a FORMAT string against values and returns the
  produced records joined by ``'\\n'`` (no trailing newline).  Repeat counts, nested
  groups, ``/``, ``:``, ``nX``, ``nH`` Hollerith, quoted strings and format reversion
  are supported, with the data edit descriptors ``Iw``, ``Iw.m``, ``I0``, ``Fw.d``,
  ``Ew.d``, ``Ew.dEe``, ``Aw``, ``A``, ``Lw``.  ``kP`` is accepted with ``0P`` meaning.
* ``lw(*values)`` is one record of list-directed output (``write(u,*) ...``).

Value types decide the conversion, as the declared Fortran type does:

========================  ======================================================
Python value              Fortran
========================  ======================================================
``int``, ``numpy.int32``  ``integer`` (list-directed width 11)
``I8(v)``                 ``integer*8`` (list-directed width 20)
``numpy.float32``         ``real`` (list-directed: 9 significant digits)
``float``                 ``real*8`` (list-directed: 17 significant digits)
``str``                   a character expression whose length is its text
``Ch(s, n)``              a ``character*n`` variable: blank-padded to n
``bool``                  ``logical``
========================  ======================================================

Character variables are stored *without* trailing blanks (see ``text.py``), so a
variable written with ``A`` editing must be passed as ``Ch(var, declared_length)``
wherever the trailing blanks reach the output.

All numeric conversion is done from the exact binary value with round-half-even, which
is what glibc ``printf`` (used by libgfortran) does; Python's ``%e``/``%f`` do the same.
"""

import math

import numpy as np


class Ch:
    """A character variable of declared length ``n`` (value stored without padding)."""

    __slots__ = ("s", "n")

    def __init__(self, s, n):
        self.s = s
        self.n = n

    def padded(self):
        s = self.s
        if len(s) >= self.n:
            return s[: self.n]
        return s + " " * (self.n - len(s))

    def __repr__(self):
        return f"Ch({self.s!r}, {self.n})"


class I8:
    """An ``integer*8`` value (list-directed width 20)."""

    __slots__ = ("v",)

    def __init__(self, v):
        self.v = int(v)


def _is_int(v):
    return isinstance(v, (int, np.integer, I8)) and not isinstance(v, bool)


def _ival(v):
    return v.v if isinstance(v, I8) else int(v)


def _fval(v):
    if isinstance(v, (bool, str, Ch)):
        raise TypeError(f"numeric edit descriptor given {v!r}")
    if isinstance(v, I8):
        return float(v.v)
    return float(v)


def _text(v):
    if isinstance(v, Ch):
        return v.padded()
    if isinstance(v, str):
        return v
    raise TypeError(f"A edit descriptor given {v!r}")


# --------------------------------------------------------------------------- edit descriptors
def _sig_digits(a, n):
    """(digits, exp10) with |a| = 0.d1..dn x 10**exp10, correctly rounded; a > 0."""
    s = "%.*e" % (max(n, 1) - 1, a)
    mant, ex = s.split("e")
    digits = mant.replace(".", "")
    return digits, int(ex) + 1


def _stars(w):
    return "*" * w


def _edit_infnan(v, w):
    if math.isnan(v):
        if w < 3:
            return _stars(w)
        return "NaN".rjust(w)
    neg = v < 0
    if w < 3 or (neg and w < 4):
        return _stars(w)
    body = "Infinity" if ((neg and w > 8) or (not neg and w > 7)) else "Inf"
    s = "-" + body if neg else body
    return s.rjust(w)


def edit_e(v, w, d, e=None):
    """Fortran ``Ew.d`` (or ``Ew.dEe``)."""
    if not math.isfinite(v):
        return _edit_infnan(v, w)
    neg = math.copysign(1.0, v) < 0
    if v == 0.0:
        digits, ex = "0" * max(d, 1), 0
    else:
        digits, ex = _sig_digits(abs(v), d)
    s = ("-" if neg else "") + "0." + digits[:d]
    aex = abs(ex)
    sgn = "-" if ex < 0 else "+"
    if e is None:
        if aex <= 99:
            s += "E%s%02d" % (sgn, aex)
        elif aex <= 999:
            s += "%s%03d" % (sgn, aex)
        else:
            return _stars(w)
    else:
        es = "%0*d" % (e, aex)
        if len(es) > e:
            return _stars(w)
        s += "E" + sgn + es
    if len(s) > w:
        alt = ("-" + s[2:]) if neg else s[1:]
        if len(alt) > w:
            return _stars(w)
        s = alt
    return s.rjust(w)


def edit_f(v, w, d):
    """Fortran ``Fw.d``."""
    if not math.isfinite(v):
        return _edit_infnan(v, w)
    neg = math.copysign(1.0, v) < 0
    body = "%.*f" % (d, abs(v))
    if d == 0:
        body += "."
    s = "-" + body if neg else body
    if len(s) > w or w == 0:
        if body.startswith("0."):
            b2 = body[1:]
            s = "-" + b2 if neg else b2
        if len(s) > w:
            return _stars(w)
    return s.rjust(w)


def edit_i(v, w, m=None):
    """Fortran ``Iw``, ``Iw.m`` and ``I0``."""
    digits = str(abs(v))
    if m is not None:
        if m == 0 and v == 0:
            digits = ""
        elif len(digits) < m:
            digits = "0" * (m - len(digits)) + digits
    s = "-" + digits if v < 0 else digits
    if w == 0:
        return s
    if len(s) > w:
        return _stars(w)
    return s.rjust(w)


# --------------------------------------------------------------------------- FORMAT parsing
# items: ('G', n, [items]) group, ('I', n, w, m), ('F', n, w, d), ('E', n, w, d, e),
#        ('A', n, w), ('L', n, w), ('X', n), ('LIT', s), ('/', n), (':',), ('P', k),
#        ('GG', n, w, d) generalised G (input only)


class _Parser:
    def __init__(self, s):
        self.s = s
        self.p = 0

    def peek(self):
        return self.s[self.p] if self.p < len(self.s) else None

    def skip_ws(self):
        while self.peek() in (" ", ","):
            self.p += 1

    def skip_spaces(self):
        while self.peek() == " ":
            self.p += 1

    def num(self):
        st = self.p
        while self.peek() is not None and self.peek().isdigit():
            self.p += 1
        if self.p > st:
            return int(self.s[st : self.p])
        return None

    def quoted(self, q):
        self.p += 1
        out = []
        while True:
            c = self.peek()
            if c is None:
                raise ValueError("unterminated string in FORMAT")
            self.p += 1
            if c == q:
                if self.peek() == q:
                    out.append(q)
                    self.p += 1
                else:
                    break
            else:
                out.append(c)
        return "".join(out)

    def items(self):
        items = []
        while True:
            self.skip_ws()
            c = self.peek()
            if c is None:
                break
            if c == ")":
                self.p += 1
                break
            if c in ("'", '"'):
                items.append(("LIT", self.quoted(c)))
                continue
            if c == "/":
                self.p += 1
                items.append(("/", 1))
                continue
            if c == ":":
                self.p += 1
                items.append((":",))
                continue
            neg = False
            if c == "-":
                neg = True
                self.p += 1
            if c == "*":  # unlimited repeat '*(...)' (Fortran 2008)
                self.p += 1
                self.skip_spaces()
                rep = 999999
            else:
                rep = self.num()
            self.skip_spaces()
            c = self.peek()
            c = c.upper() if c is not None else None
            if c == "(":
                self.p += 1
                items.append(("G", (1 if rep is None else rep), self.items()))
            elif c == "H":
                self.p += 1
                n = rep
                items.append(("LIT", self.s[self.p : self.p + n]))
                self.p += n
            elif c == "X":
                self.p += 1
                items.append(("X", rep if rep is not None else 1))
            elif c == "P":
                self.p += 1
                k = rep or 0
                items.append(("P", -k if neg else k))
            elif c == "/":
                self.p += 1
                items.append(("/", (1 if rep is None else rep)))
            elif c == "I":
                self.p += 1
                w = self.num() or 0
                m = None
                if self.peek() == ".":
                    self.p += 1
                    m = self.num()
                items.append(("I", (1 if rep is None else rep), w, m))
            elif c == "F":
                self.p += 1
                w = self.num()
                if self.peek() != ".":
                    raise ValueError(f"Fw.d expected in {self.s!r}")
                self.p += 1
                d = self.num()
                items.append(("F", (1 if rep is None else rep), w, d))
            elif c in ("E", "D"):
                self.p += 1
                w = self.num()
                if self.peek() != ".":
                    raise ValueError(f"Ew.d expected in {self.s!r}")
                self.p += 1
                d = self.num()
                e = None
                if self.peek() in ("E", "e"):
                    self.p += 1
                    e = self.num()
                items.append(("E", (1 if rep is None else rep), w, d, e))
            elif c == "G":
                self.p += 1
                w = self.num()
                d = None
                if self.peek() == ".":
                    self.p += 1
                    d = self.num()
                items.append(("GG", (1 if rep is None else rep), w, d))
            elif c == "A":
                self.p += 1
                w = self.num()
                items.append(("A", (1 if rep is None else rep), w))
            elif c == "L":
                self.p += 1
                w = self.num() or 1
                items.append(("L", (1 if rep is None else rep), w))
            elif c in ("'", '"'):
                items.append(("LIT", self.quoted(self.peek())))
            else:
                raise ValueError(f"unsupported FORMAT descriptor {c!r} in {self.s!r}")
        return items


_cache = {}


def parse_format(f):
    items = _cache.get(f)
    if items is None:
        t = f.strip()
        if not t.startswith("("):
            raise ValueError(f"FORMAT must start with '(': {f!r}")
        items = _Parser(t[1:]).items()
        _cache[f] = items
    return items


def _has_data(items):
    for it in items:
        k = it[0]
        if k == "G":
            if _has_data(it[2]):
                return True
        elif k in ("I", "F", "E", "A", "L", "GG"):
            return True
    return False


# --------------------------------------------------------------------------- FORMAT execution
class _Exec:
    def __init__(self, vals):
        self.vals = vals
        self.next = 0
        self.out = []
        self.rec = []
        self.pending = 0
        self.done = False

    def emit(self, s):
        if self.pending:
            self.rec.append(" " * self.pending)
            self.pending = 0
        self.rec.append(s)

    def newrec(self):
        self.out.append("".join(self.rec))
        self.rec = []
        self.pending = 0

    def take(self):
        if self.next < len(self.vals):
            self.next += 1
            return True, self.vals[self.next - 1]
        self.done = True
        return False, None

    def run(self, items):
        for it in items:
            k = it[0]
            if k == "G":
                for _ in range(it[1]):
                    if not self.run(it[2]):
                        return False
            elif k == "LIT":
                self.emit(it[1])
            elif k == "X":
                self.pending += it[1]
            elif k == "/":
                for _ in range(it[1]):
                    self.newrec()
            elif k == ":":
                if self.next >= len(self.vals):
                    self.done = True
                    return False
            elif k == "P":
                if it[1] != 0:
                    raise NotImplementedError("nonzero scale factor kP in output FORMAT")
            elif k == "I":
                for _ in range(it[1]):
                    ok, v = self.take()
                    if not ok:
                        return False
                    if not _is_int(v):
                        raise TypeError(f"I edit descriptor given {v!r}")
                    self.emit(edit_i(_ival(v), it[2], it[3]))
            elif k == "F":
                for _ in range(it[1]):
                    ok, v = self.take()
                    if not ok:
                        return False
                    self.emit(edit_f(_fval(v), it[2], it[3]))
            elif k == "E":
                for _ in range(it[1]):
                    ok, v = self.take()
                    if not ok:
                        return False
                    self.emit(edit_e(_fval(v), it[2], it[3], it[4]))
            elif k == "A":
                for _ in range(it[1]):
                    ok, v = self.take()
                    if not ok:
                        return False
                    full = _text(v)
                    w = it[2]
                    if w is None:
                        s = full
                    elif w >= len(full):
                        s = " " * (w - len(full)) + full
                    else:
                        s = full[:w]
                    self.emit(s)
            elif k == "L":
                for _ in range(it[1]):
                    ok, v = self.take()
                    if not ok:
                        return False
                    self.emit(("T" if v else "F").rjust(it[2]))
            elif k == "GG":
                raise NotImplementedError("G edit descriptor in output FORMAT")
        return True


def fmt(format, *vals):
    """Apply a Fortran FORMAT to ``vals``; records are joined by ``'\\n'``."""
    items = parse_format(format)
    ex = _Exec(vals)
    rev_start = 0
    for k in range(len(items) - 1, -1, -1):
        if items[k][0] == "G":
            rev_start = k
            break
    first = True
    has_data = _has_data(items)
    while True:
        cont = ex.run(items if first else items[rev_start:])
        first = False
        if not cont or ex.done:
            break
        if ex.next >= len(vals):
            break
        if not has_data:
            break
        ex.newrec()  # format reversion: new record
    ex.out.append("".join(ex.rec))
    return "\n".join(ex.out)


# --------------------------------------------------------------------------- list-directed output
def _list_real(v, prec, expw):
    width = prec + expw + 5  # 16 for real*4, 25 for real*8
    if math.isnan(v):
        return "NaN".rjust(width)
    neg = math.copysign(1.0, v) < 0
    a = abs(v)
    sgn = "-" if neg else ""
    fw = width - (expw + 2)
    if a == 0.0:
        s = sgn + "0." + "0" * (prec - 1)
        return s.rjust(fw) + " " * (expw + 2)
    if math.isinf(a):
        return (sgn + "Infinity").rjust(width)
    s = "%.*e" % (prec - 1, a)
    mant, ex = s.split("e")
    e = int(ex)
    digits = mant.replace(".", "")
    if -1 <= e < prec:
        if e >= 0:
            ip = e + 1
            body = digits[:ip] + "." + digits[ip:]
        else:
            body = "0." + digits
        return (sgn + body).rjust(fw) + " " * (expw + 2)
    es = "%s%s.%sE%s%0*d" % (sgn, digits[:1], digits[1:], "-" if e < 0 else "+", expw, abs(e))
    return es.rjust(width)


def list_real4(v):
    """gfortran list-directed ``REAL(4)`` field (16 characters, 9 significant digits)."""
    return _list_real(float(v), 9, 2)


def list_real8(v):
    """gfortran list-directed ``REAL(8)`` field (25 characters, 17 significant digits)."""
    return _list_real(float(v), 17, 3)


def lw(*vals):
    """One record of gfortran list-directed output, ``write(u,*) vals``."""
    if not vals:
        return ""
    parts = []
    first = True
    last_char = False
    for v in vals:
        is_char = isinstance(v, (str, Ch))
        if first:
            parts.append(" ")
            first = False
        elif not (is_char and last_char):
            parts.append(" ")
        last_char = is_char
        if isinstance(v, bool) or isinstance(v, np.bool_):
            parts.append("T" if v else "F")
        elif isinstance(v, I8):
            parts.append("%20d" % v.v)
        elif isinstance(v, (int, np.integer)):
            parts.append("%11d" % int(v))
        elif isinstance(v, np.float32):
            parts.append(list_real4(v))
        elif isinstance(v, (float, np.float64)):
            parts.append(list_real8(v))
        elif isinstance(v, Ch):
            parts.append(v.padded())
        elif isinstance(v, str):
            parts.append(v)
        else:
            raise TypeError(f"list-directed write of {v!r}")
    return "".join(parts)
