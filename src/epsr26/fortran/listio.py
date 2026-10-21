"""Fortran *input*: list-directed reads and the fixed formats EPSR uses.

Port of the Rust ``fortran/listread.rs``.

List-directed semantics (gfortran):

* values are separated by blanks, commas or record ends, so one READ may consume
  several records;
* when a READ statement completes, the rest of the current record is skipped;
* ``/`` terminates the statement leaving the remaining items unchanged;
* ``r*value`` repeat counts and null values (``,,``) are honoured;
* reals accept ``1.0``, ``.5``, ``5.``, ``1e5``, ``1.0D+01``, ``1.0+01`` and are converted
  with correct rounding directly to real*4 (gfortran uses ``strtof``);
* reading an integer from ``1.0`` is an error.

A list-directed READ is written as::

    r = units.list_read(10)          # read(10,*,err=..,end=..) n, x, y
    n = r.int(n)                     # each getter takes the old value, kept on null
    x = r.real(x)
    y = r.real(y)
    r.end()                          # optional; only needed for a READ with no items

The getters raise ``FortranEOF`` / ``FortranIOError`` (``end=`` / ``err=`` branches).
"""

import numpy as np

from .stop import FortranEOF, FortranIOError


def _normalise_real(tok):
    s = tok.strip()
    if not s:
        raise FortranIOError(5010, "empty real")
    low = s.lower()
    if low.lstrip("+-") in ("inf", "infinity", "nan"):
        return low
    out = []
    i = 0
    n = len(s)
    if s[0] in "+-":
        out.append(s[0])
        i = 1
    nd = 0
    seen_dot = False
    while i < n and (s[i].isdigit() or (s[i] == "." and not seen_dot)):
        if s[i] == ".":
            seen_dot = True
        else:
            nd += 1
        out.append(s[i])
        i += 1
    if nd == 0:
        raise FortranIOError(5010, f"bad real {tok!r}")
    if i < n:
        c = s[i].lower()
        j = i
        if c in "edq":
            j += 1
        elif c not in "+-":
            raise FortranIOError(5010, f"bad real {tok!r}")
        out.append("e")
        if j < n and s[j] in "+-":
            out.append(s[j])
            j += 1
        st = j
        while j < n and s[j].isdigit():
            out.append(s[j])
            j += 1
        if j != n or j == st:
            raise FortranIOError(5010, f"bad real {tok!r}")
    return "".join(out)


def parse_real(tok):
    """A Fortran real constant -> numpy.float32 (correctly rounded from the decimal)."""
    return np.float32(_normalise_real(tok))


def parse_dble(tok):
    """A Fortran real constant -> Python float (real*8)."""
    return float(_normalise_real(tok))


def parse_int(tok):
    s = tok.strip()
    s2 = s[1:] if s.startswith("+") else s
    body = s2[1:] if s2.startswith("-") else s2
    if not body or not body.isdigit():
        raise FortranIOError(5010, f"bad integer {tok!r}")
    return int(s2)


_NULL = object()
_SLASH = object()


def tokenize(line):
    """Split one record into list-directed tokens (str, _NULL or _SLASH)."""
    v = []
    i = 0
    n = len(line)
    last_was_comma = True
    pending = False
    while i < n:
        c = line[i]
        if c in " \t\r":
            i += 1
            continue
        if c == ",":
            if last_was_comma and not pending:
                v.append(_NULL)
            last_was_comma = True
            pending = False
            i += 1
            continue
        if c == "/":
            v.append(_SLASH)
            return v
        if c in ("'", '"'):
            q = c
            s = []
            i += 1
            while i < n:
                if line[i] == q:
                    if i + 1 < n and line[i + 1] == q:
                        s.append(q)
                        i += 2
                        continue
                    i += 1
                    break
                s.append(line[i])
                i += 1
            v.append("".join(s))
            last_was_comma = False
            pending = True
            continue
        st = i
        while i < n and line[i] not in " \t,/\r":
            i += 1
        tok = line[st:i]
        p = tok.find("*")
        if p > 0 and tok[:p].isdigit():
            r = int(tok[:p])
            val = tok[p + 1 :]
            for _ in range(r):
                v.append(val if val else _NULL)
            last_was_comma = False
            pending = True
            continue
        v.append(tok)
        last_was_comma = False
        pending = True
    return v


class ListRead:
    """One list-directed READ statement in progress.

    ``source`` is a callable returning the next record (raising ``FortranEOF``), or None
    for an internal read from ``text``.
    """

    def __init__(self, source=None, text=None):
        self.source = source
        self.toks = tokenize(text) if text is not None else []
        self.k = 0
        self.started = text is not None
        self.terminated = False

    def _next(self):
        if self.terminated:
            return _NULL
        while True:
            if self.k < len(self.toks):
                t = self.toks[self.k]
                self.k += 1
                if t is _SLASH:
                    self.terminated = True
                    return _NULL
                return t
            if self.source is None:
                raise FortranEOF("end of internal record")
            line = self.source()
            self.started = True
            self.toks = tokenize(line)
            self.k = 0

    def int(self, old=None):
        t = self._next()
        return old if t is _NULL else parse_int(t)

    def real(self, old=None):
        t = self._next()
        return old if t is _NULL else parse_real(t)

    def dble(self, old=None):
        t = self._next()
        return old if t is _NULL else parse_dble(t)

    def word(self, old=None):
        """A character item: an unquoted word or a quoted string."""
        t = self._next()
        return old if t is _NULL else t

    def logical(self, old=None):
        t = self._next()
        if t is _NULL:
            return old
        u = t.lstrip(".").upper()
        if u.startswith("T"):
            return True
        if u.startswith("F"):
            return False
        raise FortranIOError(5010, f"bad logical {t!r}")

    def end(self):
        """End the statement; a READ with no items still consumes one record."""
        if not self.started and self.source is not None:
            self.source()
            self.started = True


def read_internal(text):
    """``read(text,*) ...``: a list-directed read from a character variable."""
    return ListRead(text=text)


# --------------------------------------------------------------------------- formatted input
def read_formatted(fmtstr, line):
    """``read(u, '(fmt)') ...`` of one record: returns the list of values.

    Supports ``Iw``, ``Fw.d``, ``Ew.d``, ``Gw.d``, ``Aw``, ``A``, ``nX`` and literal skips.
    Blanks in numeric fields are ignored (BN, the gfortran default for files); an all-blank
    numeric field is zero.  ``A`` without width takes the rest of the record.
    """
    from .fformat import parse_format

    items = parse_format(fmtstr)
    vals = []
    pos = [0]

    def field(w):
        s = line[pos[0] : pos[0] + w]
        pos[0] += w
        return s

    def run(items):
        for it in items:
            k = it[0]
            if k == "G":
                for _ in range(it[1]):
                    run(it[2])
            elif k == "X":
                pos[0] += it[1]
            elif k == "LIT":
                pos[0] += len(it[1])
            elif k == "I":
                for _ in range(it[1]):
                    s = field(it[2]).replace(" ", "")
                    vals.append(parse_int(s) if s else 0)
            elif k in ("F", "E", "GG"):
                d = it[3]
                for _ in range(it[1]):
                    s = field(it[2]).replace(" ", "")
                    if not s:
                        vals.append(np.float32(0.0))
                    elif "." not in s and d:
                        # implied decimal point: the last d digits of the mantissa are
                        # the fraction (Fortran rule for a field without '.')
                        m = len(s)
                        for j in range(1, len(s)):
                            if s[j] in "eEdDqQ+-":
                                m = j
                                break
                        sign = s[0] if s[0] in "+-" else ""
                        mant = s[len(sign) : m].rjust(d, "0")
                        vals.append(parse_real(sign + mant[:-d] + "." + mant[-d:] + s[m:] if d else s))
                    else:
                        vals.append(parse_real(s))
            elif k == "A":
                for _ in range(it[1]):
                    w = it[2]
                    if w is None:
                        s = line[pos[0] :]
                        pos[0] = len(line)
                    else:
                        s = field(w)
                    vals.append(s.rstrip(" "))
            elif k == "/":
                raise NotImplementedError("'/' in an input FORMAT")
            elif k in (":", "P"):
                pass
            else:
                raise NotImplementedError(f"input FORMAT item {it!r}")

    run(items)
    return vals
