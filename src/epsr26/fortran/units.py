"""Fortran I/O units: OPEN, CLOSE, READ, WRITE, BACKSPACE, REWIND, INQUIRE.

A sequential formatted file is held in memory as a list of records with a record
pointer, which gives Fortran's semantics directly:

* ``open`` positions at the first record (``position='append'``: after the last one);
  with ``status='unknown'`` a missing file is created empty, as gfortran does;
* a WRITE replaces everything from the current record on (a written record becomes the
  last record of the file), so reading a file to its end and then writing appends;
* ``backspace``/``rewind`` move the pointer;
* the file is written to disk on ``close`` (and on ``close_all`` at the end of a program).

Unit 5 (and ``*`` for READ) is standard input, read on demand; units 6 (``*`` for WRITE)
and 0 are standard output and error, written immediately.  A WRITE to a unit that was
never opened goes to ``fort.N`` in the current folder, as in gfortran.

Translation patterns::

    units.open(10, fname, status='old')              # open(10,file=fname,status='old')
    line = units.read_line(10)                       # read(10,'(a)') line
    r = units.list_read(10); n = r.int(n)            # read(10,*) n
    units.write(10, fmt('(i5,f10.3)', n, x))         # write(10,'(i5,f10.3)') n,x
    units.write(6, lw('Value', x))                   # write(6,*) 'Value',x
    units.write(6, 'prompt: ', advance=False)        # write(6,'(a)',advance='no') ...
    units.close(10)

``FortranEOF`` and ``FortranIOError`` (``stop.py``) stand for ``end=``/``err=``/``iostat``.
"""

import builtins
import os
import sys

from .listio import ListRead, read_formatted
from .stop import FortranEOF, FortranIOError

ENCODING = "latin-1"


class Unit:
    def __init__(self, number, path):
        self.number = number
        self.path = path
        self.lines = []
        self.pos = 0
        self.dirty = False
        self.partial = None  # pending non-advancing output
        self.at_eof = False

    def next_record(self):
        if self.partial is not None:
            self._finish_partial()
        if self.pos >= len(self.lines):
            self.at_eof = True
            raise FortranEOF(f"end of file on unit {self.number} ({self.path})")
        self.pos += 1
        return self.lines[self.pos - 1]

    def write_record(self, rec, advance=True):
        if self.partial is not None:
            rec = self.partial + rec
            self.partial = None
        if not advance:
            self.partial = rec
            return
        if self.pos < len(self.lines):
            del self.lines[self.pos :]
        self.lines.append(rec)
        self.pos += 1
        self.dirty = True

    def _finish_partial(self):
        p = self.partial
        self.partial = None
        self.write_record(p)

    def flush(self):
        if self.partial is not None:
            self._finish_partial()
        if self.dirty:
            with builtins.open(self.path, "w", encoding=ENCODING, newline="") as fh:
                if self.lines:
                    fh.write("\n".join(self.lines))
                    fh.write("\n")
            self.dirty = False


class _Stdin:
    number = 5

    def next_record(self):
        line = sys.stdin.readline()
        if line == "":
            raise FortranEOF("end of file on standard input")
        if line.endswith("\n"):
            line = line[:-1]
        if line.endswith("\r"):
            line = line[:-1]
        return line


class _Stdout:
    def __init__(self, stream):
        self.stream = stream

    def write_record(self, rec, advance=True):
        s = self.stream()
        s.write(rec if not advance else rec + "\n")
        s.flush()

    def flush(self):
        self.stream().flush()


_units = {}
_stdin = _Stdin()
_stdout = _Stdout(lambda: sys.stdout)
_stderr = _Stdout(lambda: sys.stderr)


def _key(unit):
    if unit == "*":
        return None
    return int(unit)


def _get_read(unit):
    k = _key(unit)
    if k is None or (k == 5 and k not in _units):
        return _stdin
    u = _units.get(k)
    if u is None:
        u = _open_default(k)
    return u


def _get_write(unit):
    k = _key(unit)
    if k is None or (k == 6 and k not in _units):
        return _stdout
    if k == 0 and k not in _units:
        return _stderr
    u = _units.get(k)
    if u is None:
        u = _open_default(k)
    return u


def _open_default(k):
    return open(k, f"fort.{k}", status="unknown")


def _load(path):
    with builtins.open(path, "r", encoding=ENCODING, newline="") as fh:
        text = fh.read()
    if text == "":
        return []
    lines = text.split("\n")
    if text.endswith("\n"):
        lines.pop()
    return [ln[:-1] if ln.endswith("\r") else ln for ln in lines]


def open(unit, file, status=None, position=None, action=None, form=None):
    """``open(unit, file=file, status=status, position=position)``.

    Re-opening a unit already connected to the same file keeps the connection and its
    position when STATUS is absent or 'old'; any other STATUS is an error (iostat 5002)
    and the unit stays connected, as in gfortran.
    """
    k = int(unit)
    path = file.strip()
    if k in _units:
        u = _units[k]
        if os.path.abspath(u.path) == os.path.abspath(path):
            if status is None or status.lower() == "old":
                return u
            raise FortranIOError(5002, "Cannot change STATUS parameter in OPEN statement")
        close(k)
    st = (status or "unknown").lower()
    exists = os.path.isfile(path)
    if st == "old" and not exists:
        raise FortranIOError(2, f"Cannot open file '{path}': No such file or directory")
    if st == "new" and exists:
        raise FortranIOError(17, f"Cannot open file '{path}': File exists")
    if os.path.isdir(path):
        raise FortranIOError(21, f"Cannot open file '{path}': Is a directory")
    u = Unit(k, path)
    if exists and st != "replace" and st != "scratch":
        try:
            u.lines = _load(path)
        except OSError as e:
            raise FortranIOError(13, f"Cannot open file '{path}': {e.strerror}") from None
    else:
        if (action or "").lower() != "read":
            try:
                with builtins.open(path, "w", encoding=ENCODING):
                    pass
            except OSError as e:
                raise FortranIOError(2, f"Cannot open file '{path}': {e.strerror}") from None
    if position is not None and position.lower() == "append":
        u.pos = len(u.lines)
    _units[k] = u
    return u




def close(unit, status=None):
    k = int(unit)
    u = _units.pop(k, None)
    if u is None:
        return
    if status is not None and status.lower() == "delete":
        try:
            os.remove(u.path)
        except OSError:
            pass
        return
    u.flush()


def close_all():
    """Flush and close every open unit (end of a program)."""
    for k in list(_units):
        close(k)
    sys.stdout.flush()


def reset():
    """Forget all units without writing them (start of a program)."""
    _units.clear()


def save_state():
    """Detach the open units (the shell keeps its script unit open while a program runs)."""
    saved = dict(_units)
    _units.clear()
    return saved


def restore_state(saved):
    close_all()
    _units.update(saved)


def flush(unit):
    u = _units.get(_key(unit) if unit != "*" else 6)
    if u is not None:
        u.flush()


def is_open(unit):
    return int(unit) in _units


def unit_path(unit):
    u = _units.get(int(unit))
    return u.path if u is not None else None


def inquire_exists(file):
    """``inquire(file=file, exist=ex)``."""
    return os.path.exists(file.strip())


# --------------------------------------------------------------------------- READ
def read_line(unit):
    """``read(unit,'(a)') line``: the next record (trailing blanks kept off by callers)."""
    return _get_read(unit).next_record()


def skip_record(unit):
    """``read(unit,*)`` with no items."""
    _get_read(unit).next_record()


def list_read(unit):
    """Start a list-directed READ on ``unit``; see ``listio.ListRead``."""
    return ListRead(source=_get_read(unit).next_record)


def read_fmt(unit, fmtstr):
    """``read(unit, fmtstr) ...``: one record through an input FORMAT -> list of values."""
    return read_formatted(fmtstr, _get_read(unit).next_record())


def backspace(unit):
    u = _units.get(int(unit))
    if u is not None:
        if u.partial is not None:
            u._finish_partial()
        if u.at_eof:
            u.at_eof = False  # after an end-of-file, backspace goes back before the EOF
        elif u.pos > 0:
            u.pos -= 1


def rewind(unit):
    u = _units.get(int(unit))
    if u is not None:
        if u.partial is not None:
            u._finish_partial()
        u.pos = 0
        u.at_eof = False


def endfile(unit):
    u = _units.get(int(unit))
    if u is not None:
        del u.lines[u.pos :]
        u.dirty = True


# --------------------------------------------------------------------------- WRITE
def write(unit, text="", advance=True):
    """WRITE a record (``text`` may hold several records separated by ``'\\n'``)."""
    u = _get_write(unit)
    recs = text.split("\n")
    for r in recs[:-1]:
        u.write_record(r)
    u.write_record(recs[-1], advance)
