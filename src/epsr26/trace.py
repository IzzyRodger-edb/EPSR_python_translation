"""Call tracing: see how and when each translated routine is called.

Every translated routine is decorated with ``@routine``.  The decorator only registers
the function, so it costs nothing while tracing is off.  ``enable()`` replaces each
registered function *in its module's namespace* by a tracing wrapper, and ``disable()``
puts the originals back.  This works because translated code always calls routines
through their module (``molfit.molfit(...)``, or a bare name within the same module),
never through a name imported with ``from ... import``.

* ``EPSR_TRACE=calls`` (environment) or ``trace on [pattern]`` (shell) writes an indented
  call tree to ``trace.log``, prefixed with the EPSR iteration number when known, e.g.
  ``iter 3 > molfit > movem > en_change``.
* ``EPSR_TRACE=summary`` or ``trace summary`` gives call counts and cumulative time.
"""

import fnmatch
import functools
import os
import sys
import time

_registry = []  # functions decorated with @routine
_active = {}  # (module, name) -> original function while tracing
_stack = []
_counts = {}
_times = {}
_log = None
_logpath = "trace.log"
_pattern = "*"
_write_calls = True
_enabled = False
context = {"iter": None}


def routine(f):
    """Register a translated Fortran routine for tracing."""
    _registry.append(f)
    return f


def _matches(full, name, pattern):
    """``pattern`` is one fnmatch pattern or several separated by commas."""
    for pat in pattern.split(","):
        pat = pat.strip()
        if pat and (fnmatch.fnmatch(full, pat) or fnmatch.fnmatch(name, pat)):
            return True
    return False


def _wrap(f):
    name = f.__name__
    mod = f.__module__.rsplit(".", 1)[-1]
    full = f"{mod}.{name}"

    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        _stack.append(name)
        if _log is not None and _write_calls:
            it = context.get("iter")
            prefix = f"iter {it} > " if it is not None else ""
            _log.write(prefix + " > ".join(_stack) + "\n")
        t = time.perf_counter()
        try:
            return f(*args, **kwargs)
        finally:
            dt = time.perf_counter() - t
            _stack.pop()
            _counts[full] = _counts.get(full, 0) + 1
            _times[full] = _times.get(full, 0.0) + dt

    wrapper.__wrapped_routine__ = f
    return wrapper


def enable(pattern="*", logpath="trace.log", calls=True):
    """Start tracing the routines whose ``module.name`` matches ``pattern``."""
    global _log, _pattern, _logpath, _write_calls, _enabled
    disable()
    _enabled = True
    _pattern = pattern
    _logpath = logpath
    _write_calls = calls
    _counts.clear()
    _times.clear()
    if calls:
        _log = open(logpath, "a", encoding="utf-8")
    for f in _registry:
        mod = sys.modules.get(f.__module__)
        if mod is None:
            continue
        full = f"{f.__module__.rsplit('.', 1)[-1]}.{f.__name__}"
        if not _matches(full, f.__name__, pattern):
            continue
        if getattr(mod, f.__name__, None) is f:
            _active[(f.__module__, f.__name__)] = f
            setattr(mod, f.__name__, _wrap(f))


def refresh():
    """Wrap routines registered since ``enable()`` (modules imported later)."""
    if not _enabled:
        return
    for f in _registry:
        key = (f.__module__, f.__name__)
        if key in _active:
            continue
        mod = sys.modules.get(f.__module__)
        if mod is None:
            continue
        full = f"{f.__module__.rsplit('.', 1)[-1]}.{f.__name__}"
        if not _matches(full, f.__name__, _pattern):
            continue
        if getattr(mod, f.__name__, None) is f:
            _active[key] = f
            setattr(mod, f.__name__, _wrap(f))


def disable():
    global _log, _enabled
    _enabled = False
    for (modname, name), f in _active.items():
        mod = sys.modules.get(modname)
        if mod is not None:
            setattr(mod, name, f)
    _active.clear()
    if _log is not None:
        _log.close()
        _log = None


def is_enabled():
    return _enabled


def summary(out=None):
    """Print call counts and cumulative (inclusive) time per routine."""
    out = out or sys.stdout
    if not _counts:
        out.write("No trace data (tracing is off or nothing ran).\n")
        return
    out.write(f"{'routine':40s} {'calls':>10s} {'time [s]':>12s}\n")
    for k in sorted(_times, key=lambda k: -_times[k]):
        out.write(f"{k:40s} {_counts[k]:10d} {_times[k]:12.3f}\n")


def reset_counts():
    _counts.clear()
    _times.clear()


def enable_from_environment():
    mode = os.environ.get("EPSR_TRACE", "")
    if mode == "calls":
        enable(os.environ.get("EPSR_TRACE_PATTERN", "*"))
    elif mode == "summary":
        enable(os.environ.get("EPSR_TRACE_PATTERN", "*"), calls=False)
