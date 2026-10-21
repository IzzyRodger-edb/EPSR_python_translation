"""Program-level runtime: command-line arguments (``getarg``/``iargc``) and the clock."""

import time

_argv = []
_t0 = time.time()


def set_args(args):
    """Install the command-line arguments of the program about to run (without argv[0])."""
    global _argv
    _argv = [str(a) for a in args]


def iargc():
    """``iargc()`` / ``command_argument_count()``."""
    return len(_argv)


def getarg(i):
    """``call getarg(i, arg)`` / ``get_command_argument(i, arg)``: '' when out of range."""
    i = int(i)
    if 1 <= i <= len(_argv):
        return _argv[i - 1]
    return ""


def cpu_time():
    """``call cpu_time(t)`` (real*4 seconds)."""
    import numpy as np

    return np.float32(time.process_time())


def system_clock():
    """``call system_clock(count, count_rate)`` with gfortran's default integer*4 rate."""
    t = time.monotonic()
    return int(t * 1000) & 0x7FFFFFFF, 1000


def date_and_time():
    """``call date_and_time(date, time, zone, values)`` -> (date, time, zone, values[1..8])."""
    lt = time.localtime()
    t = time.time()
    ms = int((t - int(t)) * 1000)
    off = -time.altzone if lt.tm_isdst > 0 else -time.timezone
    sgn = "+" if off >= 0 else "-"
    zone = "%s%02d%02d" % (sgn, abs(off) // 3600, (abs(off) % 3600) // 60)
    date = "%04d%02d%02d" % (lt.tm_year, lt.tm_mon, lt.tm_mday)
    tim = "%02d%02d%02d.%03d" % (lt.tm_hour, lt.tm_min, lt.tm_sec, ms)
    values = [0, lt.tm_year, lt.tm_mon, lt.tm_mday, off // 60, lt.tm_hour, lt.tm_min, lt.tm_sec, ms]
    return date, tim, zone, values


def system(command):
    """``call system(command)`` / ``execute_command_line``: run through /bin/sh in the
    foreground, recorded in the shell's process registry; returns the exit status."""
    from ..shell import processes

    return processes.run(command.strip())


def getenv(name):
    """``call getenv(name, value)`` / ``get_environment_variable``: '' if unset."""
    import os

    return os.environ.get(name.strip(), "")
