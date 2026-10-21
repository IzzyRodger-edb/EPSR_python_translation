"""Every external process started by EPSR goes through here (PLAN.md section 4.3).

Each child is started in its own process group and recorded, so that the shell can
terminate all of them (SIGTERM, then SIGKILL after a grace period) when it shuts down.
This is what removes the residual processes the Fortran shell leaves behind after Ctrl+C.

* ``run(cmd)``: ``call system(cmd)`` -- run through ``/bin/sh`` in the foreground and
  wait; returns the exit status.  Ctrl+C while waiting stops the child's group and
  re-raises ``KeyboardInterrupt``.
* ``spawn(cmd)``: start a background viewer (gnuplot, Jmol) and return at once.
* ``terminate_all()``: stop every recorded child that is still running, including the
  process groups of foreground commands that left background members behind
  (``gnuplot ... &`` started through ``call system``).
"""

import os
import signal
import subprocess
import time

_children = []  # (Popen, command, is_viewer)
_groups = []  # (pgid, command): groups whose leader exited but which still have members
GRACE_SECONDS = 3.0


def _popen(cmd, **kw):
    # a string runs as glibc system(3) runs it, argv "sh -c cmd" (errors read "sh: 1: ...")
    if isinstance(cmd, str):
        return subprocess.Popen(["sh", "-c", cmd], executable="/bin/sh", start_new_session=True, **kw)
    return subprocess.Popen(cmd, start_new_session=True, **kw)


def _tty_fd():
    import sys

    try:
        fd = sys.stdin.fileno()
        if os.isatty(fd) and os.tcgetpgrp(fd) == os.getpgrp():
            return fd
    except (OSError, ValueError, AttributeError):
        pass
    return None


def run(cmd, cwd=None):
    """Run ``cmd`` in the foreground (as ``call system(cmd)``); returns the exit status.

    On a terminal the child gets its own process group and becomes the terminal's
    foreground group while it runs (as under bash), so editors work and Ctrl+C goes to the
    child.  Without a terminal, Ctrl+C reaches the shell, which stops the child's group.
    """
    import sys

    sys.stdout.flush()
    fd = _tty_fd()
    if fd is None:
        p = _popen(cmd, cwd=cwd)
    else:
        a = ["sh", "-c", cmd] if isinstance(cmd, str) else cmd
        p = subprocess.Popen(a, executable="/bin/sh" if isinstance(cmd, str) else None, process_group=0, cwd=cwd)
        old = signal.signal(signal.SIGTTOU, signal.SIG_IGN)
        try:
            os.tcsetpgrp(fd, p.pid)
        except OSError:
            pass
        finally:
            signal.signal(signal.SIGTTOU, old)
    _children.append((p, cmd, False))
    try:
        return p.wait()
    except KeyboardInterrupt:
        _stop_group(p)
        raise
    finally:
        if fd is not None:
            old = signal.signal(signal.SIGTTOU, signal.SIG_IGN)
            try:
                os.tcsetpgrp(fd, os.getpgrp())
            except OSError:
                pass
            finally:
                signal.signal(signal.SIGTTOU, old)
        _forget(p)
        # 'gnuplot ... &' leaves members in the group after /bin/sh has exited
        if _group_alive(p.pid):
            _groups.append((p.pid, cmd))


def spawn(cmd, viewer=True):
    """Start ``cmd`` in the background; it is stopped when the shell exits."""
    p = _popen(cmd, stdin=subprocess.DEVNULL)
    _children.append((p, cmd, viewer))
    return p


def _group_alive(pgid):
    try:
        os.killpg(pgid, 0)
        return True
    except (ProcessLookupError, PermissionError):
        return False


def _stop_pgid(pgid, grace=GRACE_SECONDS):
    try:
        os.killpg(pgid, signal.SIGTERM)
    except ProcessLookupError:
        return
    t = time.monotonic() + grace
    while time.monotonic() < t:
        if not _group_alive(pgid):
            return
        time.sleep(0.05)
    try:
        os.killpg(pgid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def _forget(p):
    for k, (q, _, _) in enumerate(_children):
        if q is p:
            del _children[k]
            return


def _stop_group(p, grace=GRACE_SECONDS):
    if p.poll() is not None:
        return
    try:
        os.killpg(p.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    t = time.monotonic() + grace
    while time.monotonic() < t:
        if p.poll() is not None:
            return
        time.sleep(0.05)
    try:
        os.killpg(p.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    try:
        p.wait(timeout=1.0)
    except subprocess.TimeoutExpired:
        pass


def live():
    """[(pid, command, is_viewer)] of the children still running."""
    out = []
    for p, cmd, viewer in list(_children):
        if p.poll() is None:
            out.append((p.pid, cmd, viewer))
        else:
            _forget(p)
    for g in list(_groups):
        if _group_alive(g[0]):
            out.append((g[0], g[1], True))
        else:
            _groups.remove(g)
    return out


def terminate_all(keep_viewers=False):
    for p, cmd, viewer in list(_children):
        if keep_viewers and viewer:
            continue
        _stop_group(p)
        _forget(p)
    if not keep_viewers:
        for g in list(_groups):
            _stop_pgid(g[0])
            _groups.remove(g)
