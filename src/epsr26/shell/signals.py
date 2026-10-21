"""Ctrl+C handling and the clean shutdown sequence (PLAN.md section 4.3).

| situation                              | behaviour                                          |
|----------------------------------------|----------------------------------------------------|
| Ctrl+C at the prompt                   | discard the line, new prompt                       |
| Ctrl+C during ``epsr``                 | the current iteration becomes the last one: its    |
|                                        | results are written and the files saved            |
| 2nd Ctrl+C within 3 s, or Ctrl+C in    | abort the program (``KeyboardInterrupt``), back to |
| any other program                      | the prompt with a warning                          |
| SIGTERM, SIGHUP (terminal closed)      | shutdown sequence, exit                            |

Shutdown sequence (also run at exit): save the history, pause a running script in
``runflag.txt``, stop every child process group, remove the lock file.
"""

import atexit
import os
import signal
import sys
import time

from . import history, processes

ABORT_WINDOW = 3.0
LOCK_NAME = ".epsrshell.lock"


class _State:
    mode = "prompt"  # prompt | program | epsr
    stop_requested = False
    last_press = 0.0
    home = None
    keep_viewers = False
    script_pauser = None  # callable that pauses a running script
    done = False


state = _State()


# --------------------------------------------------------------------------- Ctrl+C
def _on_sigint(signum, frame):
    now = time.monotonic()
    if state.mode == "epsr":
        if not state.stop_requested:
            state.stop_requested = True
            state.last_press = now
            sys.stdout.write("\n\nEPSRshell> Ctrl+C: EPSR will stop after the current iteration and save its files."
                             "\nEPSRshell> Press Ctrl+C again within 3 s to abort at once"
                             " (files may then be incomplete).\n")
            sys.stdout.flush()
            return
        if now - state.last_press > ABORT_WINDOW:
            state.last_press = now
            sys.stdout.write("\nEPSRshell> Still finishing the current iteration."
                             " Press Ctrl+C again within 3 s to abort.\n")
            sys.stdout.flush()
            return
    raise KeyboardInterrupt


def stop_requested():
    """Polled by the epsr driver at the end of every iteration."""
    return state.stop_requested


class running:
    """Context manager around a program run: ``with signals.running('epsr'): ...``."""

    def __init__(self, program):
        self.mode = "epsr" if program == "epsr" else "program"

    def __enter__(self):
        state.mode = self.mode
        state.stop_requested = False
        return self

    def __exit__(self, *exc):
        state.mode = "prompt"
        state.stop_requested = False
        return False


# --------------------------------------------------------------------------- shutdown
def shutdown():
    """The shutdown sequence; safe to call more than once."""
    if state.done:
        return
    state.done = True
    try:
        history.save_history_list()
    except Exception:
        pass
    if state.script_pauser is not None:
        try:
            state.script_pauser()
        except Exception:
            pass
    processes.terminate_all(keep_viewers=state.keep_viewers)
    remove_lock()
    try:
        sys.stdout.flush()
    except Exception:
        pass


def _on_terminate(signum, frame):
    shutdown()
    os._exit(128 + signum)


def install(home):
    state.home = home
    signal.signal(signal.SIGINT, _on_sigint)
    signal.signal(signal.SIGTERM, _on_terminate)
    if hasattr(signal, "SIGHUP"):
        signal.signal(signal.SIGHUP, _on_terminate)
    atexit.register(shutdown)


# --------------------------------------------------------------------------- lock file
def _lock_path():
    return os.path.join(state.home or os.getcwd(), LOCK_NAME)


def _alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def take_lock():
    """Create the lock file; returns the PID of another live shell in this folder, if any."""
    p = _lock_path()
    other = None
    try:
        with open(p) as fh:
            pid = int(fh.read().split()[0])
        if pid != os.getpid() and _alive(pid):
            other = pid
    except (OSError, ValueError, IndexError):
        pass
    if other is None:
        try:
            with open(p, "w") as fh:
                fh.write(f"{os.getpid()}\n")
        except OSError:
            pass
    return other


def remove_lock():
    p = _lock_path()
    try:
        with open(p) as fh:
            pid = int(fh.read().split()[0])
        if pid == os.getpid():
            os.remove(p)
    except (OSError, ValueError, IndexError):
        pass
