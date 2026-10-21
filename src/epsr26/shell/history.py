"""Line editing and history: the Linux ``FCreadline.f90`` (module jmk_readline_module) on
Python's GNU ``readline``.

* ``readline(prompt)`` reads an editable line; arrow keys walk the history.
* ``h`` on its own shows the numbered history; ``!n`` recalls line n and asks for
  confirmation (Enter accepts), as in the Fortran.  ``!`` followed by anything but a
  digit is returned as typed (the shell runs it as a system command).
* ``save_history_list()`` / ``get_history_list()`` write and merge ``FChistory.txt`` in
  the shell's start folder, de-duplicated in order of first use, exactly as the Fortran
  does (the programs the shell runs can then share the history).
"""

import os

try:
    import readline as _rl
except ImportError:  # pragma: no cover - non-GNU platforms
    _rl = None

HISTORY_FILE = "FChistory.txt"
_home = None


def set_home(folder):
    global _home
    _home = folder


def _path():
    return os.path.join(_home or os.getcwd(), HISTORY_FILE)


def _entries():
    if _rl is None:
        return []
    return [_rl.get_history_item(i) for i in range(1, _rl.get_current_history_length() + 1)]


def show_history_list():
    print("\nCurrent history list:\n")
    for i, line in enumerate(_entries(), 1):
        print("%6d %s" % (i, line))
    print("\nTo use a particular line type !<line number>")


def get_history_line(linenumber):
    e = _entries()
    if 0 < linenumber <= len(e):
        return e[linenumber - 1].rstrip(" "), True
    return "", False


def readline(prompt):
    """FCreadline.f90:readline -- read an editable input line (``EOFError`` on Ctrl+D)."""
    while True:
        line = input(prompt.rstrip() + " ")
        if line == "h":
            show_history_list()
            found = False
        elif line[:1] == "!" and line[1:2].isdigit():  # '!cmd' (not a number) is a shell command
            try:
                linenumber = int(line[1:].split()[0])
            except (ValueError, IndexError):
                linenumber = 0
            line, found = get_history_line(linenumber)
            if found:
                test = input("\n" + prompt.rstrip() + " " + line.rstrip() + " ? (Press <enter> to accept.)")
                found = len(test.strip()) == 0
        else:
            found = True
            # Don't add blank lines to history, but return them to calling program
            if _rl is not None and line.strip():
                _rl.add_history(line)
        if not found:
            print()  # Blank line to fit with EPSR format
            continue
        return line


def save_history_list():
    """FCreadline.f90:save_history_list -- write the de-duplicated history."""
    seen = []
    for line in _entries():
        t = line[:255].rstrip(" ")
        if t and t not in seen:
            seen.append(t)
    try:
        with open(_path(), "w") as fh:
            for t in seen:
                fh.write(t + "\n")
    except OSError:
        pass


def get_history_list():
    """FCreadline.f90:get_history_list -- add lines of FChistory.txt not yet in the history."""
    if _rl is None:
        return
    current = [e.rstrip(" ") for e in _entries()]
    try:
        with open(_path()) as fh:
            lines = fh.read().split("\n")
    except OSError:
        return
    for line in lines:
        t = line.rstrip(" ")
        if t and t not in current:
            _rl.add_history(t[:255])
            current.append(t)


def setup():
    """History is managed explicitly (as FCreadline does), not by input()."""
    if _rl is not None and hasattr(_rl, "set_auto_history"):
        _rl.set_auto_history(False)
