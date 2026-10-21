"""Shell built-ins that replace the ``system_commands.txt`` templates (fix 1, PLAN.md 4.2).

In the Fortran, ``ls``, ``cd``, ``md``, ``del`` and ``edit`` build a command line from
``system_commands.txt`` (``system_ls``, ``system_termination``, ...).  The template
shipped with EPSR26 holds the Windows values (``dir/b``, ``cd/d``, ``\\``), so on Linux
``ls`` runs ``dir/b`` and ``cd`` splits paths on ``\\``.  These versions work directly.

As in the Fortran, the shell process itself stays in the folder where it was started
(``runflag.txt`` and ``FChistory.txt`` live there); ``cd`` changes ``directory_name``, the
folder every program is given as its input folder.  ``ls``, ``system`` and external
commands run in that folder.
"""

import glob
import os
import shlex
import shutil

from ..common import directory_name as dn
from ..common import system_commands as sc
from . import processes

_previous = None


def current():
    """The working folder (``directory_name`` without the trailing separator)."""
    d = dn.directory_name[: dn.nchardir] if dn.nchardir else os.getcwd()
    return d.rstrip("/") or "/"


def set_current(path):
    path = os.path.abspath(path)
    d = path if path.endswith("/") else path + "/"
    dn.directory_name = d
    dn.nchardir = len(d)


def pwd(words):
    print(f"\nCurrent folder for input files: {current()}/")


def cd(words):
    """``cd``, ``cd home``, ``cd .``: back to the start folder (Fortran behaviour);
    ``cd ..``, ``cd -``, ``cd ~``, absolute and relative paths as in bash."""
    global _previous
    here = current()
    if len(words) < 2 or words[1] in ("home", "."):
        target = dn.home_directory[: dn.ncharhomedir] or os.getcwd()
    elif words[1] == "-":
        if _previous is None:
            print("\ncd> No previous folder")
            return
        target = _previous
    else:
        arg = " ".join(words[1:])
        arg = os.path.expanduser(arg)
        target = arg if os.path.isabs(arg) else os.path.join(here, arg)
    target = os.path.normpath(target)
    if not os.path.isdir(target):
        print(f"\nNo such folder! : {target}/")
    else:
        _previous = here
        set_current(target)
    print(f"\nCurrent folder for input files: {current()}/")


def ls(words):
    """``ls``/``dir``: with no arguments list the folder; otherwise ``/bin/ls args`` there."""
    if len(words) < 2:
        names = sorted(os.listdir(current()))
        if not names:
            return
        width = shutil.get_terminal_size((80, 24)).columns
        disp = [n + ("/" if os.path.isdir(os.path.join(current(), n)) else "") for n in names]
        w = max(len(n) for n in disp) + 2
        ncol = max(1, width // w)
        nrow = (len(disp) + ncol - 1) // ncol
        for r in range(nrow):
            row = [disp[c * nrow + r] for c in range(ncol) if c * nrow + r < len(disp)]
            print("".join(n.ljust(w) for n in row).rstrip())
        return
    processes.run(["ls"] + words[1:], cwd=current())


def md(words):
    if len(words) < 2:
        name = input("\nmake_directory> Type the directory name to create: ").strip()
    else:
        name = " ".join(words[1:])
    if not name:
        return
    path = name if os.path.isabs(name) else os.path.join(current(), name)
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as e:
        print(f"\nmake_directory> {e.strerror}: {path}")


def delete(words):
    if len(words) < 2:
        print("\ndelete_file> No files specified to delete. ")
        return
    for pat in words[1:]:
        path = pat if os.path.isabs(pat) else os.path.join(current(), pat)
        matches = sorted(glob.glob(path))
        if not matches:
            print(f"\ndelete_file> No such file: {pat}")
        for m in matches:
            try:
                if os.path.isdir(m):
                    print(f"\ndelete_file> {m} is a folder, not deleted")
                else:
                    os.remove(m)
            except OSError as e:
                print(f"\ndelete_file> {e.strerror}: {m}")


def edit(words):
    """``edit file``: $EDITOR, then ``system_edit`` from system_commands.txt, then nano."""
    fnamefull = ""
    if len(words) > 1:
        f = " ".join(words[1:])
        fnamefull = f if "/" in f else os.path.join(current(), f)
    print(f"\nedit_file> {fnamefull}")
    editor = os.environ.get("VISUAL") or os.environ.get("EDITOR")
    if not editor:
        se = sc.system_edit[: sc.ilensystem_edit].strip().strip('"')
        if se and shutil.which(shlex.split(se)[0]) and se.lower() != "notepad":
            editor = se
    if not editor:
        editor = "nano" if shutil.which("nano") else "vi"
    cmd = shlex.split(editor) + ([fnamefull] if fnamefull else [])
    print(f"\nedit_file> {' '.join(shlex.quote(c) for c in cmd)}")
    processes.run(cmd)


def system(words, text=None):
    """epsrshell_standalone.f:system_command -- replaced: ``system cmd ...`` or ``!cmd ...``
    runs through /bin/sh in the working folder (the Fortran ran it in the start folder)."""
    if text is None:
        text = " ".join(words[1:])
    if not text.strip():
        print("\nsystem_command> No system command specified. ")
        return
    processes.run(text, cwd=current())


def jobs(words):
    live = processes.live()
    if not live:
        print("\nNo child processes running.")
        return
    print()
    for pid, cmd, viewer in live:
        c = cmd if isinstance(cmd, str) else " ".join(cmd)
        print(f"{pid:8d} {'viewer ' if viewer else ''}{c}")


def trace_cmd(words):
    """``trace on [pattern]``, ``trace off``, ``trace summary``."""
    from .. import trace

    sub = words[1] if len(words) > 1 else ""
    if sub == "on":
        pat = words[2] if len(words) > 2 else "*"
        trace.enable(pat, logpath=os.path.join(current(), "trace.log"))
        print(f"\ntrace> call tracing on for '{pat}', writing {current()}/trace.log")
    elif sub == "off":
        trace.disable()
        print("\ntrace> call tracing off")
    elif sub == "summary":
        trace.summary()
    else:
        print("\ntrace> usage: trace on [pattern] | trace off | trace summary")
