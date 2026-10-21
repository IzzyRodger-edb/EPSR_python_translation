"""The table of EPSR programs and ``run_program``, which runs one of them in-process.

Each Fortran ``program`` is translated into ``programs/<fortran file>.py`` with a
``main()`` function.  ``run_program`` gives it what a fresh Fortran process had: new
COMMON storage (``common.reset()``), no open units, and its command-line arguments.
"""

import importlib
import sys
import time

# program name (as typed in epsrshell / the Fortran binary name) -> translated module
PROGRAMS = {
    "epsr": "epsr_standalone_rev3",
    "partials": "partials",
    "coord": "coord",
    "triangles": "triangles",
    "torangles": "torangles",
    "nwts": "nwts",
    "xwts": "xwts",
    "epsrwts": "epsrwts",
    "epsrwtsx": "epsrwtsx",
    "makeato": "makeato",
    "mixato": "mixato",
    "fmole": "fmole",
    "changeato": "changeato",
    "voids": "voids",
    "clusters": "clusters",
    "chains": "chains_sp_criterion",
    "rings": "rings_sp_criterion",
    "lsi": "lsi",
    "sdf": "sdf_rev2",
    "sharm": "sharm_rev2",
    "sdfcube": "sdfcube",
    "mapgrentropy": "mapgrentropy",
    "mapgr": "mapgrquad",
    "writexyz": "writexyz",
    "readxyz": "readxyz",
    "readpdb": "readpdb",
    "fluctuations": "fluctuations",
    "addato": "addato",
    "dockato": "dockato",
    "makelattice": "makelattice",
    "makelatticeato": "makelatticeato",
    "makemole": "makemole",
    "randomise": "introtcluster",
    "bonds": "bondsclus",
    "local_density": "local_density",
    "alignmole": "alignmole",
    "upset": "setup_generic",
    "plot": "plot",
    "plot3djmol": "plot3djmol",
    "runjmol": "runjmol",
    # binary names the shell's run_program uses after translating the command
    "mapgrquad": "mapgrquad",
}


class NotTranslated(Exception):
    pass


def load(name):
    mod = PROGRAMS.get(name)
    if mod is None:
        raise NotTranslated(f"no program {name!r}")
    try:
        return importlib.import_module(f"{__name__}.{mod}")
    except ModuleNotFoundError as e:
        if e.name and e.name.endswith(mod):
            raise NotTranslated(f"program {name!r} ({mod}.py) is not translated yet") from None
        raise


def run_program(name, args):
    """Run program ``name`` with Fortran command-line ``args``; returns the exit status.

    ``FortranStop`` ends the program normally (its code is the status); the files it
    wrote are flushed in every case, as the Fortran runtime does at exit.
    """
    from .. import common, trace
    from ..fortran import runtime, units
    from ..fortran.stop import FortranStop

    module = load(name)
    trace.refresh()
    runtime.set_args(args)
    common.reset()
    units.reset()
    trace.context["iter"] = None
    if trace.is_enabled():
        trace.reset_counts()
    status = 0
    t0 = time.perf_counter()
    try:
        module.main()
    except FortranStop as e:
        if e.message:
            sys.stderr.write(f"STOP {e.message}\n")
        elif e.code:
            sys.stderr.write(f"STOP {e.code}\n")
        status = e.code if isinstance(e.code, int) else 0
    finally:
        units.close_all()
        runtime.last_run_seconds = time.perf_counter() - t0
    return status
