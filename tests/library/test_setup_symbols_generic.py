"""library/setup_symbols_generic.py against gfortran.

data/setup_symbols_reference.txt is the output of data/setup_symbols_driver.f linked with
fortran_src/setup_symbols_generic.f and parse.f (build command in the driver header), run
in data/setup_symbols (the .symbols files used by setup_symbols_generic).  The driver
calls every keyword-table routine, then setup_symbols_generic on good and bad .symbols
files, dumping the symbol table after each call; this test does the same in Python and
compares everything written to unit 6 byte for byte.
"""

import contextlib
import io
import os
import pathlib
import shutil
import subprocess

import pytest

from epsr26 import common
from epsr26.common import commands, symbols
from epsr26.fortran import fmt, units
from epsr26.library import setup_symbols_generic as ssg

DATA = pathlib.Path(__file__).parent / "data"
ROOT = pathlib.Path(__file__).parents[2]

TABLES = ["changeato", "partials", "coord", "clusters", "rings", "chains", "lsi", "sharm", "sdf",
          "sdfcube", "mapgrentropy", "plot2d", "splot2d", "plot3d", "plot3djmol", "mapgr",
          "triangles", "torangles", "voids", "fluctuations", "nwts", "xwts", "bragginten",
          "sharm"]
CASES = ["good extra words", "nondist", "noeol", "noprompt", "badprompt", "noext", "nodataprompt",
         "nodptext", "blankinside", "badkey", "eofinside", "nodefault", "notype", "notext",
         "missing", "  good"]


def _w(s):
    units.write(6, s)


def _dump():
    s = symbols
    _w(fmt("(a,6i6)", "n ", s.nsymbol, s.nsymboldata, s.nsymboluptodata, s.nsymbolplus, s.nfiles,
           s.irefrho))
    _w(fmt("(a,i6,a)", "prompt ", s.ilenprompt, " |" + s.prompt + "|"))
    _w(fmt("(a,i6,a)", "dataprompt ", s.ilendataprompt, " |" + s.dataprompt + "|"))
    _w(fmt("(a)", "dataprompttext |" + s.dataprompttext + "|"))
    _w(fmt("(a,i6,a)", "extinp ", s.ilenextinp, " |" + s.extinp + "|"))
    for i in range(1, s.nsymbol + 1):
        _w(fmt("(i4,1x,i3,1x,a,a,a,a,a,a,a,a,a,a)", i, s.iflag[i], s.symbol[i].ljust(11), "|",
               s.symbol_type[i].ljust(1), "|", s.symbol_value[i], "|", s.symbol_default[i], "|",
               s.symbol_text[i], "|"))


def _run_python():
    common.reset()
    units.reset()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        for name in TABLES:
            _w("=== " + name)
            getattr(ssg, "setup_symbols_" + name)()
            _dump()
        for case in CASES:
            _w("=== generic " + case.rstrip())
            ierr = ssg.setup_symbols_generic(case, 99)
            _w(fmt("(a,i4)", "ssg_ierr ", ierr))
            _w(fmt("(a,3i6)", "commands ", commands.nwords, commands.ncf[1], commands.ncl[1]))
            _dump()
        units.close_all()
    units.reset()
    return buf.getvalue()


def test_matches_gfortran_reference(monkeypatch):
    monkeypatch.chdir(DATA / "setup_symbols")
    got = _run_python()
    ref = (DATA / "setup_symbols_reference.txt").read_text()
    assert got.splitlines() == ref.splitlines()
    assert got == ref


@pytest.mark.skipif(shutil.which("gfortran") is None, reason="gfortran not available")
def test_reference_is_current(tmp_path):
    """Rebuild the driver and check the stored reference is what gfortran writes now."""
    f = ROOT / "fortran_src"
    exe = tmp_path / "ssg"
    subprocess.run(["gfortran", "-O2", "-ffixed-line-length-none", "-fno-backslash", "-w", f"-I{f}",
                    "-o", str(exe), str(DATA / "setup_symbols_driver.f"),
                    str(f / "setup_symbols_generic.f"), str(f / "parse.f")],
                   check=True, cwd=tmp_path)
    out = subprocess.run([str(exe)], check=True, capture_output=True, text=True,
                         cwd=DATA / "setup_symbols", env={**os.environ, "OMP_NUM_THREADS": "1"})
    assert out.stdout == (DATA / "setup_symbols_reference.txt").read_text()
