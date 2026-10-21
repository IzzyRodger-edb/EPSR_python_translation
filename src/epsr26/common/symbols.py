"""Generated from fortran_src/symbols.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``symbols.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401
from .epsr_dimension import mfile, mpairs  # noqa: F401

msymbol = 1146  # mpairs+36+10*(mfile-1)

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of symbols.inc."""
    global symbol, symbol_value, symbol_default, symbol_save, symbol_type, symbol_text, fname_inp, fnameinp
    global fname_run, fnamerun, fname_ato, fnameato, fname_pcof, fnamepcof, fname_template, fnametemplate
    global prompt, dataprompt, datatext, dataprompttext, extinp, extdot, symbol_setup, nsymboluptodata
    global nsymboldata, nsymbol, nsymbolplus, iflag, isymboltype, ilenfname_inp, ilenfnameinp, ilenfname_run
    global ilenfnamerun, ilenfname_ato, ilenfnameato, ilenfname_pcof, ilenfnamepcof, ilenfname_template, ilenfnametemplate, nfiles
    global ilenprompt, ilendataprompt, ilendatatext, ilenextinp, ilenextinpdot, ilenextdot, ilensetup, isave
    global inotrecog, irefrho
    # common /symbol_c/
    symbol = char_array(msymbol)
    symbol_value = char_array(msymbol)
    symbol_default = char_array(msymbol)
    symbol_save = ''
    symbol_type = char_array(msymbol)
    symbol_text = char_array(msymbol)
    fname_inp = ''
    fnameinp = ''
    fname_run = ''
    fnamerun = ''
    fname_ato = ''
    fnameato = ''
    fname_pcof = ''
    fnamepcof = ''
    fname_template = ''
    fnametemplate = ''
    prompt = ''
    dataprompt = ''
    datatext = ''
    dataprompttext = ''
    extinp = ''
    extdot = ''
    symbol_setup = ''
    # common /symbol_i/
    nsymboluptodata = 0
    nsymboldata = 0
    nsymbol = 0
    nsymbolplus = 0
    iflag = int_array(msymbol)
    isymboltype = 0
    ilenfname_inp = 0
    ilenfnameinp = 0
    ilenfname_run = 0
    ilenfnamerun = 0
    ilenfname_ato = 0
    ilenfnameato = 0
    ilenfname_pcof = 0
    ilenfnamepcof = 0
    ilenfname_template = 0
    ilenfnametemplate = 0
    nfiles = 0
    ilenprompt = 0
    ilendataprompt = 0
    ilendatatext = 0
    ilenextinp = 0
    ilenextinpdot = 0
    ilenextdot = 0
    ilensetup = 0
    isave = 0
    inotrecog = 0
    irefrho = 0


reset()
