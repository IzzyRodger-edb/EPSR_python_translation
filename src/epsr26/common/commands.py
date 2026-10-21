"""Generated from fortran_src/commands.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``commands.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401

ncommand_max = 100
ncommand_line_max = 200
nworddim = 201
nstsize = 3020

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of commands.inc."""
    global command_text, runflag, runflagglobal, sometext, ncommand, nwords, ncf, ncl
    global ilencommand_text, script_fname, script_status, ncharscript, ncommand_stream
    # common /commands_c/
    command_text = ''
    runflag = ''
    runflagglobal = ''
    sometext = ''
    # common /commands_i/
    ncommand = 0
    nwords = 0
    ncf = int_array(nworddim)
    ncl = int_array(nworddim)
    ilencommand_text = 0
    # common /script_c/
    script_fname = ''
    script_status = ''
    # common /script_i/
    ncharscript = 0
    # common /command_line_i/
    ncommand_stream = 0


reset()
