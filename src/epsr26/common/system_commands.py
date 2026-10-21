"""Generated from fortran_src/system_commands.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``system_commands.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of system_commands.inc."""
    global system_ls, system_cd, system_termination, system_join, system_gnu, system_jmol, system_binaries, system_bin_ext
    global system_plot2d, system_plot3d, system_plot3djmol, system_pgout, system_md, system_edit, system_del, system_plotato
    global system_quotes, system_mopac, system_mopac_norbitals, ilensystem_ls, ilensystem_cd, ilensystem_termination, ilensystem_join, ilensystem_gnu
    global ilensystem_jmol, ilensystem_binaries, ilensystem_bin_ext, ilensystem_plot2d, ilensystem_plot3d, ilensystem_plot3djmol, ilensystem_pgout, ilensystem_md
    global ilensystem_edit, ilensystem_del, ilensystem_plotato, ilensystem_quotes, ilensystem_mopac, ilensystem_mopac_norbitals, gnu_pointtype, gnu_pointsize
    # common /system_commands_c/
    system_ls = ''
    system_cd = ''
    system_termination = ''
    system_join = ''
    system_gnu = ''
    system_jmol = ''
    system_binaries = ''
    system_bin_ext = ''
    system_plot2d = ''
    system_plot3d = ''
    system_plot3djmol = ''
    system_pgout = ''
    system_md = ''
    system_edit = ''
    system_del = ''
    system_plotato = ''
    system_quotes = ''
    system_mopac = ''
    system_mopac_norbitals = ''
    # common /system_commands_i/
    ilensystem_ls = 0
    ilensystem_cd = 0
    ilensystem_termination = 0
    ilensystem_join = 0
    ilensystem_gnu = 0
    ilensystem_jmol = 0
    ilensystem_binaries = 0
    ilensystem_bin_ext = 0
    ilensystem_plot2d = 0
    ilensystem_plot3d = 0
    ilensystem_plot3djmol = 0
    ilensystem_pgout = 0
    ilensystem_md = 0
    ilensystem_edit = 0
    ilensystem_del = 0
    ilensystem_plotato = 0
    ilensystem_quotes = 0
    ilensystem_mopac = 0
    ilensystem_mopac_norbitals = 0
    gnu_pointtype = 0
    # common /system_commands_r/
    gnu_pointsize = f32(0.0)


reset()
