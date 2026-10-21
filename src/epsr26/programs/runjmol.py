"""Translation of fortran_src/runjmol.f: the ``runjmol`` command (start Jmol in the
working folder, using ``system_cd``, ``system_join`` and ``system_jmol`` from
system_commands.txt)."""

from ..common import directory_name, system_commands
from ..fortran import fstr, lw, runtime, substr, units
from ..library import file_and_directory_commands
from ..trace import routine


@routine
def main():
    """runjmol.f:standalone -- program to perform the EPSR runjmol command."""
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_jmol()


@routine
def run_jmol():
    """runjmol.f:run_jmol -- generate the command that runs Jmol and run it."""
    sc = system_commands
    # now generate the command that will actually run this file
    plot_command = fstr(
        substr(sc.system_cd, 1, sc.ilensystem_cd)
        + " "
        + substr(sc.system_quotes, 1, sc.ilensystem_quotes)
        + substr(directory_name.directory_name, 1, directory_name.nchardir)
        + substr(sc.system_quotes, 1, sc.ilensystem_quotes)
        + " " + substr(sc.system_join, 1, sc.ilensystem_join)
        + " " + substr(sc.system_jmol, 1, sc.ilensystem_jmol),
        256,
    )
    units.write(6, lw(plot_command))
    ierr = runtime.system(plot_command)
    units.write(6, lw(ierr))
