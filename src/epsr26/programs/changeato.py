"""Translation of fortran_src/changeato.f: ``program standalone`` and
``subroutine changeato`` -- make changes to a .ato file.

``setup_symbols_changeato`` is in setup_symbols_generic.f (library/setup_symbols_generic.py,
imported inside the routine: it is translated by another work package); the interactive
loop is process_symbols.f.
"""

from ..common import clustercom, commands, epsr16, molecom, symbols  # noqa: F401
from ..fortran import Ch, fmt, lw, substr, units
from ..library import FCreadline_windows
from ..library import ato_read_and_save
from ..library import file_and_directory_commands
from ..library import process_symbols
from ..trace import routine

_F599 = "(/'Problems with .ato file ',a/'Routine will exit now'/)"
_F32 = "(/'changeato> There are ',i2,' types of atom in this file'/)"
_F33 = "(1x,'Atom type ',i2,' has label ',a3)"


def main():
    """changeato.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    changeato()


@routine
def changeato():
    """changeato.f:changeato -- make changes to a .ato file."""
    from ..library import setup_symbols_generic

    # use jmk_readline_module
    FCreadline_windows.get_history_list()
    #
    # setup the changeato symbols
    #
    setup_symbols_generic.setup_symbols_changeato()
    #
    # setup the changeato input file (i.e. the .ato file) for reading
    #
    # NB: ierrout is not initialised in the Fortran; setup_input_file sets it
    ierrout = 0
    _nout, ierrout = file_and_directory_commands.setup_input_file(10, 6, ierrout)
    if ierrout != 0:
        return
    #
    # not wanted to leave this file open
    #
    units.close(10)
    #
    # initialise number of atom and molecule types to zero
    #
    molecom.ntype = 0
    clustercom.nmoltype = 0
    #
    # read .ato file - need to do this before setting the symbol tables
    #
    molecom.sizefactor = 1.0
    units.write(6, lw(substr(symbols.fnameinp, 1, symbols.ilenfnameinp)))
    ierr = 0
    (epsr16.rmaxcell, epsr16.cellmax, epsr16.upol, epsr16.rho, epsr16.amubar, epsr16.amubar2, epsr16.ndist,
     ierr) = ato_read_and_save.minitc(symbols.fnameinp, epsr16.rmaxcell, epsr16.cellmax, epsr16.upol,
                                      epsr16.rho, epsr16.amubar, epsr16.amubar2, epsr16.ndist, ierr)
    if ierr != 0:
        units.write(6, fmt(_F599, substr(symbols.fnameinp, 1, symbols.ilenfnameinp)))
        return
    units.write(6, fmt(_F32, molecom.ntype))
    for it in range(1, molecom.ntype + 1):
        units.write(6, fmt(_F33, it, Ch(molecom.type[it], 3)))
    #
    # set the number of "files" in this input file
    #
    symbols.nfiles = molecom.ntype
    #
    # set up the initial values of symbol_value based on the values read from the .ato file
    #
    process_symbols.write_symbols_changeato(1, symbols.nsymbol)
    #
    # read commands from the keyboard and process them
    #
    _nout, ierrout = process_symbols.process_symbols(5, 6, ierrout, 0)

    FCreadline_windows.save_history_list()
    return
