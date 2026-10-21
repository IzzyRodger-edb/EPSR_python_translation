"""Translation of fortran_src/setup_generic.f: ``program setup_standalone`` (the
``upset`` program, run by the shell's ``setup`` command) and ``subroutine setup``.

Sets up and modifies the input files for EPSR and the auxiliary programs.  The
setup_symbols_* routines come from library/setup_symbols_generic.py,
library/setup_epsr_rev1.py (``setup_symbols_epsr``); ``read_plot3djmol`` from
library/plot_routines_g95.py (no active Fortran file defines it, see that module).
"""

from ..common import commands, symbols
from ..fortran import FortranIOError, FortranStop, fmt, fstr, lw, substr, units
from ..library import file_and_directory_commands
from ..trace import routine

# 2007-09-30. Modified read_pcof and write_pcof to include Gaussian values
# as well as minimum allowed distances

_F100 = (
    "(/1x,'What do you want to setup? Options are:-'//"
    "1x,'epsr <filename>  - to create/edit EPSR input filess'/"
    "1x,'epsr17 <filename> - to read data from EPSR17 input files'/"
    "1x,'epsr16 <filename> - to read data from EPSR16 input files'/"
    "1x,'bragginten <filename> - to create/edit bragginten input files'/"
    "1x,'partials <filename> - to create/edit partials input files'/"
    "1x,'clusters <filename> - to create/edit clusters input files'/"
    "1x,'voids <filename> - to create/edit voids input files'/"
    "1x,'coord <filename> - to create/edit coord input files'/"
    "1x,'chains <filename> - to create/edit chains input files'/"
    "1x,'rings <filename> - to create/edit rings input files'/"
    "1x,'lsi <filename> - to create/edit lsi input files'/"
    "1x,'triangles <filename> - to create/edit triangles input files'/"
    "1x,'torangles <filename> - to create/edit torangles input files'/"
    "1x,'sharm <filename> - to create etc. SHARM input files'/"
    "1x,'sdf <filename> - to create etc. SDF input files'/"
    "1x,'sdfcube <filename> - to create etc. SDFCUBE input files'/"
    "1x,'fluctuations <filename> - to create fluctuations input files'/"
    "1x,'nwts <filename> - to create neutron weights input files'/"
    "1x,'xwts <filename> - to create x-ray weights input files'/"
    "1x,'mapgrentropy <filename>  - to create/edit mapgrentropy input filess'/"
    "1x,'q - to quit setup'/"
    "1x,'plot2d <filename> - to create/edit input file for 2D plotting'/"
    "1x,'splot2d <filename> - to create/edit input file for 2D plotting with gnuplot'/"
    "1x,'plot3d <filename> - to create/edit input file for 3D plotting'/"
    "1x,'plot3djmol <filename> - to create/edit input file for 3D plotting with Jmol'/"
    ")"
)
_F101 = "(a)"
_F114 = "(/'setup> ')"


def main():
    """setup_generic.f:setup_standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    setup()
    raise FortranStop(0)


@routine
def setup():
    """setup_generic.f:setup -- sets up and modifies the various input files for EPSR and
    associated programs."""
    from ..library import FCreadline_windows, parse, plot_routines, plot_routines_g95
    from ..library import process_symbols, read_and_save_files, setup_epsr_rev1
    from ..library import setup_symbols_generic as ssg
    from ..library import stringutility

    # use String_Utility
    # use jmk_readline_module
    # character*20 command
    cm = commands
    sy = symbols
    # NB: ierrout is not initialised in the Fortran when setup_input_file is skipped
    # (epsr16/epsr17); 0 is used here.
    ierrout = 0
    FCreadline_windows.get_history_list()
    #
    # first determine the type of setup required
    #
    if cm.nwords < 2:
        units.write(6, fmt(_F100))
        cm.nwords = 0
        while cm.nwords == 0:
            cm.sometext = ''
            #            write(6,114,advance='no')
            #            read(5,101,iostat=ierr) sometext
            cm.sometext = FCreadline_windows.readline(cm.sometext, 'setup> ', len_line=cm.nstsize)
            cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
        nfirst = 1
    else:
        nfirst = 2
    command = fstr(substr(cm.sometext, cm.ncf[nfirst], cm.ncl[nfirst]), 20)
    command = fstr(stringutility.strlowcase(command), 20)
    sy.symbol_setup = command
    sy.ilensetup = len(sy.symbol_setup)
    #
    # first setup the symbols
    #
    if command == 'epsr':
        setup_epsr_rev1.setup_symbols_epsr()
    elif command == 'mapgr':
        ssg.setup_symbols_mapgr()
    elif command == 'mapgrentropy':
        ssg.setup_symbols_mapgrentropy()
    elif command == 'bragginten':
        # set up the clusters symbol tables and default values
        ssg.setup_symbols_bragginten()
    elif command == 'partials':
        # set up the partials symbol tables and default values
        ssg.setup_symbols_partials()
    elif command == 'clusters':
        # set up the clusters symbol tables and default values
        ssg.setup_symbols_clusters()
    elif command == 'coord':
        # set up the coord symbol tables and default values
        ssg.setup_symbols_coord()
    elif command == 'chains':
        ssg.setup_symbols_chains()
    elif command == 'rings':
        ssg.setup_symbols_rings()
    elif command == 'lsi':
        ssg.setup_symbols_lsi()
    elif command == 'sharm':
        # set up the SHARM symbol tables and default values
        ssg.setup_symbols_sharm()
    elif command == 'sdf':
        # set up the SDF symbol tables and default values
        ssg.setup_symbols_sdf()
    elif command == 'sdfcube':
        # set up the SDFCUBE symbol tables and default values
        ssg.setup_symbols_sdfcube()
    elif command == 'plot2d':
        # set up the plot2d symbol tables and default values
        ssg.setup_symbols_plot2d()
    elif command == 'splot2d':
        ssg.setup_symbols_splot2d()
    elif command == 'plot3d':
        # set up the plot3d symbol tables and default values
        ssg.setup_symbols_plot3d()
    elif command == 'plot3djmol':
        # set up the plot3djmol symbol tables and default values
        ssg.setup_symbols_plot3djmol()
    elif command == 'triangles':
        ssg.setup_symbols_triangles()
    elif command == 'torangles':
        ssg.setup_symbols_torangles()
    elif command == 'voids':
        ssg.setup_symbols_voids()
    elif command == 'fluctuations':
        ssg.setup_symbols_fluctuations()
    elif command == 'nwts' or command == 'epsrwts':
        ssg.setup_symbols_nwts()
    elif command == 'xwts':
        ssg.setup_symbols_xwts()
    elif command == 'q':
        #
        # q means quit without doing anything.
        #
        return
    else:
        # NB: ssg_ierr is implicitly REAL in the Fortran and passed to an integer dummy;
        # any nonzero integer result is a nonzero bit pattern, so the test is the same.
        ssg_ierr = 0
        ssg_ierr = ssg.setup_symbols_generic(command, ssg_ierr)
        if ssg_ierr != 0:
            units.write(6, fmt(_F114), advance=False)
            units.write(6, lw("Problem setting up symbols", " from .symbols file. Exiting setup."))
            return
    #
    # now setup the input file
    #
    nin = 11
    nout = 6
    if command != 'epsr16' and command != 'epsr17':
        nout, ierrout = file_and_directory_commands.setup_input_file(nin, nout, ierrout)
    #
    # ierrout .lt. 0 signifies the file could not be opened or on reading encountered an
    # error, but an input file has been successfully defined using the default values,
    # so these values are now available for viewing and changing.
    #
    if command == 'plot2d':
        nout, ierrout = plot_routines.read_plot2d(nin, nout, ierrout)
        if ierrout <= 0:
            _, ierrout = plot_routines.process_plot_symbols(5, 6, ierrout, 0)
    elif command == 'splot2d':
        nout, ierrout = plot_routines.read_plot_input_file(nin, nout, ierrout)
        if ierrout <= 0:
            _, ierrout = plot_routines.process_plot_symbols(5, 6, ierrout, 0)
    elif command == 'plot3d':
        nout, ierrout = plot_routines.read_plot3d(nin, nout, ierrout)
        if ierrout <= 0:
            _, ierrout = plot_routines.process_plot_symbols(5, 6, ierrout, 0)
    elif command == 'plot3djmol':
        nout, ierrout = plot_routines.read_plot_input_file(nin, nout, ierrout)
        if ierrout != 0:
            ierrout = 0
            # Re-open the file that would have been closed by read_input_file
            try:
                units.open(nin, sy.fnameinp, status='old')
            except FortranIOError as e:
                ierrout = e.iostat
            nout, ierrout = plot_routines_g95.read_plot3djmol(nin, nout, ierrout)
        _, ierrout = plot_routines.process_plot_symbols(5, 6, ierrout, 0)
    else:
        nout, ierrout = read_and_save_files.read_input_file(nin, nout, ierrout)
        if ierrout <= 0:
            _, ierrout = process_symbols.process_symbols(5, 6, ierrout, 0)
    FCreadline_windows.save_history_list()
    return
