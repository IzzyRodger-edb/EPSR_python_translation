"""Translation of fortran_src/plot_routines.f: the ``plot`` command of the EPSR shell.

Sets up the plot symbols, reads/writes ``plot_defaults.txt``, processes the plot
commands, and writes a gnuplot script ``plotNN<file>.gnu`` which is then run with
``call system(system_gnu ... ' -')`` -- in the foreground (``runtime.system``), as the
Fortran does (no ``&``).

Conventions used throughout:

* ``commands.sometext`` (character*3020) is the shared text buffer; ``_parse(nspc)`` is
  ``call parse(sometext,nwords,ncf,ncl,nspc,nstsize,nworddim)`` on the COMMON variables.
* ``FORMAT`` statements ending in ``$`` (no advance) are written with ``advance=False``
  and the ``$`` removed from the format string.
* ``read(..., iostat=ierr)`` uses the private helpers ``_read_list`` (internal
  list-directed read) and ``_read_line`` (``read(u,'(a)')``); on an error the items not
  yet read keep their old values, and on end-of-file a ``'(a)'`` read leaves the
  variable unchanged, as gfortran does.
"""

from ..common import commands, directory_name, plot_pars, symbols, system_commands
from ..fortran import (
    Ch,
    FortranEOF,
    FortranIOError,
    f32,
    findex,
    findex_padded,
    fmt,
    fstr,
    lw,
    pad,
    r4,
    read_internal,
    runtime,
    substr,
    units,
)
from ..trace import routine
from . import parse

# --------------------------------------------------------------------------- FORMATs
# process_plots
_PP100 = "(/a,1x,'Available plot types are:-'//)"
_PP101 = "(1x,i3,'  ',a)"
_PP102 = "(1x,i3)"
_PP104 = "(/a,1x,' Try again (y or n)? ')"  # ,$
# setup_plot_filenames
_SPF102 = "(/a,1x,'Specify a valid plot type between',i2,' and ',i2)"
_SPF103 = "(/a,1x,'No extensions supplied for plottype ',i3)"
_SPF104 = "(/a,1x,'Block offset and spacing not defined correctly')"
_SPF125 = "(/'setup_plot_filenames> This file is empty!')"
_SPF126 = "(/'setup_plot_filenames> There are ',i3,' blocks in the file ',a)"
_SPF127 = "(/'setup_plot_filenames> ',a)"
# find_ncolumn
_FNC125 = "(/'find_ncolumn> ',2(1x,i5))"
_FNC126 = "(/'find_ncolumn> ',6(1x,i5))"
# setup_blocknumbers
_SBN103 = "(/'setup_blocknumbers> A valid filename must be defined before ','block numbers can be set.')"
_SBN100 = ("(/'setup_blocknumbers> There are ',i3,' blocks. Specify which ones to plot'"
           "//'(Use spaces between numbers, start with a number,'/' and use n1 - n2 to define a range): ')")  # ,$
_SBN102 = "(/a,1x,'No column numbers specified - try again: ')"  # ,$
_SBN235 = "(/'setup_blocknumbers> Number of plotting columns: ',i5)"
_SBN236 = "(/'setup_blocknumbers> Need to define valid block numbers for this file')"
# gnuplot
_GP122 = "(/'GNUplot> Block numbers need be defined before ','data can be plot.')"
_GP150 = "(i2.2)"
_GP101 = "('reset',/,'set title ',1h',a,1h')"
_GP102 = "('set xlabel ',1h',a,1h')"
_GP103 = "('set ylabel ',1h',a,1h')"
_GP104 = "('set xrange [',e10.3,':',e10.3,']')"
_GP105 = "('set yrange [',e10.3,':',e10.3,']')"
_GP113 = "('set logscale ',a1)"
_GP114 = "('set nologscale ',a1)"
_GP106 = "('set label ',i3,1x,1h',a,1h',' at ',e10.3,',',e10.3)"
_GP107 = "('set style line',i3,' lt ',i3,' pt ',i3,' ps ',f7.2)"
_GP108 = "('plot \\\\')"  # gfortran -fno-backslash: two backslashes reach the file
_GP109 = "(1h',a,1h',' u ',i3,':(',f10.5,'*column(',i3,')+',f8.3,') notitle w ',a,' ls ',i3)"  # ,$
_GP110 = "(', \\\\')"
_GP111 = ("(1h',a,1h',' u ',i3,':(abs(column(',i3,'))>=1.0e-5?',f10.5,'*column(',i3,')+',f8.3"
          ",':1/0) notitle w p ls ',i3)")  # ,$
# process_plot_symbols
_PPS203 = "(/a,' ',a,' - ',a)"
_PPS202 = "(/a,': ',a,' ? ')"  # ,$
_PPS201 = "(/a,' ')"  # ,$
_PPS401 = "(/a,' Current data have not been saved.'//'Press <Enter> to save, or q to exit without saving: ')"  # ,$
# process_plot_commands
_PPC250 = "(i3)"
_PPC204 = "(//a,' ',a,' - ',a/)"
_PPC207 = "(/'List of parameters:-'//)"
_PPC209 = "(a10,' - ',a65)"
_PPC210 = "(/'Other commands:-'/)"
_PPC211 = "(a20,'-',1x,a45)"
_PPC208 = ("('ls or dir <filename> - Lists files in the current directory'"
           "/'cd <directory name>  - Changes the current directory'"
           "/'pwd                  - Show the current directory'"
           "/'save                 - Saves the current values to file'"
           "/'default              - Sets numerical values to their defaults'"
           "/'help or ?            - Show this list'"
           "/'<value>              - New value for the listed parameter'"
           "/'<parameter> <value>  - Specified parameter and its new value'"
           "/'n                    - Revert to the previous value'"
           "/'<Enter> or d            - Accept current parameter and go forwards'"
           "/'u                    - Accept current parameter and ','go backwards'"
           "/'q                    - Return to EPSR shell without saving'"
           "/'bye or e             - Return to EPSR shell'"
           "/)")
_PPC100 = "(/'process_commands> Command not recognised: ',a,': ',a)"
# check_plot_symbol
_CPS200 = "(/a,1x,'Cannot change this value!')"
_CPS201 = "(/a,1x,'Value ',a,' is not an integer - try again!')"
_CPS202 = "(/a,1x,'Value ',a,' is not a real number - try again!')"
_CPS203 = "(/a,1x,'ERROR - tried to generate too many symbol values')"
_CPS206 = "(/a,1x,'No coefficients file specified! - routine will exit'/)"
_CPS207 = "(/a,1x,'setone must be 0 or 1!')"
# find_plot_par
_FPP100 = ("(/a,' SHARM input file ',a,' does not exist'"
           "//'Will try to get number of coefficients directly',' from coefficients file')")
_FPP103 = ("(/a,' Problems with file ',a"
           "//'Will try to get number of coefficients directly',' from coefficients file')")
# save_plot_files
_SVP303 = "(/a,1x,'Current name of file is \"',a,'\"')"
_SVP302 = "(/a,1x,' Specify the name of the ',a,' input file: ')"  # ,$
_SVP304 = "(/a,1x,'Writing to input file \"',a,'\"')"
_SVP301 = "(/a,1x,'File \"',a,'\" already exists.'//' Do you want to overwrite it (y or n)? ')"  # ,$
# save_files_plots
_SFP229 = "(a,15x,'Title of this file')"
_SFP204 = "(/a,1x,i3/)"
_SFP101 = "(a11,1x,a,15x,a)"
_SFP102 = "('q')"
_SFP103 = "(/'save_files_plots> Could not open specified file for writing.')"
# save_files_plot2d / save_files_plot3d
_SF2D102 = "(a,11x)"
_SF3D102 = "(a,11x,'sphere radius, (r,theta,phi), (r,g,b colour indices)')"
_SF103 = "(a)"
_SF101 = "(a,15x,a)"


# --------------------------------------------------------------------------- helpers
def _parse(nspc):
    """``call parse(sometext,nwords,ncf,ncl,nspc,nstsize,nworddim)`` on the commands COMMON."""
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl,
                                  nspc, commands.nstsize, commands.nworddim)


def _word(i):
    """``sometext(ncf(i):ncl(i))``."""
    return substr(commands.sometext, commands.ncf[i], commands.ncl[i])


def _words(i, j):
    """``sometext(ncf(i):ncl(j))``."""
    return substr(commands.sometext, commands.ncf[i], commands.ncl[j])


def _read_list(text, kinds, olds):
    """``read(text,*,iostat=ios) v1, v2, ...``: kinds is a string of 'i'/'r'.

    Returns (values, ios); items read before an error or end of record keep the new
    value, the others their old one.
    """
    vals = list(olds)
    r = read_internal(text)
    try:
        for k, kind in enumerate(kinds):
            if kind == 'i':
                vals[k] = r.int(vals[k])
            else:
                vals[k] = r.real(vals[k])
    except FortranEOF:
        return vals, -1
    except FortranIOError as e:
        return vals, e.iostat
    return vals, 0


def _read_line(unit, old, n):
    """``read(unit,'(a)',iostat=ios) s`` for a character*n s -> (s, ios)."""
    try:
        return fstr(units.read_line(unit), n), 0
    except FortranEOF:
        return old, -1
    except FortranIOError as e:
        return old, e.iostat


def _open(unit, file, status):
    """``open(unit,file=file,status=status,iostat=ios)`` -> ios."""
    try:
        units.open(unit, file, status=status)
    except FortranIOError as e:
        return e.iostat
    return 0


def _prompt():
    """``prompt(1:ilenprompt)``."""
    return substr(symbols.prompt, 1, symbols.ilenprompt)


# --------------------------------------------------------------------------- routines
@routine
def plot_new():
    """plot_routines.f:plot_new -- the ``plot`` command: set up the plot symbols, read
    the plot defaults file and process the plot commands.
    """
    from . import file_and_directory_commands

    # NB: ierrout is not initialised in the Fortran; setup_input_file sets it.
    ierrout = 0
    #
    # set up plot symbols
    #
    setup_symbols_plot()
    #
    # now attempt to open the file containing the default plot values
    #
    commands.sometext = fstr('plot ' + substr(symbols.extinp, 1, symbols.ilenextinp), 3020)
    nin = 11
    nout = 6
    nout, ierrout = file_and_directory_commands.setup_input_file(nin, nout, ierrout)
    #
    # now read it
    #
    if ierrout == 0:
        nout, ierrout = read_plot_input_file(nin, nout, ierrout)
    #
    # ierrout .lt. 0 signifies the file could not be opened or on reading encountered an
    # error, but an input file has been successfully defined using the default values,
    # so these values are now available for viewing and changing.
    #
    _, ierrout = process_plot_symbols(5, 6, ierrout, 0)
    return


def _set_symbol(iref, sym, typ, text, flag, default):
    """The five assignments repeated for every symbol in setup_symbols_plot."""
    symbols.symbol[iref] = fstr(sym, 11)
    symbols.symbol_type[iref] = fstr(typ, 1)
    symbols.symbol_text[iref] = fstr(text, 256)
    symbols.iflag[iref] = flag
    symbols.symbol_default[iref] = fstr(default, 256)


@routine
def setup_symbols_plot():
    """plot_routines.f:setup_symbols_plot -- the plot symbols and their default values.

    iflag decides whether this parameter can be changed
    a parameter can only be changed if ireset >= iflag
    """
    iref = 0

    iref = iref + 1
    _set_symbol(iref, 'l', 'c', 'Lists available plot types', 0, '0')

    iref = iref + 1
    _set_symbol(iref, 'p', 'i', 'Plot using the current or specified plot type', 0, ' ')
    plot_pars.icountp = iref

    iref = iref + 1
    _set_symbol(iref, 'npt', 'i', 'Number of types of plot', 0, '1')

    symbols.nsymboluptodata = iref

    iref = iref + 1
    _set_symbol(iref, 'type', 'c', 'Type of plot', 1, '<undefined>')

    iref = iref + 1
    _set_symbol(iref, 'ext', 'c', 'File extension(s) used to search for plot files', 0, '<undefined>')

    iref = iref + 1
    _set_symbol(iref, 'f', 'c', 'File name to plot', 0, '<undefined>')
    plot_pars.icountf = iref

    iref = iref + 1
    _set_symbol(iref, 'b', 'c', 'Block numbers to plot (e.g. 1 2 - 5 9 - 6)', 0, '<undefined>')
    plot_pars.icountb = iref

    iref = iref + 1
    _set_symbol(iref, 'bspace', 'i', 'Spacing between blocks in plot file', 0, '2')

    iref = iref + 1
    _set_symbol(iref, 'boffset', 'i', 'Column number(s) for first block(s) in plot file(s)', 0, '2 2')

    iref = iref + 1
    _set_symbol(iref, 'xmin', 'r', 'Minimum value on x-axis', 0, '0')

    iref = iref + 1
    _set_symbol(iref, 'xmax', 'r', 'Maximum value on x-axis', 0, '1')

    iref = iref + 1
    _set_symbol(iref, 'ymin', 'r', 'Minimum value on y-axis', 0, '0')

    iref = iref + 1
    _set_symbol(iref, 'ymax', 'r', 'Maximum value on y-axis', 0, '1')

    iref = iref + 1
    _set_symbol(iref, 'ydel', 'r', 'Spacing between plots', 0, '1')

    iref = iref + 1
    _set_symbol(iref, 'yfact', 'r', 'Factor(s) outside y values', 0, '1.0 1.0')

    iref = iref + 1
    _set_symbol(iref, 'xf', 'r', 'Fractional x coordinate of labels', 0, '0.75')

    iref = iref + 1
    _set_symbol(iref, 'yf', 'r', 'Fractional y coordinate of labels', 0, '0.25')

    iref = iref + 1
    _set_symbol(iref, 'title', 'c', 'Title of plot', 0, '<undefined>')

    iref = iref + 1
    _set_symbol(iref, 'xlabel', 'c', 'x-axis label', 0, '<undefined>')

    iref = iref + 1
    _set_symbol(iref, 'ylabel', 'c', 'y-axis label', 0, '<undefined>')

    iref = iref + 1
    _set_symbol(iref, 'logx', 'c', 'log scale x (yes or no)', 0, 'no')

    iref = iref + 1
    _set_symbol(iref, 'logy', 'c', 'log scale y (yes or no)', 0, 'no')

    iref = iref + 1
    _set_symbol(iref, 'xcolumn', 'i', 'column number for x values [normally 1 or 0]', 0, '1')

    iref = iref + 1
    _set_symbol(iref, 'histogram', 'c', 'histogram data (yes or no)', 0, 'no')

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    symbols.nsymbol = iref
    #
    # initialise the symbol values
    #
    for icount in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[icount] = symbols.symbol_default[icount]
    #
    # initialise number of plot types
    #
    symbols.nfiles = 1
    #
    # setup the output prompt
    #
    symbols.prompt = fstr('plot>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')
    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('pt', 80)
    symbols.dataprompttext = fstr('Sets the specified plot type ', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1
    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('plot_defaults.txt', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1
    return


@routine
def process_plots(nin, nout, icount, ifilecount, ierr):
    """plot_routines.f:process_plots -- plots a specified set of results from EPSR or its
    associated routines (commands l, f, b, p of the plot symbols).

    Returns (nout, icount, ifilecount, ierr).
    """
    from . import file_and_directory_commands

    #
    # first determine the requested command
    #
    command = fstr(symbols.symbol[icount], 20)
    #
    # save the associated command text
    #
    iref = file_and_directory_commands.ireff(icount, ifilecount)
    input_text = fstr(symbols.symbol_value[iref], 256)  # noqa: F841
    #     write(6,*) iref,' ',input_text(1:50)
    if command == 'l':
        #
        # save the reference integer
        #
        irefl = iref
        units.write(nout, fmt(_PP100, _prompt()))
        #
        # set up the reference integer to get the symbol_values
        #
        icountwr = symbols.nsymboluptodata + 1
        ifile = 0
        # NB: iplottype is not initialised in the Fortran; it is set before it is written
        iplottype = 0
        while ifile < symbols.nfiles:
            ifile = ifile + 1
            iref = file_and_directory_commands.ireff(icountwr, ifile)
            commands.sometext = fstr(symbols.symbol_value[iref], 3020)
            _parse(5)
            #
            # attempt to read the first word as a type number
            #
            if commands.nwords > 0:
                (iplottype,), ierr = _read_list(_word(1), 'i', (iplottype,))
                if ierr != 0:
                    #
                    # if it isn't set a plot type number
                    #
                    iplottype = ifile
                    ibegin = 1
                elif commands.nwords > 1:
                    ibegin = 2
                else:
                    ibegin = 0
            else:
                ibegin = 0
            if ibegin > 0:
                units.write(nout, fmt(_PP101, iplottype, _words(ibegin, max(ibegin, commands.nwords))))
            else:
                units.write(nout, fmt(_PP102, iplottype))
        symbols.symbol_value[irefl] = ''
    elif command == 'f':
        nout, ierr = setup_plot_filenames(nout, ierr)
        if ierr != 0:
            units.write(nout, fmt(_PP104, _prompt()), advance=False)
            # NB: yorn is not initialised in the Fortran (it is read before use)
            yorn = ''
            r = units.list_read(nin)
            yorn = fstr(r.word(yorn), 1)
            if yorn == 'y':
                symbols.symbol_value[plot_pars.icountf] = 'search'
                nout, ierr = setup_plot_filenames(nout, ierr)
    elif command == 'b':
        nout, ierr = setup_blocknumbers(nin, nout, ierr)
    elif command == 'p':
        #
        # perform the plot. The plot type is determined from the value of
        # symbol_value(icountp)
        #
        plot_pars.plottype = read_internal(symbols.symbol_value[plot_pars.icountp]).int(plot_pars.plottype)
        ifilecount = plot_pars.plottype
        icount = symbols.nsymboluptodata + 1
        iref = file_and_directory_commands.ireff(icount, ifilecount)
        symbols.symbol_save = symbols.symbol_value[iref]
        nout, ierr = gnuplot(nin, nout, ierr)
    else:
        #
        # if none of the above we set the plottype to ifilecount
        #
        if ifilecount > 0 and ifilecount <= symbols.nfiles:
            plot_pars.plottype = ifilecount
        else:
            plot_pars.plottype = 1
        symbols.symbol_value[plot_pars.icountp] = fstr(lw(plot_pars.plottype), 256)
    return nout, icount, ifilecount, ierr


@routine
def setup_plot_filenames(nout, ierr):
    """plot_routines.f:setup_plot_filenames -- sets up the filenames for plotting and
    determines the number of available columns.

    Returns (nout, ierr).
    """
    from . import file_and_directory_commands

    # NB: iaccept (character*1) is not initialised in the Fortran; search_for_files sets it
    iaccept = ''
    #
    # set the number of plots to zero - this will be used to signal
    # whether valid plotting files have been defined
    #
    plot_pars.nplot = 0
    #
    # get the active plot type
    #
    (plot_pars.plottype,), ierr = _read_list(symbols.symbol_value[plot_pars.icountp], 'i', (plot_pars.plottype,))
    if ierr != 0 or plot_pars.plottype < 1 or plot_pars.plottype > symbols.nfiles:
        units.write(nout, fmt(_SPF102, 'setup_plot_filenames>', 1, symbols.nfiles))
        ierr = 1
        return nout, ierr
    #
    # get the extensions associated with this plottype
    #
    itest = symbols.nsymboluptodata + 2
    iref = file_and_directory_commands.ireff(itest, plot_pars.plottype)
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    _parse(5)
    if commands.nwords > 0:
        plot_pars.ext1 = fstr(_word(1), 20)
        plot_pars.ilenext1 = int(commands.ncl[1] - commands.ncf[1] + 1)
        plot_pars.nplot = 1
        if commands.nwords > 1:
            plot_pars.ext2 = fstr(_word(2), 20)
            plot_pars.ilenext2 = int(commands.ncl[2] - commands.ncf[2] + 1)
            plot_pars.nplot = 2
    else:
        units.write(nout, fmt(_SPF103, 'setup_plot_filenames>', plot_pars.plottype))
        ierr = 1
        plot_pars.nplot = 0
        return nout, ierr
    #
    # get the values of bspace and boffset - spacing between blocks in file
    # and offset for first block
    #
    iref = iref + 1
    (plot_pars.bspace,), ierr1 = _read_list(symbols.symbol_value[iref], 'i', (plot_pars.bspace,))
    iref = iref + 1

    # With this version we introduce the possibility of having different offsets
    # for each dataset when two are plotted together. We have to ensure this does
    # not crash older versions of plot_defaults.txt which will only have one entry
    # for boffset

    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    _parse(5)
    if commands.nwords > 1:
        (plot_pars.boffset, plot_pars.boffset2), ierr2 = _read_list(
            symbols.symbol_value[iref], 'ii', (plot_pars.boffset, plot_pars.boffset2))
    else:
        (plot_pars.boffset,), ierr2 = _read_list(symbols.symbol_value[iref], 'i', (plot_pars.boffset,))
        plot_pars.boffset2 = plot_pars.boffset
    if ierr1 != 0 or ierr2 != 0:
        units.write(nout, fmt(_SPF104, 'setup_plot_filenames>'))
        plot_pars.nplot = 0
        ierr = 1
        return nout, ierr
    #
    # identify a file to be plotted, ignoring any characters after the first
    # dot.
    #
    commands.sometext = fstr(symbols.symbol_value[plot_pars.icountf], 3020)
    ilenuptodot = findex(commands.sometext, '.') - 1
    if ilenuptodot > 0:
        commands.sometext = fstr(substr(commands.sometext, 1, ilenuptodot), 3020)
    _parse(5)
    if (commands.nwords < 1 or symbols.symbol_value[plot_pars.icountf] == 'search'
            or symbols.symbol_value[plot_pars.icountf] == '<undefined>'):
        iaccept = file_and_directory_commands.search_for_files(plot_pars.ext1, plot_pars.ilenext1, iaccept)
        #
        # remove any leading or trailing spaces in the filename
        #
        _parse(5)
    else:
        #
        # Add the correct extension back in.
        #
        commands.sometext = fstr(_words(1, commands.nwords) + pad(plot_pars.ext1, 20), 3020)
        _parse(5)
    #
    # regenerate the filename into a single word
    #
    _, commands.nwords, plot_pars.fnameplot, plot_pars.ilenfnameplot = parse.concatenate(
        1, commands.nwords, plot_pars.fnameplot, plot_pars.ilenfnameplot,
        commands.nwords, commands.ncf, commands.ncl, commands.sometext)
    #
    # build the full path plot filename(s)
    #
    plot_pars.fname1, plot_pars.ilenfname1 = file_and_directory_commands.get_full_filename(
        plot_pars.fnameplot, plot_pars.ilenfnameplot, plot_pars.fname1, plot_pars.ilenfname1)

    #
    # strip out anything beyond the extension that has been supplied
    #
    # (index with the blank-padded character*20 ext1, as in the Fortran)
    ilenuptoext = findex_padded(plot_pars.fnameplot, pad(plot_pars.ext1, 20), 256) - 1
    #
    # open the requested files and determine the number of columns
    #
    find_ncolumn(plot_pars.fname1)
    if plot_pars.nblock == 0:
        units.write(nout, fmt(_SPF125))
        ierr = 1
        plot_pars.nplot = 0
        return nout, ierr
    else:
        units.write(nout, fmt(_SPF126, plot_pars.nblock, substr(plot_pars.fname1, 1, plot_pars.ilenfname1)))
    ierr = 0
    #
    # save the good value of fnameplot
    #
    symbols.symbol_value[plot_pars.icountf] = fstr(substr(plot_pars.fnameplot, 1, ilenuptoext), 256)
    #
    # generate extension for corresponding data file
    #
    if plot_pars.nplot > 1:
        #
        # add on the new extension
        #
        plot_pars.fnameplot = fstr(substr(plot_pars.fnameplot, 1, ilenuptoext)
                                   + substr(plot_pars.ext2, 1, plot_pars.ilenext2), 256)
        plot_pars.ilenfnameplot = ilenuptoext + plot_pars.ilenext2
        #
        # create the full file name
        #
        plot_pars.fname2, plot_pars.ilenfname2 = file_and_directory_commands.get_full_filename(
            plot_pars.fnameplot, plot_pars.ilenfnameplot, plot_pars.fname2, plot_pars.ilenfname2)
        units.write(nout, fmt(_SPF127, 'setup_plot_filenames> ', substr(plot_pars.fname2, 1, plot_pars.ilenfname2)))
    return nout, ierr


@routine
def find_ncolumn(fname):
    """plot_routines.f:find_ncolumn -- number of blocks (and their labels) in the plot file.

    this routine applies strictly to EPSR generated files and is used by the plotting
    command. A more general purpose routine is get_nterms.
    """
    from . import file_and_directory_commands

    #
    # get the value of xcolumn
    #
    text = fstr('xcolumn', 20)
    # NB: ifound is not initialised in the Fortran; search_for_symbol only uses an
    # in-range value that already names the symbol, so any start value gives the same result
    ifound = 0
    ifound = file_and_directory_commands.search_for_symbol(text, ifound)
    iref = file_and_directory_commands.ireff(ifound, plot_pars.plottype)
    plot_pars.xcolumn = read_internal(symbols.symbol_value[iref]).int(plot_pars.xcolumn)

    #     write(6,*) fname
    ilenfname = findex_padded(fname, ' ', 256) - 1
    ifilefind = 0  # noqa: F841
    ierr = _open(10, fname, 'old')
    if ierr != 0:
        units.write(6, lw())
        units.write(6, lw('Problems with file ', substr(fname, 1, ilenfname)))
        plot_pars.nblock = 0
        return
    # read a line of text
    commands.sometext, ierr = _read_line(10, commands.sometext, 3020)
    if ierr == 0:
        # find out how many words it contains
        _parse(30)
        ntitle = 0
        #
        # if the first word has a # in it then the following words are
        # the x axis title, and n column titles
        #
        ilent = 0
        ilent = findex(_word(1), '#')
        if ilent > 0:
            #
            # A label, if it exists at all, will appear at the head of each column
            # or pair of columns
            #
            #
            # ignore the first word which contains the # and second word which contains the
            # x-axis label (if xcolumn.gt.0)
            #
            if plot_pars.xcolumn > 0:
                boff = 3
            else:
                boff = 3
            iw = boff - 1
            while iw < commands.nwords:
                iw = iw + 1
                ntitle = ntitle + 1
                plot_pars.plabel[ntitle] = fstr(_word(iw), 29)
            #
            # now read another line of text and use this to determine the number of blocks
            #
            commands.sometext, ierr = _read_line(10, commands.sometext, 3020)
            _parse(30)
        #
        # set up the block numbers
        #
        iw = plot_pars.boffset
        plot_pars.nblock = 0
        while iw <= commands.nwords:
            plot_pars.nblock = plot_pars.nblock + 1
            #
            # if ntitle is still zero the block numbers are used to invent suitable titles
            #
            if ntitle < plot_pars.nblock:
                plot_pars.plabel[plot_pars.nblock] = fstr(lw('Block', plot_pars.nblock), 29)
                ntitle = plot_pars.nblock
            iw = iw + plot_pars.bspace
        #
        # the ratio ntitle/nblock gives the number of titles per block. If this ratio
        # is greater than 1, it means there are more titles than we need, so we need
        # to pick out the appropriate ones for the present value of boffset
        #
        units.write(6, fmt(_FNC125, ntitle, plot_pars.nblock))
        if plot_pars.nblock > 0:
            bdiv = ntitle // plot_pars.nblock
            if bdiv > 1:
                boff = plot_pars.boffset - plot_pars.xcolumn
                it = boff
                ib = 1
                while it <= ntitle:
                    plot_pars.plabel[ib] = plot_pars.plabel[it]
                    it = it + bdiv
                    ib = ib + 1
    else:
        plot_pars.nblock = 0
    units.close(10)
    units.write(6, fmt(_FNC126, plot_pars.plottype, commands.nwords, plot_pars.boffset, plot_pars.boffset2,
                       plot_pars.bspace, plot_pars.nblock))
    return


@routine
def setup_blocknumbers(nin, nout, ierr):
    """plot_routines.f:setup_blocknumbers -- the block numbers to plot (icolumn, npcolumn).

    Returns (nout, ierr).
    """
    # NB: nfirst and nlast are not initialised in the Fortran; they are only used
    # after a read (an unsuccessful read leaves 0 here)
    nfirst = 0
    nlast = 0
    #
    # first check that a valid filename has been defined - any error and the
    # routine exits
    #
    nout, ierr = setup_plot_filenames(nout, ierr)
    if ierr != 0:
        units.write(nout, fmt(_SBN103))
        return nout, ierr
    #
    # set up the plot columns
    #
    #     write(6,*) icountb,' ',symbol_value(icountb)(1:50)
    commands.sometext = fstr(symbols.symbol_value[plot_pars.icountb], 3020)
    _parse(5)
    if commands.nwords <= 0:
        units.write(nout, fmt(_SBN100, plot_pars.nblock), advance=False)
        commands.sometext = fstr(units.read_line(nin), 3020)
        _parse(5)
        while commands.nwords == 0:
            units.write(nout, fmt(_SBN102, _prompt()), advance=False)
            commands.sometext = fstr(units.read_line(nin), 3020)
            _parse(5)
    #
    # From the input list determine the block numbers to be plotted.
    # The first word should always be a number - therefore
    # read the words until a number is found
    #
    iw = 1
    plot_pars.npcolumn = 0
    (nfirst,), ierr = _read_list(_word(iw), 'i', (nfirst,))
    iw = iw + 1
    #
    # now step through the remaining words and find out the other column numbers
    #
    while ierr == 0 and iw <= commands.nwords:
        if _word(iw) == '-':
            #
            # a range has been specified
            #
            iw = iw + 1
            (nlast,), ierr = _read_list(_word(iw), 'i', (nlast,))
            if nlast > nfirst:
                idel = 1
            elif nlast < nfirst:
                idel = -1
            else:
                idel = 0
            while nfirst != nlast:
                if nfirst > 0 and nfirst <= plot_pars.nblock:
                    plot_pars.npcolumn = plot_pars.npcolumn + 1
                    plot_pars.icolumn[plot_pars.npcolumn] = nfirst
                nfirst = nfirst + idel
        else:
            if nfirst > 0 and nfirst <= plot_pars.nblock:
                plot_pars.npcolumn = plot_pars.npcolumn + 1
                plot_pars.icolumn[plot_pars.npcolumn] = nfirst
            (nfirst,), ierr = _read_list(_word(iw), 'i', (nfirst,))
        iw = iw + 1
    #
    # save the last value
    #
    if nfirst > 0 and nfirst <= plot_pars.nblock:
        plot_pars.npcolumn = plot_pars.npcolumn + 1
        plot_pars.icolumn[plot_pars.npcolumn] = nfirst

    if plot_pars.npcolumn > 0:

        units.write(6, fmt(_SBN235, plot_pars.npcolumn))

    else:

        ierr = 1
        units.write(6, fmt(_SBN236))

    #
    # if the command is successful, save it
    #
    if ierr == 0:
        _, commands.nwords, symbols.symbol_value[plot_pars.icountb], ilen = parse.concatenate(
            1, commands.nwords, symbols.symbol_value[plot_pars.icountb], 0,
            commands.nwords, commands.ncf, commands.ncl, commands.sometext)

    return nout, ierr


@routine
def gnuplot(nin, nout, ierr):
    """plot_routines.f:gnuplot -- sets up a GNUplot input file and runs GNU plot.

    gnuplot is started with ``call system(... ' -')`` in the foreground
    (``runtime.system``), exactly as the Fortran does.
    Returns (nout, ierr).
    """
    from . import file_and_directory_commands

    #
    # now check that block numbers have been defined - any error and the
    # routine exits
    #
    nout, ierr = setup_blocknumbers(nin, nout, ierr)
    if ierr != 0:
        units.write(nout, fmt(_GP122))
        return nout, ierr
    #
    # get the other information needed to form the plot
    #
    #
    # get the extensions associated with this plottype
    #
    itest = symbols.nsymboluptodata + 5
    iref = file_and_directory_commands.ireff(itest, plot_pars.plottype)
    plot_pars.xmin = read_internal(symbols.symbol_value[iref]).real(plot_pars.xmin)
    iref = iref + 1
    plot_pars.xmax = read_internal(symbols.symbol_value[iref]).real(plot_pars.xmax)
    iref = iref + 1
    plot_pars.ymin = read_internal(symbols.symbol_value[iref]).real(plot_pars.ymin)
    iref = iref + 1
    plot_pars.ymax = read_internal(symbols.symbol_value[iref]).real(plot_pars.ymax)
    iref = iref + 1
    plot_pars.ydel = read_internal(symbols.symbol_value[iref]).real(plot_pars.ydel)

    # factor for y values

    iref = iref + 1
    plot_pars.yfact1 = f32(1.0)
    plot_pars.yfact2 = f32(1.0)
    (plot_pars.yfact1, plot_pars.yfact2), ierr = _read_list(
        symbols.symbol_value[iref], 'rr', (plot_pars.yfact1, plot_pars.yfact2))
    units.write(6, lw(substr(symbols.symbol_value[iref], 1, 20), plot_pars.yfact1, plot_pars.yfact2))
    #
    # get fractional positions of labels
    #
    iref = iref + 1
    plot_pars.xfrac = read_internal(symbols.symbol_value[iref]).real(plot_pars.xfrac)
    iref = iref + 1
    plot_pars.yfrac = read_internal(symbols.symbol_value[iref]).real(plot_pars.yfrac)
    #
    # plot title
    #
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    _parse(5)
    if commands.nwords > 0:
        _, commands.nwords, plot_pars.title, plot_pars.ilentitle = parse.concatenate(
            1, commands.nwords, plot_pars.title, plot_pars.ilentitle,
            commands.nwords, commands.ncf, commands.ncl, commands.sometext)
    else:
        plot_pars.title = ''
        plot_pars.ilentitle = 1
    #
    # set up a default title in case it is undefined
    #
    if substr(plot_pars.title, 1, plot_pars.ilentitle) == '<undefined>':
        ilenuptodot = findex(plot_pars.fnameplot, '.') - 1
        plot_pars.title = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                               + substr(plot_pars.fnameplot, 1, ilenuptodot), 256)
        plot_pars.ilentitle = directory_name.nchardir + ilenuptodot
    #
    # x-axis label
    #
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    _parse(5)
    if commands.nwords > 0:
        _, commands.nwords, plot_pars.xlabel, plot_pars.ilenxlabel = parse.concatenate(
            1, commands.nwords, plot_pars.xlabel, plot_pars.ilenxlabel,
            commands.nwords, commands.ncf, commands.ncl, commands.sometext)
    else:
        plot_pars.xlabel = ''
        plot_pars.ilenxlabel = 1
    #
    # y-axis label
    #
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    _parse(5)
    if commands.nwords > 0:
        _, commands.nwords, plot_pars.ylabel, plot_pars.ilenylabel = parse.concatenate(
            1, commands.nwords, plot_pars.ylabel, plot_pars.ilenylabel,
            commands.nwords, commands.ncf, commands.ncl, commands.sometext)
    else:
        plot_pars.ylabel = ''
        plot_pars.ilenylabel = 1
    #
    # log x-axis
    #
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    _parse(5)
    if commands.nwords > 0:
        _, _, plot_pars.logx, plot_pars.ilenlogx = parse.concatenate(
            1, 1, plot_pars.logx, plot_pars.ilenlogx,
            commands.nwords, commands.ncf, commands.ncl, commands.sometext)
    else:
        plot_pars.logx = 'no'
        plot_pars.ilenlogx = 2
    #
    # if this axis is log x, check that the minimum x value is > 0
    #
    if substr(plot_pars.logx, 1, 1) == 'y' and plot_pars.xmin <= 0.0:
        plot_pars.xmin = r4(1.0e-12)

    #
    # log y-axis
    #
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    _parse(5)
    if commands.nwords > 0:
        _, _, plot_pars.logy, plot_pars.ilenlogy = parse.concatenate(
            1, 1, plot_pars.logy, plot_pars.ilenlogy,
            commands.nwords, commands.ncf, commands.ncl, commands.sometext)
    else:
        plot_pars.logy = 'no'
        plot_pars.ilenlogy = 2
    #
    # if this axis is log y, check that the minimum y value is > 0
    #
    if substr(plot_pars.logy, 1, 1) == 'y' and plot_pars.ymin <= 0.0:
        plot_pars.ymin = r4(1.0e-12)
    units.write(6, lw(substr(plot_pars.logx, 1, 1), substr(plot_pars.logy, 1, 1)))
    #
    # column number for x values. 0 means use the item number as x value
    #
    iref = iref + 1
    plot_pars.xcolumn = read_internal(symbols.symbol_value[iref]).int(plot_pars.xcolumn)
    #
    # decide if this is a histogram or not
    #
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    _parse(5)
    if commands.nwords > 0:
        _, _, plot_pars.histogram, plot_pars.ilenhistogram = parse.concatenate(
            1, 1, plot_pars.histogram, plot_pars.ilenhistogram,
            commands.nwords, commands.ncf, commands.ncl, commands.sometext, len_return_string=20)
    else:
        plot_pars.histogram = 'no'
        plot_pars.ilenhistogram = 2
    #
    # set the dimensions of plot
    #
    plot_pars.yup = plot_pars.ymax + f32(plot_pars.npcolumn - 1) * plot_pars.ydel
    #
    # generate the type label
    #
    typelabel = fstr(fmt(_GP150, plot_pars.plottype), 2)
    #
    # generate the name of the output file
    #
    #
    # first strip out the extension in fnameplot
    #
    #
    # strip out anything beyond the . in fnameplot
    #
    ilenuptodot = findex(plot_pars.fnameplot, '.') - 1
    plot_output = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir) + 'plot'
                       + substr(typelabel, 1, 2) + substr(plot_pars.fnameplot, 1, ilenuptodot) + '.gnu', 256)
    ilenplot_output = directory_name.nchardir + ilenuptodot + 10
    ierr = _open(11, plot_output, 'unknown')
    #
    # write the title and x and y labels
    #
    units.write(11, fmt(_GP101, substr(plot_pars.title, 1, plot_pars.ilentitle)))
    units.write(11, fmt(_GP102, substr(plot_pars.xlabel, 1, plot_pars.ilenxlabel)))
    units.write(11, fmt(_GP103, substr(plot_pars.ylabel, 1, plot_pars.ilenylabel)))
    units.write(11, fmt(_GP104, plot_pars.xmin, plot_pars.xmax))
    units.write(11, fmt(_GP105, plot_pars.ymin, plot_pars.yup))
    #
    # log or linear axes
    #
    if substr(plot_pars.logx, 1, 1) == 'y':
        units.write(11, fmt(_GP113, 'x'))
    else:
        units.write(11, fmt(_GP114, 'x'))
    if substr(plot_pars.logy, 1, 1) == 'y':
        units.write(11, fmt(_GP113, 'y'))
    else:
        units.write(11, fmt(_GP114, 'y'))

    #
    # set the labels and linestyles
    #
    # x coordinate of data label
    #
    # NB: iplot is written (write(6,*) below) before it is first set when npcolumn = 0,
    # which setup_blocknumbers excludes; it then holds icolumn(npcolumn) from this loop.
    iplot = 0
    labelx = plot_pars.xmin + plot_pars.xfrac * (plot_pars.xmax - plot_pars.xmin)
    for ip in range(1, plot_pars.npcolumn + 1):
        #
        # y coordinate of data label
        #
        labely = plot_pars.yfrac + f32(ip - 1) * plot_pars.ydel
        iplot = int(plot_pars.icolumn[ip])
        ilent = findex_padded(plot_pars.plabel[iplot], ' ', 29) - 1  # noqa: F841
        units.write(11, fmt(_GP106, ip, substr(plot_pars.plabel[iplot], 1, 29), labelx, labely))
        units.write(11, fmt(_GP107, ip, ip, system_commands.gnu_pointtype, system_commands.gnu_pointsize))

    # decide whether a histogram is needed

    if substr(plot_pars.histogram, 1, 1) == 'y':
        plot_pars.style = 'histeps'
        plot_pars.ilenstyle = 7
    else:
        plot_pars.style = 'lines'
        plot_pars.ilenstyle = 5
    #
    # now write the plot command
    #
    units.write(11, fmt(_GP108))
    units.write(6, lw(plot_pars.boffset, plot_pars.boffset2))
    for ip in range(1, plot_pars.npcolumn + 1):
        yadd = f32(ip - 1) * plot_pars.ydel
        units.write(6, lw(ip, plot_pars.boffset, plot_pars.bspace, plot_pars.icolumn[ip], iplot))
        iplot = plot_pars.boffset + plot_pars.bspace * (int(plot_pars.icolumn[ip]) - 1)
        units.write(11, fmt(_GP109, substr(plot_pars.fname1, 1, plot_pars.ilenfname1), plot_pars.xcolumn,
                            plot_pars.yfact1, iplot, yadd, substr(plot_pars.style, 1, plot_pars.ilenstyle), ip),
                    advance=False)
        #
        # add the continuation line if necessary
        #
        if ip != plot_pars.npcolumn or plot_pars.nplot != 1:
            units.write(11, fmt(_GP110))
        if plot_pars.nplot == 2:
            iplot = plot_pars.boffset2 + plot_pars.bspace * (int(plot_pars.icolumn[ip]) - 1)
            units.write(11, fmt(_GP111, substr(plot_pars.fname2, 1, plot_pars.ilenfname2), plot_pars.xcolumn,
                                iplot, plot_pars.yfact2, iplot, yadd, ip), advance=False)
            if ip != plot_pars.npcolumn:
                units.write(11, fmt(_GP110))
    #
    # finish off the last line
    #
    units.write(11, lw())
    #     write(11,*) 'pause -1'
    units.close(11)
    plot_command = fstr(substr(system_commands.system_gnu, 1, system_commands.ilensystem_gnu)
                        + ' '
                        + substr(system_commands.system_quotes, 1, system_commands.ilensystem_quotes)
                        + substr(plot_output, 1, ilenplot_output)
                        + substr(system_commands.system_quotes, 1, system_commands.ilensystem_quotes)
                        + ' -', 256)
    ierr = runtime.system(plot_command)
    return nout, ierr


@routine
def read_plot_input_file(nin, nout, ierrout):
    """plot_routines.f:read_plot_input_file -- reads the PLOT input file.

    Returns (nout, ierrout).
    """
    #
    # don't do anything if there is an error in the input file - will start out with
    # default values in this case
    #
    if ierrout != 0:
        units.close(nin)
        return nout, ierrout
    #
    # Read in commands from the already opened input file, and process each command
    # in turn. At the end close this file.
    #
    command = 'next'
    icount = 1
    ifilecount = 1
    #
    # we assume that the symbols are specified in order, but the
    # process_commands routine will check each command in turn to ensure the correct
    # values are updated. The last one should be a 'q' to signify the end of
    # command input
    #
    ierr = 0
    while command != 'e' and ierr == 0:
        #
        # check the status of the file counter
        #
        if icount > symbols.nsymboldata and ifilecount < symbols.nfiles:
            ifilecount = ifilecount + 1
            icount = symbols.nsymboluptodata + 1
        elif icount > symbols.nsymbol:
            icount = 1
        #
        # read the commands
        #
        commands.nwords = 0
        # (on end of file sometext keeps the previous line, which is processed again)
        commands.sometext, ierr = _read_line(nin, commands.sometext, 3020)
        _parse(12)
        #     write(6,*) sometext(ncf(1):ncl(nwords))

        # Need to ignore blank lines read from an input file

        if commands.nwords > 0:
            nout, command, icount, ifilecount, ierrout = process_plot_commands(
                nin, nout, command, icount, ifilecount, ierrout, 1)

    #
    # finished input - close the input file
    #
    units.close(nin)
    return nout, ierrout


@routine
def process_plot_symbols(nin, nout, ierrout, ignore):
    """plot_routines.f:process_plot_symbols -- processes commands coming in on the input
    stream, nin. Diagnostics are written to the output stream nout. An 'e' on the input
    stream signifies the end of input.

    Returns (nout, ierrout).
    """
    from . import file_and_directory_commands

    #
    # set an error flag so that calling routine can tell if an error has occurred
    #
    ierrout = 0
    command = 'next'
    icount = 0
    #
    # set up a counter to say which dataset is being looked at
    #
    ifilecount = 1
    #
    # signal whether current values have been saved or not
    #
    symbols.isave = 0
    #
    # setup an integer for the number of not recognised commands
    #
    symbols.inotrecog = 0
    ierr = 0
    while command != 'e' and ierr == 0:
        if icount > 0 and icount <= symbols.nsymbol:
            #
            # get the reference value for the symbol
            #
            iref = file_and_directory_commands.ireff(icount, ifilecount)
            #     write(nout,*) icount,iref,ifilecount,nfiles
            #
            # sort out the words of symbol_value
            #
            commands.sometext = fstr(symbols.symbol_value[iref], 3020)
            _parse(5)
            ilensymbol = findex_padded(symbols.symbol[icount], ' ', 11) - 1
            st = symbols.symbol_text[icount]
            units.write(nout, fmt(_PPS203, _prompt(), substr(symbols.symbol[icount], 1, ilensymbol),
                                  substr(st, 1, max(1, findex_padded(st, '  ', 256)))))
            units.write(nout, fmt(_PPS202, substr(symbols.symbol[icount], 1, ilensymbol),
                                  _words(1, max(1, commands.nwords))), advance=False)
        else:
            #
            # reset the file counter if you have got here
            #
            #           ifilecount=1
            units.write(nout, fmt(_PPS201, _prompt()), advance=False)

        commands.sometext = ''
        commands.sometext, ierr = _read_line(nin, commands.sometext, 3020)
        _parse(11)
        nout, command, icount, ifilecount, ierrout = process_plot_commands(
            nin, nout, command, icount, ifilecount, ierrout, ignore)
        #        write(6,*) command,isave

    if symbols.isave != 1:
        #
        # plot files are always saved without question
        #
        if symbols.symbol_setup != 'plot':
            commands.nwords = 0
            commands.sometext = ''
            units.write(nout, fmt(_PPS401, _prompt()), advance=False)
            commands.sometext, ierr = _read_line(nin, commands.sometext, 3020)
            _parse(5)
            if commands.nwords > 0:
                command = fstr(_word(1), 20)
            else:
                command = ''
            if command != 'q':
                nout = save_plot_files(nin, nout)
        else:
            nout = save_plot_files(nin, nout)
    return nout, ierrout


@routine
def process_plot_commands(nin, nout, command, icount, ifilecount, ierrout, ignore):
    """plot_routines.f:process_plot_commands -- processes command in the input stream, nin.

    Diagnostics are written to the output stream nout. An 'e' on the input stream
    signifies the end of input.  'data' <number> signifies the following entries relate
    to dataset <number>. If no number is specified, then dataset 1 is assumed.
    Returns (nout, command, icount, ifilecount, ierrout).
    """
    from . import file_and_directory_commands

    # NB: ierr is a local that is not initialised in the Fortran; it is only tested
    # after check_plot_symbol has set it
    ierr = 0
    ndata = 0
    #
    # set an error flag so that calling routine can tell if an error has occurred
    #
    ierrout = 0
    #
    # if nwords .eq. 0 it means the previous change (if any) has been accepted.
    # if however symbol_save is already blank, it means this was a blank line
    # in which case the symbol counter should not be incremented
    #
    if commands.nwords > 0:
        command = fstr(_word(1), 20)
    else:
        command = 'd'
    if command == 'd' or command == 'u':
        if symbols.symbol_save != '':
            symbols.symbol_save = ''
        ifilecountold = ifilecount

        if command == 'd':

            icount = icount + 1
            #
            # check the status of the file counter
            #
            if icount > symbols.nsymboldata and ifilecount < symbols.nfiles:
                ifilecount = ifilecount + 1
                icount = symbols.nsymboluptodata + 1
            elif icount > symbols.nsymbol:
                icount = 0
                ifilecount = 1
            #
            # signal that you want to write out the file counter
            #
            if icount == symbols.nsymboluptodata + 1:
                ifilecountold = 0
            else:
                ifilecountold = ifilecount
            #     write(6,*) icount,ifilecount,ifilecountold,nsymboluptodata,nsymboldata,nsymbol
        elif command == 'u':

            icount = icount - 1
            #
            # check the status of the file counter
            #
            if icount < 0:
                icount = symbols.nsymbol
                ifilecount = symbols.nfiles
            if icount <= symbols.nsymboluptodata and ifilecount > 1:
                icount = symbols.nsymboldata
                ifilecount = ifilecount - 1
            #
            # signal that you want to write out the file counter
            #
            if icount == symbols.nsymboldata:
                ifilecountold = 0
            else:
                ifilecountold = ifilecount
        #
        # set up the data text if ifilecount has changed
        #
        if ifilecount != ifilecountold:
            commands.sometext = fstr(fmt(_PPC250, ifilecount), 3020)
            _parse(5)
            if commands.nwords > 0:
                symbols.datatext = fstr(_word(1), 80)
                symbols.ilendatatext = int(commands.ncl[1] - commands.ncf[1] + 1)
            else:
                symbols.datatext = ''
                symbols.ilendatatext = 1
            #
            # and write it!
            #
            if nin == 5:
                units.write(nout, fmt(_PPC204, substr(symbols.dataprompt, 1, symbols.ilendataprompt),
                                      substr(symbols.datatext, 1, symbols.ilendatatext),
                                      substr(symbols.dataprompttext, 1,
                                             max(1, findex_padded(symbols.dataprompttext, '  ', 80)))))
        command = 'next'
        return nout, command, icount, ifilecount, ierrout
    #
    # perform any necessary translations
    #
    if command == 'bye':
        command = 'e'
    #
    # reflect the input text into the output stream
    #
    #     write(nout,*) sometext(ncf(1):ncl(nwords))
    #
    # first check whether the value of command is one of the special values
    #
    if command == 'default':
        symbols.inotrecog = 0
        is_ = 0
        ifilecount = 1
        while is_ < symbols.nsymboluptodata:
            is_ = is_ + 1
            symbols.symbol_value[is_] = symbols.symbol_default[is_]
        while ifilecount < symbols.nfiles:
            is_ = is_ + 1
            iref = is_ + symbols.nsymbolplus * (ifilecount - 1)
            symbols.symbol_value[iref] = symbols.symbol_default[is_]
            if is_ >= symbols.nsymboldata and ifilecount < symbols.nfiles:
                is_ = symbols.nsymboluptodata
                ifilecount = ifilecount + 1
        command = 'next'
    elif command == 'ls' or command == 'dir':
        symbols.inotrecog = 0
        file_and_directory_commands.list_directory()
        command = 'next'
    elif command == 'cd' or command == 'pwd':
        symbols.inotrecog = 0
        file_and_directory_commands.set_directory_name()
        command = 'next'
    elif command == symbols.dataprompt:
        symbols.inotrecog = 0
        #
        # try to decide which data set you wish to access
        #
        ifound = 0
        if commands.nwords > 1:
            (ndata,), ierr = _read_list(_word(2), 'i', (ndata,))
            if ierr == 0:
                ifound = 1
        if ifound == 0:
            ndata = 1
        if ndata < 1:
            ndata = 1
        if ndata > symbols.nfiles:
            ndata = symbols.nfiles
        ifilecount = ndata
        icount = symbols.nsymboluptodata + 1
        commands.sometext = fstr(fmt(_PPC250, ifilecount), 3020)
        _parse(5)
        if commands.nwords > 0:
            symbols.datatext = fstr(_word(1), 80)
            symbols.ilendatatext = int(commands.ncl[1] - commands.ncf[1] + 1)
        else:
            symbols.datatext = ''
            symbols.ilendatatext = 1
        #
        # and write it!
        #
        if nin == 5:
            units.write(nout, fmt(_PPC204, substr(symbols.dataprompt, 1, symbols.ilendataprompt),
                                  substr(symbols.datatext, 1, symbols.ilendatatext),
                                  substr(symbols.dataprompttext, 1,
                                         max(1, findex_padded(symbols.dataprompttext, '  ', 80)))))
        #
        # for plotting we need to process this command
        #
        if symbols.symbol_setup == 'plot' and ignore == 0:
            iref = file_and_directory_commands.ireff(icount, ifilecount)
            # allow this command through the nochanges filter in check_symbol
            symbols.symbol_save = symbols.symbol_value[iref]
            command = 'store'
        else:
            command = 'next'
    elif command == 'help' or command == '?':
        symbols.inotrecog = 0
        #
        # write out the list of commands
        #
        units.write(nout, fmt(_PPC207))
        for is_ in range(1, symbols.nsymbol + 1):
            ilensymbol = findex_padded(symbols.symbol[is_], ' ', 11) - 1  # noqa: F841
            units.write(nout, fmt(_PPC209, Ch(symbols.symbol[is_], 11), substr(symbols.symbol_text[is_], 1, 65)))
        units.write(nout, fmt(_PPC210))
        if symbols.ilendataprompt > 0:
            text = fstr(substr(symbols.dataprompt, 1, symbols.ilendataprompt) + ' <number> ', 20)
            units.write(nout, fmt(_PPC211, substr(text, 1, 20), Ch(symbols.dataprompttext, 80)))
        units.write(nout, fmt(_PPC208))
        command = 'next'
    elif command == 'save':
        symbols.inotrecog = 0
        nout = save_plot_files(nin, nout)
        symbols.isave = 1
        command = 'next'
    elif command == 'q':
        symbols.inotrecog = 0
        symbols.isave = 1
        command = 'e'
    elif command != 'e' and command != 'next' and command != 'n':
        #
        # have not recognised this command so far
        #
        #
        # first check if a valid symbol has been called
        #
        # (we add 1 to the value of icount since that is the most likely
        # value if commands are being read from a file - this is to save search
        # time)
        #
        ifound = icount + 1
        ifound = file_and_directory_commands.search_for_symbol(command, ifound)
        if ifound > 0:
            #
            # set up the correct reference integer depending on the value of ifilecount
            #
            iref = file_and_directory_commands.ireff(ifound, ifilecount)
            #
            # save the current value of this symbol
            #
            symbols.symbol_save = symbols.symbol_value[iref]
            #
            # can only store the new value if it is defined! (nwords.gt.1)
            #
            # concatenate the text from the words, starting from the second, and place it
            # in the appropriate symbol_value
            #
            #     write(6,*) ifound,ifilecount,iref,command,' ',symbol_value(iref)(1:30)
            if commands.nwords > 1:
                _, commands.nwords, symbols.symbol_value[iref], ilenfound = parse.concatenate(
                    2, commands.nwords, symbols.symbol_value[iref], 0,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext)
                icount = ifound
                command = 'store'
            else:

                # otherwise go to this parameter and ask for the new value
                # if no new value was specified we go to this new value. Normally we do nothing
                # at that point, unless these are plot symbols, when action with the current value
                # may need to be processed

                icount = ifound
                if symbols.symbol_setup == 'plot':
                    command = 'store'
                else:
                    command = 'next'
            if ifound > symbols.nsymboldata:
                ifilecount = symbols.nfiles
        else:
            #
            # if icount is greater than zero
            # save the text in the appropriate symbol value
            # however if the input channel is not 5, it probably means the value is
            # spurious since values from a file would be preceded by a recognisable symbol
            #
            if icount > 0 and nin == 5:
                #
                # set up the correct reference integer depending on the value of ifilecount
                #
                iref = file_and_directory_commands.ireff(icount, ifilecount)
                #
                # save the current value of the symbol
                #
                symbols.symbol_save = symbols.symbol_value[iref]
                _, commands.nwords, symbols.symbol_value[iref], ilenfound = parse.concatenate(
                    1, commands.nwords, symbols.symbol_value[iref], 0,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext)
                command = 'store'
            else:
                units.write(nout, fmt(_PPC100, Ch(command, 20), ' ', substr(commands.sometext, 1, 40)))
                symbols.inotrecog = symbols.inotrecog + 1
                ierrout = 1
                command = 'next'
    #     write(6,*) command,isave
    #
    # if command is 'store' it means we need to check the new symbol_value
    # and set it back to the saved value if there is an error, or make any other changes
    # to the symbol values that are required. See check_symbol for details.
    #
    if command == 'store':

        #
        # signal that the new value will not have been written to the input file
        #
        symbols.inotrecog = 0
        symbols.isave = 0

        nout, icount, ifilecount, ierr = check_plot_symbol(nin, nout, icount, ifilecount, ierr, ignore)

        #
        # set the value of ifilecount if icount does not refer to a dataset
        #
        #           if(icount.le.nsymboluptodata) ifilecount=1
        #           if(icount.gt.nsymboldata) ifilecount=nfiles
        if ierr == 0:
            command = 'next'
        else:
            command = 'n'
        #     write(6,*) iref,ierr,symbol_value(iref)(1:30),command

    if command == 'n':

        #
        # change of value was not accepted
        #

        symbols.inotrecog = 0
        #
        # set up the correct reference integer depending on the value of ifilecount
        #
        if symbols.symbol_save != '' and icount > 0:
            #
            # the new value was not accepted so revert to the previous value
            #
            iref = file_and_directory_commands.ireff(icount, ifilecount)
            symbols.symbol_value[iref] = symbols.symbol_save

        else:
            #
            # previous saved value was accepted - nothing to change
            #
            icount = icount + 1
        command = 'next'

    return nout, command, icount, ifilecount, ierrout


@routine
def check_plot_symbol(nin, nout, icount, ifilecount, ierr, ignore):
    """plot_routines.f:check_plot_symbol -- check a new symbol value (type, permission,
    number of plot types) and process it (plots are made from here).

    Returns (nout, icount, ifilecount, ierr).
    """
    from . import file_and_directory_commands

    # NB: icheck, rcheck, nfilet, nvary, nvaryt and iaccept are not initialised in the
    # Fortran; each is set by a read (or search_for_files) before it is used
    icheck = 0
    rcheck = f32(0.0)
    nfilet = 0
    nvary = 0
    nvaryt = 0
    iaccept = ''
    #     write(6,*) symbol_value(iinitref)(1:10),iinitref

    #
    # make sure icount does not exceed the maximum allowed value
    #
    if icount < 1:
        icount = 1
    if icount > symbols.nsymbol:
        icount = symbols.nsymbol
    #     if(icount.gt.nsymboldata) ifilecount=nfiles
    iref = file_and_directory_commands.ireff(icount, ifilecount)
    #
    # set error flag to zero
    #
    ierr = 0
    #
    # check that it is allowed to change this value - if a change has been requested
    #
    if symbols.iflag[icount] > 0 and ignore == 0 and symbols.symbol_value[iref] != symbols.symbol_save:
        units.write(nout, fmt(_CPS200, _prompt()))
        symbols.symbol_value[iref] = symbols.symbol_save
        return nout, icount, ifilecount, ierr
    #
    # now check that the specified value(s) is(are) of the correct type if it is not
    # a character variable
    #
    if symbols.symbol_type[icount] != 'c':
        commands.sometext = fstr(symbols.symbol_value[iref], 3020)
        #     write(6,*) symbol_value(iref),sometext(1:50)
        _parse(11)
        #     write(6,*) symbol_value(iref),sometext(1:50),nwords,ncf(1),ncl(1)
        iw = 0
        if symbols.symbol_type[icount] == 'i':
            while iw < commands.nwords and ierr == 0:
                iw = iw + 1
                (icheck,), ierr = _read_list(_word(iw), 'i', (icheck,))
                if ierr != 0:
                    units.write(6, fmt(_CPS201, _prompt(), _word(iw)))
                    symbols.symbol_value[iref] = symbols.symbol_save
                    ierr = 0
                    return nout, icount, ifilecount, ierr
        else:
            #     write(6,*) nwords,iw,sometext(ncf(1):ncl(1))
            while iw < commands.nwords and ierr == 0:
                iw = iw + 1
                (rcheck,), ierr = _read_list(_word(iw), 'r', (rcheck,))
                #     write(6,*) ierr,rcheck
                if ierr != 0:
                    units.write(nout, fmt(_CPS202, _prompt(), _word(iw)))
                    symbols.symbol_value[iref] = symbols.symbol_save
                    ierr = 0
                    return nout, icount, ifilecount, ierr
    #
    # if icount .eq. nsymboluptodata the number of files may have been
    # changed. This means the values of symbol_value above nsymboluptodata will
    # need to be resequenced
    #
    if icount == symbols.nsymboluptodata:
        #
        # read the new value of nfiles into a temporary location
        #
        (nfilet,), ierr = _read_list(symbols.symbol_value[iref], 'i', (nfilet,))

        if ierr == 0:
            #
            # First reassign the symbol_values above nsymboldata
            #
            nfadd = (nfilet - symbols.nfiles) * symbols.nsymbolplus
            #
            # If new number of files is greater than old this must be done in reverse order
            # to avoid overwriting existing values
            #
            if nfadd > 0:
                isymbol = symbols.nsymbol + 1
                ivalue = file_and_directory_commands.ireff(symbols.nsymbol, symbols.nfiles) + 1
                while isymbol > symbols.nsymboldata + 1:
                    isymbol = isymbol - 1
                    ivalue = ivalue - 1
                    ivaluenew = ivalue + nfadd
                    #     write(6,*) nsymboluptodata,nsymboldata,nsymbol,isymbol,ivalue,ivaluenew
                    if ivaluenew <= symbols.msymbol:
                        symbols.symbol_value[ivaluenew] = symbols.symbol_value[ivalue]
                    else:
                        units.write(nout, fmt(_CPS203, _prompt()))
                        ierr = 1
                        return nout, icount, ifilecount, ierr
            elif nfadd < 0:
                isymbol = symbols.nsymboldata
                ivalue = file_and_directory_commands.ireff(symbols.nsymboldata, symbols.nfiles)
                while isymbol < symbols.nsymbol:
                    isymbol = isymbol + 1
                    ivalue = ivalue + 1
                    ivaluenew = ivalue + nfadd
                    #     write(6,*) nsymboluptodata,nsymboldata,nsymbol,isymbol,ivalue,ivaluenew
                    if ivaluenew <= symbols.msymbol:
                        symbols.symbol_value[ivaluenew] = symbols.symbol_value[ivalue]
                    else:
                        units.write(nout, fmt(_CPS203, _prompt()))
                        ierr = 1
                        return nout, icount, ifilecount, ierr
                #
                # no need to do anything if nfadd = 0
                #
            ifiletemp = symbols.nfiles
            symbols.nfiles = nfilet
            #
            # if new number of files .gt. old number of files
            # set the new values of symbol_value to the default values
            #
            if nfadd > 0:
                while ifiletemp < symbols.nfiles:
                    ifiletemp = ifiletemp + 1
                    isymbol = symbols.nsymboluptodata
                    while isymbol < symbols.nsymboldata:
                        isymbol = isymbol + 1
                        iref = file_and_directory_commands.ireff(isymbol, ifiletemp)
                        symbols.symbol_value[iref] = symbols.symbol_default[isymbol]

    #
    # plotting commands are processed by their own routine, but not during reading them in
    #
    if symbols.symbol_setup == 'plot' and ignore == 0:
        nout, icount, ifilecount, ierr = process_plots(nin, nout, icount, ifilecount, ierr)
        return nout, icount, ifilecount, ierr
    #
    # changes specific to 2D and 3D plotting
    #

    if symbols.symbol_setup == 'plot2d' or symbols.symbol_setup == 'plot3d':

        if symbols.symbol[icount] == 'shcoeffs':
            #
            # changing the name of the coefficients file - read the corresponding .dat file
            #
            text = fstr(symbols.symbol_value[file_and_directory_commands.ireff(symbols.nsymbol, symbols.nfiles)], 256)
            commands.nwords = 0
            commands.nwords = parse.parse(text, commands.nwords, commands.ncf, commands.ncl, 5, 256, commands.nworddim)
            ext = fstr(substr(text, commands.ncf[1], commands.ncl[1]), 20)
            ilenext = int(commands.ncl[1] - commands.ncf[1] + 1)
            commands.nwords = 0
            #
            # if symbol_value(1) is 'search' we need to find a valid SHARM coefficients file
            #
            if symbols.symbol_value[iref] == 'search' or symbols.symbol_value[iref] == '<undefined>':
                iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)
                _parse(5)
            else:
                commands.sometext = fstr(symbols.symbol_value[iref], 3020)
                _parse(5)
                #
                # add the extension if not supplied
                #
                ilendot = findex(_word(1), '.')
                if ilendot <= 0:
                    commands.sometext = fstr(_word(1) + pad(ext, 20), 3020)
                    commands.ncl[1] = commands.ncl[1] + ilenext
            if commands.nwords <= 0:
                units.write(nout, fmt(_CPS206, _prompt()))
                symbols.symbol_value[iref] = symbols.symbol_save
                ierr = 0
                return nout, icount, ifilecount, ierr
            else:
                symbols.symbol_value[iref] = fstr(_word(1), 256)
            nout, ierr = find_plot_par(12, nout, ierr)
        #
        # if setone has been changed, make sure new value is either 0 or 1
        #
        if symbols.symbol[icount] == 'setone':
            (icheck,), ierr = _read_list(symbols.symbol_value[iref], 'i', (icheck,))
            # NB: the test is always true (icheck cannot be both 0 and 1), as in the Fortran
            if ierr != 0 or icheck != 0 or icheck != 1:
                units.write(nout, fmt(_CPS207, _prompt()))
                symbols.symbol_value[iref] = symbols.symbol_save
                ierr = 0
                return nout, icount, ifilecount, ierr
        #
        # if nvary has changed (plot2d, plot3d) need to change symbol_text and symbol_value
        #
        if symbols.symbol[icount] == 'nvary':
            #
            # read the new value of nvary
            #
            iref = file_and_directory_commands.ireff(icount, ifilecount)
            (nvary,), ierr1 = _read_list(symbols.symbol_save, 'i', (nvary,))
            (nvaryt,), ierr2 = _read_list(symbols.symbol_value[iref], 'i', (nvaryt,))
            if ierr1 == 0 and ierr2 == 0 and nvaryt != nvary:
                iref = iref + 1
                if symbols.symbol_setup == 'plot3d':
                    if nvaryt == 1:
                        symbols.symbol_text[iref] = 'phim, thetam, chim ' + '(degrees - separated by spaces)'
                    elif nvaryt == 2:
                        symbols.symbol_text[iref] = 'thetal, phil, chim ' + '(degrees - separated by spaces)'
                    elif nvaryt == 3:
                        symbols.symbol_text[iref] = 'thetal, phil, phim ' + '(degrees - separated by spaces)'
                elif symbols.symbol_setup == 'plot2d':
                    if nvaryt == 1:
                        symbols.symbol_text[iref] = 'phil, phim, thetam, chim ' + '(degrees - separated by spaces)'
                    elif nvaryt == 2:
                        symbols.symbol_text[iref] = 'thetal, phim, thetam, chim ' + '(degrees - separated by spaces)'
                    elif nvaryt == 3:
                        symbols.symbol_text[iref] = 'thetal, phil, thetam, chim ' + '(degrees - separated by spaces)'
                    elif nvaryt == 4:
                        symbols.symbol_text[iref] = 'thetal, phil, phim, chim ' + '(degrees - separated by spaces)'
                    elif nvaryt == 5:
                        symbols.symbol_text[iref] = 'thetal, phil, phim, thetam ' + '(degrees - separated by spaces)'
                symbols.symbol_value[iref] = symbols.symbol_default[icount + 1]
            else:
                ierr = 1

    return nout, icount, ifilecount, ierr


@routine
def read_plot2d(nin, nout, ierrout):
    """plot_routines.f:read_plot2d -- sets up the plot2d input files.

    Returns (nout, ierrout).
    """
    from . import file_and_directory_commands

    # NB: nvaryt is not initialised in the Fortran (set by the read before use)
    nvaryt = 0
    #
    # return if an error occurred on trying to open the input file
    #
    if ierrout != 0:
        units.close(nin)
        return nout, ierrout
    units.write(6, lw(substr(symbols.fname_run, 1, symbols.ilenfname_run)))
    #
    #  If O.K. read in the data from the specified input file
    #
    symbols.symbol_value[1] = fstr(substr(symbols.fname_run, 1, symbols.ilenfname_run), 256)
    #
    # commands to read the rest of the input file
    #
    icount = 1
    ifilecount = 1
    ierr = 0
    while icount < symbols.nsymbol and ierr == 0:
        icount = icount + 1
        iref = file_and_directory_commands.ireff(icount, ifilecount)
        symbols.symbol_save = symbols.symbol_value[iref]
        commands.sometext, ierr = _read_line(nin, commands.sometext, 3020)
        if ierr == 0:
            _parse(5)
            if commands.nwords > 0:
                if icount == symbols.nsymboluptodata + 1:
                    _, _, symbols.symbol_value[iref], ilen = parse.concatenate(
                        1, 3, symbols.symbol_value[iref], 0,
                        commands.nwords, commands.ncf, commands.ncl, commands.sometext)
                    iref1 = iref + 1
                    _, _, symbols.symbol_value[iref1], ilen1 = parse.concatenate(
                        4, 4, symbols.symbol_value[iref1], 0,
                        commands.nwords, commands.ncf, commands.ncl, commands.sometext)
                    icount = icount + 1
                else:
                    _, commands.nwords, symbols.symbol_value[iref], ilen = parse.concatenate(
                        1, commands.nwords, symbols.symbol_value[iref], 0,
                        commands.nwords, commands.ncf, commands.ncl, commands.sometext)
                    if icount == symbols.nsymboluptodata:
                        symbols.nfiles = read_internal(symbols.symbol_value[iref]).int(symbols.nfiles)

                    # if symbol_value .eq. nvary has changed need to change symbol_text for the
                    # next line.

                    if symbols.symbol[icount] == 'nvary':
                        (nvaryt,), ierr2 = _read_list(symbols.symbol_value[iref], 'i', (nvaryt,))
                        if nvaryt == 1:
                            symbols.symbol_text[iref + 1] = 'phim, thetam, chim ' + '(degrees - separated by spaces)'
                        elif nvaryt == 2:
                            symbols.symbol_text[iref + 1] = 'thetal, phil, chim ' + '(degrees - separated by spaces)'
                        elif nvaryt == 3:
                            symbols.symbol_text[iref + 1] = 'thetal, phil, phim ' + '(degrees - separated by spaces)'

            else:
                symbols.symbol_value[iref] = ''
            if icount >= symbols.nsymboldata and ifilecount < symbols.nfiles:
                ifilecount = ifilecount + 1
                icount = symbols.nsymboluptodata
    units.close(nin)
    #
    # determine the number of coefficients in the coefficients file and whether the
    # molecules are identical or distinct
    #
    nout, ierrout = find_plot_par(nin, nout, ierrout)
    return nout, ierrout


@routine
def read_plot3d(nin, nout, ierrout):
    """plot_routines.f:read_plot3d -- sets up the plot3d input files.

    Returns (nout, ierrout).
    """
    from . import file_and_directory_commands

    # NB: nvaryt is not initialised in the Fortran (set by the read before use)
    nvaryt = 0
    #
    # return if an error occurred on trying to open the input file
    #
    if ierrout != 0:
        units.close(nin)
        return nout, ierrout
    units.write(6, lw(substr(symbols.fname_run, 1, symbols.ilenfname_run)))
    #
    #  If O.K. read in the data from the specified input file
    #
    symbols.symbol_value[1] = fstr(substr(symbols.fname_run, 1, symbols.ilenfname_run), 256)
    #
    # commands to read the rest of the input file
    #
    icount = 1
    ifilecount = 1
    ierr = 0
    while icount < symbols.nsymbol and ierr == 0:
        icount = icount + 1
        iref = file_and_directory_commands.ireff(icount, ifilecount)
        symbols.symbol_save = symbols.symbol_value[iref]
        commands.sometext, ierr = _read_line(nin, commands.sometext, 3020)
        if ierr == 0:
            _parse(5)
            if commands.nwords > 0:
                if icount == symbols.nsymboluptodata + 1:
                    _, _, symbols.symbol_value[iref], ilen = parse.concatenate(
                        1, 1, symbols.symbol_value[iref], 0,
                        commands.nwords, commands.ncf, commands.ncl, commands.sometext)
                    iref1 = iref + 1
                    _, _, symbols.symbol_value[iref1], ilen1 = parse.concatenate(
                        2, 4, symbols.symbol_value[iref1], 0,
                        commands.nwords, commands.ncf, commands.ncl, commands.sometext)
                    iref2 = iref + 2
                    _, _, symbols.symbol_value[iref2], ilen2 = parse.concatenate(
                        5, 7, symbols.symbol_value[iref2], 0,
                        commands.nwords, commands.ncf, commands.ncl, commands.sometext)
                    icount = icount + 2
                else:
                    _, commands.nwords, symbols.symbol_value[iref], ilen = parse.concatenate(
                        1, commands.nwords, symbols.symbol_value[iref], 0,
                        commands.nwords, commands.ncf, commands.ncl, commands.sometext)
                    if icount == symbols.nsymboluptodata:
                        symbols.nfiles = read_internal(symbols.symbol_value[iref]).int(symbols.nfiles)

                    # if symbol_value .eq. nvary has changed need to change symbol_text for the
                    # next line.

                    if symbols.symbol[icount] == 'nvary':
                        (nvaryt,), ierr2 = _read_list(symbols.symbol_value[iref], 'i', (nvaryt,))
                        if nvaryt == 1:
                            symbols.symbol_text[iref + 1] = 'phim, thetam, chim ' + '(degrees - separated by spaces)'
                        elif nvaryt == 2:
                            symbols.symbol_text[iref + 1] = 'thetal, phil, chim ' + '(degrees - separated by spaces)'
                        elif nvaryt == 3:
                            symbols.symbol_text[iref + 1] = 'thetal, phil, phim ' + '(degrees - separated by spaces)'

            else:
                symbols.symbol_value[iref] = ''
            if icount >= symbols.nsymboldata and ifilecount < symbols.nfiles:
                ifilecount = ifilecount + 1
                icount = symbols.nsymboluptodata
    units.close(nin)
    #
    # determine the number of coefficients in the coefficients file and whether the
    # molecules are identical or distinct
    #
    nout, ierrout = find_plot_par(nin, nout, ierrout)
    return nout, ierrout


@routine
def find_plot_par(nin, nout, ierr):
    """plot_routines.f:find_plot_par -- finds certain plotting parameters from the SHARM
    input file.

    Returns (nout, ierr).
    """
    from . import auxiliary_routines, file_and_directory_commands

    # search for the symbol "shcoeffs"

    moretext = fstr('shcoeffs', 20)
    ifound = 0
    ifound = file_and_directory_commands.search_for_symbol(moretext, ifound)
    iref = file_and_directory_commands.ireff(ifound, symbols.nfiles)
    #
    # generate the SHARM input filename
    #
    ilensymbol_value = auxiliary_routines.indexlastdot(symbols.symbol_value[iref]) - 1
    if ilensymbol_value <= 0:
        ilensymbol_value = findex_padded(symbols.symbol_value[iref], ' ', 256) - 1
    fnametemp = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                     + substr(symbols.symbol_value[iref], 1, ilensymbol_value) + '.dat', 256)
    units.write(nout, lw(substr(fnametemp, 1, directory_name.nchardir + ilensymbol_value + 4)))
    ierr = _open(nin, fnametemp, 'old')
    if ierr != 0:
        units.write(nout, fmt(_FPP100, _prompt(), substr(symbols.symbol_value[iref], 1, ilensymbol_value) + '.dat'))
        fnametemp = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                         + substr(symbols.symbol_value[iref], 1, ilensymbol_value) + '.h01', 256)
        text = fstr('.h01', 10)
        nin1 = nin + 1
        # NB: the Fortran calls get_nterms with 4 of its 5 arguments (no ierr); the
        # missing status is discarded here.  ncoeffs is not initialised (0 here).
        ncoeffs = 0
        fnametemp, ncoeffs, _ = file_and_directory_commands.get_nterms(nin1, fnametemp, text, ncoeffs, 0)
        #
        # write this to symbol_value(iref+1)
        #
        symbols.symbol_value[iref + 1] = fstr(lw(ncoeffs), 256)
        return nout, ierr
    #
    # read lines of text from this file, looking for specific keywords
    #
    ierr = 0
    centrec = ''
    centres = ''
    while ierr == 0 and (centrec == '' or centres == ''):
        commands.sometext, ierr = _read_line(nin, commands.sometext, 3020)
        _parse(11)
        #
        # check for keyword
        #
        if _word(1) == 'ncoeffs':
            symbols.symbol_value[iref + 1] = fstr(_word(2), 256)
        elif _word(1) == 'atom-c':
            centrec = fstr(_words(2, commands.nwords), 256)
        elif _word(1) == 'atom-s':
            centres = fstr(_words(2, commands.nwords), 256)
    units.close(nin)
    #
    # if centrec or centres have not been set it probably means there is an error
    # with the .dat file
    #
    if centrec == '' or centres == '':
        units.write(nout, fmt(_FPP103, _prompt(), substr(symbols.symbol_value[iref], 1, ilensymbol_value) + '.dat'))
        fnametemp = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                         + substr(symbols.symbol_value[iref], 1, ilensymbol_value) + '.h01', 256)
        text = fstr('.h01', 10)
        nin1 = nin + 1
        # NB: get_nterms called with 4 of its 5 arguments, as above
        ncoeffs = 0
        fnametemp, ncoeffs, _ = file_and_directory_commands.get_nterms(nin1, fnametemp, text, ncoeffs, 0)
        #
        # write this to symbol_value(iref+1)
        #
        symbols.symbol_value[iref + 1] = fstr(lw(ncoeffs), 256)
        ierr = 0
    else:
        ierr = 0
    #
    # decide whether molecules were the same or different
    #
    if centrec != centres:
        symbols.symbol_value[iref + 2] = '1'
    else:
        symbols.symbol_value[iref + 2] = '0'

    return nout, ierr


@routine
def save_plot_files(nin, nout):
    """plot_routines.f:save_plot_files -- choose the input file name and save the
    current values (save_files_plots / save_files_plot2d / save_files_plot3d).

    Returns nout.
    """
    from . import file_and_directory_commands

    #
    # first check whether the existing file name is O.K.
    #
    ierrfile = 1
    while ierrfile != 0:

        if symbols.ilenfname_inp > 0:

            units.write(nout, fmt(_SVP303, _prompt(), substr(symbols.fname_inp, 1, symbols.ilenfname_inp)))
            #     *//'Specify the new filename or press <Enter> to accept this name: '
            #     *,$)
            #           read(nin,100,iostat=ierr) sometext
            #           call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            #           if(nwords.gt.0) then
            #              call concatenate(1,nwords,fname_inp,ilenfname_inp
            #     *,nwords,ncf,ncl,sometext)
            #           endif

        else:

            commands.nwords = 0
            while commands.nwords == 0:
                units.write(nout, fmt(_SVP302, _prompt(), substr(symbols.symbol_setup, 1, symbols.ilensetup)),
                            advance=False)
                commands.sometext, ierr = _read_line(nin, commands.sometext, 3020)
                _parse(5)
                if commands.nwords > 0:
                    _, commands.nwords, symbols.fname_inp, symbols.ilenfname_inp = parse.concatenate(
                        1, commands.nwords, symbols.fname_inp, symbols.ilenfname_inp,
                        commands.nwords, commands.ncf, commands.ncl, commands.sometext)

        #
        # check the name for the expected extension - if not present it is added on
        #
        # first strip out any extension that may have already been supplied
        #
        ilenuptodot = findex(symbols.fname_inp, '.') - 1
        if ilenuptodot > 0:
            symbols.ilenfname_inp = ilenuptodot
        else:
            symbols.ilenfname_inp = findex_padded(symbols.fname_inp, '  ', 256) - 1
        #
        # if the expected extension has characters before the . ,then check if these are present
        # - anything beyond these characters if found will be ignored
        #
        if symbols.ilenextinpdot > 1:
            itest = findex_padded(symbols.fname_inp, substr(symbols.extinp, 1, symbols.ilenextinpdot - 1), 256)
            if itest > 0:
                symbols.ilenfname_inp = itest - 1
        #
        # reconstruct the input file name with the correct extension
        #
        if symbols.ilenfname_inp > 0:
            symbols.fname_inp = fstr(substr(symbols.fname_inp, 1, symbols.ilenfname_inp)
                                     + substr(symbols.extinp, 1, symbols.ilenextinp), 256)
        else:
            symbols.fname_inp = fstr(substr(symbols.extinp, 1, symbols.ilenextinp), 256)
        symbols.ilenfname_inp = symbols.ilenfname_inp + symbols.ilenextinp
        units.write(nout, fmt(_SVP304, _prompt(), substr(symbols.fname_inp, 1, symbols.ilenfname_inp)))
        #
        # create the full name of the input file by adding on the directory name if not present

        symbols.fnameinp, symbols.ilenfnameinp = file_and_directory_commands.get_full_filename(
            symbols.fname_inp, symbols.ilenfname_inp, symbols.fnameinp, symbols.ilenfnameinp)
        #
        # open input file. For plot files we force it to accept the current filename
        #
        if symbols.symbol_setup == 'plot':
            ierrfile = 0
        else:
            ierrfile = _open(14, symbols.fnameinp, 'new')
        if ierrfile != 0:
            commands.nwords = 0
            while commands.nwords == 0:
                commands.sometext = ''
                units.write(6, fmt(_SVP301, _prompt(), substr(symbols.fname_inp, 1, symbols.ilenfname_inp)),
                            advance=False)
                commands.sometext, ierr = _read_line(5, commands.sometext, 3020)
                _parse(5)
            if _word(1) == 'y':
                ierrfile = 0
            else:
                symbols.ilenfname_inp = 0
    #
    # so now go on and save the data!
    #
    if symbols.symbol_setup == 'plot2d':
        save_files_plot2d(14)
    elif symbols.symbol_setup == 'plot3d':
        save_files_plot3d(14)
    else:
        save_files_plots(14)

    units.close(14)
    return nout


@routine
def save_files_plots(nout):
    """plot_routines.f:save_files_plots -- writes the general input file.

    Returns nout.
    """
    from . import auxiliary_routines, file_and_directory_commands

    #
    # open the new input file
    #
    ierr = _open(nout, symbols.fnameinp, 'unknown')
    if ierr == 0:
        #
        # generate the run name - this should be the same as the input file but without
        # the last .extension
        #
        symbols.ilenfname_run = auxiliary_routines.indexlastdot(symbols.fname_inp) - 1
        symbols.fname_run = fstr(substr(symbols.fname_inp, 1, symbols.ilenfname_run), 256)
        #
        # write new input file
        #
        units.write(nout, fmt(_SFP229, substr(symbols.fname_run, 1, symbols.ilenfname_run)))
        ierr = 0
        icount = 0
        ifilecount = 1
        while icount < symbols.nsymboldata:
            icount = icount + 1
            iref = file_and_directory_commands.ireff(icount, ifilecount)
            if icount == symbols.nsymboluptodata + 1:
                units.write(nout, fmt(_SFP204, substr(symbols.dataprompt, 1, symbols.ilendataprompt), ifilecount))
                ierr = 0
            commands.sometext = fstr(symbols.symbol_value[iref], 3020)
            _parse(5)
            units.write(nout, fmt(_SFP101, Ch(symbols.symbol[icount], 11), _words(1, max(1, commands.nwords)),
                                  Ch(symbols.symbol_text[icount], 256)))
            ierr = 0
            if icount >= symbols.nsymboldata and ifilecount < symbols.nfiles:
                ifilecount = ifilecount + 1
                icount = symbols.nsymboluptodata
        #
        # write a 'q' to signify the end of input (but save is not needed)
        #
        units.write(nout, fmt(_SFP102))
        ierr = 0

    if ierr != 0:

        units.write(6, fmt(_SFP103))
        return nout

    return nout


@routine
def save_files_plot2d(nout):
    """plot_routines.f:save_files_plot2d -- writes the plot2d input file.

    Returns nout.
    """
    from . import file_and_directory_commands

    #
    # open the new input file
    #
    ierr = _open(nout, symbols.fnameinp, 'unknown')  # noqa: F841
    #
    # write .txt file
    #
    icount = 0
    ifilecount = 1
    while icount < symbols.nsymbol:
        icount = icount + 1
        iref = file_and_directory_commands.ireff(icount, ifilecount)
        #
        # lines containing sphere parameters has to be printed as a single line
        #
        if icount == symbols.nsymboluptodata + 1:
            commands.sometext = fstr(symbols.symbol_value[iref], 3020)
            _parse(5)
            _, commands.nwords, symbols.symbol_value[iref], ilen = parse.concatenate(
                1, commands.nwords, symbols.symbol_value[iref], 0,
                commands.nwords, commands.ncf, commands.ncl, commands.sometext)
            iref1 = iref + 1
            commands.sometext = fstr(symbols.symbol_value[iref1], 3020)
            _parse(5)
            _, commands.nwords, symbols.symbol_value[iref1], ilen1 = parse.concatenate(
                1, commands.nwords, symbols.symbol_value[iref1], 0,
                commands.nwords, commands.ncf, commands.ncl, commands.sometext)
            commands.sometext = fstr(substr(symbols.symbol_value[iref], 1, ilen)
                                     + ' ' + substr(symbols.symbol_value[iref1], 1, ilen1), 3020)
            units.write(nout, fmt(_SF2D102, substr(commands.sometext, 1, ilen + ilen1 + 1)))
            # c102  format(a,11x,'(x,y,z) of circle, and circle radius')
            icount = icount + 1
        else:
            commands.sometext = fstr(symbols.symbol_value[iref], 3020)
            _parse(5)
            #
            # character values cannot be printed with the corresponding text
            #
            if symbols.symbol_type[icount] == 'c':
                units.write(nout, fmt(_SF103, _words(1, max(1, commands.nwords))))
            else:
                units.write(nout, fmt(_SF101, _words(1, max(1, commands.nwords)),
                                      Ch(symbols.symbol_text[icount], 256)))
        if icount >= symbols.nsymboldata and ifilecount < symbols.nfiles:
            ifilecount = ifilecount + 1
            icount = symbols.nsymboluptodata

    return nout


@routine
def save_files_plot3d(nout):
    """plot_routines.f:save_files_plot3d -- writes the EPSR .inp and .pcof files
    (the plot3d input file).

    Returns nout.
    """
    from . import file_and_directory_commands

    #
    # open the new input file
    #
    ierr = _open(nout, symbols.fnameinp, 'unknown')  # noqa: F841
    #
    # write .txt file
    #
    icount = 0
    ifilecount = 1
    while icount < symbols.nsymbol:
        icount = icount + 1
        iref = file_and_directory_commands.ireff(icount, ifilecount)
        #
        # lines containing sphere parameters has to be printed as a single line
        #
        if icount == symbols.nsymboluptodata + 1:
            commands.sometext = fstr(symbols.symbol_value[iref], 3020)
            _parse(5)
            _, commands.nwords, symbols.symbol_value[iref], ilen = parse.concatenate(
                1, commands.nwords, symbols.symbol_value[iref], 0,
                commands.nwords, commands.ncf, commands.ncl, commands.sometext)
            iref1 = iref + 1
            commands.sometext = fstr(symbols.symbol_value[iref1], 3020)
            _parse(5)
            _, commands.nwords, symbols.symbol_value[iref1], ilen1 = parse.concatenate(
                1, commands.nwords, symbols.symbol_value[iref1], 0,
                commands.nwords, commands.ncf, commands.ncl, commands.sometext)
            iref2 = iref + 2
            commands.sometext = fstr(symbols.symbol_value[iref2], 3020)
            _parse(5)
            _, commands.nwords, symbols.symbol_value[iref2], ilen2 = parse.concatenate(
                1, commands.nwords, symbols.symbol_value[iref2], 0,
                commands.nwords, commands.ncf, commands.ncl, commands.sometext)
            commands.sometext = fstr(substr(symbols.symbol_value[iref], 1, ilen)
                                     + ' ' + substr(symbols.symbol_value[iref1], 1, ilen1)
                                     + ' ' + substr(symbols.symbol_value[iref2], 1, ilen2), 3020)
            units.write(nout, fmt(_SF3D102, substr(commands.sometext, 1, ilen + ilen1 + ilen2 + 2)))
            icount = icount + 2
        else:
            commands.sometext = fstr(symbols.symbol_value[iref], 3020)
            _parse(5)
            #
            # character values cannot be printed with the corresponding text
            #
            if symbols.symbol_type[icount] == 'c':
                units.write(nout, fmt(_SF103, _words(1, max(1, commands.nwords))))
            else:
                units.write(nout, fmt(_SF101, _words(1, max(1, commands.nwords)),
                                      Ch(symbols.symbol_text[icount], 256)))
        if icount >= symbols.nsymboldata and ifilecount < symbols.nfiles:
            ifilecount = ifilecount + 1
            icount = symbols.nsymboluptodata

    return nout
