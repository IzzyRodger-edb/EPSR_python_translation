"""Translation of fortran_src/read_and_save_files.f: reading and saving EPSR input files.

read_input_file (the .inp command file), histogram helpers (initgofz, initgofzs, outgofz,
outgofz_sharm, outgofzs_sharm), saving input files (save_files, force_save_files,
save_files_others, which also writes the .pcof and .ato files for EPSR), and the readers
of the .pcof (read_pcof) and .ato (read_ato) files named in the input.

Routines of other packages (process_symbols, file_and_directory_commands,
ato_read_and_save, FCreadline_windows, reallocation_routines, coordconstraint_mod) and
the Fortran-module data (threebody_mod, coordconstraints_mod) are imported inside the
routines that use them, to avoid import cycles.
"""

from ..common import (
    clustercom,
    commands,
    epsr16,
    iresetnumbers,
    molecom,
    potential_par,
    symbols,
)
from ..fortran import (
    PI4ATAN,
    Ch,
    FortranEOF,
    FortranIOError,
    cosf,
    dsqrt,
    f32,
    findex,
    fmt,
    fstr,
    idiv,
    lw,
    pad,
    r4,
    read_internal,
    real_array,
    sign,
    substr,
    units,
)
from ..trace import routine
from . import auxiliary_routines, parse


# --------------------------------------------------------------------------- helpers
def _read_real(text, old):
    """``read(text,*,iostat=ierr) x`` for a real: returns (x, ierr); x unchanged on error."""
    try:
        return read_internal(text).real(old), 0
    except FortranEOF:
        return old, -1
    except FortranIOError as e:
        return old, e.iostat


def _read_int(text, old):
    """``read(text,*,iostat=ierr) i`` for an integer: returns (i, ierr)."""
    try:
        return read_internal(text).int(old), 0
    except FortranEOF:
        return old, -1
    except FortranIOError as e:
        return old, e.iostat


def _star_repeat(formatline):
    """A FORMAT built with ``i1.1`` gets ``*`` as repeat count when the count is >= 10;
    gfortran then treats ``*(...)`` as an unlimited-repeat group, which formats the
    remaining items exactly as a large repeat count does (fformat has no ``*(``)."""
    return formatline.replace("*(", "999(")


def _w(word_no):
    """sometext(ncf(i):ncl(i)) of the COMMON command text."""
    return substr(commands.sometext, commands.ncf[word_no], commands.ncl[word_no])


def _parse(nspcword):
    """call parse(sometext,nwords,ncf,ncl,nspcword,nstsize,nworddim) on COMMON /commands/."""
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl,
                                  nspcword, commands.nstsize, commands.nworddim)


def _concatenate(nfirst, nlast, return_string, ncharstring):
    """call concatenate(nfirst,nlast,return_string,ncharstring,nwords,ncf,ncl,sometext)
    on COMMON /commands/ (return_string is character*256 in every call of this file).
    Returns (return_string, ncharstring); nfirst/nlast are literals or nwords itself
    here, which concatenate leaves unchanged."""
    _, _, return_string, ncharstring = parse.concatenate(
        nfirst, nlast, return_string, ncharstring, commands.nwords, commands.ncf, commands.ncl,
        commands.sometext)
    return return_string, ncharstring


# --------------------------------------------------------------------------- read_input_file
_F100 = "(a)"


@routine
def read_input_file(nin, nout, ierrout):
    """read_and_save_files.f:read_input_file -- reads the EPSR input files.

    Commands are read from the already opened unit nin and processed in turn
    (process_commands); at the end the file is closed.
    Returns (nout, ierrout).
    """
    from . import file_and_directory_commands, process_symbols  # noqa: F401  (packages A, C)

    # character*20 command
    #
    # don't do anything if there is an error in the input file - will start out with
    # default values in this case
    #
    if ierrout != 0:
        units.close(nin)
        return nout, ierrout

    #
    # set up the symbol numbers for the fnameato and fnamepcof symbols
    #
    ifnameato = 0  # noqa: F841
    ifnamepcof = 0  # noqa: F841
    #
    # set the number of potential coefficient pairs to zero - this is used to
    # signal whether a valid .pcof file has been read in
    #
    potential_par.npairspcof = 0
    #
    # signal whether any atoms have been read in
    #
    molecom.ntype = 0
    #
    # Read in commands from the already opened input file, and process each command
    # in turn. At the end close this file.
    #
    command = "next"
    icount = 1
    ifilecount = 1
    #
    # we assume that the symbols are specified in order, but the
    # process_commands routine will check each command in turn to ensure the correct
    # values are updated. The last one should be a 'q' to signify the end of
    # command input
    #
    ierr = 0
    while command != "e" and ierr == 0:
        #
        # check the status of the file counter
        #
        if icount > symbols.nsymboldata and ifilecount < symbols.nfiles:
            ifilecount = ifilecount + 1
            icount = symbols.nsymboluptodata + 1
        elif icount > symbols.nsymbol:
            icount = 1
        #
        # save the symbol number for fnameato and fnamepcof in case these files need to
        # be read as well
        #
        if symbols.symbol[icount] == "fnameato":
            ifnameato = icount  # noqa: F841
        if symbols.symbol[icount] == "fnamepcof":
            ifnamepcof = icount  # noqa: F841
        #
        # read the commands
        #
        commands.nwords = 0
        # read(nin,100,iostat=ierr) sometext -- at end of file sometext keeps its value
        try:
            commands.sometext = fstr(units.read_line(nin), commands.nstsize)
            ierr = 0
        except FortranEOF:
            ierr = -1
        except FortranIOError as e:
            ierr = e.iostat
        _parse(12)

        # Need to ignore blank lines read from an input file

        if commands.nwords > 0:
            nout, command, icount, ifilecount, ierrout = process_symbols.process_commands(
                nin, nout, command, icount, ifilecount, ierrout, 1)
        #         write(6,*) 'read_input_file> ',icount,ifilecount
    #
    # finished input - close the input file
    #
    units.close(nin)
    return nout, ierrout


# --------------------------------------------------------------------------- histograms
@routine
def initgofz(nz, ndist, gofz, nxdim, nfdim):
    """read_and_save_files.f:initgofz -- initialisation routine for histograms (zero gofz)."""
    #
    # zero gofz
    #
    for j in range(1, nz + 1):
        for id_ in range(1, ndist + 1):
            gofz[j, id_] = 0.0


@routine
def initgofzs(nz, ndist, gofz, gofzs, nxdim, nfdim):
    """read_and_save_files.f:initgofzs -- initialisation routine for histograms
    (zero gofz and gofzs)."""
    #
    # zero gofz
    #
    for j in range(1, nz + 1):
        for id_ in range(1, ndist + 1):
            gofz[j, id_] = 0.0
            gofzs[j, id_] = 0.0


def _outgofz_body(fnameinp, ext, nz, ndist, nsum, z, gofz, gzt, gzt2, titler, titledat, nxdim, nfdim):
    """Body shared by outgofz and outgofz_sharm (the two routines are identical)."""
    # real*8 fc1,fc2,a,b,a1,a2,b1,b2,error
    #
    # set up the name of the output file with extension .n
    #
    # strip off the last . extension and add a new one
    #
    ilen = auxiliary_routines.indexlastdot(fnameinp)
    fname = fstr(substr(fnameinp, 1, ilen) + substr(ext, 1, 1), 256)
    #
    # read in data to be fitted if nsum is greater than zero
    #
    # NB: lptin is not initialised in the Fortran when r_int fails before reading;
    # it is then not used (nsum is set to -1).  0 is passed.
    lptin = 0
    ierr = 0
    if nsum > 0:
        lptin, fname, ierr = auxiliary_routines.r_int(lptin, fname, ndist, z, gzt, gzt2, nxdim, nfdim, ierr)
    else:
        lptin = nz
    #
    # if an error occurs on opening this file, set nsum to -1 to start a new
    # accumulation
    #
    if ierr != 0:
        nsum = -1
    if nsum > 0:
        for id_ in range(1, ndist + 1):
            for j in range(1, nz + 1):
                if gzt2[j, id_] < r4(1.0e-20):
                    gzt2[j, id_] = 0.0
                if abs(gzt[j, id_]) < r4(1.0e-20):
                    gzt[j, id_] = 0.0
                # check that number of points read in is not less than number expected
                if j <= lptin:
                    gzt2[j, id_] = gzt2[j, id_] * gzt2[j, id_] + gzt[j, id_] * gzt[j, id_]
                else:
                    gzt[j, id_] = 0.0
                    gzt2[j, id_] = 0.0
    else:
        for id_ in range(1, ndist + 1):
            for j in range(1, nz + 1):
                gzt[j, id_] = 0.0
                gzt2[j, id_] = 0.0
    if nsum >= 0:
        nsum = nsum + 1
        fc1 = float(r4(1.0) / f32(nsum))  # 1.0/float(nsum) is single precision
        fc2 = 1.0 - fc1
    else:
        fc1 = 1.0
        fc2 = 0.0
    lpt = nz  # #points
    for id_ in range(1, ndist + 1):
        for i in range(1, nz + 1):
            a2 = fc2 * float(gzt[i, id_])
            a1 = fc1 * float(gofz[i, id_])
            a = a1 + a2
            b2 = fc2 * float(gzt2[i, id_])
            b1 = fc1 * float(gofz[i, id_] * gofz[i, id_])  # (gofz(i,id))**2 in real*4
            b = b1 + b2
            error = b - a * a
            if error > 0.0:
                error = dsqrt(error)
            else:
                error = 0.0
            gzt[i, id_] = a
            gzt2[i, id_] = error
    fname = auxiliary_routines.w_int(lpt, fname, ndist, z, gzt, gzt2, titler, titledat, nxdim, nfdim)  # read data
    return nsum


@routine
def outgofz(fnameinp, ext, nz, ndist, nsum, z, gofz, gzt, gzt2, titler, titledat, nxdim, nfdim):
    """read_and_save_files.f:outgofz -- summing and output routine for histograms.

    The running average (and its error) of gofz over nsum calls is kept in the file
    fnameinp-without-extension.<ext>01... (read back with r_int, written with w_int).
    Returns nsum.
    """
    return _outgofz_body(fnameinp, ext, nz, ndist, nsum, z, gofz, gzt, gzt2, titler, titledat, nxdim, nfdim)


@routine
def outgofz_sharm(fnameinp, ext, nz, ndist, nsum, z, gofz, gzt, gzt2, titler, titledat, nxdim, nfdim):
    """read_and_save_files.f:outgofz_sharm -- summing and output routine for histograms
    (identical to outgofz; its local arrays ll1, ll2, ll, n1, n2 are never used).

    Returns nsum.
    """
    return _outgofz_body(fnameinp, ext, nz, ndist, nsum, z, gofz, gzt, gzt2, titler, titledat, nxdim, nfdim)


@routine
def outgofzs_sharm(fnameinp, ext, nz, ndist, nsum, z, gofz, gofzs, gzt, gzt2, titler, titledat, nxdim, nfdim):
    """read_and_save_files.f:outgofzs_sharm -- summing and output routine for histograms,
    gofzs being accumulated alongside gofz (no error estimate).

    (The local arrays ll1, ll2, ll, n1, n2 of the Fortran are never used.)
    Returns nsum.
    """
    #
    # set up the name of the output file with extension .n
    #
    # strip off the last . extension and add a new one
    #
    ilen = auxiliary_routines.indexlastdot(fnameinp)
    fname = fstr(substr(fnameinp, 1, ilen) + substr(ext, 1, 1), 256)
    #
    # read in data to be fitted if nsum is greater than zero
    #
    # NB: lptin is not initialised in the Fortran when r_int fails before reading;
    # it is then not used (nsum is set to -1).  0 is passed.
    lptin = 0
    ierr = 0
    if nsum > 0:
        lptin, fname, ierr = auxiliary_routines.r_int(lptin, fname, ndist, z, gzt, gzt2, nxdim, nfdim, ierr)
    else:
        lptin = nz
    #
    # if an error occurs on opening this file, set nsum to -1 to start a new
    # accumulation
    #
    if ierr != 0:
        nsum = -1
    if nsum > 0:
        for id_ in range(1, ndist + 1):
            for j in range(1, nz + 1):
                # check that number of points read in is not less than number expected
                if j > lptin:
                    gzt[j, id_] = 0.0
                    gzt2[j, id_] = 0.0
    else:
        for id_ in range(1, ndist + 1):
            for j in range(1, nz + 1):
                gzt[j, id_] = 0.0
                gzt2[j, id_] = 0.0
    if nsum >= 0:
        nsum = nsum + 1
        fc1 = float(r4(1.0) / f32(nsum))  # 1.0/real(nsum) is single precision
        fc2 = 1.0 - fc1
    else:
        fc1 = 1.0
        fc2 = 0.0
    lpt = nz  # #points
    for id_ in range(1, ndist + 1):
        for i in range(1, nz + 1):
            a2 = fc2 * float(gzt[i, id_])
            a1 = fc1 * float(gofz[i, id_])
            b2 = fc2 * float(gzt2[i, id_])
            b1 = fc1 * float(gofzs[i, id_])
            gzt[i, id_] = a1 + a2
            gzt2[i, id_] = b1 + b2
    fname = auxiliary_routines.w_int(lpt, fname, ndist, z, gzt, gzt2, titler, titledat, nxdim, nfdim)  # read data
    return nsum


# --------------------------------------------------------------------------- save_files
_F303 = "(/a,1x,'Current name of file is \"',a,'\"')"
_F302 = "(a,1x,' Specify the name of the ',a,' input file: ')"
_F304 = "(/a,1x,'Writing to input file \"',a,'\"')"
_F301 = "(/a,1x,'File \"',a,'\" already exists.'//)"
_F305 = "(' Do you want to overwrite it (y or n)? ')"


def _ask_input_filename(nout):
    """The 'else' part of the file-name loop of save_files / force_save_files: prompt
    until a non-blank file name is given."""
    from . import FCreadline_windows  # package A (module jmk_readline_module)

    commands.nwords = 0
    while commands.nwords == 0:
        units.write(nout, lw())
        promptstring = fstr(fmt(_F302, substr(symbols.prompt, 1, symbols.ilenprompt),
                                substr(symbols.symbol_setup, 1, symbols.ilensetup)), 256)
        lenpromptstring = len(promptstring)
        commands.sometext = fstr(FCreadline_windows.readline(commands.sometext, substr(promptstring, 1, lenpromptstring),
                                                             len_line=commands.nstsize), commands.nstsize)
        _parse(5)
        if commands.nwords > 0:
            symbols.fname_inp, symbols.ilenfname_inp = _concatenate(1, commands.nwords, symbols.fname_inp,
                                                                    symbols.ilenfname_inp)


def _set_input_extension(nout):
    """Common middle part of save_files / force_save_files: check the name for the
    expected extension, write the name and build the full name."""
    from . import file_and_directory_commands  # package A

    #
    # check the name for the expected extension - if not present it is added on
    #
    # first strip out any extension that may have already been supplied
    #
    ilenuptodot = findex(symbols.fname_inp, ".") - 1
    if ilenuptodot > 0:
        symbols.ilenfname_inp = ilenuptodot
    else:
        symbols.ilenfname_inp = len(symbols.fname_inp)
    #
    # if the expected extension has characters before the . ,then check if these are present
    # - anything beyond these characters if found will be ignored
    #
    if symbols.ilenextinpdot > 1:
        itest = findex(symbols.fname_inp, substr(symbols.extinp, 1, symbols.ilenextinpdot - 1))
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
    units.write(nout, fmt(_F304, substr(symbols.prompt, 1, symbols.ilenprompt),
                          substr(symbols.fname_inp, 1, symbols.ilenfname_inp)))
    #
    # create the full name of the input file by adding on the directory name if not present

    symbols.fnameinp, symbols.ilenfnameinp = file_and_directory_commands.get_full_filename(
        symbols.fname_inp, symbols.ilenfname_inp, symbols.fnameinp, symbols.ilenfnameinp)


@routine
def save_files(nin, nout):
    """read_and_save_files.f:save_files -- ask for (or confirm) the input file name, then
    save the data (save_files_changeato or save_files_others on unit 14).

    Returns nout.
    """
    from . import FCreadline_windows, process_symbols  # packages A, C

    # character*256 promptstring
    #
    # first check whether the existing file name is O.K.
    #
    ierrfile = 1
    while ierrfile != 0:

        if symbols.ilenfname_inp > 0:

            units.write(nout, fmt(_F303, substr(symbols.prompt, 1, symbols.ilenprompt),
                                  substr(symbols.fname_inp, 1, symbols.ilenfname_inp)))

        else:

            _ask_input_filename(nout)

        _set_input_extension(nout)
        #
        # open input file. For plot files we force it to accept the current filename
        #
        if symbols.symbol_setup == "plot":
            ierrfile = 0
        else:
            try:
                units.open(14, symbols.fnameinp, status="new")
                ierrfile = 0
            except FortranIOError as e:
                ierrfile = e.iostat
        if ierrfile != 0:
            units.write(nout, fmt(_F301, substr(symbols.prompt, 1, symbols.ilenprompt),
                                  substr(symbols.fname_inp, 1, symbols.ilenfname_inp)))
            commands.nwords = 0
            while commands.nwords == 0:
                commands.sometext = ""
                promptstring = fstr(fmt(_F305), 256)
                lenpromptstring = len(promptstring)
                commands.sometext = fstr(FCreadline_windows.readline(commands.sometext,
                                                                     substr(promptstring, 1, lenpromptstring),
                                                                     len_line=commands.nstsize), commands.nstsize)
                _parse(5)
            if _w(1) == "y":
                ierrfile = 0
            else:
                symbols.ilenfname_inp = 0
    #
    # so now go on and save the data!
    #
    if symbols.symbol_setup == "changeato":
        process_symbols.save_files_changeato(14)
    else:
        save_files_others(14)

    units.close(14)
    return nout


@routine
def force_save_files(nin, nout):
    """read_and_save_files.f:force_save_files -- as save_files, but the file is always
    accepted (overwritten without asking).

    Returns nout.
    """
    from . import process_symbols  # package C

    # character*256 promptstring
    #
    # first check whether the existing file name is O.K.
    #
    ierrfile = 1
    while ierrfile != 0:

        if symbols.ilenfname_inp > 0:

            units.write(nout, fmt(_F303, substr(symbols.prompt, 1, symbols.ilenprompt),
                                  substr(symbols.fname_inp, 1, symbols.ilenfname_inp)))

        else:

            _ask_input_filename(nout)

        _set_input_extension(nout)
        #
        # open input file. For plot files we force it to accept the current filename
        #
        ierrfile = 0
    #
    # so now go on and save the data!
    #
    if symbols.symbol_setup == "changeato":
        process_symbols.save_files_changeato(14)
    else:
        save_files_others(14)

    units.close(14)
    return nout


_F229 = "(a,15x,'Title of this file')"
_F204 = "(/a,1x,i3/)"
_F101 = "(a11,1x,a,15x,a)"
_F102 = "('q')"
_F103 = "(/'save_files> Could not open specified file for writing.')"


@routine
def save_files_others(nout):
    """read_and_save_files.f:save_files_others -- writes the general input file
    (fnameinp) and, for EPSR, the .pcof file and the .ato file.

    Returns nout.
    """
    from ..common import coordconstraints_mod, threebody_mod  # Fortran module data (packages H, G)
    from . import ato_read_and_save, file_and_directory_commands, reallocation_routines  # packages D, A

    ireff = file_and_directory_commands.ireff
    tb = threebody_mod
    # integer*4 iqmin,npgausstest
    # character*256 formatline,fullline,shortline

    piconv = PI4ATAN / r4(180.0)  # 4.0*atan(1.0)/180.0
    #
    # find the symbol reference for qmin
    #
    if symbols.symbol_setup == "epsr":
        iqmin = ireff(iresetnumbers.iqminref, symbols.nfiles)
    else:
        #
        # ensure this value will not be used if it has not been defined
        #
        iqmin = symbols.nsymbol
    #
    # open the new input file
    #
    try:
        units.open(nout, symbols.fnameinp, status="old")
        ierr = 0
    except FortranIOError as e:
        ierr = e.iostat
    if ierr != 0:  # Try opening as an new file
        try:
            units.open(nout, symbols.fnameinp, status="new")
            ierr = 0
        except FortranIOError as e:
            ierr = e.iostat
    # Because gfortran gives an error when a file opened with status unknown does not exist!
    #      ierr=0

    ifilecount = 1
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
        units.write(nout, fmt(_F229, substr(symbols.fname_run, 1, symbols.ilenfname_run)))
        icount = 0
        ifilecount = 1
        while icount < symbols.nsymboldata:
            icount = icount + 1
            # For EPSR files skip over the ereq parameters which are now stored in the .pcof file
            if symbols.symbol_setup == "epsr":
                if icount > iresetnumbers.iereqfirst and icount <= iresetnumbers.iereqlast:
                    icount = iresetnumbers.iereqlast + 1
                elif icount > iqmin and icount < symbols.nsymboluptodata:
                    # for EPSR files skip over the bit related to the .pcof file - this is written
                    # to a separate file
                    icount = symbols.nsymboluptodata
            iref = ireff(icount, ifilecount)
            if icount == symbols.nsymboluptodata + 1:
                units.write(nout, fmt(_F204, substr(symbols.dataprompt, 1, symbols.ilendataprompt), ifilecount))
            commands.sometext = fstr(symbols.symbol_value[iref], commands.nstsize)
            _parse(5)
            units.write(nout, fmt(_F101, Ch(symbols.symbol[icount], 11),
                                  substr(commands.sometext, commands.ncf[1], commands.ncl[max(1, commands.nwords)]),
                                  Ch(symbols.symbol_text[icount], 256)))
            if icount >= symbols.nsymboldata and ifilecount < symbols.nfiles:
                ifilecount = ifilecount + 1
                icount = symbols.nsymboluptodata
        #
        # write a 'q' to signify the end of input (but save is not needed)
        #
        units.write(nout, fmt(_F102))

    if ierr != 0:

        units.write(6, fmt(_F103))
        return nout

    #
    # if symbol_setup .eq. epsr then save the pcof and ato files
    #
    if symbols.symbol_setup == "epsr":

        units.close(nout)
        try:
            units.open(nout, symbols.fnamepcof, status="unknown")
            ierr = 0
        except FortranIOError as e:
            ierr = e.iostat
        icount = iresetnumbers.iereqfirst + 1  # Start from ereq
        while icount < symbols.nsymboluptodata:
            iref = ireff(icount, ifilecount)
            commands.sometext = fstr(symbols.symbol_value[iref], commands.nstsize)
            _parse(5)
            units.write(nout, fmt(_F101, Ch(symbols.symbol[icount], 11),
                                  substr(commands.sometext, commands.ncf[1], commands.ncl[max(1, commands.nwords)]),
                                  Ch(symbols.symbol_text[icount], 256)))
            if icount == iresetnumbers.iereqlast:
                icount = iqmin  # Skip to bit before data files
            icount = icount + 1
        units.write(nout, fmt(_F102))
        #
        # now write the potential parameters
        #
        pp = potential_par
        icount = symbols.nsymboldata
        units.write(nout, lw(pp.npairspcof))
        for ip in range(1, pp.npairspcof + 1):
            #
            # write the atom symbols - these will be used to recognise this atom pair in other
            # simulations
            #
            units.write(nout, fmt("(a3,1x,a3,1x,5x,i0,2x,f10.5)", Ch(pp.tlab1[ip], 3), Ch(pp.tlab2[ip], 3), ip,
                                  pp.coulfac[ip]))
            #
            # save the current values of rminex in case they have been changed
            #
            icount = icount + 1
            iref = ireff(icount, symbols.nfiles)
            commands.sometext = ""
            commands.sometext = fstr(symbols.symbol_value[iref], commands.nstsize)
            _parse(11)
            nwords = commands.nwords

            # No. of Gaussians to write. We need three parameters for each Gaussian, its position
            # its sigma and its height or depth. Any error in reading these and the data are ignored.

            iword = 1
            pp.rminex[ip], ierr = _read_real(_w(iword), pp.rminex[ip])
            if iword < nwords:
                iword = iword + 1
                pp.coordex[ip], ierr = _read_real(_w(iword), pp.coordex[ip])
            else:
                pp.coordex[ip] = 0.0
            # Determine if any rmaxex has been defined.
            npgausstest = idiv(nwords - iword, 3)
            iresidual = nwords - iword - 3 * npgausstest
            if iresidual == 1:
                iword = iword + 1
                pp.rmaxex[ip], ierr = _read_real(_w(iword), pp.rmaxex[ip])
            else:
                pp.rmaxex[ip] = 0.0
                pp.pcoefexmax[ip] = 0.0
            igauss = 0
            while iword < nwords and ierr == 0 and igauss < potential_par.mpgauss:

                igauss = igauss + 1
                iword = iword + 1
                pp.rpgauss[igauss, ip], ierr = _read_real(_w(iword), pp.rpgauss[igauss, ip])
                iword = iword + 1
                pp.wpgauss[igauss, ip], ierr = _read_real(_w(iword), pp.wpgauss[igauss, ip])
                iword = iword + 1
                pp.hpgauss[igauss, ip], ierr = _read_real(_w(iword), pp.hpgauss[igauss, ip])
                if ierr != 0:
                    # FIX (FORTRAN_BUGS.md 1.4): the Fortran increments igauss here,
                    # double-counting the bad Gaussian; read_pcof decrements (Rust pcof.rs).
                    igauss = igauss - 1
            pp.npgauss[ip] = igauss
            #            write(6,*) ip,rminex(ip),pcoefex(ip)
            #
            # write the repulsive potential coefficients and Gaussian values if present
            #
            if pp.npgauss[ip] == 0:
                if pp.rmaxex[ip] <= 0.0:
                    units.write(nout, lw(pp.rminex[ip], pp.pcoefex[ip], pp.coordex[ip]))
                else:
                    units.write(nout, lw(pp.rminex[ip], pp.pcoefex[ip], pp.coordex[ip], pp.rmaxex[ip],
                                         pp.pcoefexmax[ip]))
            else:
                gvals = []
                for i in range(1, pp.npgauss[ip] + 1):
                    gvals += [pp.rpgauss[i, ip], pp.wpgauss[i, ip], pp.hpgauss[i, ip]]
                if pp.rmaxex[ip] <= 0.0:
                    formatline = fstr(fmt('("(3(1x,f10.5),",i1.1,"(1x,f9.4,1x,f9.4,1x,e11.4))")', pp.npgauss[ip]), 256)
                    units.write(nout, fmt(formatline, pp.rminex[ip], pp.pcoefex[ip], pp.coordex[ip], *gvals))
                else:
                    formatline = fstr(fmt('("(5(1x,f10.5),",i1.1,"(1x,f9.4,1x,f9.4,1x,e11.4))")', pp.npgauss[ip]), 256)
                    units.write(nout, fmt(formatline, pp.rminex[ip], pp.pcoefex[ip], pp.coordex[ip], pp.rmaxex[ip],
                                          pp.pcoefexmax[ip], *gvals))
            #
            # write the potential generation coefficients
            #
            units.write(nout, lw(*[pp.pcoefp[ic, ip] for ic in range(1, pp.ncoeffp + 1)]))
        # Write out any tb potential groups
        if tb.nthreeb > 0:
            for ithreeb in range(1, tb.nthreeb + 1):
                #
                # save the current values of the three-body potential in case they have been changed
                #
                icount = icount + 1
                iref = ireff(icount, symbols.nfiles)
                commands.sometext = ""
                commands.sometext = fstr(symbols.symbol_value[iref], commands.nstsize)
                _parse(11)
                nwords = commands.nwords
                # FIX (FORTRAN_BUGS.md 1.1): the Fortran re-parses the tb line into the
                # in-memory arrays, storing the angles and the angle width in radians
                # (read_pcof stores cosines and width 1.0) and the amplitudes without abs(),
                # which corrupts the three-body potential from the next iteration on.  As in
                # Rust (threebody.rs:write_pcof_threebody) the re-parsed values are kept in
                # locals, with the same arithmetic, so the written line is byte-identical.
                grp = [None] + [tb.threebgrp[j, ithreeb] for j in range(1, 4)]
                amp = tb.threebamp[ithreeb]
                awidth = tb.threebawidth[ithreeb]
                rwidth = tb.threebrwidth[ithreeb]
                rmax = [None, tb.threebrmax[1, ithreeb], tb.threebrmax[2, ithreeb]]
                samemol = [None] + [tb.threebsamemol[j, ithreeb] for j in range(1, 4)]
                nangs = int(tb.nthreebangs[ithreeb])
                angs = real_array(max(nangs, 9))
                amps = real_array(max(nangs, 9))
                for j in range(1, nangs + 1):
                    angs[j] = tb.threebangs[j, ithreeb]
                    amps[j] = tb.threebamps[j, ithreeb]
                if nwords > 12:
                    for i in range(1, 4):
                        # read(sometext(ncf(i):ncl(i)),'(a3)') threebgrp(i,ithreeb)
                        grp[i] = fstr(_w(i), 3)
                    i = 4
                    amp = read_internal(_w(i)).real(amp)
                    i = i + 1
                    awidth = read_internal(_w(i)).real(awidth)
                    awidth = piconv * awidth
                    i = i + 1
                    rwidth = read_internal(_w(i)).real(rwidth)
                    for j in range(1, 3):
                        i = i + 1
                        rmax[j] = read_internal(_w(i)).real(rmax[j])
                    for j in range(1, 4):
                        i = i + 1
                        samemol[j] = read_internal(_w(i)).int(samemol[j])
                    nangs = idiv(nwords - i, 2)
                    if nangs > 9:
                        nangs = 9
                    if nangs > tb.nthreebangsmax:
                        # (the Fortran takes nthreebangs(nthreeb), not (ithreeb); never reached
                        # since read_pcof already took the maximum of the same values)
                        tb.nthreebangsmax = int(tb.nthreebangs[tb.nthreeb])
                    tb.threebangs = reallocation_routines.reallocate2d_r(tb.threebangs, tb.nthreebangsmax, tb.nthreeb)
                    tb.threebamps = reallocation_routines.reallocate2d_r(tb.threebamps, tb.nthreebangsmax, tb.nthreeb)
                    for j in range(1, nangs + 1):
                        i = i + 1
                        angs[j] = read_internal(_w(i)).real(angs[j])
                        angs[j] = piconv * angs[j]
                        i = i + 1
                        amps[j] = read_internal(_w(i)).real(amps[j])
                formatline = fstr(fmt('("(a,3(1x,a3),1x,f10.2,2(1x,f7.2),2(1x,f10.5),3(1x,i1),",i1.1,"(1x,f7.2))")',
                                      2 * nangs), 256)
                units.write(6, fmt("(a)", formatline))
                vals = ["tb"] + [Ch(grp[j], 3) for j in range(1, 4)] + [amp, awidth / piconv, rwidth]
                vals += [rmax[j] for j in range(1, 3)] + [samemol[j] for j in range(1, 4)]
                for j in range(1, nangs + 1):
                    vals += [angs[j] / piconv, amps[j]]
                units.write(nout, fmt(_star_repeat(formatline), *vals))
        if coordconstraints_mod.ncoordconstraints > 0:
            for icoord in range(1, coordconstraints_mod.ncoordconstraints + 1):
                cc = coordconstraints_mod.coordconstraints[icoord]
                formatline = fstr(fmt('("(a,1x,a,1x,i3,2(1x,f12.3,1x,f8.3),",i2.2,"(1x,a,2(1x,f8.3)))")', cc.ns), 256)
                units.write(6, fmt("(a,1x,i5)", formatline, cc.ns))
                vals = ["coordcon", Ch(cc.atc, 3), cc.nc, cc.cctemp, cc.ccdistort, cc.ccgaussweight, cc.ccgausswidth]
                for i in range(1, cc.ns + 1):
                    vals += [Ch(cc.ats[i], 3), cc.ds[i], cc.dsmax[i]]
                shortline = fstr(fmt(_star_repeat(formatline), *vals), 256)
                fullline = shortline
                if cc.hastb:
                    # (this format and the next end with one ')' too many; gfortran, and
                    # fformat, stop at the matching parenthesis)
                    formatline = fstr(fmt('("(1x,a,2(1x,f12.3,1x,f8.3),",i2.2,"(1x,f8.3)))")', cc.ntbangles), 256)
                    units.write(6, fmt("(a,1x,i5)", formatline, cc.ntbangles))
                    vals = ["angles", cc.tbweight, cc.tbdistort, cc.tbgaussweight, cc.tbgausswidth]
                    vals += [cc.tbangles[i] for i in range(1, cc.ntbangles + 1)]
                    shortline = fstr(fmt(_star_repeat(formatline), *vals), 256)
                    fullline = fstr(fullline + pad(shortline, 256), 256)
                if cc.hastbintra:
                    formatline = fstr(fmt('("(1x,a,2(1x,f12.3,1x,f8.3),",i2.2,"(1x,a3,1x,f8.3)))")',
                                          cc.ntbintraangles), 256)
                    units.write(6, fmt("(a,i12)", formatline, cc.ntbintraangles))
                    vals = ["anglesintra", cc.tbintraweight, cc.tbintradistort, cc.tbintragaussweight,
                            cc.tbintragausswidth]
                    for i in range(1, cc.ntbintraangles + 1):
                        vals += [Ch(cc.tbintraats[i], 3), cc.tbintraangles[i]]
                    shortline = fstr(fmt(_star_repeat(formatline), *vals), 256)
                    fullline = fstr(fullline + pad(shortline, 256), 256)
                units.write(nout, fmt("(a)", Ch(fullline, 256)))
        units.close(nout)
        #
        # now write the specified .ato file
        # but only write the .ato file if it contains some molecules!
        #
        if clustercom.nmol > 0:
            ato_read_and_save.moloutc(symbols.fnameato)
    else:
        units.close(nout)

    return nout


# --------------------------------------------------------------------------- read_pcof
_F103_PCOF = "(/a,1x,'No .pcof file specified! - routine will exit'/)"
_F104 = "(/'.pcof file \"',a,'\"')"
_F101_PCOF = "(/a,1x,'Problems with specified .pcof file: ',a/' - will use defaults or existing values')"


def _list_read_reals(unit, arr, col, n):
    """``read(unit,*,iostat=ierr) (arr(ic,col),ic=1,n)``: returns ierr; the items read
    before an error keep their new values."""
    try:
        r = units.list_read(unit)
        for ic in range(1, n + 1):
            arr[ic, col] = r.real(arr[ic, col])
        r.end()
    except FortranEOF:
        return -1
    except FortranIOError as e:
        return e.iostat
    return 0


@routine
def read_pcof(nin, nout, icount, ifilecount, ierrout):
    """read_and_save_files.f:read_pcof -- sets up the EPSR pcof values and symbols.

    This should always be done AFTER the .ato file has been read.  Reads the .pcof file
    named by symbol icount (the parameters before 'q', the pair potentials, and any
    three-body (tb) and coordination-constraint (coordcon) lines), or sets up defaults.
    Returns (nout, icount, ierrout).
    """
    from ..common import coordconstraints_mod, threebody_mod  # Fortran module data (packages H, G)
    from . import coordconstraint_mod, file_and_directory_commands, reallocation_routines  # packages H, A

    ireff = file_and_directory_commands.ireff
    tb = threebody_mod
    pp = potential_par
    rr = reallocation_routines
    # integer icount,ifilecount,ifound,ierrout,ierrread
    # integer irefireset,irefiinit
    # character*256 formatline,text
    # character*20 ext,command
    # character*1 iaccept
    # integer ilenext,ilenlastdot,npgausstest,npcoeffread,ilendash
    # real sumamp
    # NB: ireset, iinit, roverlap, rminfac are implicit locals of read_pcof; roverlap and
    # rminfac are uninitialised when their symbol cannot be read (0 is used), and so are
    # the lengths returned by concatenate (outputs only).
    roverlap = f32(0.0)
    rminfac = f32(0.0)
    ilenvalue = 0
    ilensymbol = 0
    ilensymbol_value = 0
    iaccept = ""
    ierrread = 0

    pi = PI4ATAN  # 4.0*atan(1.0)
    piconv = pi / r4(180.0)
    #
    # ensure the input stream number is different from the stream number where data are
    # currently being read
    #
    nin1 = nin + 2
    # Initialise the three-body potential
    tb.nthreeb = 0
    tb.nthreebangsmax = 0
    # Initialise the coordination constraints
    coordconstraints_mod.ncoordconstraints = 0
    #
    # set up the integer corresponding to fnamepcof - it is assumed that symbol(icount)
    # has already been checked before calling this routine
    #
    ifound = ireff(icount, ifilecount)
    #      write(6,*) 'read_pcof> ',icount,ifilecount,ifound
    #
    # get the symbol number for ireset and iinit, and read ireset. If not found
    # then assume ireset = 1
    #
    if iresetnumbers.iresetref > 0 and iresetnumbers.iresetref <= symbols.nsymbol:
        irefireset = ireff(iresetnumbers.iresetref, ifilecount)
        irefiinit = ireff(iresetnumbers.iinitref, ifilecount)
        ireset = read_internal(symbols.symbol_value[irefireset]).int(0)
        iinit = read_internal(symbols.symbol_value[irefiinit]).int(0)
    else:
        ireset = 1
        irefireset = 0
        iinit = 1
        irefiinit = 0
    #
    # setup the expected extension for .pcof files
    #
    ext = ".pcof"
    ilenext = 5
    #
    # open the .pcof file and read in the relevant parameters from there
    #
    if symbols.symbol_value[ifound] == "search" or symbols.symbol_value[ifound] == "<undefined>":
        #
        # search for a suitable pcof file
        #
        commands.nwords = 0
        while commands.nwords <= 0:
            iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)
            _parse(5)
            if commands.nwords <= 0:
                units.write(nout, fmt(_F103_PCOF, substr(symbols.prompt, 1, symbols.ilenprompt)))
                if symbols.symbol_save != "":
                    symbols.symbol_value[ifound] = symbols.symbol_save
                else:
                    symbols.symbol_value[ifound] = symbols.symbol_default[ifound]
                ierrout = 0
                return nout, icount, ierrout
        symbols.fname_pcof, symbols.ilenfname_pcof = _concatenate(1, commands.nwords, symbols.fname_pcof,
                                                                  symbols.ilenfname_pcof)
    else:
        commands.sometext = fstr(symbols.symbol_value[ifound], commands.nstsize)
        _parse(5)
        symbols.fname_pcof, symbols.ilenfname_pcof = _concatenate(1, commands.nwords, symbols.fname_pcof,
                                                                  symbols.ilenfname_pcof)

    # ensure this file as the correct extension

    ilenlastdot = auxiliary_routines.indexlastdot(symbols.fname_pcof) - 1
    if ilenlastdot > 0:
        symbols.fname_pcof = fstr(substr(symbols.fname_pcof, 1, ilenlastdot) + pad(ext, 20), 256)
        symbols.ilenfname_pcof = ilenlastdot + ilenext
    else:
        symbols.fname_pcof = fstr(substr(symbols.fname_pcof, 1, symbols.ilenfname_pcof) + pad(ext, 20), 256)
        symbols.ilenfname_pcof = symbols.ilenfname_pcof + ilenext
    #
    # convert this to a full path name if necessary
    #
    symbols.fnamepcof, symbols.ilenfnamepcof = file_and_directory_commands.get_full_filename(
        symbols.fname_pcof, symbols.ilenfname_pcof, symbols.fnamepcof, symbols.ilenfnamepcof)
    #
    # open the file if possible and read the commands
    #
    try:
        units.open(nin1, symbols.fnamepcof, status="old")
        ierr = 0
    except FortranIOError as e:
        ierr = e.iostat
    if ierr == 0:
        #
        # reset the number of symbols to the number of data symbols
        #
        symbols.nsymbol = symbols.nsymboldata
        #
        # get the symbol number for 'ereq' - this is now read in from pcof file, if present.
        #
        icount = iresetnumbers.iereqfirst
        #
        # sum of all the read errors
        #
        ierrread = 0
        #
        # set up a store for the symbol number of ncoeffp
        #
        incoeffpref = 0  # noqa: F841
        #
        # initialise command
        #
        command = ""
        ierr = 0
        while command != "q" and ierr == 0:
            icount = icount + 1
            #
            # ignore lines with nothing in them
            #
            commands.nwords = 0
            while commands.nwords == 0 and ierr == 0:
                # read(nin1,100,iostat=ierr) sometext -- at end of file sometext keeps its value
                try:
                    commands.sometext = fstr(units.read_line(nin1), commands.nstsize)
                    ierr = 0
                except FortranEOF:
                    ierr = -1
                except FortranIOError as e:
                    ierr = e.iostat
                ierrread = ierrread + abs(ierr)
                _parse(11)
            command = fstr(_w(1), 20)
            #
            # reflect the input text into the output stream
            #
            #           if(nwords.gt.0) write(nout,*) (sometext(ncf(i):ncl(i)),' ',i=1,nwords)

            #
            # don't waste time searching for a symbol which is not there for sure
            #
            if commands.nwords > 1:
                #
                # otherwise check the symbol number and save the new value
                #
                icount = file_and_directory_commands.search_for_symbol(command, icount)
                if icount != 0 and icount <= symbols.nsymbol:
                    iref = ireff(icount, ifilecount)
                    #      write(6,*) 'read_pcof> ',icount,ifilecount,iref
                    symbols.symbol_value[iref], ilenvalue = _concatenate(2, commands.nwords,
                                                                         symbols.symbol_value[iref], ilenvalue)
                    #      write(6,*) command,' ',iref,icount,' ',symbol_value(iref)
                #
                # save the symbol number corresponding to ncoeffp
                #
                if command == "ncoeffp":
                    iresetnumbers.incoeffp = icount
                if command == "roverlap":
                    iresetnumbers.iroverlap = icount
                if command == "rminfac":
                    iresetnumbers.irminfac = icount

        #
        # Set up the array size for energy pair values
        #
        #         command='menerg'
        #         call search_for_symbol(command,icount)
        #         iref=ireff(icount,ifilecount)
        #         read(symbol_value(iref),*) menerg

        #         command='menergadd'
        #         call search_for_symbol(command,icount)
        #         iref=ireff(icount,ifilecount)
        #         read(symbol_value(iref),*) menergadd
        #
        # get the number of empirical potential coefficients
        #
        iref = ireff(iresetnumbers.incoeffp, ifilecount)
        pp.ncoeffp, ierr = _read_int(symbols.symbol_value[iref], pp.ncoeffp)
        #         write(6,*) ncoeffp
        # Make sure this number does not exceed the allowed dimension
        if pp.ncoeffp > potential_par.mcoeff:
            pp.ncoeffp = potential_par.mcoeff
        # Write the new value to the symbol table
        commands.sometext = fstr(lw(pp.ncoeffp), commands.nstsize)
        _parse(5)
        symbols.symbol_value[iref], ilensymbol_value = _concatenate(1, 1, symbols.symbol_value[iref], ilensymbol_value)
        ierrread = ierrread + abs(ierr)
        iref = ireff(iresetnumbers.iroverlap, ifilecount)
        roverlap, ierr = _read_real(symbols.symbol_value[iref], roverlap)
        ierrread = ierrread + abs(ierr)
        iref = ireff(iresetnumbers.irminfac, ifilecount)
        rminfac, ierr = _read_real(symbols.symbol_value[iref], rminfac)
        ierrread = ierrread + abs(ierr)
        #
        # now get from this file, or create, the rminex values
        #
        pp.npairspcof = 0
        try:
            r = units.list_read(nin1)
            pp.npairspcof = r.int(pp.npairspcof)
            ierr = 0
        except FortranEOF:
            ierr = -1
        except FortranIOError as e:
            ierr = e.iostat
        #     write(6,*) npairspcof,ncoeffp
        ierrread = ierrread + abs(ierr)

        units.write(nout, fmt(_F104, substr(symbols.fname_pcof, 1, symbols.ilenfname_pcof)))

        for irefp in range(1, pp.npairspcof + 1):
            pp.efacp[irefp] = 1.0  # Set the factor which controls whether this potential will be refined
            try:
                commands.sometext = fstr(units.read_line(nin1), commands.nstsize)
                ierr = 0
            except FortranEOF:
                ierr = -1
            except FortranIOError as e:
                ierr = e.iostat
            ierrread = ierrread + abs(ierr)
            _parse(11)
            pp.tlab1[irefp] = fstr(_w(1), 3)
            pp.tlab2[irefp] = fstr(_w(2), 3)
            pp.coulfac[irefp] = 1.0
            if commands.nwords > 3 and ireset != 1:
                pp.coulfac[irefp], ierr = _read_real(_w(4), pp.coulfac[irefp])
            #      write(6,*) irefp,npairspcof,nfiles,nwords,nsymbol
            #
            # generate the new symbol name
            #
            symbols.nsymbol = symbols.nsymbol + 1
            nsymbol = symbols.nsymbol
            iref = ireff(nsymbol, symbols.nfiles)
            symbols.symbol[nsymbol] = fstr(_w(1) + "-" + _w(2), 11)
            symbols.iflag[nsymbol] = 0
            symbols.symbol_default[nsymbol] = ""
            symbols.symbol_type[nsymbol] = "r"
            try:
                commands.sometext = fstr(units.read_line(nin1), commands.nstsize)
                ierr = 0
            except FortranEOF:
                ierr = -1
            except FortranIOError as e:
                ierr = e.iostat
            _parse(11)
            nwords = commands.nwords

            # No. of Gaussians to read. We need three parameters for each Gaussian, its position
            # its sigma and its height or depth. Any error in reading these and the data are ignored.

            iword = 1
            # 2021-02-28 We  have added the expected coordination number, coordex(irefp) after pcoefex here
            pp.rminex[irefp], ierr = _read_real(_w(iword), pp.rminex[irefp])
            iword = iword + 1
            pp.pcoefex[irefp], ierr = _read_real(_w(iword), pp.pcoefex[irefp])
            # Make sure if the specified radius is zero or ireset = 2, the corresponding coefficient is zero
            if pp.rminex[irefp] == 0.0 or ireset == 2:
                pp.pcoefex[irefp] = 0.0
            if iword < nwords:
                iword = iword + 1
                pp.coordex[irefp], ierr = _read_real(_w(iword), pp.coordex[irefp])
            else:
                pp.coordex[irefp] = 0.0
            pp.rmaxex[irefp] = 0.0
            pp.pcoefexmax[irefp] = 0.0
            # Check whether rmaxex has been specified. If it has there will be 2+npgauss*3 words left.
            # If not there will be only npgauss*3 words left.
            npgausstest = idiv(nwords - iword, 3)
            iwordresidual = nwords - iword - npgausstest * 3
            if iwordresidual == 2:
                iword = iword + 1
                pp.rmaxex[irefp], ierr = _read_real(_w(iword), pp.rmaxex[irefp])
                iword = iword + 1
                pp.pcoefexmax[irefp], ierr = _read_real(_w(iword), pp.pcoefexmax[irefp])
                if pp.rmaxex[irefp] == 0.0:
                    pp.pcoefexmax[irefp] = 0.0
            # Read Gaussians if present
            igauss = 0
            while iword < nwords and ierr == 0 and igauss < potential_par.mpgauss:

                igauss = igauss + 1
                iword1 = iword + 1
                iword2 = iword + 2
                iword3 = iword + 3
                pp.rpgauss[igauss, irefp], ierr = _read_real(_w(iword1), pp.rpgauss[igauss, irefp])
                pp.wpgauss[igauss, irefp], ierr = _read_real(_w(iword2), pp.wpgauss[igauss, irefp])
                pp.hpgauss[igauss, irefp], ierr = _read_real(_w(iword3), pp.hpgauss[igauss, irefp])

                if ierr != 0:
                    igauss = igauss - 1

                iword = iword3

            pp.npgauss[irefp] = igauss
            #
            # remove absurdly small values of pcoefex
            #
            if abs(pp.pcoefex[irefp]) < r4(1.0e-10):
                pp.pcoefex[irefp] = 0.0
            ierrread = ierrread + abs(ierr)
            #
            # check that the minimum allowed radius for this pair is correct
            # However for Q or q atoms the minimum radius can be zero
            # NB: sometext is the numbers line here, not the labels line, so the q/Q
            # exemption never applies (reproduced; Rust pcof.rs:258).
            #
            c1 = substr(commands.sometext, commands.ncf[1], commands.ncf[1])
            c2 = substr(commands.sometext, commands.ncf[2], commands.ncf[2])
            if (abs(pp.rminex[irefp]) < roverlap and c1 != "q" and c1 != "Q" and c2 != "q" and c2 != "Q"):
                pp.rminex[irefp] = sign(roverlap, pp.rminex[irefp])

            #
            # generate the appropriate symbol value
            #
            commands.sometext = ""
            if pp.npgauss[irefp] == 0:
                if pp.rmaxex[irefp] <= 0.0:
                    commands.sometext = fstr(fmt("(1x,f10.5,a,f10.5)", pp.rminex[irefp], " ", pp.coordex[irefp]),
                                             commands.nstsize)
                else:
                    commands.sometext = fstr(fmt("((1x,f10.5,a,f10.5,1x,f10.5))", pp.rminex[irefp], " ",
                                                 pp.coordex[irefp], pp.rmaxex[irefp]), commands.nstsize)
            else:
                gvals = []
                for igauss in range(1, pp.npgauss[irefp] + 1):
                    gvals += [pp.rpgauss[igauss, irefp], pp.wpgauss[igauss, irefp], pp.hpgauss[igauss, irefp]]
                if pp.rmaxex[irefp] <= 0.0:
                    formatline = fstr(fmt('("((1x,f10.5,a,f10.5),",i1.1,"(1x,f9.4,1x,f9.4,1x,e11.4))")',
                                          pp.npgauss[irefp]), 256)
                    commands.sometext = fstr(fmt(formatline, pp.rminex[irefp], " ", pp.coordex[irefp], *gvals),
                                             commands.nstsize)
                else:
                    formatline = fstr(fmt('("((1x,f10.5,a,f10.5,1x,f10.5),",i1.1,"(1x,f9.4,1x,f9.4,1x,e11.4))")',
                                          pp.npgauss[irefp]), 256)
                    commands.sometext = fstr(fmt(formatline, pp.rminex[irefp], " ", pp.coordex[irefp],
                                                 pp.rmaxex[irefp], *gvals), commands.nstsize)
            ierr = 0  # (the internal writes have iostat=ierr)
            ierrread = ierrread + abs(ierr)
            _parse(11)
            symbols.symbol_value[iref], ilensymbol = _concatenate(1, commands.nwords, symbols.symbol_value[iref],
                                                                  ilensymbol)
            symbols.symbol_text[nsymbol] = "Minimum (and maximum) separation and P,W,H of Gaussians for this atom pair"
            units.write(nout, lw(substr(symbols.symbol[nsymbol], 1, 10), substr(symbols.symbol_value[iref], 1, ilensymbol)))
            #
            # read the potential coefficients or reset them, depending on the value of ireset
            #
            # The first one is not used by the program and will be used to signal a factor for the remainder
            ierr = _list_read_reals(nin1, pp.pcoefp, irefp, pp.ncoeffp)
            ierrread = ierrread + abs(ierr)
            if ireset != 0:
                # Reset them to zero if requested. The first one is not used by the program and will be used to
                # signal a factor for the remainder
                if ireset < 2:
                    pp.pcoefp[1, irefp] = 1.0
                for ic in range(2, pp.ncoeffp + 1):
                    pp.pcoefp[ic, irefp] = 0.0
            pp.efacp[irefp] = pp.pcoefp[1, irefp]
            # Check whether any of the additional terms requires the EP to be held fixed: to achieve this the
            # position or width needs to be set < 0
            for ig in range(1, pp.npgauss[irefp] + 1):
                if pp.rpgauss[ig, irefp] < 0.0 or pp.wpgauss[ig, irefp] < 0.0:
                    pp.efacp[irefp] = 0.0
            pp.pcoefp[1, irefp] = pp.efacp[irefp]
        #
        # try to match each pair of atoms in the current .ato file with one of the
        # sets of coefficients read in. If not then create a new set of coefficients.
        #
        for ip in range(1, molecom.npairs + 1):
            it = molecom.n1extr[ip]
            jt = molecom.n2extr[ip]
            irefpfound = 0
            irefp = 0
            while irefp < pp.npairspcof and irefpfound == 0:
                irefp = irefp + 1
                if ((pp.tlab1[irefp] == molecom.type[it] and pp.tlab2[irefp] == molecom.type[jt])
                        or (pp.tlab1[irefp] == molecom.type[jt] and pp.tlab2[irefp] == molecom.type[it])):
                    irefpfound = irefp
            #
            # if the correct pair was found, then save it
            #
            if irefpfound > 0:
                pp.irefpcof[ip] = irefpfound
                units.write(6, fmt("(a,2(1x,i5),2(1x,a))", "read_pcof> Pair found: ", ip, pp.irefpcof[ip],
                                   Ch(pp.tlab1[irefpfound], 3), Ch(pp.tlab2[irefpfound], 3)))
            else:
                #
                # otherwise create a new set of coefficients
                #
                pp.npairspcof = pp.npairspcof + 1
                n = pp.npairspcof
                pp.irefpcof[ip] = n
                pp.tlab1[n] = molecom.type[it]
                pp.tlab2[n] = molecom.type[jt]
                pp.pcoefp[1, n] = 1.0
                for ic in range(2, pp.ncoeffp + 1):
                    pp.pcoefp[ic, n] = 0.0

                # Do not set the minimum distance for qtype atoms

                ti1 = substr(molecom.type[it], 1, 1)
                tj1 = substr(molecom.type[jt], 1, 1)
                if ti1 != "q" and ti1 != "Q" and tj1 != "q" and tj1 != "Q":

                    pp.rminex[n] = rminfac * molecom.dcor[ip]
                    if pp.rminex[n] < roverlap:
                        pp.rminex[n] = roverlap

                else:

                    pp.rminex[n] = 0.0

                pp.pcoefex[n] = 0.0
                pp.rmaxex[n] = 0.0
                pp.pcoefexmax[n] = 0.0
                pp.npgauss[n] = 0
                # FIX (FORTRAN_BUGS.md 1.3): the Fortran leaves coulfac, coordex and efacp of
                # the new pair unset (COMMON data, i.e. 0 or stale), so the pair had no Coulomb
                # term and its EP was not refined.  Give it the defaults a pair read from file
                # gets (Rust pcof.rs).
                pp.coulfac[n] = 1.0
                pp.efacp[n] = 1.0
                pp.coordex[n] = 0.0
                #
                # generate the new symbol name
                #
                commands.sometext = fstr(lw(Ch(molecom.type[it], 3), " ", Ch(molecom.type[jt], 3)), commands.nstsize)
                _parse(11)
                symbols.nsymbol = symbols.nsymbol + 1
                nsymbol = symbols.nsymbol
                iref = ireff(nsymbol, symbols.nfiles)
                symbols.symbol[nsymbol] = fstr(_w(1) + "-" + _w(2), 11)
                symbols.iflag[nsymbol] = 0
                symbols.symbol_default[nsymbol] = ""
                symbols.symbol_type[nsymbol] = "r"
                #
                # generate the appropriate symbol value
                #
                commands.sometext = ""
                commands.sometext = fstr(lw(pp.rminex[n]), commands.nstsize)
                ierr = 0
                ierrread = ierrread + abs(ierr)
                _parse(11)
                symbols.symbol_value[iref], ilensymbol = _concatenate(1, 1, symbols.symbol_value[iref], ilensymbol)
                symbols.symbol_text[nsymbol] = "Minimum separation and P,W,H of Gaussians for this atom pair"
                units.write(nout, lw(substr(symbols.symbol[nsymbol], 1, 10),
                                     substr(symbols.symbol_value[iref], 1, ilensymbol)))

        # Get any three body terms or coordination constraints that have been defined
        while ierr == 0:
            try:
                commands.sometext = fstr(units.read_line(nin1), commands.nstsize)
                ierr = 0
            except FortranEOF:
                ierr = -1
            except FortranIOError as e:
                ierr = e.iostat
            #
            # parse this text into words
            #
            if ierr == 0:
                _parse(11)
                nwords = commands.nwords
                if nwords > 0:
                    #                  write(6,'(a,i3,1x,a)') 'minitc> ',nwords,sometext(ncf(1):ncl(nwords))
                    # Trap the program here in case a three body potential group has been read in
                    if _w(1) == "tb" and nwords > 13:
                        # Three body group to be read in
                        tb.nthreeb = tb.nthreeb + 1
                        nthreeb = tb.nthreeb
                        tb.threebgrp = rr.reallocate2d_c(tb.threebgrp, 3, 3, nthreeb)
                        tb.threebamp = rr.reallocate1d_r(tb.threebamp, nthreeb)
                        tb.threebrmax = rr.reallocate2d_r(tb.threebrmax, 2, nthreeb)
                        tb.threebawidth = rr.reallocate1d_r(tb.threebawidth, nthreeb)
                        tb.threebrwidth = rr.reallocate1d_r(tb.threebrwidth, nthreeb)
                        tb.nthreebangs = rr.reallocate1d_i(tb.nthreebangs, nthreeb)
                        tb.threebsamemol = rr.reallocate2d_i(tb.threebsamemol, 3, nthreeb)
                        iref = 1
                        for i in range(1, 4):
                            iref = iref + 1
                            # read(sometext(ncf(iref):ncl(iref)),'(a3)') threebgrp(i,nthreeb)
                            tb.threebgrp[i, nthreeb] = fstr(_w(iref), 3)
                        # Ensure third atom is same as first (for time being), otherwise this three body
                        # group is ignored
                        #                     if(threebgrp(1,nthreeb).eq.threebgrp(3,nthreeb)) then
                        iref = iref + 1
                        tb.threebamp[nthreeb] = read_internal(_w(iref)).real(tb.threebamp[nthreeb])
                        iref = iref + 1
                        tb.threebawidth[nthreeb] = read_internal(_w(iref)).real(tb.threebawidth[nthreeb])
                        tb.threebawidth[nthreeb] = 1.0
                        iref = iref + 1
                        tb.threebrwidth[nthreeb] = read_internal(_w(iref)).real(tb.threebrwidth[nthreeb])
                        for i in range(1, 3):
                            iref = iref + 1
                            tb.threebrmax[i, nthreeb] = read_internal(_w(iref)).real(tb.threebrmax[i, nthreeb])
                        for i in range(1, 4):
                            iref = iref + 1
                            tb.threebsamemol[i, nthreeb] = read_internal(_w(iref)).int(tb.threebsamemol[i, nthreeb])
                        tb.nthreebangs[nthreeb] = idiv(nwords - iref, 2)
                        units.write(6, lw("read_pcof> No. of tb angles: ", nthreeb, tb.nthreebangs[nthreeb]))
                        if tb.nthreebangs[nthreeb] > 9:
                            tb.nthreebangs[nthreeb] = 9
                        if tb.nthreebangs[nthreeb] > tb.nthreebangsmax:
                            tb.nthreebangsmax = int(tb.nthreebangs[nthreeb])
                        tb.threebangs = rr.reallocate2d_r(tb.threebangs, tb.nthreebangsmax, nthreeb)
                        tb.threebamps = rr.reallocate2d_r(tb.threebamps, tb.nthreebangsmax, nthreeb)
                        sumamp = f32(0.0)
                        for i in range(1, tb.nthreebangs[nthreeb] + 1):
                            iref = iref + 1
                            tb.threebangs[i, nthreeb] = read_internal(_w(iref)).real(tb.threebangs[i, nthreeb])
                            tb.threebangs[i, nthreeb] = cosf(piconv * tb.threebangs[i, nthreeb])
                            iref = iref + 1
                            tb.threebamps[i, nthreeb] = read_internal(_w(iref)).real(tb.threebamps[i, nthreeb])
                            tb.threebamps[i, nthreeb] = abs(tb.threebamps[i, nthreeb])
                            sumamp = sumamp + tb.threebamps[i, nthreeb]
                        # Normalise individual terms to total sum of amplitudes
                        if sumamp > 0.0:
                            pass
                            #                            do i=1,nthreebangs(nthreeb)
                            # 2016-02-09 Cut this out for compatibility with earlier simulations
                            #                               threebamps(i,nthreeb)=threebamps(i,nthreeb)/sumamp
                            #                            end do
                        else:
                            tb.threebamp[nthreeb] = 0.0
                        #                         threebawidth(nthreeb)=threebawidth(nthreeb)*threebangs(1,nthreeb)/pi
                        #                         threebangs(1,nthreeb)=cos(threebangs(1,nthreeb))+threebawidth(nthreeb)*log(2.0)
                        #                     else
                        #                         nthreeb=nthreeb-1
                        #                     end if
                        #
                        # generate the new symbol name
                        #
                        symbols.nsymbol = symbols.nsymbol + 1
                        nsymbol = symbols.nsymbol
                        iref = ireff(nsymbol, symbols.nfiles)
                        symbols.symbol[nsymbol] = fstr(_w(1), 11)
                        symbols.iflag[nsymbol] = 0
                        symbols.symbol_default[nsymbol] = ""
                        symbols.symbol_type[nsymbol] = "c"
                        symbols.symbol_value[iref], ilensymbol = _concatenate(2, commands.nwords,
                                                                              symbols.symbol_value[iref], ilensymbol)
                        symbols.symbol_text[nsymbol] = ("Three atom labels, amplitude, width, angle width, two dist.,"
                                                        " and angles")
                        units.write(nout, lw(substr(symbols.symbol[nsymbol], 1, 10),
                                             substr(symbols.symbol_value[iref], 1, ilensymbol)))
                    elif _w(1) == "coordcon" and nwords > 8:
                        coordconstraints_mod.ncoordconstraints = coordconstraints_mod.ncoordconstraints + 1
                        ncc = coordconstraints_mod.ncoordconstraints
                        coordconstraints_mod.coordconstraints = coordconstraint_mod.reallocate1d_cc(
                            coordconstraints_mod.coordconstraints, ncc)
                        text = fstr(substr(commands.sometext, commands.ncf[2], commands.ncl[nwords]), 256)
                        units.write(6, lw("read_pcof> No. coord. constraints: ", ncc, " ", text))
                        # text is character*256: the callee sees len(sometext) = 256
                        coordconstraints_mod.coordconstraints[ncc] = coordconstraint_mod.init_constraint(
                            coordconstraints_mod.coordconstraints[ncc], text, len_sometext=256)
                        #                  call find_constraints(coordconstraints(ncoordconstraints),nmol,natom,xadd,yadd,zadd,xa,ya,za,mref
                        #     1,itype,ntype,type)
        #
        # if any error occurred during input, then ireset is set to 1
        #
        if ierrread != 0:
            ireset = 1
            iinit = 1
            ierrout = ierrread
        else:
            ierrout = 0
        units.close(nin1)
    else:
        units.write(6, fmt(_F101_PCOF, substr(symbols.prompt, 1, symbols.ilenprompt),
                           substr(symbols.fname_pcof, 1, symbols.ilenfname_pcof)))
        units.write(6, lw(pp.npairspcof))
        #
        # get the default parameters needed to set up the potential coefficients
        #
        iref = ireff(iresetnumbers.incoeffp, ifilecount)
        pp.ncoeffp, ierr = _read_int(symbols.symbol_value[iref], pp.ncoeffp)

        # Make sure this number does not exceed the allowed dimension

        if pp.ncoeffp > potential_par.mcoeff:
            pp.ncoeffp = potential_par.mcoeff
        iref = ireff(iresetnumbers.iroverlap, ifilecount)
        roverlap, ierr = _read_real(symbols.symbol_value[iref], roverlap)
        iref = ireff(iresetnumbers.irminfac, ifilecount)
        rminfac, ierr = _read_real(symbols.symbol_value[iref], rminfac)
        #
        # the value of npairspcof will be used to signify whether a set of coefficients
        # has already been defined. If so we leave well alone, but keep the new name
        # for the .pcof file so this will be written when the save is done.
        #
        if pp.npairspcof <= 0:
            #
            # set ireset to 1 to ensure the .pcof file and other parameters are initialised correctly
            # on saving the input files and running the simulation
            #
            ireset = 1
            iinit = 1
            ierrread = 1
        else:
            ierrout = 0
            ierrread = 0
            iref = ireff(symbols.nsymboldata, symbols.nfiles)
            for ip in range(1, pp.npairspcof + 1):
                iref = iref + 1
                units.write(6, lw(pp.rminex[ip], " ", substr(symbols.symbol_value[iref], 1, 20)))

    #
    # if there was an error on reading the .pcof file we generate default values
    #
    if ierrread != 0 or ireset == 1:
        #
        # set the number of symbols to the current maximum
        #
        symbols.nsymbol = symbols.nsymboldata
        #
        # set up minimum overlap distances for each atom pair. These are based on the values of
        # the L-J sigma for each atom, but can be altered as required.
        #
        pp.npairspcof = molecom.npairs
        units.write(nout, fmt(_F104, substr(symbols.fname_pcof, 1, symbols.ilenfname_pcof)))
        for ip in range(1, molecom.npairs + 1):
            it = molecom.n1extr[ip]
            jt = molecom.n2extr[ip]

            # Do not set the minimum distance for qtype atoms

            ti1 = substr(molecom.type[it], 1, 1)
            tj1 = substr(molecom.type[jt], 1, 1)
            if ti1 != "q" and ti1 != "Q" and tj1 != "q" and tj1 != "Q":

                pp.rminex[ip] = rminfac * molecom.dcor[ip]
                if pp.rminex[ip] < roverlap:
                    pp.rminex[ip] = roverlap

            else:

                pp.rminex[ip] = 0.0

            pp.pcoefex[ip] = 0.0
            pp.coordex[ip] = 0.0
            pp.rmaxex[ip] = 0.0
            pp.pcoefexmax[ip] = 0.0
            pp.irefpcof[ip] = ip
            pp.coulfac[pp.irefpcof[ip]] = 1.0
            pp.tlab1[ip] = molecom.type[it]
            pp.tlab2[ip] = molecom.type[jt]
            pp.npgauss[ip] = 0
            #
            # generate the new symbol name
            #
            commands.sometext = fstr(lw(Ch(pp.tlab1[ip], 3), " ", Ch(pp.tlab2[ip], 3)), commands.nstsize)
            _parse(11)
            symbols.nsymbol = symbols.nsymbol + 1
            nsymbol = symbols.nsymbol
            iref = ireff(nsymbol, symbols.nfiles)
            symbols.symbol[nsymbol] = fstr(_w(1) + "-" + _w(2), 11)
            symbols.iflag[nsymbol] = 0
            symbols.symbol_default[nsymbol] = ""
            #
            # generate the appropriate symbol value
            #
            commands.sometext = ""
            commands.sometext = fstr(lw(pp.rminex[ip], " ", pp.coordex[ip]), commands.nstsize)
            _parse(11)
            symbols.symbol_value[iref], ilensymbol = _concatenate(1, 1, symbols.symbol_value[iref], ilensymbol)
            symbols.symbol_text[nsymbol] = ("Minimum (coordination number and maximum) separation and P,W,H of"
                                            " Gaussians for this atom pair")
            units.write(nout, lw(substr(symbols.symbol[nsymbol], 1, 10), substr(symbols.symbol_value[iref], 1, ilensymbol)))
            #
            # Initialize the potential coefficients
            #
            pp.pcoefp[1, ip] = 1.0
            for ic in range(2, pp.ncoeffp + 1):
                pp.pcoefp[ic, ip] = 0.0

        ierrout = 0

    symbols.symbol_value[ifound] = fstr(substr(symbols.fname_pcof, 1, symbols.ilenfname_pcof), 256)
    #
    # write the new values of ireset and iinit in case they have been changed
    #
    if iresetnumbers.iresetref != 0:
        commands.sometext = fstr(lw(ireset), commands.nstsize)
        _parse(5)
        symbols.symbol_value[irefireset], ilensymbol_value = _concatenate(1, 1, symbols.symbol_value[irefireset],
                                                                          ilensymbol_value)
    if iresetnumbers.iinitref != 0:
        commands.sometext = fstr(lw(iinit), commands.nstsize)
        _parse(5)
        symbols.symbol_value[irefiinit], ilensymbol_value = _concatenate(1, 1, symbols.symbol_value[irefiinit],
                                                                         ilensymbol_value)

    return nout, icount, ierrout


# --------------------------------------------------------------------------- read_ato
_F103_ATO = "(/a,1x,'No .ato file specified! - routine will exit'/)"
_F202 = "(/'Attempting to read file: ',a/)"
_F31 = "('No atoms were read from this file')"


@routine
def read_ato(nin, nout, icount, ifilecount, ierrout):
    """read_and_save_files.f:read_ato -- reads the EPSR .ato file named by symbol icount
    (minitc), and updates the rho symbol.

    Returns (nout, ierrout).
    """
    from . import ato_read_and_save, file_and_directory_commands  # packages D, A

    ireff = file_and_directory_commands.ireff
    # integer*4 isizefactor,irho
    # character*20 ext,text
    # character*1 iaccept
    iaccept = ""  # NB: uninitialised in the Fortran; intent of search_for_files is output
    #
    # ensure the input stream number is different from the stream number where data are
    # currently being read
    #
    nin1 = nin + 2  # noqa: F841
    #
    # set up the integer corresponding to fnameato - it is assumed that symbol(icount)
    # has already been checked before calling this routine
    #
    ifound = ireff(icount, ifilecount)
    #
    # get the symbol number for sizefactor and read sizefactor. If not found then assume
    # sizefactor is 1.0
    #
    text = "sizefactor"
    isizefactor = 0
    isizefactor = file_and_directory_commands.search_for_symbol(text, isizefactor)
    if isizefactor > 0:
        iref = ireff(isizefactor, ifilecount)
        molecom.sizefactor = read_internal(symbols.symbol_value[iref]).real(molecom.sizefactor)
    else:
        molecom.sizefactor = f32(1.0)

    # set up the expected extension of the .ato file

    ext = ".ato"
    ilenext = 4

    # open the .ato file and read in the relevant parameters from there

    if symbols.symbol_value[ifound] == "search" or symbols.symbol_value[ifound] == "<undefined>":

        # search for a suitable ato file

        commands.nwords = 0
        while commands.nwords <= 0:
            iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)
            _parse(5)
            if commands.nwords <= 0:
                units.write(nout, fmt(_F103_ATO, substr(symbols.prompt, 1, symbols.ilenprompt)))
                if symbols.symbol_save != "":
                    symbols.symbol_value[ifound] = symbols.symbol_save
                else:
                    symbols.symbol_value[ifound] = symbols.symbol_default[ifound]
                ierrout = 0
                return nout, ierrout
        symbols.fname_ato, symbols.ilenfname_ato = _concatenate(1, commands.nwords, symbols.fname_ato,
                                                                symbols.ilenfname_ato)

    else:

        commands.sometext = fstr(symbols.symbol_value[ifound], commands.nstsize)
        _parse(5)
        symbols.fname_ato, symbols.ilenfname_ato = _concatenate(1, commands.nwords, symbols.fname_ato,
                                                                symbols.ilenfname_ato)

    # ensure this file as the correct extension

    ilenlastdot = auxiliary_routines.indexlastdot(symbols.fname_ato) - 1
    if ilenlastdot > 0:
        symbols.fname_ato = fstr(substr(symbols.fname_ato, 1, ilenlastdot) + pad(ext, 20), 256)
        symbols.ilenfname_ato = ilenlastdot + ilenext
    else:
        symbols.fname_ato = fstr(substr(symbols.fname_ato, 1, symbols.ilenfname_ato) + pad(ext, 20), 256)
        symbols.ilenfname_ato = symbols.ilenfname_ato + ilenext

    #
    # convert this to a full path name if necessary
    #
    symbols.fnameato, symbols.ilenfnameato = file_and_directory_commands.get_full_filename(
        symbols.fname_ato, symbols.ilenfname_ato, symbols.fnameato, symbols.ilenfnameato)

    units.write(nout, fmt(_F202, substr(symbols.fnameato, 1, symbols.ilenfnameato)))
    #
    # now attempt to read the specified .ato file
    #
    #
    # initialise number of atom and molecule types to zero, after saving them
    #
    ntypesave = molecom.ntype
    nmoltypesave = clustercom.nmoltype
    molecom.ntype = 0
    clustercom.nmoltype = 0
    (epsr16.rmaxcell, epsr16.cellmax, epsr16.upol, epsr16.rho, epsr16.amubar, epsr16.amubar2, epsr16.ndist,
     ierrout) = ato_read_and_save.minitc(symbols.fnameato, epsr16.rmaxcell, epsr16.cellmax, epsr16.upol, epsr16.rho,
                                         epsr16.amubar, epsr16.amubar2, epsr16.ndist, ierrout)

    if ierrout != 0:
        #
        # the file did not exist or had errors, but we will leave any existing
        # atom parameters as they are and write them on saving to the new file
        #
        units.write(nout, fmt(_F31))
        #
        # restore the number of atom and molecule types
        #
        molecom.ntype = ntypesave
        clustercom.nmoltype = nmoltypesave
        ierrout = 0
    else:
        pass
        #
        # write out the list of atom types in this file
        #
        #         write(nout,32) ntype
        # 32    format(/'There are ',i2,' types of atom in this file'/)
        #         do it=1,ntype
        #            write(nout,33) it,type(it)
        # 33    format(1x,'Atom type ',i2,' has label ',a3)
        #         end do

    # save the current fnameato symbol_value

    symbols.symbol_value[ifound] = fstr(substr(symbols.fname_ato, 1, symbols.ilenfname_ato), 256)
    #
    # update the value of rho and save it
    #
    irho = 0
    text = "rho"
    irho = file_and_directory_commands.search_for_symbol(text, irho)
    if irho != 0:
        iref = ireff(irho, ifilecount)
        symbols.symbol_value[iref] = fstr(lw(epsr16.rho), 256)

    return nout, ierrout
