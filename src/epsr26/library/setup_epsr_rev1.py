"""Translation of fortran_src/setup_epsr_rev1.f: read_epsr16, read_epsr17, setup_symbols_epsr.

``setup_symbols_epsr`` defines the ``.EPSR.inp`` keyword table and its defaults.  (This
is the revision linked into ``epsr``; ``setup_epsr.f`` is an older, incompatible one.)
"""

from ..common import commands, epsr_dimension, iresetnumbers, symbols
from ..fortran import Ch, FortranEOF, FortranIOError, findex, findex_padded, fmt, fstr, lw, read_internal, substr, units
from ..trace import routine
from . import parse

_F499 = "(/'changeato> No .inp file specified! - routine will exit'/)"
_F102 = "(/'setup_input_file> Reading input file: \"',a,'\"')"
_F104 = "(/'setup_input_file> Run name in input file is different from ','filename specified')"
_F103 = "(/'setup_input_file> Problems with specified input file: \"',a,'\"'//'   - will use default values')"


@routine
def read_epsr16(nin, nout, ierrout):
    """setup_epsr_rev1.f:read_epsr16 -- read an old EPSR16 .inp file.

    sets up the EPSR input files
    Returns (nout, ierrout).
    """
    from . import auxiliary_routines, file_and_directory_commands, read_and_save_files

    ext = ""  # character*20
    text = ""  # character*256
    iaccept = ""  # character*1
    #
    # set up the symbol numbers for the fnameato and fnamepcof symbols
    #
    ifnameato = 0  # noqa: F841
    ifnamepcof = 0  # noqa: F841
    #
    # set the number of potential coefficient pairs to zero - this is used to
    # signal whether a valid .pcof file has been read in
    #
    npairspcof = 0  # noqa: F841  (a local here: potential_par.inc is not included)
    #
    # signal whether any atoms have been read in
    #
    ntype = 0  # noqa: F841  (a local here: molecom.inc is not included)
    #
    # generate the name of the input file
    #
    # search for the word 'EPSR16' in the input stream, and use the next word to define the
    # input filename, if present
    #
    iw = 1
    while iw <= commands.nwords and substr(commands.sometext, commands.ncf[iw], commands.ncl[iw]) != "epsr16":
        iw = iw + 1
    if iw >= commands.nwords:
        #
        # no filename found so prompt for it from the user
        #
        ext = fstr(".inp", 20)
        ilenext = 4
        commands.nwords = 0
        while commands.nwords <= 0:
            iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)
            commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                          commands.nstsize, commands.nworddim)
            if commands.nwords <= 0:
                units.write(6, fmt(_F499))
                return nout, ierrout
        # c100  format(a)
        iw = 0
    iw = iw + 1
    iw, commands.nwords, symbols.fname_inp, symbols.ilenfname_inp = parse.concatenate(
        iw, commands.nwords, symbols.fname_inp, symbols.ilenfname_inp,
        commands.nwords, commands.ncf, commands.ncl, commands.sometext, len_return_string=256)
    #
    # if this filename has a .ext we do nothing and will attempt to read it as is,
    # warts and all!
    #
    text = fstr(symbols.fname_inp, 256)
    ilenuptodot = findex(symbols.fname_inp, ".") - 1
    #     write(6,*) fname_inp(1:ilenfname_inp),ilenfname_inp,ilenuptodot
    if ilenuptodot < 0:
        #
        # .ext not supplied - therefore we need to add on the correct one
        #
        ilenuptodot = symbols.ilenfname_inp
        symbols.fname_inp = fstr(substr(symbols.fname_inp, 1, ilenuptodot) + ".inp", 256)
        symbols.ilenfname_inp = ilenuptodot + 4
    #
    # create the full name of the input file by adding on the directory name if not present
    #
    symbols.fnameinp, symbols.ilenfnameinp = file_and_directory_commands.get_full_filename(
        symbols.fname_inp, symbols.ilenfname_inp, symbols.fnameinp, symbols.ilenfnameinp)
    #
    # open this file if possible
    #
    try:
        units.open(nin, symbols.fnameinp, status="old")
        ierrout = 0
    except FortranIOError as e:
        ierrout = e.iostat
    if ierrout == 0:
        #
        #  If O.K. read in the data from the specified input file
        #
        units.write(nout, fmt(_F102, substr(symbols.fname_inp, 1, symbols.ilenfname_inp)))
        #
        # find the number of characters in this filename up to the last dot
        #
        text = fstr(symbols.fname_inp, 256)
        ilenuptodot = auxiliary_routines.indexlastdot(text) - 1
        #
        # read the name of the run - this should correspond to fname_inp but without the
        # last .extension. It will be used to construct the names of input
        # and output files (except .ato and .pcof files which may be named independently
        # from the name of the run)
        #
        ierrout = 0
        try:
            r = units.list_read(nin)
            symbols.fname_run = fstr(r.word(symbols.fname_run), 256)
        except FortranEOF:
            ierrout = -1
        except FortranIOError as e:
            ierrout = e.iostat
        symbols.ilenfname_run = findex_padded(symbols.fname_run, " ", 256) - 1
        #
        # remove the .extension if supplied
        #
        #
        # check that this name corresponds to filename (without the .extension) specified.
        # If not signal the error but carry on, using the correct name when writing the file
        #
        if (substr(symbols.fname_run, 1, symbols.ilenfname_run).rstrip(" ")
                != substr(symbols.fname_inp, 1, ilenuptodot).rstrip(" ")):
            units.write(nout, fmt(_F104))
            symbols.fname_run = fstr(substr(symbols.fname_inp, 1, ilenuptodot), 256)
            symbols.ilenfname_run = ilenuptodot
        #
        # generate the long form of fname_run
        #
        symbols.fnamerun, symbols.ilenfnamerun = file_and_directory_commands.get_full_filename(
            symbols.fname_run, symbols.ilenfname_run, symbols.fnamerun, symbols.ilenfnamerun)
        #
        # signal that if an error occurred the attempt to open the input file
        # was unsuccessful, but we can still read commands from the terminal
        #
        if ierrout != 0:
            ierrout = -1
            return nout, ierrout
        # NFILE = number of data files to fit (9 maximum)
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 36, 36, 1)
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
        units.write(6, lw(symbols.nfiles))
        for in_ in range(1, symbols.nfiles + 1):
            # INFILE = name of data file for S(Q) data
            file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 37, 37, in_)
            # FNAMEW = name of file containing weighting factors for this dataset
            file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 38, 38, in_)
            # type of input file (1,2, 3 or 4, 5)
            #     NRTYPE=1, Q, S(Q) ASCII values, 1 set per line
            #     NRTYPE=2, Q, S(Q), E(Q) ASCII values, 1 set per line
            #     NRTYPE=3,4 N, Q, S(Q), E(Q) ASCII values, 1 set per line, 4 blank lines
            # (GENIE SHOW DATA format, NRTYPE=3 corresponds to histogram format, 4 = point
            # mode)
            #     NRTYPE=5, Q, S(Q), E(Q) ASCII values, 1 set per line, using GUDRUN .DCS format (histogram)
            # RSHMIN = number of nodes in potential
            # SZEROS is the Q=0 limit for this file (currently not used. In the past it has been used to
            # normalise g(r) to give the correct zero limit).
            # TWEAK is a simple factor to multiply the data by prior to any processing. Normally 1.0.
            # EFILEREQ is the requested energy amplitude for this file - increase it to improve fit.
            # However these values will only work if PRESSREQ is .le. zero. If PRESSREQ is .gt. 0 then
            # EFILEREQ is set to PRESSREQ.
            # ** However if EFILEREQ is set to zero, then it will not be set by PRESSREQ. Instead the
            # simulation will calculate the structure factor for this dataset but will not refine the
            # empirical potential against these data.
            # ** Note that it may not be possible to achieve the specified energies for particular
            # datasets, particularly if a dataset contains more than one partial structure factor. This
            # is because the other datasets may be also making demands on each partial, so that the
            # different datasets may be in conflict.
            file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 39, 43, in_)
        # FEEDBACK = overall feedback factor. Must be =< 1.0. Should only be set to 1.0 if the
        # matrix of weight coefficients is well conditioned.
        # Can be increased or decreased to emphasize or de-emphasize the data.
        # POTFAC decides whether potential will be refined. POTFAC=0.0 means no refinement,
        # (but simulation will proceed as usual) otherwise refinement will take place
        # PRESSREQ is the requested absolute energy (kJ/mole) of empirical potential.
        # This can be used to control degree of fit and energy/pressure for all datasets together.
        # SIZEFACTOR (should be >=1) is used to shrink molecules during initial refinement when
        # molecular entanglement may occur. Advisable to set potfac to 0 if sizefactor > 1
        #
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 1, 1, 1)
        # NQ = No. of Q values (Maximum = 1000)
        # QSTEP = requested step in  Q
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 5, 6, 1)
        # NB: nq is an uninitialised local in the Fortran; 0 is kept if the value is null
        nq = read_internal(symbols.symbol_value[5]).int(0)
        if nq > epsr_dimension.mpts:
            nq = epsr_dimension.mpts
        symbols.symbol_value[5] = fstr(lw(nq), 256)
        # IRESET = 1 to set EP to zero, else IRESET = 0.
        # IINIT = 1 to zero all accumulators, else IINIT = 0. IINIT is automatically set whenever
        #     IRESET = 1
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 7, 7, 1)
        # will have to reset since EPSR17 requires different files and format to EPSR16
        ireset = 1
        symbols.symbol_value[7] = fstr(lw(ireset), 256)
        symbols.symbol_value[8] = fstr(lw(ireset), 256)
        # write these values to the appropriate symbol values
        # NTIMES = no. of trials per iteration between potential calculations (typically 5)
        # NITER = no. of iterations (typically 1)
        # NSUMS is the number of previous iterations that have already been accumulated
        #     If NSUMS is -1, then a running average is used (determined from INTER).
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 9, 11, 1)
        # INTRA is the number of intermolecular moves that are made before an INTRA-
        # molecular move (molecular vibration) is made - typically INTRA=50 or 100
        # INTER is the number of iterations used to form running averages - after calling
        # subroutine INITA, the value of NESUMT will be set to INTER
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 12, 12, 1)
        # number of iterations to occur after intramolecular disordering before
        # accumulating energy and pressure averages.
        # RHO no. density in scattering units per A**3 (this is determined from
        # the .ATO file - this value is written here for information only)
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 14, 14, 1)
        # CELLST = step size in linear region of r - controls number of r values
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 15, 15, 1)
        # resolution function of data (Gaussian) FWHM(Q) = FWHM + FWHMQ*Q
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 16, 17, 1)
        # NSMOOp = defines whether a background based on Rmin will be subtracted
        # from the data prior to fitting the potential
        # (0 for no background subtraction)
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 18, 18, 1)
        # read filename for molecular coordinates
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 19, 19, 1)
        nout, ierrout = read_and_save_files.read_ato(nin, nout, 19, 1, ierrout)
        # read filename for potential coefficients - this is not defined for EPSR16
        # so invent a new name
        symbols.symbol_value[20] = fstr(lw(substr(symbols.fname_run, 1, symbols.ilenfname_run) + ".pcof"), 256)
        commands.sometext = fstr(symbols.symbol_value[20], commands.nstsize)
        commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                      commands.nstsize, commands.nworddim)
        _one, commands.nwords, symbols.fname_pcof, symbols.ilenfname_pcof = parse.concatenate(
            1, commands.nwords, symbols.fname_pcof, symbols.ilenfname_pcof,
            commands.nwords, commands.ncf, commands.ncl, commands.sometext, len_return_string=256)
        #
        # convert this to a full path name if necessary
        #
        symbols.fnamepcof, symbols.ilenfnamepcof = file_and_directory_commands.get_full_filename(
            symbols.fname_pcof, symbols.ilenfname_pcof, symbols.fnamepcof, symbols.ilenfnamepcof)

        # write these values to the appropriate symbol values
        # read minimum Q value for molecular fit
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 21, 21, 1)
        units.close(nin)
    else:
        symbols.fname_run = fstr(substr(symbols.fname_inp, 1, ilenuptodot), 256)
        symbols.ilenfname_run = ilenuptodot
        units.write(nout, fmt(_F103, substr(symbols.fname_inp, 1, symbols.ilenfname_inp)))
        #
        # generate the long form of fname_run
        #
        symbols.fnamerun, symbols.ilenfnamerun = file_and_directory_commands.get_full_filename(
            symbols.fname_run, symbols.ilenfname_run, symbols.fnamerun, symbols.ilenfnamerun)
        #
        # signify that the input file could not be opened or read, but a run
        # name has been successfully defined with the default values
        #
        ierrout = -1
        return nout, ierrout
    return nout, ierrout


@routine
def read_epsr17(nin, nout, ierrout):
    """setup_epsr_rev1.f:read_epsr17 -- read an EPSR17 .inp file.

    sets up the EPSR input files
    Returns (nout, ierrout).
    """
    from . import auxiliary_routines, file_and_directory_commands, read_and_save_files

    ext = ""  # character*20
    iaccept = ""  # character*1
    text = ""  # character*256
    #
    # set up the symbol numbers for the fnameato and fnamepcof symbols
    #
    ifnameato = 0  # noqa: F841
    ifnamepcof = 0  # noqa: F841
    #
    # set the number of potential coefficient pairs to zero - this is used to
    # signal whether a valid .pcof file has been read in
    #
    npairspcof = 0  # noqa: F841  (a local here: potential_par.inc is not included)
    #
    # signal whether any atoms have been read in
    #
    ntype = 0  # noqa: F841  (a local here: molecom.inc is not included)
    #
    # generate the name of the input file
    #
    # search for the word 'EPSR17' in the input stream, and use the next word to define the
    # input filename, if present
    #
    iw = 1
    while iw <= commands.nwords and substr(commands.sometext, commands.ncf[iw], commands.ncl[iw]) != "epsr17":
        iw = iw + 1
    if iw >= commands.nwords:
        #
        # no filename found so prompt for it from the user
        #
        ext = fstr(".inp", 20)
        ilenext = 4
        commands.nwords = 0
        while commands.nwords <= 0:
            iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)
            commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                          commands.nstsize, commands.nworddim)
            if commands.nwords <= 0:
                units.write(6, fmt(_F499))
                return nout, ierrout
        # c100  format(a)
        iw = 0
    iw = iw + 1
    iw, commands.nwords, symbols.fname_inp, symbols.ilenfname_inp = parse.concatenate(
        iw, commands.nwords, symbols.fname_inp, symbols.ilenfname_inp,
        commands.nwords, commands.ncf, commands.ncl, commands.sometext, len_return_string=256)
    #
    # if this filename has a .ext we do nothing and will attempt to read it as is,
    # warts and all!
    #
    text = fstr(symbols.fname_inp, 256)
    ilenuptodot = findex(symbols.fname_inp, ".") - 1
    #     write(6,*) fname_inp(1:ilenfname_inp),ilenfname_inp,ilenuptodot
    if ilenuptodot < 0:
        #
        # .ext not supplied - therefore we need to add on the correct one
        #
        ilenuptodot = symbols.ilenfname_inp
        symbols.fname_inp = fstr(substr(symbols.fname_inp, 1, ilenuptodot) + ".inp", 256)
        symbols.ilenfname_inp = ilenuptodot + 4
    #
    # create the full name of the input file by adding on the directory name if not present
    #
    symbols.fnameinp, symbols.ilenfnameinp = file_and_directory_commands.get_full_filename(
        symbols.fname_inp, symbols.ilenfname_inp, symbols.fnameinp, symbols.ilenfnameinp)
    #
    # open this file if possible
    #
    try:
        units.open(nin, symbols.fnameinp, status="old")
        ierrout = 0
    except FortranIOError as e:
        ierrout = e.iostat
    if ierrout == 0:
        #
        #  If O.K. read in the data from the specified input file
        #
        units.write(nout, fmt(_F102, substr(symbols.fname_inp, 1, symbols.ilenfname_inp)))
        #
        # find the number of characters in this filename up to the last dot
        #
        text = fstr(symbols.fname_inp, 256)
        ilenuptodot = auxiliary_routines.indexlastdot(text) - 1
        #
        # read the name of the run - this should correspond to fname_inp but without the
        # last .extension. It will be used to construct the names of input
        # and output files (except .ato and .pcof files which may be named independently
        # from the name of the run)
        #
        ierrout = 0
        try:
            r = units.list_read(nin)
            symbols.fname_run = fstr(r.word(symbols.fname_run), 256)
        except FortranEOF:
            ierrout = -1
        except FortranIOError as e:
            ierrout = e.iostat
        symbols.ilenfname_run = findex_padded(symbols.fname_run, " ", 256) - 1
        #
        # remove the .extension if supplied
        #
        #
        # check that this name corresponds to filename (without the .extension) specified.
        # If not signal the error but carry on, using the correct name when writing the file
        #
        if (substr(symbols.fname_run, 1, symbols.ilenfname_run).rstrip(" ")
                != substr(symbols.fname_inp, 1, ilenuptodot).rstrip(" ")):
            units.write(nout, fmt(_F104))
            symbols.fname_run = fstr(substr(symbols.fname_inp, 1, ilenuptodot), 256)
            symbols.ilenfname_run = ilenuptodot
        #
        # generate the long form of fname_run
        #
        symbols.fnamerun, symbols.ilenfnamerun = file_and_directory_commands.get_full_filename(
            symbols.fname_run, symbols.ilenfname_run, symbols.fnamerun, symbols.ilenfnamerun)
        #
        # signal that if an error occurred the attempt to open the input file
        # was unsuccessful, but we can still read commands from the terminal
        #
        units.write(6, lw(ierrout))
        if ierrout != 0:
            ierrout = -1
            return nout, ierrout
        # NFILE = number of data files to fit (9 maximum)
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 36, 36, 1)
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
        units.write(6, lw(symbols.nfiles))
        for in_ in range(1, symbols.nfiles + 1):
            # INFILE = name of data file for S(Q) data
            file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 37, 37, in_)
            # FNAMEW = name of file containing weighting factors for this dataset
            file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 38, 38, in_)
            # type of input file (1,2, 3 or 4, 5)
            #     NRTYPE=1, Q, S(Q) ASCII values, 1 set per line
            #     NRTYPE=2, Q, S(Q), E(Q) ASCII values, 1 set per line
            #     NRTYPE=3,4 N, Q, S(Q), E(Q) ASCII values, 1 set per line, 4 blank lines
            # (GENIE SHOW DATA format, NRTYPE=3 corresponds to histogram format, 4 = point
            # mode)
            #     NRTYPE=5, Q, S(Q), E(Q) ASCII values, 1 set per line, using GUDRUN .DCS format (histogram)
            # RSHMIN = number of nodes in potential
            # SZEROS is the Q=0 limit for this file (currently not used. In the past it has been used to
            # normalise g(r) to give the correct zero limit).
            # TWEAK is a simple factor to multiply the data by prior to any processing. Normally 1.0.
            # EFILEREQ is the requested energy amplitude for this file - increase it to improve fit.
            # However these values will only work if PRESSREQ is .le. zero. If PRESSREQ is .gt. 0 then
            # EFILEREQ is set to PRESSREQ.
            # ** However if EFILEREQ is set to zero, then it will not be set by PRESSREQ. Instead the
            # simulation will calculate the structure factor for this dataset but will not refine the
            # empirical potential against these data.
            # ** Note that it may not be possible to achieve the specified energies for particular
            # datasets, particularly if a dataset contains more than one partial structure factor. This
            # is because the other datasets may be also making demands on each partial, so that the
            # different datasets may be in conflict.
            file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 39, 43, in_)
        # FEEDBACK = overall feedback factor. Must be =< 1.0. Should only be set to 1.0 if the
        # matrix of weight coefficients is well conditioned.
        # Can be increased or decreased to emphasize or de-emphasize the data.
        # POTFAC decides whether potential will be refined. POTFAC=0.0 means no refinement,
        # (but simulation will proceed as usual) otherwise refinement will take place
        # PRESSREQ is the requested absolute energy (kJ/mole) of empirical potential.
        # This can be used to control degree of fit and energy/pressure for all datasets together.
        # SIZEFACTOR (should be >=1) is used to shrink molecules during initial refinement when
        # molecular entanglement may occur. Advisable to set potfac to 0 if sizefactor > 1
        #
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 1, 4, 1)
        # NQ = No. of Q values (Maximum = 1000)
        # QSTEP = requested step in  Q
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 5, 6, 1)
        # NB: nq is an uninitialised local in the Fortran; 0 is kept if the value is null
        nq = read_internal(symbols.symbol_value[5]).int(0)

        if nq > epsr_dimension.mpts:
            nq = epsr_dimension.mpts
        symbols.symbol_value[5] = fstr(lw(nq), 256)

        # IRESET = 1 to set EP to zero, else IRESET = 0.
        # IINIT = 1 to zero all accumulators, else IINIT = 0. IINIT is automatically set whenever
        #     IRESET = 1
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 7, 8, 1)
        ireset = 1
        symbols.symbol_value[7] = fstr(lw(ireset), 256)
        symbols.symbol_value[8] = fstr(lw(ireset), 256)
        # NTIMES = no. of trials per iteration between potential calculations (typically 5)
        # NITER = no. of iterations (typically 1)
        # NSUMS is the number of previous iterations that have already been accumulated
        #     If NSUMS is -1, then a running average is used (determined from INTER).
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 9, 11, 1)
        # INTRA is the number of intermolecular moves that are made before an INTRA-
        # molecular move (molecular vibration) is made - typically INTRA=50 or 100
        # INTER is the number of iterations used to form running averages - after calling
        # subroutine INITA, the value of NESUMT will be set to INTER
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 12, 13, 1)
        # number of iterations to occur after intramolecular disordering before
        # accumulating energy and pressure averages.
        # RHO no. density in scattering units per A**3 (this is determined from
        # the .ATO file - this value is written here for information only)
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 14, 14, 1)
        # CELLST = step size in linear region of r - controls number of r values
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 15, 15, 1)
        # resolution function of data (Gaussian) FWHM(Q) = FWHM + FWHMQ*Q
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 16, 17, 1)
        # NSMOOp = defines whether a background based on Rmin will be subtracted
        # from the data prior to fitting the potential
        # (0 for no background subtraction)
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 18, 18, 1)
        # read filename for molecular coordinates
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 19, 19, 1)
        nout, ierrout = read_and_save_files.read_ato(nin, nout, 19, 1, ierrout)
        # read filename for potential coefficients
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 20, 20, 1)
        units.write(nout, lw(Ch(symbols.symbol_value[20], 256)))
        commands.sometext = fstr(symbols.symbol_value[20], commands.nstsize)
        commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                      commands.nstsize, commands.nworddim)
        _one, commands.nwords, symbols.fname_pcof, symbols.ilenfname_pcof = parse.concatenate(
            1, commands.nwords, symbols.fname_pcof, symbols.ilenfname_pcof,
            commands.nwords, commands.ncf, commands.ncl, commands.sometext, len_return_string=256)
        #
        # convert this to a full path name if necessary
        #
        symbols.fnamepcof, symbols.ilenfnamepcof = file_and_directory_commands.get_full_filename(
            symbols.fname_pcof, symbols.ilenfname_pcof, symbols.fnamepcof, symbols.ilenfnamepcof)

        # read minimum Q value for molecular fit
        file_and_directory_commands.get_a_line(nin, 10, commands.nstsize, 21, 21, 1)
        units.close(nin)
    else:
        symbols.fname_run = fstr(substr(symbols.fname_inp, 1, ilenuptodot), 256)
        symbols.ilenfname_run = ilenuptodot
        units.write(nout, fmt(_F103, substr(symbols.fname_inp, 1, symbols.ilenfname_inp)))
        #
        # generate the long form of fname_run
        #
        symbols.fnamerun, symbols.ilenfnamerun = file_and_directory_commands.get_full_filename(
            symbols.fname_run, symbols.ilenfname_run, symbols.fnamerun, symbols.ilenfnamerun)
        #
        # signify that the input file could not be opened or read, but a run
        # name has been successfully defined with the default values
        #
        ierrout = -1
        return nout, ierrout
    return nout, ierrout


@routine
def setup_symbols_epsr():
    """setup_epsr_rev1.f:setup_symbols_epsr -- the .EPSR.inp keyword table and defaults.

    sets the symbols and their default values for input to EPSR.

    set up the symbol table for the EPSR .inp file

    iflag decides whether this parameter can be changed
    a parameter can only be changed if ireset >= iflag
    iflag=1 means all accumulators are reset
    iflag=2 means potential coefficients are zerod as well
    Returns nothing.
    """

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('feedback', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Confidence factor - should be < 1. [0.8]', 256)
    symbols.iflag[iref] = -1
    symbols.symbol_default[iref] = fstr('0.9', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('potfac', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('>0.0 to enable potential refinement, 0.0 to inhibit', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    iresetnumbers.iereqfirst = iref  # ereq and parameters will be read and written to .pcof file

    iref = iref + 1
    symbols.symbol[iref] = fstr('ereq', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('EP amplitude[0.0].', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ereqmin', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum value for ereq [0.0]. Set 0 to ignore.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ereqmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum value for ereq [0.0]. Set 0 to ignore.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ereqstep', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Step size and unit size for ereq control [0] (Try 1 to initiate control)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('thresh', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Control baseline [0.2]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('bias', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Controls the bias on the steps in ereq. 0 means unbiassed steps. [2]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('sfreq', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Sampling frequency for trend line, and no. of times kT for automatic ereq to begin [50 1]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('50 1', 256)

    iresetnumbers.iereqlast = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('rspcrmin', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum distance for calculating the R-space coefficient [1.00]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.00', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rspcfrac', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Fraction of R-space coefficient in control level [0.2]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.5', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('num_threds', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('No. of parallel threads to be used, (0 to let program choose) [0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('porodpwr', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Q-value and power for Porod scattering [0 0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nmolcell', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Average number of molecules in a cell [10]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('10', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('control-p', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Control pressure [a, b and/or c axes control plus value plus step, none 0 0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('none 0 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ref_intra', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Weighting on EP and Coulomb terms for intra-molecular structure [0 0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('sizefactor', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Multiplying factor for box dimension, decline rate, and threshold. [1.0 0.9 0.0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0 0.9 0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nq', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of Q values. [600]', 256)
    symbols.iflag[iref] = -1
    symbols.symbol_default[iref] = fstr('600', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('qstep', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Size of Q step [1/A]. [0.05]', 256)
    symbols.iflag[iref] = -1
    symbols.symbol_default[iref] = fstr('0.05', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('qinscalefc', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Scale factor on input Q values and broadening power for self terms. [1.0 0.0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0 0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ireset', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('1: complete reset; 2: sets the Empirical Potential to zero', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iresetnumbers.iresetref = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('iinit', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Sets accumulators to zero. Recalculates r and Q. [1]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iresetnumbers.iinitref = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('ntimes', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of MC cycles between potential refinements. [5]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('5', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('niter', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of potential refinements before exitting. [1]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsumt', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of iterations already accumulated. [-1 with reset]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('-1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('intra', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of iterations between molecule shakes. [100]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('100', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rotfreq', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of iterations between internal rotation moves. [5]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('5', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('inter', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of iterations in running averages. [5]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('5', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rho', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Atomic number density - will be derived from .ato file', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('0.1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('cellst', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Size of r step [A]. [0.03]', 256)
    symbols.iflag[iref] = -1
    symbols.symbol_default[iref] = fstr('0.03', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmaxgr', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Range of g(r) and F.T. (0.0 will use half the cell box) [0.0]', 256)
    symbols.iflag[iref] = -1
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ngrsamples', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Requested no. of origin molecules to sample g(r). (0 will use 1000 molecules) [0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('fwhm', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Resolution width - Q independent term. [0.0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('fwhmq', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Resolution width - Q dependent term. [0.02 for SLS]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.02', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsmoop', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('1 means background subtraction is ON, 0 means OFF', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnameato', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of .ato file', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnamepcof', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of potential coefficients file.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('revlorch', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Broadening factor in Q space. [0.0 0.0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0 0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('qwidthqmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Broadening and maximum Q for Bragg peak calculation', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.01 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('mplicities', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('No. of unit cells along a, b and c for Bragg peak calculation', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1 1 1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('hklqmin', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum value of qhkl to be used, minimum radius for Bragg g(r), and Debye-Waller factor', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0 0.0 0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('diffuse', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('No. of unit cells along a, b and c for diffuse scattering calculation, maximum l, steps in l and m [0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nlatrefine', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('No. of lattice refinement steps [0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('latrefine', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Refinement steps for each of a, b, c, alpha, beta, gamma [0 0 0 0 0 0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0 0 0 0 0 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rejrate', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Rejection rate [0.75]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.75', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('qmin', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Qmin for Fourier transforms and for potential fits. [0.05 0.0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.05 0.0', 256)

    iresetnumbers.iqminref = iref  # The following parameters are stored in the .pcof file

    iref = iref + 1
    symbols.symbol[iref] = fstr('roverlap', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum allowed intermolecular separation between two atoms', 256)
    symbols.iflag[iref] = -3
    symbols.symbol_default[iref] = fstr('0.0', 256)
    iresetnumbers.iroverlap = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('rminfac', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Factor to set the minimum separation between pairs. [0.0]', 256)
    symbols.iflag[iref] = -3
    symbols.symbol_default[iref] = fstr('0.0', 256)
    iresetnumbers.irminfac = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('rminpt', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Radius at which potential truncation begins.', 256)
    symbols.iflag[iref] = -3
    symbols.symbol_default[iref] = fstr('9', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmaxpt', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Radius at which potential truncation goes to 0.0.', 256)
    symbols.iflag[iref] = -3
    symbols.symbol_default[iref] = fstr('12', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rbroad', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Controls potential decay. [0.0]', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('expecf', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Potential decay for short distance repulsive term and gradient limit. [0.3 0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.3 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('reppottype', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Repulsive potential type: exponential or harmonic. [exponential]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('exponential', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('addpottype', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Additional potential type: Gaussian or modmorse. [Gaussian]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('Gaussian', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rcharge', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Calculates energy due to molecular polarisation. [0.0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('power', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Repulsive power in Lennard-Jones potential. [12]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('12', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ncoeffp', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of coefficients used to define the EP. Set by program.', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('120', 256)
    iresetnumbers.incoeffp = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('pdmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum distance of Empirical Potential. Set by program.', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('12.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('pdstep', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Spacing between coefficients in r. Set by program.', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('0.1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('npitss', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of steps for refining the potential. [1000]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1000', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('paccept', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Acceptance factor for potential refinement. [0.0005]', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('0.0005', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('Gaussian', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Select "Gaussian" for Gaussian representation of EP. Otherwise Poisson. [Poisson]', 256)
    symbols.iflag[iref] = -3
    symbols.symbol_default[iref] = fstr('Poisson', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('psigma2', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Width for empirical potential functions. [0.01]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.01', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('refpotfac', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Factor to apply to reference potential. [1.0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ndata', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of data files to be fit by EPSR', 256)
    symbols.iflag[iref] = -1
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('datafile', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of data file to be fit', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('wtsfile', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of weights file for this data set', 256)
    symbols.iflag[iref] = -1
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nrtype', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Data type - see User Manual for more details', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('5', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rshmin', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum radius [A] - used for background subtraction', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.7', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('szeros', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Zero limit - 0 means use first data point for Q=0', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('tweak', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Scaling factor for this data set. [1.0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('efilereq', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Requested energy amplitude for this data set [1.0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('porodadd', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Additional amplitude of Porod scattering for this dataset [0.0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    symbols.nsymboldata = iref
    symbols.nsymbol = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    #
    # write the defaults to symbol_value
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
    symbols.prompt = fstr("Setup epsr>", 80)
    symbols.ilenprompt = findex(symbols.prompt, ">")
    #
    # set data prompt text
    #
    symbols.dataprompt = fstr("data", 80)
    symbols.ilendataprompt = 4
    symbols.dataprompttext = fstr("jumps to the specified dataset", 80)
    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr(".EPSR.inp", 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, " ", 20) - 1
    #
    # default number of files
    #
    symbols.nfiles = 1
    return
