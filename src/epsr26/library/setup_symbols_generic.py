"""Translation of fortran_src/setup_symbols_generic.f: the keyword tables of the auxiliary programs.

RH: setup_symbols.f
for each of the analysis rountines sets up arrays/tables which are
used to read the input files (the *.dat files).

``setup_symbols_generic`` (hand-written below) reads a ``<program>.symbols`` file; every
other routine (``setup_symbols_changeato``, ``setup_symbols_partials``, ...) is a
hard-coded table and is generated, statement by statement, by
``tools/gen/gen_setup_symbols_generic.py`` (do not edit those by hand: change the
generator and rerun it).
"""

from ..common import chains, commands, lsi, symbols
from ..fortran import (
    Ch,
    FortranEOF,
    FortranIOError,
    findex,
    findex_padded,
    fmt,
    fstr,
    iostat_of,
    lw,
    pad,
    read_internal,
    substr,
    units,
)
from ..trace import routine
from . import parse

_F102 = "(/'setup_symbols_generic> ',a)"
_F107 = "(1x,i3,1x,a11,1x,a30,1x,a30)"


def _mem(sometext2, symbols_fname):
    """The memory gfortran's parse(sometext2,...,nstsize,...) reads: sometext2, then symbols_fname.

    Not a Fortran routine: see the NB in setup_symbols_generic.  Beyond symbols_fname
    the scan never goes (it stops after the 13 blanks that follow the file name).
    """
    return pad(sometext2, 256) + pad(symbols_fname, 256)


@routine
def setup_symbols_generic(auxprogname, ssg_ierr):
    """setup_symbols_generic.f:setup_symbols_generic -- table from a <program>.symbols file.

    Routine to parse a .symbols file to allow external programs to use setup.

    In the .symbols file required parameters are:
    - extension   (the extension for the .dat file)
    - prompt      (the prompt used in the setup program)

    And if using the analysis routines which calculate for more than one set of input
    parameters (ndist > 0), then
    - dataprompt   (the generic name of these different parameter sets)
    - dataprompttext (a description of this for the  .dat file)
    are required.

    Other than this all other symbols are declared by prefixing a line with the word
    "symbol" followed by the value for that symbol. Then on following lines need to
    specify default, type and text description for each symbol (these are indicated by
    the key words symbol_default, symbol_type, and symbol_text. There are three types a
    symbol can have: real -- specified by an r character, integer -- specified by an i,
    and character -- indicated by a c character.

    Blank lines are ignored.  (Not between the symbol_default/_type/_text lines: there a
    blank line is "Unrecognised", as in the Fortran.)

    Example .symbols file:

    extension .NDIST_EXAMPLE.dat
    prompt setup example>
    dataprompt ndist_example
    dataprompttext   ndist_example text, nice

    symbol fnameato
    symbol_default <undefined>
    symbol_type  c
    symbol_text  A nice description

    symbol ndist
    symbol_default 1
    symbol_type i
    symbol_text Number of distributions/variable sets to calculate over

    auxprogname is character*256.  Returns ssg_ierr.
    """
    sometext2 = ""  # character*256
    symbols_fname = ""  # character*256
    val = ""  # character*256
    ext = ""  # character*8
    cmd = ""  # character*20

    # 102  format(/'setup_symbols_generic> ',a)
    # 100  format(a)

    commands.nwords = parse.parse(auxprogname, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    # Need to parse up to space...
    units.write(6, fmt(_F102, " command passed in: "
                       + substr(auxprogname, commands.ncf[1], commands.ncl[1])))
    units.write(6, fmt(_F102, " command passed in: " + pad(auxprogname, 256)))

    # extention for template file name
    ext = fstr(".symbols", 8)
    ilenext = 8

    # Get file name to read:
    symbols_fname = fstr(substr(auxprogname, commands.ncf[1], commands.ncl[1]) + pad(ext, 8), 256)
    ilen_sf = commands.ncl[1] - commands.ncf[1] + ilenext + 1

    # initialise variables which need to be checked that they have been
    # parsed:
    symbols.extinp = ""
    symbols.prompt = ""
    symbols.dataprompt = ""
    symbols.dataprompttext = ""
    symbols.nsymbol = 0
    symbols.nsymboldata = 0
    symbols.nsymboluptodata = 0

    #     open this file if possible
    nin = 14
    ierr = 0
    try:
        units.open(nin, substr(symbols_fname, 1, ilen_sf), status="old")
    except FortranIOError as e:
        ierr = iostat_of(e)
    if ierr == 0:
        units.write(6, fmt(_F102, "Processing " + substr(symbols_fname, 1, ilen_sf)))
        # process the .symbols file now
        ierr = 0
        lnum = 0
        iref = 0
        while ierr == 0:
            commands.nwords = 0
            # Have to use sometext2 rather than sometext as setup is using
            # sometext
            sometext2 = ""
            try:
                sometext2 = fstr(units.read_line(nin), 256)  # read next line
            except (FortranEOF, FortranIOError) as e:
                ierr = iostat_of(e)
            # NB: parse is told the text is nstsize (3020) characters long, but sometext2 is
            # character*256: on a blank line (or at end of file) it scans on into the next
            # local in memory, symbols_fname, and finds its first word (gfortran 15, -O2
            # -fopenmp; checked by tests/library/test_setup_symbols_generic.py).  _mem
            # reproduces that; on a line with text the scan stops long before column 256.
            mem = _mem(sometext2, symbols_fname)
            commands.nwords = parse.parse(mem, commands.nwords, commands.ncf, commands.ncl, 12,
                                          commands.nstsize, commands.nworddim)
            # If there are words then process them; ignore blank lines:
            cmd = ""
            if commands.nwords > 0:
                # cmd is first word on each line; val is the rest of the line
                cmd = fstr(substr(mem, commands.ncf[1], commands.ncl[1]), 20)
                val = fstr(substr(mem, commands.ncf[2], commands.ncl[commands.nwords]), 256)
                if cmd == "extension":
                    units.write(6, lw("Found extension, extension = ", Ch(val, 256)))
                    symbols.extinp = fstr(val, 20)
                    symbols.ilenextinp = findex_padded(symbols.extinp, " ", 20) - 1

                elif cmd == "prompt":
                    units.write(6, lw("Found prompt, prompt = ", Ch(val, 256)))
                    symbols.prompt = fstr(val, 80)
                    symbols.ilenprompt = -1
                    symbols.ilenprompt = findex(symbols.prompt, ">")
                    if symbols.ilenprompt <= 0:
                        units.write(6, fmt(_F102, "Prompt must end in a > character"))
                        units.write(6, fmt(_F102, "Exiting"))
                        ssg_ierr = -1
                        return ssg_ierr

                elif cmd == "dataprompt":
                    units.write(6, lw("Found dataprompt, dataprompt = ", Ch(val, 256)))
                    symbols.dataprompt = fstr(val, 80)
                    symbols.ilendataprompt = findex_padded(symbols.dataprompt, " ", 80) - 1

                elif cmd == "dataprompttext":
                    units.write(6, lw("Found dataprompttext, dataprompttext = ", Ch(val, 256)))
                    symbols.dataprompttext = fstr(val, 80)

                elif cmd == "symbol":
                    units.write(6, lw("Found new symbol, symbol = ", Ch(val, 256)))
                    iref = iref + 1  # found another new symbol
                    symbols.symbol[iref] = fstr(val, 11)
                    symbols.iflag[iref] = 0

                    if val == "ndist":
                        symbols.nsymboluptodata = iref
                        units.write(6, fmt(_F102, "Found ndist symbol, therefore all"
                                           + "symbols beyond this are part of each "
                                           + "distribution."))

                    # Check that get default, type, and text -- so read next 3 lines:
                    il = 0
                    # variables to check that parsed def:
                    parsed_def = -1
                    parsed_type = -1
                    parsed_text = -1
                    while il < 3 and ierr == 0:
                        sometext2 = ""
                        commands.nwords = 0
                        try:
                            sometext2 = fstr(units.read_line(nin), 256)  # read next line
                        except (FortranEOF, FortranIOError) as e:
                            ierr = iostat_of(e)
                        mem = _mem(sometext2, symbols_fname)  # NB: see the outer loop
                        commands.nwords = parse.parse(mem, commands.nwords, commands.ncf, commands.ncl,
                                                      12, commands.nstsize, commands.nworddim)
                        if commands.nwords >= 0:  # ignore blank lines
                            # NB: always true, so a blank line (or end of file) here is not
                            # ignored: parse finds the first word of symbols_fname (see above),
                            # and the line is "Unrecognised value ... : <name>.symbols".
                            cmd = fstr(substr(mem, commands.ncf[1], commands.ncl[1]), 20)
                            val = fstr(substr(mem, commands.ncf[2], commands.ncl[commands.nwords]), 256)

                            if cmd == "symbol_default":
                                parsed_def = 1
                                symbols.symbol_default[iref] = fstr(val, 256)
                            elif cmd == "symbol_type":
                                parsed_type = 1
                                symbols.symbol_type[iref] = fstr(val, 1)
                            elif cmd == "symbol_text":
                                parsed_text = 1
                                symbols.symbol_text[iref] = fstr(val, 256)
                            else:
                                units.write(6, fmt(_F102, " "))
                                units.write(6, lw("Unrecognised value for symbol ",
                                                  Ch(symbols.symbol[iref], 11), " : ", Ch(cmd, 20),
                                                  " in line ", lnum))
                                units.write(6, fmt(_F102, "Exiting setup_symbol_generic"))
                                ssg_ierr = -1
                                return ssg_ierr

                            il = il + 1
                        lnum = lnum + 1
                    # end of loop thru' def, type and text lines
                    if parsed_def <= 0:
                        units.write(6, fmt(_F102, " "))
                        units.write(6, lw("No default found for symbol: ", Ch(symbols.symbol[iref], 11)))
                        units.write(6, fmt(_F102, "Exiting setup_symbol_generic"))
                        ssg_ierr = -1
                        return ssg_ierr
                    if parsed_type <= 0:
                        units.write(6, fmt(_F102, " "))
                        units.write(6, lw("No type found for symbol: ", Ch(symbols.symbol[iref], 11)))
                        units.write(6, fmt(_F102, "Exiting setup_symbol_generic"))
                        ssg_ierr = -1
                        return ssg_ierr
                    if parsed_text <= 0:
                        units.write(6, fmt(_F102, " "))
                        units.write(6, lw("No text found for symbol: ", Ch(symbols.symbol[iref], 11)))
                        units.write(6, fmt(_F102, "Exiting setup_symbol_generic"))
                        ssg_ierr = -1
                        return ssg_ierr

            lnum = lnum + 1

        units.write(6, fmt(_F102, " "))
        units.write(6, lw("Parsed ", lnum, " lines in " + substr(symbols_fname, 1, ilen_sf)))
        units.write(6, lw("Found ", iref, " symbols "))

        units.close(nin)

    else:
        # If there is a problem we exit setup:
        units.write(6, fmt(_F102, "Problem reading file: " + substr(symbols_fname, 1, ilen_sf)))
        units.write(6, fmt(_F102, "Exiting setup_symbol_generic"))
        ssg_ierr = -1
        return ssg_ierr

    # Check that have parsed the important variables: prompt,
    # extension...
    if symbols.extinp == "":
        units.write(6, fmt(_F102, " "))
        units.write(6, lw("No extension specified in file : ", substr(symbols_fname, 1, ilen_sf)))
        units.write(6, fmt(_F102, "Extension is required."))
        units.write(6, fmt(_F102, "Exiting setup_symbol_generic"))
        ssg_ierr = -1

        return ssg_ierr
    if symbols.prompt == "":
        units.write(6, fmt(_F102, " "))
        units.write(6, lw("No prompt specified in file : ", substr(symbols_fname, 1, ilen_sf)))
        units.write(6, fmt(_F102, "Prompt is required."))
        units.write(6, fmt(_F102, "Exiting setup_symbol_generic"))
        ssg_ierr = -1
        return ssg_ierr

    if symbols.dataprompt == "" and symbols.nsymboluptodata != 0:
        units.write(6, fmt(_F102, " "))
        units.write(6, lw("No dataprompt specified in file : ", substr(symbols_fname, 1, ilen_sf)))
        units.write(6, fmt(_F102, "If have set ndist then dataprompt is required."))
        units.write(6, fmt(_F102, "Exiting setup_symbol_generic"))
        ssg_ierr = -1
        return ssg_ierr

    if symbols.dataprompttext == "" and symbols.nsymboluptodata != 0:
        units.write(6, fmt(_F102, " "))
        units.write(6, lw("No dataprompttext specified in file : ", substr(symbols_fname, 1, ilen_sf)))
        units.write(6, fmt(_F102, "If have set ndist then dataprompttext"
                           + " is required."))
        units.write(6, fmt(_F102, "Exiting setup_symbol_generic"))
        ssg_ierr = -1
        return ssg_ierr

    # nsymboldata specifies the end of the file -- the number of symbols
    # to read from the file
    symbols.nsymboldata = iref

    # If not set any distributions set nsymboluptodata to be beyond the
    # number of symbols nsymbol:
    if symbols.nsymboluptodata == 0:
        symbols.nsymboluptodata = iref + 1
        # the character constant is continued onto the next source line: the
        # two halves join with nothing in between ("...therefore nodistributions...")
        units.write(6, fmt(_F102, "ndist not found, therefore no"
                           "distributions specified"))
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata
    symbols.nsymbol = iref

    # Initialise symbol values and print symbols:
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F107, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256),
                           Ch(symbols.symbol_type[is_], 1)))
        # 107     format(1x,i3,1x,a11,1x,a30,1x,a30)

    #     get the number of files in this input file
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1

    units.write(6, fmt(_F102, "Complete setup_symbol_generic",
                       "-- symbols() have been setup"))
    ssg_ierr = 0
    return ssg_ierr


# --- GENERATED BELOW by tools/gen/gen_setup_symbols_generic.py; do not edit by hand ---
_F100 = '(1x,i3,1x,a11,1x,a30)'
_F101_rings = "('Maximum ring length (max ',i2,')')"  # setup_symbols_rings
_F101_chains = "('Maximum chain length (max ',i3,')')"  # setup_symbols_chains
_F101_lsi = "('Number of LSI values (max ',i3,')')"  # setup_symbols_lsi


@routine
def setup_symbols_changeato():
    """setup_symbols_generic.f:setup_symbols_changeato -- keyword table and defaults.

    sets the symbols and their default values for input to EPSR.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to EPSR.
    #
    #
    # set up the symbol table for the EPSR .inp file
    #

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('bond', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Intra-molecular bond lengths. Type y to action.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('n', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('label', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Change atom labels. Type y to action.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('n', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('tether', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('To change tethered status of specified molecule types. Type y to action.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('n', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('remove', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Remove specified molecule types. Type y to action.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('n', 256)

    iref = iref + 1
    symbols.irefrho = iref
    symbols.symbol[iref] = fstr('rho', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Density of this .ato file. [atoms/A**3]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('temp', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Temperature of this .ato file.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('300', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('tol', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Weighting for tethered atoms', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('stepmi', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Intramolecular translation step.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('stepri', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Intramolecular rotation step.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('stepr', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Intermolecular rotation step.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('stepm', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Intermolecular translation step.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('vibtemp', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Vibrational energy weight [65].', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('65', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('angtemp', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Requested bond angle energy weight [3.0].', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('dihtemp', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Requested dihedral energy weight [10.0].', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('10.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ecore', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Energy parameter for non-bonded intramolecular correlations', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('dcore', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Distance parameter for non-bonded intramolecular correlations', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('box', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Shape of simulation box - [cubic] or [noncubic] or [change]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('cubic', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('findmols', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Generate molecules from atoms? Type y.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('n', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('getmolfs', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Generate .mol files from atoms? Type y.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('n', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('addbond', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Create a new intra-molecular bond? Type y.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('n', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('removebond', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Remove an intra-molecular bond? Type y.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('n', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('add_w_rot', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Add a WHOLE molecule rotational group? Type y.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('n', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('axisshift', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Shift molecules of one type along specified axis? Type y.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('n', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('axisrotate', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Rotate molecules of one type around specified axis? Type y.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('n', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('alignmole', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Align molecules along specified axis? Type y.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('n', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('setimage', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Set an image atom for specified molecule? Type y.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('n', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('unsetimage', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Unset an image atom for specified molecule? Type y.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('n', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('removerot', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Remove a rotational group? Type y.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('n', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('natomtypes', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of atom types in this file - set by .ato file', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('molscalefc', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Generate molecular scale factors. Type y to action.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('n', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('atlabel', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom label in .ato file', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr(' ', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('epsilon', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Lennard-Jones epsilon [kJ/mol].', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.65', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('sigma', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Lennard-Jones sigma [A].', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('3.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atmass', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Atomic mass [amu].', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('12', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('charge', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Coulomb charge [e].', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('qradius', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Radius of charge distribution [A].', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atsymbol', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Chemical symbol for this component.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('iexchange', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('= 1 if atom exchanges with other atoms, else 0.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('-1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('natomsph', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('No. of atoms represented by sphere (1.0).', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('radminsph', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum radius for coarse graining [A] (0.0).', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('radmaxsph', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum radius for coarse graining [A] (0.0).', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('cgmolfac', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Coarse graining factor for this type of atom (1.0).', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    symbols.nsymbol = iref
    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    symbols.nfiles = 1

    #
    # initialise the symbol_values to their default values
    #
    for icount in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[icount] = symbols.symbol_default[icount]

    #
    # signal the type of symbols that have been set up
    #
    symbols.prompt = fstr('changeato>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')
    #
    # set data prompt text
    #
    symbols.dataprompt = fstr('atom', 80)
    symbols.ilendataprompt = 4
    symbols.dataprompttext = fstr('jumps to the specified atom type number', 80)

    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.ato', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_partials():
    """setup_symbols_generic.f:setup_symbols_partials -- keyword table and defaults.

    sets the symbols and their default values for input to clusters.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to clusters.
    #
    #
    # set up the symbol table for the PARTIALS input file
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if ireset >= iflag
    #

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnameato', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of .ato file ', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nz', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of r values [max: 1000]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('300', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rnormalise', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Controls the r- or z-range of distribution [1.0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsumt', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of configurations already accumulated', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('-1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom-c', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Centre atom type [*]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('*', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom-z', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Centre atom type for rho-z [*]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('*', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom-r', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Centre atom type for rho-r [*]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('*', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ndist', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of site-site distributions', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref
    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    symbols.nsymbol = iref

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)

    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup partials>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')
    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('partial', 80)
    symbols.dataprompttext = fstr('jumps to the specified cluster ', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1
    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.PARTIALS.dat', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_coord():
    """setup_symbols_generic.f:setup_symbols_coord -- keyword table and defaults.

    sets the symbols and their default values for input to clusters.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to clusters.
    #
    #
    # set up the symbol table for the COORD input file
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if ireset >= iflag
    #

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnameato', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of .ato file ', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsize', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Maximum coordination number (max 200)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('200', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsumt', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of configurations already accumulated', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('-1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom-c', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Centre atom type [*]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('*', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom-z', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Centre atom type for rho-z [*]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('*', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom-r', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Centre atom type for rho-r [*]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('*', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ndist', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of coordination number distributions', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom1', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom type 1 to define a coordination number', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom2', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom type 2 to define a coordination number', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmin', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum distance for this bond', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum distance for this bond', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rminrmaxc', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum and maximum distance for atom 2 from atom c', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1 2', 256)

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    symbols.nsymbol = iref

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)

    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup coord>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')
    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('distribution', 80)
    symbols.dataprompttext = fstr('jumps to the specified distribution ', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1
    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.COORD.dat', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_clusters():
    """setup_symbols_generic.f:setup_symbols_clusters -- keyword table and defaults.

    sets the symbols and their default values for input to clusters.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to clusters.
    #
    #
    # set up the symbol table for the CLUSTERS input file
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if ireset >= iflag
    #

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnameato', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of .ato file ', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsize', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Maximum cluster size', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('5000', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsumt', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of configurations already accumulated', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('-1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ndist', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of cluster distributions', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom1', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom type 1 to define a bond', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom2', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom type 2 to define a bond', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmin', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum distance for this bond', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum distance for this bond', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2', 256)

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    symbols.nsymbol = iref

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)

    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup clusters>', 80)

    symbols.ilenprompt = findex(symbols.prompt, '>')
    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('cluster', 80)
    symbols.dataprompttext = fstr('Specifies a cluster ', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1
    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.CLUSTERS.dat', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_rings():
    """setup_symbols_generic.f:setup_symbols_rings -- keyword table and defaults.

    sets the symbols and their default values for input to clusters.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to clusters.
    #
    #
    # set up the symbol table for the COORD input file
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if ireset >= iflag
    #

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnameato', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of .ato file ', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsize', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr(fmt(_F101_rings, chains.mring), 256)
    # 101 format('Maximum ring length (max ',i2,')')
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr(lw(chains.mring), 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsumt', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of configurations already accumulated', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('-1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ndist', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of ring distributions', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom1', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom type 1 for first atom in first molecule', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom2', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom type 2 for second atom in first molecule', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom3', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom type 3 for atom in second molecule', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmin12', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum distance for atom1-atom2', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax12', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum distance for atom1-atom2', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmin13', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum distance for atom1-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax13', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum distance for atom1-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmin23', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum distance for atom2-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax23', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum distance for atom2-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ang213min', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum included angle for atom2-atom1-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ang213max', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum included angle for atom2-atom1-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('180.0', 256)

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    symbols.nsymbol = iref

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)

    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup rings>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')
    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('ring', 80)
    symbols.dataprompttext = fstr('jumps to the specified ring ', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1
    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.RINGS.dat', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_chains():
    """setup_symbols_generic.f:setup_symbols_chains -- keyword table and defaults.

    sets the symbols and their default values for input to clusters.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to clusters.
    #
    #
    # set up the symbol table for the COORD input file
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if ireset >= iflag
    #

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnameato', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of .ato file ', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsize', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr(fmt(_F101_chains, chains.mchain), 256)
    # 101 format('Maximum chain length (max ',i3,')')
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr(lw(chains.mchain), 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsumt', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of configurations already accumulated', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('-1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ndist', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of chain distributions', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom1', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom type 1 for first atom in first molecule', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom2', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom type 2 for second atom in first molecule', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom3', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom type 3 for atom in second molecule', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmin12', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum distance for atom1-atom2', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax12', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum distance for atom1-atom2', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmin13', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum distance for atom1-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax13', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum distance for atom1-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmin23', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum distance for atom2-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax23', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum distance for atom2-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ang213min', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum included angle for atom2-atom1-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ang213max', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum included angle for atom2-atom1-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('180.0', 256)

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    symbols.nsymbol = iref

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)
    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1
    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup chains>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')
    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('chain', 80)
    symbols.dataprompttext = fstr('jumps to the specified chain ', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1
    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.CHAINS.dat', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_lsi():
    """setup_symbols_generic.f:setup_symbols_lsi -- keyword table and defaults.

    sets the symbols and their default values for input to chains.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to chains.
    #
    #
    # set up the symbol table for the LSI input file
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if ireset >= iflag
    #

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnameato', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of .ato file ', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsize', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr(fmt(_F101_lsi, lsi.mchain), 256)
    # 101 format('Number of LSI values (max ',i3,')')
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr(lw(lsi.mchain), 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('lsiminmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum and maximum value of LSI [-1 2]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('-1 2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsumt', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of configurations already accumulated', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('-1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ndist', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of chain distributions', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom1', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom type 1 for first atom in first molecule', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom2', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom type 2 for second atom in first molecule', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom3', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom type 3 for atom in second molecule', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmin12', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum distance for atom1-atom2', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax12', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum distance for atom1-atom2', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmin13', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum distance for atom1-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax13', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum distance for atom1-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmin23', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum distance for atom2-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax23', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum distance for atom2-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ang213min', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum included angle for atom2-atom1-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ang213max', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum included angle for atom2-atom1-atom3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('180.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rtrimax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum atomic separation for triangle/q calculation', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nbonded', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Average number of bonded molecules per atom1', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0 0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nnonbonded', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Average number of non-bonded molecules per atom1', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0 0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('avbond', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Average bonded atom1-atom3 distance', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0 0.0', 256)

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    symbols.nsymbol = iref

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)
    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1
    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup lsi>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')
    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('lsi', 80)
    symbols.dataprompttext = fstr('jumps to the specified lsi ', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1
    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.LSI.dat', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_sharm():
    """setup_symbols_generic.f:setup_symbols_sharm -- keyword table and defaults.

    sets the symbols and their default values for input to SHARM.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to SHARM.
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if iflag=0
    #

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnameato', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of .ato file ', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nr', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Number of radius values (max 200)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('100', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum radius for spherical harmonic coefficients', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('10', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsumt', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of configurations already accumulated', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ncoeffs', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of coefficients (program calculates this)', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('l1values', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('L1 values (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('l2values', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('L2 values (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('lvalues', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('L values (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('n1step', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Step in N1 values', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('n2step', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Step in N2 values', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom-c', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Central molecule - list of centre atom types', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('axisc1', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('First axis definition for central molecule', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('axisc2', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Second axis definition for central molecule', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom-s', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Second molecule - list of centre atom types', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('axiss1', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('First axis definition for second molecule', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('axiss2', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Second axis definition for second molecule', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('test', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Set phi1,theta1,chi1,phi2,theta2,phi2,thetal,phil, to perform the orientation test [notest]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('notest', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('averagemol', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Output average molecular coordinates for central molecule (y or no) [y]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('y', 256)

    symbols.nsymbol = iref
    symbols.nsymboldata = iref
    symbols.nsymboluptodata = iref + 1
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata
    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1
    #
    # now write the defaults to symbol_value
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]

    #
    # set the number of files - not used in SHARM
    #
    symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('Setup sharm>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')

    #
    # set up the data prompt - not used in SHARM
    #
    symbols.dataprompt = fstr(' ', 80)
    symbols.dataprompttext = fstr(' ', 80)
    symbols.ilendataprompt = 0

    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.SHARM.dat', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_sdf():
    """setup_symbols_generic.f:setup_symbols_sdf -- keyword table and defaults.

    sets the symbols and their default values for input to SDF.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to SDF.
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if iflag=0
    #

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnameato', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of .ato file ', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nr', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Number of radius values (max 200)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('130', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum radius for spherical harmonic coefficients', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('13', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsumt', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of configurations already accumulated', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ncoeffs', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of coefficients (program calculates this)', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('l1values', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('L1 values (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('l2values', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('L2 values (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('lvalues', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('L values (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('n1step', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Step in N1 values', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('n2step', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Step in N2 values', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom-c', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('List of centre atom types. First will be origin.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('pairdist-c', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('List of pair distances for centre molecule.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('fracdist-c', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Acceptance half-width (fraction of pair distance).', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('axisc1', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('First axis definition for centre molecule.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('axisc2', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Second axis definition for centre molecule.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom-s', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('List of second atom types. First will be origin.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('pairdist-s', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('List of pair distances for second molecule.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('fracdist-s', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Acceptance half-width (fraction of pair distance).', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('axiss1', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('First axis definition for second molecule.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('axiss2', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Second axis definition for second molecule.', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('test', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Set phi1,theta1,chi1,phi2,theta2,phi2,thetal,phil, to perform the orientation test [notest]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('notest', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('averagemol', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Output average molecular coordinates for central molecule (y or no) [y]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('y', 256)

    symbols.nsymbol = iref
    symbols.nsymboldata = iref
    symbols.nsymboluptodata = iref + 1
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata
    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1
    #
    # now write the defaults to symbol_value
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]

    #
    # set the number of files - not used in SHARM
    #
    symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('Setup sdf>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')

    #
    # set up the data prompt - not used in SDF
    #
    symbols.dataprompt = fstr(' ', 80)
    symbols.dataprompttext = fstr(' ', 80)
    symbols.ilendataprompt = 0
    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.SDF.dat', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_sdfcube():
    """setup_symbols_generic.f:setup_symbols_sdfcube -- keyword table and defaults.

    sets the symbols and their default values for input to SHARM.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to SHARM.
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if iflag=0
    #

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnameato', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of .ato file ', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('mpixel', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of pixels along (half of) one axis', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('20', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum distance for spatial density function', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('10', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsumt', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of configurations already accumulated', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom-c', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Central molecule - list of centre atom types', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('axisc1', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('First axis definition for central molecule', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('axisc2', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Second axis definition for central molecule', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom-s', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Second molecule - list of centre atom types', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('axiss1', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('First axis definition for second molecule', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('axiss2', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Second axis definition for second molecule', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('averagemol', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Output average molecular coordinates for central molecule (y or no) [y]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('y', 256)

    symbols.nsymbol = iref
    symbols.nsymboldata = iref
    symbols.nsymboluptodata = iref + 1
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata
    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1
    #
    # now write the defaults to symbol_value
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]

    #
    # set the number of files - not used in SHARM
    #
    symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('Setup sdfcube>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')

    #
    # set up the data prompt - not used in SHARM
    #
    symbols.dataprompt = fstr(' ', 80)
    symbols.dataprompttext = fstr(' ', 80)
    symbols.ilendataprompt = 0

    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.SDFCUBE.dat', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_mapgrentropy():
    """setup_symbols_generic.f:setup_symbols_mapgrentropy -- keyword table and defaults.

    sets the symbols and their default values for input to clusters.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to clusters.
    #
    #
    # set up the symbol table for the plot3d input file
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if iflag=0
    #
    # define symbols and default values
    #
    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('shcoeffs', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of file containing spherical harmonic coefficients', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ncoeffs', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('no. of coefficients - ' + 'determined from coefficients file', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr(' ', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ident', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('= 0 for identical molecules,' + ' else 1 if different', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr(' ', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsmoo', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('rmin and number of smoothings on coefficients on radial, spatial and orientational coeffs. [0 0 0 0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0 0 0 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('mirrorsym', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Apply mirror symmetry (y or no)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('no', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ntrials', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of random trials [10000]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('10000', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nfreq', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Write frequency', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('100', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('molndens', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Molecule no. density (mols per A**3) [0.0334]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0334', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('molmass', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Molecule mass (Da) [18]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('18', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('moltemp', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Temperature of simulation (K) [298]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('298', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('molrot', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Molecule rotation factor [water = 44921.624]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('44921.624', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('molvib', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Molecule vibration factor [water = 0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    # Last symbol is the extension of the file containing the spherical harmonic coefficients
    # Used only for setup, not in the plotting program itself.

    iref = iref + 1
    symbols.symbol[iref] = fstr('extcoeff', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('extension of coefficients file', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('.SHARM.h01', 256)

    symbols.nsymbol = iref
    symbols.nsymboluptodata = iref
    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)

    #
    # get the number of files in this input file
    #
    symbols.nfiles = 0
    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup mapgrentropy>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')
    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('none', 80)
    symbols.dataprompttext = fstr('jumps to nothing', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1

    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.mapgrentropy.dat', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_plot2d():
    """setup_symbols_generic.f:setup_symbols_plot2d -- keyword table and defaults.

    sets the symbols and their default values for input to clusters.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to clusters.
    #
    #
    # set up the symbol table for the plot2d input file
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if iflag=0
    #
    # define symbols and default values
    #
    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('bakcolour', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Background colour index [1]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('shcoeffs', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of file containing SHARM coefficients', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ncoeffs', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('no. of SHARM coefficients - ' + 'determined from coefficients file', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr(' ', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ident', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('= 0 for identical molecules,' + ' else 1 if different', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr(' ', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('setone', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('0 sets first coefficient to zero - normally 1', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsmoo', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('number of smoothings on coefficients', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('4', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('xdim', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('x dimension of plotting square', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('10', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ydim', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('y dimension of plotting square', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('10', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmin', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('minimum radius of plot', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ncircles', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of circles at centre of plot (max 25)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('xyzcircle', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('(x,y,z) coordinates for circle (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0 0 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('radcircle', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Radius of circle', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    iref = iref + 1
    symbols.symbol[iref] = fstr('nplotxy', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('no. of plots along x- and y-axis [set at 1 1]', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('1 1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('aspect', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('aspect ratio of plot [1.0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('use_l1_l2', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('use l1 and l2 (1 or 0)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('use_n1_n2', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('use n1 and n2 (1 or 0)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('use_m', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('use m2 (1 or 0)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nvary', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('vary thetal (1), phil (2), phim(3), thetam(4), chim(5)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ph_th_ch', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('phil, phim, thetam, chim (degrees - separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0 0 0 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('axespar', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('axes character size, line width and colour' + ' (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.5 2 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('plottitle', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('plot title', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr(' ', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('titlecoord', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('(x,y) coords. of title, and character size' + ' (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.1 -1.3 2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('colour', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Colour scale [2 - 8]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('3', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('zmin_zmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('range of intensites to be plotted' + ' (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1 5', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('contrast', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Contrast factor [0.0 - 1.0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('contover', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Contours overlaid? - set to N - cannot be changed', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('N', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('extraline', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('extra lines (0) - cannot be set', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('extratext', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('extra text (0) - cannot be set', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('0', 256)

    # Last symbol is the extension of the file containing the spherical harmonic coefficients
    # Used only for setup, not in the plotting program itself.

    iref = iref + 1
    symbols.symbol[iref] = fstr('extcoeff', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('extension of coefficients file', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('.SHARM.h01', 256)

    symbols.nsymbol = iref

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)

    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup plot2d>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')
    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('circle', 80)
    symbols.dataprompttext = fstr('jumps to the specified circle ', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1

    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.plot2d.txt', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_splot2d():
    """setup_symbols_generic.f:setup_symbols_splot2d -- keyword table and defaults.

    sets the symbols and their default values for input to clusters.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to clusters.
    #
    #
    # set up the symbol table for the plot2d input file
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if iflag=0
    #
    # define symbols and default values
    #
    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('shcoeffs', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of file containing SHARM coefficients', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ncoeffs', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('no. of SHARM coefficients - ' + 'determined from coefficients file', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr(' ', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ident', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('= 0 for identical molecules,' + ' else 1 if different', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsmoo', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('number of smoothings on coefficients', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('xdim', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('x dimension of plotting square', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('10', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ydim', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('y dimension of plotting square', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('10', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('mpixel', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of pixels in each direction [50]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('50', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmin_rmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('minimum radius of plot', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2.0 10', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('use_l1_l2', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('use l1, l2 and l (1 or 0) [1 0 1]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1 0 1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('use_n1_n2', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('use n1 and n2 (1 or 0)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('use_m', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('use m2 (1 or 0)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('mirrorsym', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Apply mirror symmetry (y or no) [no]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('no', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nvary', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('vary thetal (1), phil (2), phim(3), thetam(4), chim(5)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ph_th_ch', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('phil, phim, thetam, chim (degrees - separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0 0 0 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ncolour', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of colours needed for palette [1]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('colourfrac', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('colour fraction (0 - 1) [0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rgbvalues', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('RGB values for this colour fraction (separated by spaces) [1 1 1]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1 1 1', 256)

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    iref = iref + 1
    symbols.symbol[iref] = fstr('ncontour', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('no of contours [0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('wherecont', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Where to put contours (surface, base, both) [base]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('base', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('addcons', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Constant to add to final reconstruction [1.0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('surfonly', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Surface only? (y or no) [no]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('no', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('title', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Title of plot', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('xlabel', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('x-axis label', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('x [A]', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ylabel', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('y-axis label', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('y [A]', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('setgrid', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('set grid> (y or no) [no]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('no', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('smoothness', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('degree of interpolation along x and y (1 1 = no interpolation, 0 0 = optimal [0 0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('plotrotate', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Rotation of plot (deg) [0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    # Last symbol is the extension of the file containing the spherical harmonic coefficients
    # Used only for setup, not in the plotting program itself.

    iref = iref + 1
    symbols.symbol[iref] = fstr('extcoeff', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('extension of coefficients file', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('.SHARM.h01', 256)

    symbols.nsymbol = iref
    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)

    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1
    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup splot2d>', 80)
    symbols.ilenprompt = len(symbols.prompt)
    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('colour', 80)
    symbols.dataprompttext = fstr('jumps to the specified colour ', 80)
    symbols.ilendataprompt = len(symbols.dataprompt)

    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.splot2d.txt', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_plot3d():
    """setup_symbols_generic.f:setup_symbols_plot3d -- keyword table and defaults.

    sets the symbols and their default values for input to clusters.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to clusters.
    #
    #
    # set up the symbol table for the plot3d input file
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if iflag=0
    #
    # define symbols and default values
    #
    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('shcoeffs', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of file containing spherical harmonic coefficients', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ncoeffs', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('no. of coefficients - ' + 'determined from coefficients file', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr(' ', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ident', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('= 0 for identical molecules,' + ' else 1 if different', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr(' ', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('setone', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('0 sets first coefficient to zero - normally 1', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsmoo', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('number of smoothings on coefficients', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('4', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('radmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('maximum radius of plotting box', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('5', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nplotxy', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('no. of plots along x- and y-axis [set at 1 1]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1 1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('aspect', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('aspect ratio of plot [1.0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmin_rmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('minimum and maximum radius of plot', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2.0 5.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('surfra', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('fractional isosurface level (-ve for absolute)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.15', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('use_l1_l2', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('use l1 and l2 (1 or 0)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('use_n1_n2', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('use n1 and n2 (1 or 0)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('use_m', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('use m2 (1 or 0)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nvary', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('vary (thetal, phil) (1), (thetam, phim) (2), ' + '(thetam, chim) (3)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ph_th_ch', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('phim, thetam, chim (degrees - separated by spaces', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0 0 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsphere', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('number of spheres at centre of plot (max 25)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('radsphere', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('radius of sphere', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rthetaphi', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('(r,theta,phi) coordinates for sphere' + ' (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0 0 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rgbsph', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('red green blue fractions for sphere' + ' (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.7 0.7 0.7', 256)

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    iref = iref + 1
    symbols.symbol[iref] = fstr('axespar', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('axes character size, line width and colour' + ' (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.5 2 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('plottitle', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('plot title', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr(' ', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('titlecoord', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('(x,y) coords. of title, and character size' + ' (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.1 -1.3 2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('blank', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('A blank line required by plotting program ', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr(' ', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rgbbak', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('red green blue fractions for background' + ' (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.8 0.8 1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ishade', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('ishade (1-8): 0 means no shading,' + ' -ve means inverted shading', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('-1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rgbobj', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('red green blue fractions for object' + ' (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1 1 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('lightcoord', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('(x,y,z) coordinates for light source' + ' (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2 2 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('fadefc', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('fade factor (0 = no fading, 1=full fading)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('itrans', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('transparency of object (0=0%,1=25%,2=50%,3=75%)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('appearance', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('difuse, shine, polish and contrast', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0 0.0 1.0 1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rotelev', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('rotation and elevation of viewing point (deg.)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('15 35', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('extraline', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('extra lines (0) - cannot be set', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('extratext', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('extra text (0) - cannot be set', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('0', 256)

    # Last symbol is the extension of the file containing the spherical harmonic coefficients
    # Used only for setup, not in the plotting program itself.

    iref = iref + 1
    symbols.symbol[iref] = fstr('extcoeff', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('extension of coefficients file', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('.SHARM.h01', 256)

    symbols.nsymbol = iref

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)

    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup plot3d>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')
    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('sphere', 80)
    symbols.dataprompttext = fstr('jumps to the specified sphere ', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1

    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.plot3d.txt', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_plot3djmol():
    """setup_symbols_generic.f:setup_symbols_plot3djmol -- keyword table and defaults.

    sets the symbols and their default values for input to clusters.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to clusters.
    #
    #
    # set up the symbol table for the plot3djmol input file
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if iflag=0
    #
    # define symbols and default values
    #
    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('shcoeffs', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of file containing spherical harmonic coefficients', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ncoeffs', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('no. of coefficients - determined from coefficients file', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr(' ', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ident', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('0 for identical molecules, else 1 if different' + ' - determined from SHARM/SDF.dat file', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('setone', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Central molecule: show all central molecule atoms (0), SHARM/SDF construct (1), or spheres (2) [1]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsmoo', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('number of smoothings on coefficients [4]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('4', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('radmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('maximum radius of plotting box [5]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('5', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('mpixel', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('no. of pixels per axis [20]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('20', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmin_rmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('minimum and maximum radius of plot', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2.0 5.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('surfra', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('fractional isosurface level (-ve for absolute)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.15', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('use_l1_l2', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('use l1, l2 and l (1 or 0) [1 0 1]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1 0 1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('use_n1_n2', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('use n1 and n2 (1 or 0)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('use_m', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('use m2 (1 or 0)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('mirrorsym', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('apply mirror symmetry (y or no)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('no', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nvary', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('vary (thetal, phil) (1), (thetam, phim) (2), (thetam, chim) (3)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ph_th_ch', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('phim, thetam, chim (degrees - separated by spaces', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0 0 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsphere', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('number of spheres at centre of plot', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('chsymbol', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('chemical symbol', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('C', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rthetaphi', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('atomic position (r,theta,phi, separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0 0.0 0.0', 256)

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    iref = iref + 1
    symbols.symbol[iref] = fstr('rgbbak', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('red green blue fractions for background' + ' (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.8 0.8 1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rgbobj', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('red green blue fractions for surface' + ' (separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1 1 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('itrans', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('transparency of surface (0=opaque, 1=translucent, 2=transparent)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('imesh', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('mesh (0=nomesh,1=mesh,2=dots)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rotelev', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('rotation and elevation of viewing point (deg.)[0 0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('scalefact', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('scale factor for SDF or OCF [1.0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('tlucency', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('translucency of surface (if itrans=1) (0 - 1) [0.5]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.5', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rotatemol', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('rotate central molecule if nvary > 1 (y or no) [y]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('y', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('spin', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('spin picture? (x, y, z, or off) [off]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('off', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nextfile', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Filename to create another surface [<none>]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<none>', 256)

    # Last symbol is the extension of the file containing the spherical harmonic coefficients
    # Used only for setup, not in the plotting program itself.

    iref = iref + 1
    symbols.symbol[iref] = fstr('extcoeff', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('extension of coefficients file', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('.SHARM.h01', 256)

    symbols.nsymbol = iref

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)

    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup plot3djmol>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')
    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('sphere', 80)
    symbols.dataprompttext = fstr('jumps to the specified sphere ', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1

    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.CUBE.txt', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_mapgr():
    """setup_symbols_generic.f:setup_symbols_mapgr -- keyword table and defaults.

    sets the symbols and their default values for input to mapgr
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to mapgr
    #
    #
    # set up the symbol table for the mapgr input file
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if iflag=0
    #
    # define symbols and default values
    #
    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('shcoeffs', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of file containing spherical harmonic coefficients', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ncoeffs', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('no. of coefficients - determined from coefficients file', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr(' ', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ident', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('0 for identical molecules, else 1 if different' + ' - determined from SHARM/SDF.dat file', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsmoo', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('rmin and number of smoothings on coefficients on radial, spatial and orientational coeffs. [0 0 0 0]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0 0 0 0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('mirrorsym', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('apply mirror symmetry (y or no)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('no', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nviews', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('number of views (directions) to calculate g(r)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('grtype', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('1 = radial, 2 = spatial, 3 = full orientational', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('theta_phi', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('thetal, phil, phim, thetam, chim - separated by spaces)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0 0 0 0 0', 256)

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    # Last symbol is the extension of the file containing the spherical harmonic coefficients
    # Used only for setup, not in the plotting program itself.

    iref = iref + 1
    symbols.symbol[iref] = fstr('extcoeff', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('extension of coefficients file', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('.SHARM.h01', 256)

    symbols.nsymbol = iref

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)

    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup mapgr>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')
    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('view', 80)
    symbols.dataprompttext = fstr('jumps to the specified view ', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1

    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.mapgr.txt', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_triangles():
    """setup_symbols_generic.f:setup_symbols_triangles -- keyword table and defaults.

    sets the symbols and their default values for input to tri-angles.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to tri-angles.
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if ireset >= iflag
    #

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnameato', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of .ato file ', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsize', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of theta bins [180]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('180', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsumt', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of configurations already accumulated', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('-1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom-c', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Centre atom type [*]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('*', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom-z', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Centre atom type for rho-z [*]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('*', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom-r', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Centre atom type for rho-r [*]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('*', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ndist', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of tri-angle distributions', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom1', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom type 1 to define a triangle', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom2', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom type 2 to define a triangle', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom3', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom type 3 to define a triangle', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ltype', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Single atom3 types [1], or multiple atom3 types [2]?', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmin12', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum distance for atoms 1 and 2', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax12', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum distance for atoms 1 and 2', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmin23', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum distance for atoms 2 and 3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax23', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum distance for atoms 2 and 3', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rminrmaxc', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum and maximum distance for atom 2 from atom c', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1 2', 256)

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    symbols.nsymbol = iref

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)

    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup triangles>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')
    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('triangles', 80)
    symbols.dataprompttext = fstr('jumps to the specified triangle ', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1
    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.TRI.dat', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_torangles():
    """setup_symbols_generic.f:setup_symbols_torangles -- keyword table and defaults.

    sets the symbols and their default values for input to torangles.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to torangles.
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if ireset >= iflag
    #

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnameato', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of .ato file ', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsize', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of cos(theta) bins [180]', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('180', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsumt', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of configurations already accumulated', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('-1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ndist', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of tri-angle distributions', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom1', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom type to define which molecule type', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('XXX', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ltyp1', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('First atom number to define the axis', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ltyp2', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Second atom number to define the axis', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ltyp3', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Two atom numbers for vector from atom 1', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('3 3', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ltyp4', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Two atom numbers for vector from atom 2', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('4 4', 256)

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    symbols.nsymbol = iref

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)

    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup torangles>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')
    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('torangles', 80)
    symbols.dataprompttext = fstr('jumps to the specified torangle ', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1
    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.TOR.dat', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_voids():
    """setup_symbols_generic.f:setup_symbols_voids -- keyword table and defaults.

    sets the symbols and their default values for input to voids.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to voids.
    #
    #
    # set up the symbol table for the VOIDS input file
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if ireset >= iflag
    #

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnameato', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of .ato file ', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsize', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('No. of histogram points (max 1000)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1000', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('npix', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('No. of pixels for void calculation (max 201)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('100', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('voidmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum void volume (in multiples of atomic volume)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('500.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsumt', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of configurations already accumulated', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('-1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ndist', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of void distributions', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('radius', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Average distance from atoms to define a void', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('include', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('List of atom types to include', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr(' ', 256)

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    symbols.nsymbol = iref

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)

    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup voids>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')
    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('void', 80)
    symbols.dataprompttext = fstr('Specifies a void ', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1
    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.VOIDS.dat', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_fluctuations():
    """setup_symbols_generic.f:setup_symbols_fluctuations -- keyword table and defaults.

    sets the symbols and their default values for input to fluctuations.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to fluctuations.
    #
    #
    # set up the symbol table for the FLUCTUATIONS input file
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if ireset >= iflag
    #

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnameato', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of .ato file ', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsize', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('No. of histogram points (max 1000)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1000', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsumt', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of samples already accumulated', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('-1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsamples', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of samples to take from this configuration', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1000', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ndist', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of volume sizes', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Radius of volume to probe', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('10', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ntotal', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Sum of N', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsqtotal', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Sum of N^2', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('szero', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('S(0)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    symbols.nsymbol = iref

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)

    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup fluctuations>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')

    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('size', 80)
    symbols.dataprompttext = fstr('Specifies a size to measure fluctuations ', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1
    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.FLUCTUATIONS.dat', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_nwts():
    """setup_symbols_generic.f:setup_symbols_nwts -- keyword table and defaults.

    sets the symbols and their default values for input to nwts.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to nwts.
    #
    # set up the symbol table for the NWTS input file
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if ireset >= iflag
    #

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnameato', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of .ato file ', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    # iref=iref+1
    # symbol(iref)='filename'
    # symbol_type(iref)='c'
    # symbol_text(iref)='Base name for wts file (no spaces)'
    # iflag(iref)=0
    # symbol_default(iref)='<undefined>'

    iref = iref + 1
    symbols.symbol[iref] = fstr('output', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Is the output to be per atom (1) or per molecule (2)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('nsamples', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of samples (1,2, or 3)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('xfrac', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Fraction of first sample in mixture', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('normtot', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Normalise totals to (0 = nothing, 1 = <b>^2, 2 = <b^2>', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('normdif', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Normalise first order difference (0 = no, 1 = yes)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ncomponent', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of components (determined from ato file)', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom label - set from .ato file', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atsymbol', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom symbol - set from .ato file', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('iexchange', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('1 if this atom exchanges with other atoms, 0 otherwise', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('abundances', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('List of mass numbers and abundances. 0 for natural isotope', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0 1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('isubs', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('This atom isotope substituted? (1 = yes, 0 = no)', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('abund2', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('List of mass numbers and abundances for second sample', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('0 1', 256)

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    symbols.nsymbol = iref

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)

    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup nwts>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')

    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('component', 80)
    symbols.dataprompttext = fstr('Specifies isotopic abundances of a component', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1
    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.NWTS.dat', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_xwts():
    """setup_symbols_generic.f:setup_symbols_xwts -- keyword table and defaults.

    sets the symbols and their default values for input to xwts.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to xwts.
    #
    # set up the symbol table for the XWTS input file
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if ireset >= iflag
    #

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnameato', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of .ato file ', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    # iref=iref+1
    # symbol(iref)='filename'
    # symbol_type(iref)='c'
    # symbol_text(iref)='Base name for wts file (no spaces)'
    # iflag(iref)=0
    # symbol_default(iref)='<undefined>'

    iref = iref + 1
    symbols.symbol[iref] = fstr('output', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Is the output to be per atom (1) or per molecule (2)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('normtot', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Normalise totals to (0 = nothing, 1 = <b>^2, 2 = <b^2>', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('2', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('iafformaff', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Use the IAFFs (1) or MAFFs (2) for normalisation?', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ncomponent', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of components (determined from ato file)', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom label - set from .ato file', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('atsymbol', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom symbol to specify form factor', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('maffq', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Modified atomic form factors alpha values (0)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('maffdelta', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Modified atomic form factors delta value (0)', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0', 256)

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    symbols.nsymbol = iref

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)

    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup xwts>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')

    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('component', 80)
    symbols.dataprompttext = fstr('Specifies form factors for a component', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1
    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.XWTS.dat', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return


@routine
def setup_symbols_bragginten():
    """setup_symbols_generic.f:setup_symbols_bragginten -- keyword table and defaults.

    sets the symbols and their default values for input to braggintensities.
    Returns nothing.
    """
    #
    # sets the symbols and their default values for input to braggintensities.
    #
    # set up the symbol table for the braggintensities input file
    #
    # iflag decides whether this parameter can be changed
    # a parameter can only be changed if ireset >= iflag
    #

    iref = 0

    iref = iref + 1
    symbols.symbol[iref] = fstr('fnameato', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Name of .ato file ', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('<undefined>', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('mplicities', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('No. of unit cells along a, b and c', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1 1 1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('hklqmin', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Minimum value of qhkl to be used', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('rmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Maximum value of r - used to create convolution function', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('qstepqmax', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Step in Q and maximum Q for output', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.001 1', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('qwidth', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Width of broadening function', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('0.005', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('ncomponent', 11)
    symbols.symbol_type[iref] = fstr('i', 1)
    symbols.symbol_text[iref] = fstr('Number of components (determined from ato file)', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('1', 256)

    symbols.nsymboluptodata = iref

    iref = iref + 1
    symbols.symbol[iref] = fstr('atom', 11)
    symbols.symbol_type[iref] = fstr('c', 1)
    symbols.symbol_text[iref] = fstr('Atom label - set from .ato file', 256)
    symbols.iflag[iref] = 1
    symbols.symbol_default[iref] = fstr('0', 256)

    iref = iref + 1
    symbols.symbol[iref] = fstr('scatlen', 11)
    symbols.symbol_type[iref] = fstr('r', 1)
    symbols.symbol_text[iref] = fstr('Neutron scattering length', 256)
    symbols.iflag[iref] = 0
    symbols.symbol_default[iref] = fstr('1.0', 256)

    symbols.nsymboldata = iref
    symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata

    symbols.nsymbol = iref

    #
    # initialise symbol_value to the default values
    #
    for is_ in range(1, symbols.nsymbol + 1):
        symbols.symbol_value[is_] = symbols.symbol_default[is_]
        units.write(6, fmt(_F100, is_, Ch(symbols.symbol[is_], 11), Ch(symbols.symbol_value[is_], 256)))
        # 100 format(1x,i3,1x,a11,1x,a30)

    #
    # get the number of files in this input file
    #
    ierr = 0
    try:
        symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:
        symbols.nfiles = 1

    #
    # setup the output prompt
    #
    symbols.prompt = fstr('setup bragginten>', 80)
    symbols.ilenprompt = findex(symbols.prompt, '>')

    #
    # set up the data prompt
    #
    symbols.dataprompt = fstr('component', 80)
    symbols.dataprompttext = fstr('Specifies neutron scattering length for a component', 80)
    symbols.ilendataprompt = findex_padded(symbols.dataprompt, ' ', 80) - 1
    #
    # setup the extension for the input file
    #
    symbols.extinp = fstr('.BRAGGINTEN.dat', 20)
    symbols.ilenextinp = findex_padded(symbols.extinp, ' ', 20) - 1

    return
