"""Translation of fortran_src/process_symbols.f: the keyword/command engine that reads,
checks and edits every symbol value (``.EPSR.inp`` and the other setup files), plus the
``changeato`` editing routines.

Routines (Fortran order): process_symbols, process_commands, check_symbol,
save_files_changeato, write_symbols_changeato, change_symbol_changeato, change_bond,
change_label, change_tether, shift_molecules, rotate_molecules, align_molecules,
rotate_mol, add_bond, remove_bond, remove_molecule_types, remove_molecule,
cleanup_types, add_whole_rot, add_image_atom, make_mol_files,
set_molecule_scale_factors.

Other library files are reached through ``_LazyModule`` proxies (see below): this file
and read_and_save_files.f call each other, and several of the callees are translated by
other work packages.  A proxy resolves the attribute on every access, so tracing and
``common.reset()`` (which replace module attributes) are honoured exactly as with a
plain module import.
"""

import importlib

from ..common import clustercom, commands, epsr16, epsr_dimension, iresetnumbers, lattice, molecom, potential_par, symbols
from ..fortran import (
    PI4ATAN,
    Ch,
    FortranEOF,
    FortranIOError,
    acosf,
    atan2f,
    cosd_intrinsic,
    cosf,
    f32,
    findex,
    findex_padded,
    fmt,
    fstr,
    ifix,
    int_array,
    lw,
    pad,
    powf,
    r4,
    read_internal,
    real_array,
    sind_intrinsic,
    sinf,
    sqrtf,
    substr,
    units,
)
from ..trace import routine
from . import parse


# --------------------------------------------------------------------------- helpers
class _LazyModule:
    """A module imported on first use (private helper, see the module docstring).

    ``names`` are candidate module names, tried in order: the attribute is taken from the
    first importable candidate that defines it (this covers the ato_manipulation_rev4
    part-file and the not yet fixed Python name of FCreadline_windows.f90).  Attribute
    assignment is forwarded to the first importable candidate (for data modules).
    """

    def __init__(self, *names):
        object.__setattr__(self, "_names", names)

    def _modules(self):
        mods = []
        for n in self._names:
            try:
                mods.append(importlib.import_module(n))
            except ModuleNotFoundError as e:
                if e.name != n:
                    raise
        return mods

    def __getattr__(self, attr):
        mods = self._modules()
        for m in mods:
            if hasattr(m, attr):
                return getattr(m, attr)
        raise AttributeError(f"{attr!r} not found in any of {self._names}")

    def __setattr__(self, attr, value):
        mods = self._modules()
        if not mods:
            raise ModuleNotFoundError(f"none of {self._names} can be imported")
        setattr(mods[0], attr, value)


_LIB = __name__.rsplit(".", 1)[0]  # epsr26.library
_COMMON = _LIB.rsplit(".", 1)[0] + ".common"

auxiliary_routines = _LazyModule(_LIB + ".auxiliary_routines")
ato_manipulation_rev4 = _LazyModule(_LIB + ".ato_manipulation_rev4", _LIB + ".ato_manipulation_rev4_part2")
ato_read_and_save = _LazyModule(_LIB + ".ato_read_and_save")
file_and_directory_commands = _LazyModule(_LIB + ".file_and_directory_commands")
fcreadline_windows = _LazyModule(_LIB + ".fcreadline_windows", _LIB + ".FCreadline_windows",
                                 _LIB + ".jmk_readline_module")
get_molecules = _LazyModule(_LIB + ".get_molecules")
lib_molecom_mod = _LazyModule(_LIB + ".molecom_mod")  # routines of module molecom_mod
pfroutines_mod = _LazyModule(_LIB + ".pfroutines_mod")  # routines of module pfroutines
plot_routines = _LazyModule(_LIB + ".plot_routines")
read_and_save_files = _LazyModule(_LIB + ".read_and_save_files")
molecom_mod = _LazyModule(_COMMON + ".molecom_mod")  # data of module molecom_mod


def _read_int(text, old):
    """``read(text,*,iostat=ierr) i`` -> (i, ierr); i keeps its old value on error."""
    try:
        return read_internal(text).int(old), 0
    except FortranEOF:
        return old, -1
    except FortranIOError as e:
        return old, e.iostat


def _read_real(text, old):
    """``read(text,*,iostat=ierr) x`` (real*4) -> (x, ierr)."""
    try:
        return read_internal(text).real(old), 0
    except FortranEOF:
        return old, -1
    except FortranIOError as e:
        return old, e.iostat


def _read_word(text, old, n):
    """``read(text,*,iostat=ierr) c`` (character*n) -> (c, ierr)."""
    try:
        w = read_internal(text).word(None)
        return (old if w is None else fstr(w, n)), 0
    except FortranEOF:
        return old, -1
    except FortranIOError as e:
        return old, e.iostat


def _read_line5(old, n):
    """``read(5,'(a)',iostat=ierr) s`` (character*n) -> (s, ierr); s unchanged on error."""
    try:
        return fstr(units.read_line(5), n), 0
    except FortranEOF:
        return old, -1
    except FortranIOError as e:
        return old, e.iostat


def _elem_or0(a, i):
    """``a(i)`` of a molecom_mod allocatable that may be unallocated or too short:
    0.0 then (the Fortran reads unallocated memory; Rust changeato.rs get_or0)."""
    if a is None or i < 1 or i >= a.shape[0]:
        return f32(0.0)
    return a[i]


def _parse_sometext(nspc, nchar=None, text=None):
    """``call parse(sometext,nwords,ncf,ncl,nspc,nstsize,nworddim)`` on the COMMON words."""
    if text is None:
        text = commands.sometext
    if nchar is None:
        nchar = commands.nstsize
    commands.nwords = parse.parse(text, commands.nwords, commands.ncf, commands.ncl, nspc, nchar,
                                  commands.nworddim)


def _word(k):
    """``sometext(ncf(k):ncl(k))``."""
    return substr(commands.sometext, commands.ncf[k], commands.ncl[k])


def _pr():
    """``prompt(1:ilenprompt)``."""
    return substr(symbols.prompt, 1, symbols.ilenprompt)


# --------------------------------------------------------------------------- process_symbols
_PS_F203 = "(/a,' ',a,' - ',a/)"
_PS_F401 = ("(/a,' Current data have not been saved.'"
            "//'Press <Enter> to save, or q to exit without saving: ')")


@routine
def process_symbols(nin, nout, ierrout, ignore):
    """process_symbols.f:process_symbols -- processes commands coming in on the input
    stream, nin.

    Diagnostics are written to the output stream nout. An 'e' on the input stream
    signifies the end of input. 'data' <number> signifies the following entries relate
    to dataset <number>. If no number is specified, then dataset 1 is assumed.
    Returns (nout, ierrout).
    """
    promptstring = ""  # character*256
    #
    # set an error flag so that calling routine can tell if an error has occurred
    #
    ierrout = 0
    command = "next"  # character*20
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
    while command != "e" and ierr == 0:
        if icount > 0 and icount <= symbols.nsymbol:
            #
            # get the reference value for the symbol
            #
            iref = file_and_directory_commands.ireff(icount, ifilecount)
            #
            # sort out the words of symbol_value
            #
            commands.sometext = fstr(symbols.symbol_value[iref], 3020)
            _parse_sometext(5)
            ilensymbol = findex_padded(symbols.symbol[icount], " ", 11) - 1
            units.write(nout, fmt(_PS_F203, _pr(), substr(symbols.symbol[icount], 1, ilensymbol),
                                  substr(symbols.symbol_text[icount], 1,
                                         max(1, findex_padded(symbols.symbol_text[icount], "  ", 256)))))
            promptstring = fstr(fmt('(a,": ",a," ? ")', substr(symbols.symbol[icount], 1, ilensymbol),
                                    substr(commands.sometext, commands.ncf[1],
                                           commands.ncl[max(1, commands.nwords)])), 256)
            # 202   format(/a,': ',a,' ? ')
        else:
            #
            # reset the file counter if you have got here
            #
            #           ifilecount=1
            units.write(6, lw())  # Blank line
            promptstring = fstr(fmt("(a,1x)", _pr()), 256)
            # 201   format(/a,' ')

        lenpromptstring = len(promptstring)
        commands.sometext = ""
        commands.sometext = fstr(fcreadline_windows.readline(commands.sometext,
                                                             substr(promptstring, 1, lenpromptstring)), 3020)
        #         read(nin,100,iostat=ierr) sometext
        _parse_sometext(11)
        #         write(nout,*) nwords,' ',sometext(ncf(1):ncl(nwords)
        nout, command, icount, ifilecount, ierrout = process_commands(nin, nout, command, icount, ifilecount,
                                                                      ierrout, ignore)
        #        write(6,*) command,isave

    if symbols.isave != 1:
        #
        # plot files are always saved without question
        #
        if symbols.symbol_setup != "plot":
            commands.nwords = 0
            commands.sometext = ""
            units.write(nout, fmt(_PS_F401, _pr()), advance=False)
            try:
                commands.sometext = fstr(units.read_line(nin), 3020)
                ierr = 0
            except FortranEOF:
                ierr = -1
            except FortranIOError as e:
                ierr = e.iostat
            _parse_sometext(5)
            if commands.nwords > 0:
                command = fstr(_word(1), 20)
            else:
                command = ""
            if command != "q":
                nout = read_and_save_files.save_files(nin, nout)
        else:
            nout = read_and_save_files.save_files(nin, nout)
    fcreadline_windows.save_history_list()
    return nout, ierrout


# --------------------------------------------------------------------------- process_commands
_PC_F204 = "(//a,' ',a,' - ',a/)"
_PC_F207 = "(/'List of parameters:-'//)"
_PC_F209 = "(a10,' - ',a65)"
_PC_F210 = "(/'Other commands:-'/)"
_PC_F211 = "(a20,'-',1x,a45)"
_PC_F208 = (
    "("
    "'ls or dir <filename> - Lists files in the current directory'"
    "/'cd <directory name>  - Changes the current directory'"
    "/'pwd                  - Show the current directory'"
    "/'save                 - Saves the current values to file'"
    "/'forcesave            - Saves the current values to file directly without checking if OK to do so'"
    "/'default              - Sets numerical values to their defaults'"
    "/'help or ?            - Show this list'"
    "/'<value>              - New value for the listed parameter'"
    "/'<parameter> <value>  - Specified parameter and its new value'"
    "/'n                    - Revert to the previous value'"
    "/'<Enter> or d            - Accept current parameter and go forwards'"
    "/'u                    - Accept current parameter and '"
    ",'go backwards'"
    "/'q                    - Return to EPSR shell without saving'"
    "/'bye or e             - Return to EPSR shell'"
    "/)"
)
_PC_F100 = "(/'process_commands> Command not recognised: ',a,': ',a)"


@routine
def process_commands(nin, nout, command, icount, ifilecount, ierrout, ignore):
    """process_symbols.f:process_commands -- processes the command in sometext (COMMON).

    Diagnostics are written to the output stream nout. An 'e' on the input stream
    signifies the end of input. 'data' <number> signifies the following entries relate
    to dataset <number>. If no number is specified, then dataset 1 is assumed.
    Returns (nout, command, icount, ifilecount, ierrout).
    """
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
        command = "d"
    if command == "d" or command == "u":
        if symbols.symbol_save != "":
            symbols.symbol_save = ""
        ifilecountold = ifilecount

        if command == "d":

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
        elif command == "u":

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
            commands.sometext = fstr(fmt("(i3)", ifilecount), 3020)  # 250 format(i3)
            _parse_sometext(5)
            if commands.nwords > 0:
                symbols.datatext = fstr(_word(1), 80)
                symbols.ilendatatext = int(commands.ncl[1] - commands.ncf[1] + 1)
            else:
                symbols.datatext = ""
                symbols.ilendatatext = 1
            #
            # and write it!
            #
            if nin == 5:
                units.write(nout, fmt(_PC_F204, substr(symbols.dataprompt, 1, symbols.ilendataprompt),
                                      substr(symbols.datatext, 1, symbols.ilendatatext),
                                      substr(symbols.dataprompttext, 1,
                                             max(1, findex_padded(symbols.dataprompttext, "  ", 80)))))
        command = "next"
        return nout, command, icount, ifilecount, ierrout
    #
    # perform any necessary translations
    #
    if command == "bye":
        command = "e"
    #
    # reflect the input text into the output stream
    #
    #     write(nout,*) sometext(ncf(1):ncl(nwords))
    #
    # first check whether the value of command is one of the special values
    #
    if command == "default":
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
        command = "next"
    elif command == "ls" or command == "dir":
        symbols.inotrecog = 0
        file_and_directory_commands.list_directory()
        command = "next"
    elif command == "cd" or command == "pwd":
        symbols.inotrecog = 0
        file_and_directory_commands.set_directory_name()
        command = "next"
    elif command == symbols.dataprompt:
        symbols.inotrecog = 0
        #
        # try to decide which data set you wish to access
        #
        ifound = 0
        ndata = 0  # NB: only read when set below
        if commands.nwords > 1:
            ndata, ierr = _read_int(_word(2), ndata)
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
        commands.sometext = fstr(fmt("(i3)", ifilecount), 3020)
        _parse_sometext(5)
        if commands.nwords > 0:
            symbols.datatext = fstr(_word(1), 80)
            symbols.ilendatatext = int(commands.ncl[1] - commands.ncf[1] + 1)
        else:
            symbols.datatext = ""
            symbols.ilendatatext = 1
        #
        # and write it!
        #
        if nin == 5:
            units.write(nout, fmt(_PC_F204, substr(symbols.dataprompt, 1, symbols.ilendataprompt),
                                  substr(symbols.datatext, 1, symbols.ilendatatext),
                                  substr(symbols.dataprompttext, 1,
                                         max(1, findex_padded(symbols.dataprompttext, "  ", 80)))))
        #
        # for plotting we need to process this command
        #
        if symbols.symbol_setup == "plot" and ignore == 0:
            iref = file_and_directory_commands.ireff(icount, ifilecount)
            # allow this command through the nochanges filter in check_symbol
            symbols.symbol_save = symbols.symbol_value[iref]
            command = "store"
        else:
            command = "next"
    elif command == "help" or command == "?":
        symbols.inotrecog = 0
        #
        # write out the list of commands
        #
        units.write(nout, fmt(_PC_F207))
        for is_ in range(1, symbols.nsymbol + 1):
            ilensymbol = findex_padded(symbols.symbol[is_], " ", 11) - 1  # noqa: F841
            units.write(nout, fmt(_PC_F209, Ch(symbols.symbol[is_], 11), substr(symbols.symbol_text[is_], 1, 65)))
        units.write(nout, fmt(_PC_F210))
        if symbols.ilendataprompt > 0:
            text = fstr(substr(symbols.dataprompt, 1, symbols.ilendataprompt) + " <number> ", 20)
            units.write(nout, fmt(_PC_F211, substr(text, 1, 20), Ch(symbols.dataprompttext, 80)))
        units.write(nout, fmt(_PC_F208))
        command = "next"
    elif command == "save":
        symbols.inotrecog = 0
        nout = read_and_save_files.save_files(nin, nout)
        symbols.isave = 1
        command = "next"
    elif command == "forcesave":
        symbols.inotrecog = 0
        # NB: the Fortran passes a third actual argument (.true.) that force_save_files
        # does not declare; it is ignored (see docs/SIGNATURE_EXCEPTIONS.md)
        nout = read_and_save_files.force_save_files(nin, nout)
        symbols.isave = 1
        command = "next"
    elif command == "q":
        symbols.inotrecog = 0
        symbols.isave = 1
        command = "e"
    elif command != "e" and command != "next" and command != "n":
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
        ilenfound = 0
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
                    2, commands.nwords, symbols.symbol_value[iref], ilenfound, commands.nwords,
                    commands.ncf, commands.ncl, commands.sometext)
                icount = ifound
                command = "store"
            else:

                # otherwise go to this parameter and ask for the new value
                # if no new value was specified we go to this new value. Normally we do nothing
                # at that point, unless these are plot symbols, when action with the current value
                # may need to be processed

                icount = ifound
                if symbols.symbol_setup == "plot":
                    command = "store"
                else:
                    command = "next"
            if ifound > symbols.nsymboldata:
                ifilecount = symbols.nfiles
            if ifound < symbols.nsymboluptodata + 1:
                ifilecount = 1
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
                    1, commands.nwords, symbols.symbol_value[iref], ilenfound, commands.nwords,
                    commands.ncf, commands.ncl, commands.sometext)
                command = "store"
            else:
                units.write(nout, fmt(_PC_F100, Ch(command, 20), " ", substr(commands.sometext, 1, 40)))
                symbols.inotrecog = symbols.inotrecog + 1
                ierrout = 1
                command = "next"
    #     write(6,*) command,isave
    #
    # if command is 'store' it means we need to check the new symbol_value
    # and set it back to the saved value if there is an error, or make any other changes
    # to the symbol values that are required. See check_symbol for details.
    #
    ierr = 0  # NB: always assigned by check_symbol before use
    if command == "store":

        #
        # signal that the new value will not have been written to the input file
        #
        symbols.inotrecog = 0
        symbols.isave = 0

        nout, icount, ierr = check_symbol(nin, nout, icount, ifilecount, ierr, ignore)

        #
        # set the value of ifilecount if icount does not refer to a dataset
        #
        #           if(icount.le.nsymboluptodata) ifilecount=1
        #           if(icount.gt.nsymboldata) ifilecount=nfiles
        if ierr == 0:
            command = "next"
        else:
            command = "n"
        #     write(6,*) iref,ierr,symbol_value(iref)(1:30),command

    if command == "n":

        #
        # change of value was not accepted
        #

        symbols.inotrecog = 0
        #
        # set up the correct reference integer depending on the value of ifilecount
        #
        if symbols.symbol_save != "" and icount > 0:
            #
            # the new value was not accepted so revert to the previous value
            #
            iref = file_and_directory_commands.ireff(icount, ifilecount)
            symbols.symbol_value[iref] = symbols.symbol_save
            #
            # for some modules this may mean we need to reset more than just the symbol_value
            #
            if symbols.symbol_setup == "changeato":
                nout, icount, ierr = check_symbol(nin, nout, icount, ifilecount, ierr, ignore)

        else:
            #
            # previous saved value was accepted - nothing to change
            #
            icount = icount + 1
        command = "next"

    return nout, command, icount, ifilecount, ierrout


# --------------------------------------------------------------------------- check_symbol
_CS_F200 = "(/a,1x,'Cannot change this value!')"
_CS_F201 = "(/a,1x,'Value ',a,' is not an integer - try again!')"
_CS_F202 = "(/a,1x,'Value ',a,' is not a real number - try again!')"
_CS_F203 = "(/a,1x,'ERROR - tried to generate too many symbol values')"
_CS_F205 = "(/a,1x,'WARNING! Atom type ',a,' not found in .ato file ',a)"
_CS_F501 = "(/a,1x,'WARNING! Number of bonds specified is smaller ','than expected. Should be: ',i3)"
_CS_F502 = "(/a,1x,'WARNING! No need for axis atoms in this case')"
_CS_F120 = "(/'Specify the extension to search on: ')"
_CS_F121 = "(/a,1x,'No file specified - value will not be changed'/)"
_CS_F206 = "(/a,1x,'No coefficients file specified! - routine will exit'/)"
_CS_F300 = ("(/'check_symbol> Too many potential coefficients defined.'"
            ",' Will reduce to maximum allowed (',i3,')')")
_CS_F252 = "(/'check_symbol> Need to specify mass numbers and abundances in pairs: ')"


@routine
def check_symbol(nin, nout, icount, ifilecount, ierr, ignore):
    """process_symbols.f:check_symbol -- checks a new symbol_value and makes any other
    changes to the symbol values (and program state) that the change requires.

    Returns (nout, icount, ierr).
    """
    ext = ""  # character*20
    text = ""  # character*256
    iaccept = ""  # character*1
    icheck = 0
    rcheck = f32(0.0)
    ierrout = 0  # NB: local of check_symbol, set by read_ato/read_pcof
    ilenfound = 0
    ilentext = 0
    ilen = 0
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
        units.write(nout, fmt(_CS_F200, _pr()))
        symbols.symbol_value[iref] = symbols.symbol_save
        return nout, icount, ierr
    #
    # now check that the specified value(s) is(are) of the correct type if it is not
    # a character variable
    #
    if symbols.symbol_type[icount] != "c":
        commands.sometext = fstr(symbols.symbol_value[iref], 3020)
        #     write(6,*) symbol_value(iref),sometext(1:50)
        _parse_sometext(11)
        #     write(6,*) symbol_value(iref),sometext(1:50),nwords,ncf(1),ncl(1)
        iw = 0
        if symbols.symbol_type[icount] == "i":
            while iw < commands.nwords and ierr == 0:
                iw = iw + 1
                icheck, ierr = _read_int(_word(iw), icheck)
                if ierr != 0:
                    units.write(6, fmt(_CS_F201, _pr(), _word(iw)))
                    symbols.symbol_value[iref] = symbols.symbol_save
                    ierr = 0
                    return nout, icount, ierr
        else:
            #     write(6,*) nwords,iw,sometext(ncf(1):ncl(1))
            while iw < commands.nwords and ierr == 0:
                iw = iw + 1
                rcheck, ierr = _read_real(_word(iw), rcheck)
                #     write(6,*) ierr,rcheck
                if ierr != 0:
                    units.write(nout, fmt(_CS_F202, _pr(), _word(iw)))
                    symbols.symbol_value[iref] = symbols.symbol_save
                    ierr = 0
                    return nout, icount, ierr
    #
    # if icount .eq. nsymboluptodata the number of files may have been
    # changed. This means the values of symbol_value above nsymboluptodata will
    # need to be resequenced
    #
    if icount == symbols.nsymboluptodata:
        #
        # read the new value of nfiles into a temporary location
        #
        nfilet, ierr = _read_int(symbols.symbol_value[iref], 0)

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
                        units.write(nout, fmt(_CS_F203, _pr()))
                        ierr = 1
                        return nout, icount, ierr
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
                        units.write(nout, fmt(_CS_F203, _pr()))
                        ierr = 1
                        return nout, icount, ierr
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
    # if the name of the ato file has changed, we need to read in the
    # new values
    #
    if symbols.symbol[icount] == "fnameato":
        #     write(6,*) 'Got to fnameato check_symbol'
        nout, ierrout = read_and_save_files.read_ato(nin, nout, icount, ifilecount, ierrout)

    #
    # if the name of the pcof file has changed, we need to read in the
    # new values
    #
    if symbols.symbol[icount] == "fnamepcof":
        #     write(6,*) 'Got to fnamepcof check_symbol'
        icountsave = icount
        nout, icount, ierrout = read_and_save_files.read_pcof(nin, nout, icount, ifilecount, ierrout)
        icount = icountsave

    #
    # ATOM CHECK
    #
    # If the first four characters of symbol are 'atom', it means the associated variable
    # is a list of one or more atom types. We need to check that these atom types actually
    # exist in the current .ato file.
    #
    if substr(symbols.symbol[icount], 1, 4) == "atom":
        #
        # First, set up the list of atoms
        #
        commands.sometext = fstr(symbols.symbol_value[iref], 3020)
        _parse_sometext(5)

        # now check that the atom types actually exist in the current .ato file

        iw = 0
        while iw < commands.nwords:
            iw = iw + 1

            # Search for a * in the atom type: atom types will only be matched up to the *

            ncstar = findex(_word(iw), "*")
            if ncstar > 0:
                ncsearch = ncstar - 1
                ncluse = commands.ncf[iw] + ncsearch - 1
            else:
                ncsearch = commands.ncl[iw] - commands.ncf[iw] + 1
                ncluse = commands.ncl[iw]
            #
            # step through the atom types and check that one exists with the same type
            # as that specified
            #
            it = 0
            ifound = 0
            while it < molecom.ntype and ifound == 0:
                it = it + 1

                # If number of characters being searched is .gt. 0, then check the atom type

                if ncsearch > 0:

                    # NB: type(it)(1:ncsearch) with ncsearch > 3 reads past the character*3
                    # element in the Fortran; blank padding is used here (as in the Rust)
                    if substr(molecom.type[it], 1, ncsearch) == substr(commands.sometext, commands.ncf[iw], ncluse):
                        ifound = it

                else:

                    # Otherwise always accept this atom type

                    ifound = it

            if ifound == 0:
                units.write(nout, fmt(_CS_F205, _pr(), _word(iw), substr(symbols.fname_ato, 1, symbols.ilenfname_ato)))
                symbols.symbol_value[iref] = symbols.symbol_save

    #
    # changes specific to SHARM
    #

    if symbols.symbol_setup == "sharm":

        # reformat the symbol_values (11) and (14) so there is only 1 space between words

        if iref == 11 or iref == 14:
            commands.sometext = fstr(symbols.symbol_value[iref], 3020)
            _parse_sometext(5)
            _, commands.nwords, symbols.symbol_value[iref], ilenfound = parse.concatenate(
                1, commands.nwords, symbols.symbol_value[iref], ilenfound, commands.nwords,
                commands.ncf, commands.ncl, commands.sometext)

        #
        # if centres is the same as centrec set the axis definition for the second molecule
        # to be the same as the central molecule
        #
        if iref >= 14 and symbols.symbol_value[11] == symbols.symbol_value[14]:
            # set axiss1 same as axisc1
            symbols.symbol_value[15] = symbols.symbol_value[12]
            # set axiss2 same as axisc2
            symbols.symbol_value[16] = symbols.symbol_value[13]

    #
    # changes specific to SDF
    #

    if symbols.symbol_setup == "sdf":

        # reformat the symbol_valueS (11) and (16) so there is only 1 space between words

        if iref == 11 or iref == 16:
            commands.sometext = fstr(symbols.symbol_value[iref], 3020)
            _parse_sometext(5)
            _, commands.nwords, symbols.symbol_value[iref], ilenfound = parse.concatenate(
                1, commands.nwords, symbols.symbol_value[iref], ilenfound, commands.nwords,
                commands.ncf, commands.ncl, commands.sometext)

        elif iref == 12:

            # Check that the correct number of bonds have been read in

            _parse_sometext(5, 256, symbols.symbol_value[11])
            nbondexpect = commands.nwords * (commands.nwords - 1) // 2
            _parse_sometext(5, 256, symbols.symbol_value[12])
            if commands.nwords < nbondexpect:
                units.write(6, fmt(_CS_F501, _pr(), nbondexpect))

        elif iref >= 16:

            # For second molecule if it only has one atom, ensure that no atoms are specified for
            # the coordinate axes

            _parse_sometext(5, 256, symbols.symbol_value[16])
            naxisexpect = commands.nwords
            _parse_sometext(5, 256, symbols.symbol_value[19])
            if naxisexpect == 1 and commands.nwords > 1:
                units.write(6, fmt(_CS_F502, _pr()))
                symbols.symbol_value[19] = fstr(substr(symbols.symbol_value[19], commands.ncf[1], commands.ncl[1]), 256)
            _parse_sometext(5, 256, symbols.symbol_value[20])
            if naxisexpect == 1 and commands.nwords > 1:
                units.write(6, fmt(_CS_F502, _pr()))
                symbols.symbol_value[20] = fstr(substr(symbols.symbol_value[20], commands.ncf[1], commands.ncl[1]), 256)

        elif iref == 17:

            # Check that the correct number of bonds have been read in

            _parse_sometext(5, 256, symbols.symbol_value[16])
            nbondexpect = commands.nwords * (commands.nwords - 1) // 2
            _parse_sometext(5, 256, symbols.symbol_value[17])
            if commands.nwords < nbondexpect:
                units.write(6, fmt(_CS_F501, _pr(), nbondexpect))

        if symbols.symbol_value[11] == symbols.symbol_value[16]:

            # if atom-s is the same as atom-c set the axis definition for the second molecule
            # to be the same as the central molecule

            # set bond distances to be the same for both molecules
            symbols.symbol_value[17] = symbols.symbol_value[12]
            # set acceptance fractions to be the same
            symbols.symbol_value[18] = symbols.symbol_value[13]
            # set axiss1 same as axisc1
            symbols.symbol_value[19] = symbols.symbol_value[14]
            # set axiss2 same as axisc2
            symbols.symbol_value[20] = symbols.symbol_value[15]

    # if requested search for a suitable data or wts file

    if ((symbols.symbol[icount] == "datafile" or symbols.symbol[icount] == "wtsfile")
            and symbols.symbol_value[iref] == "search"):
        if symbols.symbol[icount] == "datafile":
            commands.nwords = 0
            while commands.nwords <= 0:
                #
                # specify the extension to search on
                #
                units.write(nout, fmt(_CS_F120), advance=False)
                commands.sometext = fstr(units.read_line(nin), 3020)  # 100 format(a)
                _parse_sometext(5)
            _, _, text, ilentext = parse.concatenate(1, 1, text, ilentext, commands.nwords, commands.ncf,
                                                     commands.ncl, commands.sometext)
            if ilentext <= 20:
                ext = fstr(substr(text, 1, ilentext), 20)
                ilenext = ilentext
            else:
                ext = ""
                ilenext = 1
        else:
            ext = ".wts"
            ilenext = 4
        iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)
        if iaccept != "y":
            units.write(nout, fmt(_CS_F121, _pr()))
            symbols.symbol_value[iref] = symbols.symbol_save
            ierr = 0
        else:
            _parse_sometext(5)
            _, commands.nwords, symbols.symbol_value[iref], ilen = parse.concatenate(
                1, commands.nwords, symbols.symbol_value[iref], ilen, commands.nwords,
                commands.ncf, commands.ncl, commands.sometext)
            ierr = 0
    # Changes specific to mapgrentropy
    if symbols.symbol_setup == "mapgrentropy" or symbols.symbol_setup == "mapgr":

        if symbols.symbol[icount] == "shcoeffs":
            #
            # changing the name of the coefficients file - read the corresponding .dat file
            #
            text = fstr(symbols.symbol_value[file_and_directory_commands.ireff(symbols.nsymbol, symbols.nfiles)], 256)
            commands.nwords = 0
            _parse_sometext(5, 256, text)
            ext = fstr(substr(text, commands.ncf[1], commands.ncl[1]), 20)
            ilenext = int(commands.ncl[1] - commands.ncf[1] + 1)
            commands.nwords = 0
            #
            # if symbol_value(1) is 'search' we need to find a valid SHARM coefficients file
            #
            if symbols.symbol_value[iref] == "search" or symbols.symbol_value[iref] == "<undefined>":
                units.write(6, fmt("(a)", substr(ext, 1, ilenext)))
                iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)
                _parse_sometext(5)
            else:
                commands.sometext = fstr(symbols.symbol_value[iref], 3020)
                _parse_sometext(5)
                #
                # add the extension if not supplied
                #
                ilendot = findex(_word(1), ".")
                if ilendot <= 0:
                    commands.sometext = fstr(_word(1) + pad(ext, 20), 3020)
                    commands.ncl[1] = commands.ncl[1] + ilenext
            if commands.nwords <= 0:
                units.write(nout, fmt(_CS_F206, _pr()))
                symbols.symbol_value[iref] = symbols.symbol_save
                ierr = 0
                return nout, icount, ierr
            else:
                symbols.symbol_value[iref] = fstr(_word(1), 256)
            nout, ierr = plot_routines.find_plot_par(12, nout, ierr)
    # Other changes specific to mapgr
    if symbols.symbol_setup == "mapgr":
        #
        # if nvary has changed mapgr needs to change symbol_text for all the files.
        #
        if symbols.symbol[icount] == "nvary":
            #
            # read the new value of nvary
            #
            iref = file_and_directory_commands.ireff(icount, ifilecount)
            nvary, ierr1 = _read_int(symbols.symbol_save, 0)
            nvaryt, ierr2 = _read_int(symbols.symbol_value[iref], 0)
            if ierr1 == 0 and ierr2 == 0 and nvaryt != nvary:
                for i in range(1, symbols.nfiles + 1):
                    iref = file_and_directory_commands.ireff(symbols.nsymboluptodata, i)
                    if nvaryt == 1:
                        iref = iref + 1
                        symbols.symbol_text[iref] = "thetal and phil of view (degrees - separated by spaces)"
                        iref = iref + 1
                        symbols.symbol_text[iref] = "phim, thetam, chim " + "(degrees - separated by spaces)"
                    elif nvaryt == 2:
                        iref = iref + 1
                        symbols.symbol_text[iref] = "thetam and phim of view (degrees - separated by spaces)"
                        iref = iref + 1
                        symbols.symbol_text[iref] = "thetal, phil, chim " + "(degrees - separated by spaces)"
                    elif nvaryt == 3:
                        iref = iref + 1
                        symbols.symbol_text[iref] = "thetam and chim of view (degrees - separated by spaces)"
                        iref = iref + 1
                        symbols.symbol_text[iref] = "thetal, phil, phim " + "(degrees - separated by spaces)"
    #
    # if nvary has changed (plot2d, splot2d, plot3d) need to change symbol_text and symbol_value
    #
    #
    # changes specific to EPSR
    #
    if symbols.symbol_setup == "epsr":

        #
        # if a parameter with a negative iflag is changed, we need to set iinit or ireset to 1
        # (but not when reading from .pcof file - otherwise ireset will always be set to 1!)
        #
        irefiinit = file_and_directory_commands.ireff(iresetnumbers.iinitref, ifilecount)
        irefireset = file_and_directory_commands.ireff(iresetnumbers.iresetref, ifilecount)
        if ignore == 0:
            if symbols.iflag[iref] < 0:
                symbols.symbol_value[irefiinit] = "1"
            if symbols.iflag[iref] == -2:
                symbols.symbol_value[irefireset] = "1"
            elif symbols.iflag[iref] < -2:
                symbols.symbol_value[irefireset] = "2"
        #
        # if ireset has been set, need to set iinit to 1 as well
        #
        if symbols.symbol_value[irefireset] != "0":
            symbols.symbol_value[irefiinit] = "1"

        # Ensure the number of potential coefficients does not exceed allowed dimension

        if symbols.symbol[icount] == "ncoeffp":

            irefncoeffp = file_and_directory_commands.ireff(iresetnumbers.incoeffp, ifilecount)
            potential_par.ncoeffp, ierr = _read_int(symbols.symbol_value[irefncoeffp], potential_par.ncoeffp)

            # Make sure this number does not exceed the allowed dimension

            if potential_par.ncoeffp > potential_par.mcoeff:

                units.write(6, fmt(_CS_F300, potential_par.mcoeff))
                potential_par.ncoeffp = potential_par.mcoeff
                commands.sometext = fstr(lw(potential_par.ncoeffp), 3020)
                ierr = 0
                _parse_sometext(5)
                _, commands.nwords, symbols.symbol_value[irefncoeffp], ilenfound = parse.concatenate(
                    1, commands.nwords, symbols.symbol_value[irefncoeffp], ilenfound, commands.nwords,
                    commands.ncf, commands.ncl, commands.sometext)

    # Changes specific to nwts

    if symbols.symbol_setup == "nwts":

        #         write(6,*) iref,icount,ifilecount

        # If fnameato has been set then get the number of components

        if symbols.symbol[icount] == "fnameato" and symbols.symbol_value[iref] != symbols.symbol_default[iref]:

            nfilet = molecom.ntype
            if nfilet < 10:
                symbols.symbol_value[symbols.nsymboluptodata] = fstr(fmt("(i1)", nfilet), 256)
            elif nfilet < 100:
                symbols.symbol_value[symbols.nsymboluptodata] = fstr(fmt("(i2)", nfilet), 256)
            else:
                symbols.symbol_value[symbols.nsymboluptodata] = fstr(fmt("(i3)", nfilet), 256)

            ifiletemp = symbols.nfiles
            symbols.nfiles = nfilet
            units.write(6, lw(nfilet, ifiletemp, symbols.nfiles))
            #
            # if new number of files .gt. old number of files
            # set the new values of symbol_value to the default values
            #
            while ifiletemp < symbols.nfiles:
                ifiletemp = ifiletemp + 1
                isymbol = symbols.nsymboluptodata
                units.write(6, lw(ifiletemp, isymbol))
                while isymbol < symbols.nsymboldata:
                    isymbol = isymbol + 1
                    iref = file_and_directory_commands.ireff(isymbol, ifiletemp)
                    symbols.symbol_value[iref] = symbols.symbol_default[isymbol]

            # Now set the specific values given by the .ato file

            il = symbols.nsymboluptodata
            if1 = 0
            #      write(6,*) il,if1
            while if1 < symbols.nfiles:
                if1 = if1 + 1
                #      write(6,*) if1,il,type(if1),' ',atsymbol(if1)
                symbols.symbol_value[il + 1] = molecom.type[if1]
                symbols.symbol_value[il + 2] = molecom.atsymbol[if1]
                if symbols.symbol_value[il + 2] != "H":
                    symbols.symbol_value[il + 3] = fstr(fmt("(i1)", 0), 256)
                il = il + symbols.nsymbolplus

        # Check that output is 1 or 2

        elif symbols.symbol[icount] == "output":

            itest = read_internal(symbols.symbol_value[iref]).int(0)
            if itest < 1:
                itest = 1
            if itest > 2:
                itest = 2
            symbols.symbol_value[iref] = fstr(fmt("(i1)", itest), 256)

        # Get the number of samples

        elif symbols.symbol[icount] == "nsamples":

            nsamples = read_internal(symbols.symbol_value[iref]).int(0)
            if nsamples > 3:
                nsamples = 3
            if nsamples < 1:
                nsamples = 1
            symbols.symbol_value[iref] = fstr(fmt("(i1)", nsamples), 256)

            # xfrac can only be set if number of samples is 3

            if nsamples == 3:

                symbols.iflag[icount + 1] = 0

            else:

                symbols.iflag[icount + 1] = 1
                symbols.symbol_value[iref + 1] = symbols.symbol_default[icount + 1]

            # nsamples is 2 or 3, then isotope substitution is allowed, otherwise not

            if nsamples < 2:

                # Set the flags on the individual components to prevent isotope substitution flags being set

                symbols.iflag[symbols.nsymboluptodata + 5] = 1
                symbols.iflag[symbols.nsymboluptodata + 6] = 1
                il = symbols.nsymboluptodata
                for if1 in range(1, symbols.nfiles + 1):

                    symbols.symbol_value[il + 5] = symbols.symbol_default[symbols.nsymboluptodata + 5]
                    symbols.symbol_value[il + 6] = symbols.symbol_default[symbols.nsymboluptodata + 6]
                    il = il + symbols.nsymbolplus

            else:

                symbols.iflag[symbols.nsymboluptodata + 5] = 0
                symbols.iflag[symbols.nsymboluptodata + 6] = 0

        # Make sure normtot has an allowed value

        elif symbols.symbol[icount] == "normtot":

            itest = read_internal(symbols.symbol_value[iref]).int(0)
            if itest < 0:
                itest = 0
            if itest > 2:
                itest = 2
            symbols.symbol_value[iref] = fstr(lw(itest), 256)

        # Make sure normdif has an allowed value

        elif symbols.symbol[icount] == "normdif":

            itest = read_internal(symbols.symbol_value[iref]).int(0)
            if itest < 0:
                itest = 0
            if itest > 1:
                itest = 1
            symbols.symbol_value[iref] = fstr(lw(itest), 256)

        # Only 'H' atoms are allowed to exchange

        elif symbols.symbol[icount] == "iexchange":

            if symbols.symbol_value[iref - 1] != "H":

                symbols.symbol_value[iref] = fstr(fmt("(i1)", 0), 256)

            itest = read_internal(symbols.symbol_value[iref]).int(0)
            if itest < 0:
                itest = 0
            if itest > 1:
                itest = 1
            symbols.symbol_value[iref] = fstr(lw(itest), 256)

        # For abundances need to ensure there are an even number entries

        elif symbols.symbol[icount] == "abundances" or symbols.symbol[icount] == "abund2":

            commands.sometext = fstr(symbols.symbol_value[iref], 3020)
            _parse_sometext(5)
            while commands.nwords != 2 * (commands.nwords // 2):
                units.write(6, fmt(_CS_F252), advance=False)
                commands.sometext = fstr(units.read_line(5), 3020)
                _parse_sometext(5)
            # NB: with nwords = 0 the Fortran reads ncl(0) (past the array start); index 0
            # (= 0) gives the same empty value
            symbols.symbol_value[iref] = fstr(substr(commands.sometext, commands.ncf[1],
                                                     commands.ncl[commands.nwords]), 256)

    # Changes specific to xwts

    if symbols.symbol_setup == "xwts":

        #         write(6,*) iref,icount,ifilecount

        # If fnameato has been set then get the number of components

        if symbols.symbol[icount] == "fnameato" and symbols.symbol_value[iref] != symbols.symbol_default[iref]:

            nfilet = molecom.ntype
            if nfilet < 10:
                symbols.symbol_value[symbols.nsymboluptodata] = fstr(fmt("(i1)", nfilet), 256)
            elif nfilet < 100:
                symbols.symbol_value[symbols.nsymboluptodata] = fstr(fmt("(i2)", nfilet), 256)
            else:
                symbols.symbol_value[symbols.nsymboluptodata] = fstr(fmt("(i3)", nfilet), 256)

            ifiletemp = symbols.nfiles
            symbols.nfiles = nfilet
            #      write(6,*) nfilet,ifiletemp,nfiles
            #
            # if new number of files .gt. old number of files
            # set the new values of symbol_value to the default values
            #
            while ifiletemp < symbols.nfiles:
                ifiletemp = ifiletemp + 1
                isymbol = symbols.nsymboluptodata
                #     write(6,*) ifiletemp,isymbol
                while isymbol < symbols.nsymboldata:
                    isymbol = isymbol + 1
                    iref = file_and_directory_commands.ireff(isymbol, ifiletemp)
                    symbols.symbol_value[iref] = symbols.symbol_default[isymbol]

            # Now set the specific values given by the .ato file

            il = symbols.nsymboluptodata
            if1 = 0
            #      write(6,*) il,if1
            while if1 < symbols.nfiles:
                if1 = if1 + 1
                #      write(6,*) if1,il,type(if1),' ',atsymbol(if1)
                symbols.symbol_value[il + 1] = molecom.type[if1]
                if symbols.symbol_value[il + 2] == symbols.symbol_default[symbols.nsymboluptodata + 2]:
                    symbols.symbol_value[il + 2] = molecom.atsymbol[if1]
                il = il + symbols.nsymbolplus

        # Check that output is 1 or 2

        elif symbols.symbol[icount] == "output":

            itest = read_internal(symbols.symbol_value[iref]).int(0)
            if itest < 1:
                itest = 1
            if itest > 2:
                itest = 2
            symbols.symbol_value[iref] = fstr(fmt("(i1)", itest), 256)

        # Normtot cannot be less than 0 or greater than 2

        elif symbols.symbol[icount] == "normtot":

            itest = read_internal(symbols.symbol_value[iref]).int(0)
            units.write(6, lw(itest))
            if itest < 0:
                itest = 0
            if itest > 2:
                itest = 2
            symbols.symbol_value[iref] = fstr(fmt("(i1)", itest), 256)

    #
    # changes specific to changeato
    #

    if symbols.symbol_setup == "changeato":

        ierr = change_symbol_changeato(icount, ifilecount, ierr)

    # Changes specific to bragginten

    # Changes specific to xwts

    if symbols.symbol_setup == "bragginten":

        #         write(6,*) iref,icount,ifilecount

        # If fnameato has been set then get the number of components

        if symbols.symbol[icount] == "fnameato" and symbols.symbol_value[iref] != symbols.symbol_default[iref]:

            nfilet = molecom.ntype
            if nfilet < 10:
                symbols.symbol_value[symbols.nsymboluptodata] = fstr(fmt("(i1)", nfilet), 256)
            elif nfilet < 100:
                symbols.symbol_value[symbols.nsymboluptodata] = fstr(fmt("(i2)", nfilet), 256)
            else:
                symbols.symbol_value[symbols.nsymboluptodata] = fstr(fmt("(i3)", nfilet), 256)

            ifiletemp = symbols.nfiles
            symbols.nfiles = nfilet
            #      write(6,*) nfilet,ifiletemp,nfiles
            #
            # if new number of files .gt. old number of files
            # set the new values of symbol_value to the default values
            #
            while ifiletemp < symbols.nfiles:
                ifiletemp = ifiletemp + 1
                isymbol = symbols.nsymboluptodata
                #     write(6,*) ifiletemp,isymbol
                while isymbol < symbols.nsymboldata:
                    isymbol = isymbol + 1
                    iref = file_and_directory_commands.ireff(isymbol, ifiletemp)
                    symbols.symbol_value[iref] = symbols.symbol_default[isymbol]

            # Now set the specific values given by the .ato file

            il = symbols.nsymboluptodata
            if1 = 0
            #      write(6,*) il,if1
            while if1 < symbols.nfiles:
                if1 = if1 + 1
                #      write(6,*) if1,il,type(if1),' ',atsymbol(if1)
                symbols.symbol_value[il + 1] = molecom.type[if1]
                # Q type atoms cannot scatter, therefore set scattering length for these atoms to zero
                if molecom.qtypea[if1] > 0:
                    symbols.symbol_value[il + 2] = "0.0"
                il = il + symbols.nsymbolplus

    return nout, icount, ierr


# --------------------------------------------------------------------------- changeato
@routine
def save_files_changeato(nout):
    """process_symbols.f:save_files_changeato -- closes unit nout and writes out the
    revised .ato file (moloutc)."""
    #
    # close the opened file as it is not needed
    #
    units.close(nout)
    #
    # write out the revised .ato file
    #
    ato_read_and_save.moloutc(symbols.fnameinp)
    return


@routine
def write_symbols_changeato(nfirst, nlast):
    """process_symbols.f:write_symbols_changeato -- writes the values read in from the
    .ato file to the appropriate symbol_value.

    Returns (nfirst, nlast).
    """
    #
    # ensure nfirst and nlast are within the allowed range
    #
    if nfirst < 1:
        nfirst = 1
    if nfirst <= symbols.nsymboldata:
        ifilecount = 1
    else:
        ifilecount = symbols.nfiles
    if nlast > symbols.nsymbol:
        nlast = symbols.nsymbol
    #
    # do all the symbols that can be defined from the .ato file - the others
    # are set to their default values
    #
    sv = symbols.symbol_value
    icount = nfirst - 1
    while icount <= nlast:
        if icount >= symbols.nsymboldata and ifilecount < symbols.nfiles:
            ifilecount = ifilecount + 1
            icount = symbols.nsymboluptodata
        icount = icount + 1
        iref = file_and_directory_commands.ireff(icount, ifilecount)
        s = symbols.symbol[icount]
        if s == "bond":
            sv[iref] = symbols.symbol_default[iref]
        elif s == "label":
            sv[iref] = symbols.symbol_default[iref]
        elif s == "rho":
            sv[iref] = fstr(lw(epsr16.rho), 256)
        elif s == "temp":
            sv[iref] = fstr(lw(molecom.simtk), 256)
        elif s == "tol":
            sv[iref] = fstr(lw(molecom.tolout), 256)
        elif s == "stepmi":
            sv[iref] = fstr(lw(molecom.stepmi), 256)
        elif s == "stepri":
            sv[iref] = fstr(lw(molecom.stepri), 256)
        elif s == "stepr":
            sv[iref] = fstr(lw(molecom.stepr), 256)
        elif s == "stepm":
            sv[iref] = fstr(lw(molecom.stepm), 256)
        elif s == "vibtemp":
            sv[iref] = fstr(lw(molecom.tempb), 256)
        elif s == "angtemp":
            sv[iref] = fstr(lw(molecom.tempangreq), 256)
        elif s == "dihtemp":
            sv[iref] = fstr(lw(molecom.tempdihreq), 256)
        elif s == "ecore":
            sv[iref] = fstr(lw(molecom.ecore), 256)
        elif s == "dcore":
            sv[iref] = fstr(lw(molecom.dcore), 256)
        elif s == "box":
            if lattice.cubic:
                sv[iref] = fstr(lw("cubic"), 256)
            else:
                sv[iref] = fstr(lw("noncubic"), 256)
        elif s == "addbond":
            sv[iref] = symbols.symbol_default[iref]
        elif s == "removebond":
            sv[iref] = symbols.symbol_default[iref]
        elif s == "natomtypes":
            sv[iref] = fstr(lw(molecom.ntype), 256)
            symbols.nfiles = molecom.ntype
        elif s == "atlabel":
            sv[iref] = fstr(lw(Ch(molecom.type[ifilecount], 3)), 256)
        elif s == "epsilon":
            sv[iref] = fstr(lw(molecom.epscor[ifilecount]), 256)
        elif s == "sigma":
            sv[iref] = fstr(lw(molecom.sigcor[ifilecount]), 256)
        elif s == "atmass":
            sv[iref] = fstr(lw(molecom.amassa[ifilecount]), 256)
        elif s == "charge":
            sv[iref] = fstr(lw(molecom.charge[ifilecount]), 256)
        elif s == "qradius":
            sv[iref] = fstr(lw(molecom.qradius[ifilecount]), 256)
        elif s == "atsymbol":
            sv[iref] = fstr(lw(Ch(molecom.atsymbol[ifilecount], 3)), 256)
        elif s == "iexchange":
            sv[iref] = fstr(lw(molecom.iexchange[ifilecount]), 256)
        elif s == "natomsph":
            sv[iref] = fstr(lw(_elem_or0(molecom_mod.natomsphere, ifilecount)), 256)
        elif s == "radminsph":
            sv[iref] = fstr(lw(_elem_or0(molecom_mod.radminsphere, ifilecount)), 256)
        elif s == "radmaxsph":
            sv[iref] = fstr(lw(_elem_or0(molecom_mod.radmaxsphere, ifilecount)), 256)
        elif s == "cgmolfac":
            sv[iref] = fstr(lw(_elem_or0(molecom_mod.pseudofacmol, ifilecount)), 256)

    return nfirst, nlast


_CC_F102 = "(/'changeato> Current a, b, and c axes: ',3(1x,f10.5))"
_CC_F100 = ("(/'changeato> Type a factor for each axis, or one factor for all axes,'"
            "/' (Use -ve for actual values, or press <Enter> to leave unchanged: ')")
_CC_F101 = "(/'changeato> Type new alpha, beta, gamma values or press <Enter> to leave unchanged: ')"
_CC_F209 = "(/'changeato> Did not recognise chemical symbol ',a3,' for atom type ',a3)"


@routine
def change_symbol_changeato(icount, ifilecount, ierr):
    """process_symbols.f:change_symbol_changeato -- resets the value of a specified
    parameter in the .ato file if requested by calling program.

    Returns ierr.
    """
    alpha = f32(0.0)
    beta = f32(0.0)
    gam = f32(0.0)
    newalpha = f32(0.0)
    newbeta = f32(0.0)
    newgam = f32(0.0)
    rtemp = f32(0.0)
    afac = f32(0.0)
    bfac = f32(0.0)
    cfac = f32(0.0)
    newphib = f32(0.0)
    newthetac = f32(0.0)
    newphic = f32(0.0)
    itemp = 0
    ctemp = ""  # character*3
    tokens = ""  # character(len=256)
    ntf = int_array(commands.nworddim)
    ntl = int_array(commands.nworddim)
    ntokens = 0

    def first_token():  # tokens(ntf(1):ntl(1))
        return substr(tokens, ntf[1], ntl[1])

    def parse_tokens():  # call parse(tokens,ntokens,ntf,ntl,11,nstsize,nworddim)
        return parse.parse(tokens, ntokens, ntf, ntl, 11, commands.nstsize, commands.nworddim)

    #
    # initial error signal
    #
    ierr = 0

    iref = file_and_directory_commands.ireff(icount, ifilecount)
    s = symbols.symbol[icount]

    if s == "bond":

        if symbols.symbol_value[iref] == "y":
            change_bond()
        symbols.symbol_value[iref] = symbols.symbol_save

    elif s == "label":

        if symbols.symbol_value[iref] == "y":
            change_label()
        symbols.symbol_value[iref] = symbols.symbol_save

    elif s == "tether":

        # Get the command and break into tokens
        tokens = fstr(symbols.symbol_value[iref], 256)
        ntokens = parse_tokens()
        if first_token() == "y":
            change_tether(tokens, ntokens, ntf, ntl)
        symbols.symbol_value[iref] = symbols.symbol_save

    elif s == "remove":

        if symbols.symbol_value[iref] == "y":
            remove_molecule_types()
        symbols.symbol_value[iref] = symbols.symbol_save

    elif s == "rho":

        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0 and rtemp > 0.0:
            rhoold = molecom_mod.natomrcg / (lattice.alatx * lattice.blaty * lattice.clatz)
            epsr16.rho = rtemp
            vol = rhoold / epsr16.rho
            dmaxfac = powf(vol, r4(1.0) / r4(3.0))
            afac = dmaxfac
            bfac = dmaxfac
            cfac = dmaxfac
            newphib = lattice.phib
            newthetac = lattice.thetac
            newphic = lattice.phic
            ato_manipulation_rev4.change_boxsize(afac, bfac, cfac, newphib, newthetac, newphic)
            #      write(6,*) rhoold,rho,dmaxfac,alat,blat,clat

    elif s == "temp":

        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom.simtk = rtemp

    elif s == "tol":

        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom.tolout = rtemp

    elif s == "stepmi":

        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom.stepmi = rtemp

    elif s == "stepri":

        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom.stepri = rtemp

    elif s == "stepr":

        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom.stepr = rtemp

    elif s == "stepm":

        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom.stepm = rtemp

    elif s == "vibtemp":

        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom.tempb = rtemp

    elif s == "angtemp":

        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom.tempangreq = rtemp

    elif s == "dihtemp":

        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom.tempdihreq = rtemp

    elif s == "ecore":

        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom.ecore = rtemp

    elif s == "dcore":

        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom.dcore = rtemp

    elif s == "addbond":

        if symbols.symbol_value[iref] == "y":
            add_bond()
        symbols.symbol_value[iref] = symbols.symbol_save

    elif s == "removebond":

        if symbols.symbol_value[iref] == "y":
            remove_bond()
        symbols.symbol_value[iref] = symbols.symbol_save

    elif s == "box":

        # Get the command and break into tokens
        tokens = fstr(symbols.symbol_value[iref], 256)
        ntokens = parse_tokens()
        # Only need to do anything if symbol_value is 'noncubic' or if cubic = .false. and symbol_value is 'cubic'
        # Note that to change shape of a cubic box you need to type 'noncubic'.
        # Otherwise there is no change in box shape, so nothing to do.
        if (lattice.cubic and first_token() == "noncubic") or first_token() == "change":
            # First deal with case of going from cubic box to non-cubic box, or wanting to change a non-cubic box.
            units.write(6, fmt(_CC_F102, lattice.alat, lattice.blat, lattice.clat))
            # Get the new axes
            if ntokens > 3:
                commands.sometext = fstr(substr(tokens, ntf[2], ntl[4]), 3020)
            else:
                units.write(6, fmt(_CC_F100), advance=False)
                commands.sometext = ""
                commands.sometext, ierr = _read_line5(commands.sometext, 3020)
            _parse_sometext(11)
            if commands.nwords == 1:
                afac = read_internal(_word(1)).real(afac)
                bfac = afac
                cfac = afac
            elif commands.nwords == 3:
                r = read_internal(substr(commands.sometext, commands.ncf[1], commands.ncl[3]))
                afac = r.real(afac)
                bfac = r.real(bfac)
                cfac = r.real(cfac)
            else:
                afac = f32(1.0)
                bfac = f32(1.0)
                cfac = f32(1.0)
            if afac < 0.0:
                afac = abs(afac) / lattice.alat
            if bfac < 0.0:
                bfac = abs(bfac) / lattice.blat
            if cfac < 0.0:
                cfac = abs(cfac) / lattice.clat

            #            write(6,'(a,3(1x,f10.5))') 'changeato> Current polar angles phib, thetac and phic: ',phib,thetac,phic
            alpha, beta, gam = auxiliary_routines.convertphibthetacphic(lattice.phib, lattice.thetac, lattice.phic,
                                                                         alpha, beta, gam)
            units.write(6, fmt("(a,3(1x,f10.5))", "changeato> Current crystallographic angles alpha, beta, gamma: ",
                               alpha, beta, gam))
            if ntokens > 6:
                commands.sometext = fstr(substr(tokens, ntf[5], ntl[7]), 3020)
            else:
                units.write(6, fmt(_CC_F101), advance=False)
                commands.sometext = ""
                commands.sometext, ierr = _read_line5(commands.sometext, 3020)
            _parse_sometext(11)
            if commands.nwords == 3:
                r = read_internal(substr(commands.sometext, commands.ncf[1], commands.ncl[3]))
                newalpha = r.real(newalpha)
                newbeta = r.real(newbeta)
                newgam = r.real(newgam)
                newphib, newthetac, newphic = auxiliary_routines.convertalphabetagamma(
                    newalpha, newbeta, newgam, newphib, newthetac, newphic)
                units.write(6, fmt("(a,3(1x,f10.5))", "changeato> New polar angles phib, thetac and phic: ",
                                   newphib, newthetac, newphic))
            else:
                newphib = lattice.phib
                newthetac = lattice.thetac
                newphic = lattice.phic

            ato_manipulation_rev4.change_boxsize(afac, bfac, cfac, newphib, newthetac, newphic)

        elif (not lattice.cubic) and first_token() == "cubic":

            # Now for the case when going from non-cubic to cubic

            # Get the volume of the current simulation box

            volold = lattice.alatx * lattice.blaty * lattice.clatz

            # Get the new axes

            lattice.alat = f32(1.0)
            lattice.blat = f32(1.0)
            lattice.clat = f32(1.0)
            lattice.phib = f32(90.0)
            lattice.thetac = f32(0.0)
            lattice.phic = f32(0.0)

            # set up lattice parameters for applying the minimum image convention
            # NB: cosd/sind are not declared EXTERNAL here, so gfortran's intrinsics are used

            lattice.alatx = lattice.alat
            lattice.blatx = lattice.blat * cosd_intrinsic(lattice.phib)
            lattice.blaty = lattice.blat * sind_intrinsic(lattice.phib)
            lattice.clatx = lattice.clat * sind_intrinsic(lattice.thetac) * cosd_intrinsic(lattice.phic)
            lattice.claty = lattice.clat * sind_intrinsic(lattice.thetac) * sind_intrinsic(lattice.phic)
            lattice.clatz = lattice.clat * cosd_intrinsic(lattice.thetac)

            # Get the volume of this new box

            vol = volold / (lattice.alatx * lattice.blaty * lattice.clatz)

            # Scale the box dimensions to give the save volume as previously

            dmaxfac = powf(vol, r4(1.0) / r4(3.0))
            lattice.alat = lattice.alat * dmaxfac
            lattice.alatx = lattice.alatx * dmaxfac
            lattice.blat = lattice.blat * dmaxfac
            lattice.blatx = lattice.blatx * dmaxfac
            lattice.blaty = lattice.blaty * dmaxfac
            lattice.clat = lattice.clat * dmaxfac
            lattice.clatx = lattice.clatx * dmaxfac
            lattice.claty = lattice.claty * dmaxfac
            lattice.clatz = lattice.clatz * dmaxfac
            #
            # define the Cartesian coordinates of the lattice vectors
            #
            lattice.alatc[1] = lattice.alatx
            lattice.alatc[2] = 0.0
            lattice.alatc[3] = 0.0
            lattice.blatc[1] = lattice.blatx
            lattice.blatc[2] = lattice.blaty
            lattice.blatc[3] = 0.0
            lattice.clatc[1] = lattice.clatx
            lattice.clatc[2] = lattice.claty
            lattice.clatc[3] = lattice.clatz
            lattice.bxax = lattice.blatx / lattice.alatx
            lattice.cxax = lattice.clatx / lattice.alatx
            lattice.cyby = lattice.claty / lattice.blaty

            #
            # calculate the perpendicular distance between the faces of the lattice
            #
            perpab = auxiliary_routines.perpd(lattice.alatc, lattice.blatc, lattice.clatc)
            perpbc = auxiliary_routines.perpd(lattice.blatc, lattice.clatc, lattice.alatc)
            perpca = auxiliary_routines.perpd(lattice.clatc, lattice.alatc, lattice.blatc)
            units.write(6, lw(perpab, perpbc, perpca))
            #
            # set up maximum radius to be used for Fourier transforms and to truncate
            # the potential. This is defined as the minimum perpendicular distance between
            # planes.
            #
            # G(r) and potentials will be calculated to rmax, or maximum allowed distance
            # based on number of points and step size, whichever is smaller.
            #
            clustercom.dmax = min(perpab, perpbc, perpca)
            units.write(6, lw(clustercom.nmol, clustercom.dmax))
            epsr16.cell = clustercom.dmax
            lattice.cubic = True

        # Calculate the new density and assign the value to the appropriate symbol

        epsr16.rho = f32(molecom.natomr) / (lattice.alatx * lattice.blaty * lattice.clatz)
        symbols.symbol_value[symbols.irefrho] = fstr(lw(epsr16.rho), 256)

        # Make sure a sensible value is left, in case the typed value was nonsense

        if lattice.cubic:
            symbols.symbol_value[iref] = fstr(lw("cubic"), 256)
        else:
            symbols.symbol_value[iref] = fstr(lw("noncubic"), 256)

    elif s == "findmols":

        if symbols.symbol_value[iref] == "y":
            get_molecules.get_molecules()
        symbols.symbol_value[iref] = symbols.symbol_save

    elif s == "getmolfs":

        if symbols.symbol_value[iref] == "y":
            make_mol_files()
        symbols.symbol_value[iref] = symbols.symbol_save

    elif s == "add_w_rot":

        # Get the command and break into tokens
        tokens = fstr(symbols.symbol_value[iref], 256)
        ntokens = parse_tokens()
        if first_token() == "y":
            add_whole_rot(tokens, ntokens, ntf, ntl)
        symbols.symbol_value[iref] = symbols.symbol_save

    elif s == "axisshift":

        # Get the command and break into tokens
        tokens = fstr(symbols.symbol_value[iref], 256)
        ntokens = parse_tokens()
        if first_token() == "y":
            shift_molecules(tokens, ntokens, ntf, ntl)
        symbols.symbol_value[iref] = symbols.symbol_save

    elif s == "axisrotate":

        # Get the command and break into tokens
        tokens = fstr(symbols.symbol_value[iref], 256)
        ntokens = parse_tokens()
        if first_token() == "y":
            rotate_molecules(tokens, ntokens, ntf, ntl)
        symbols.symbol_value[iref] = symbols.symbol_save

    elif s == "alignmole":

        # Get the command and break into tokens
        tokens = fstr(symbols.symbol_value[iref], 256)
        ntokens = parse_tokens()
        if first_token() == "y":
            align_molecules(tokens, ntokens, ntf, ntl)
        symbols.symbol_value[iref] = symbols.symbol_save

    elif s == "setimage":

        # Get the command and break into tokens
        tokens = fstr(symbols.symbol_value[iref], 256)
        ntokens = parse_tokens()
        if first_token() == "y":
            add_image_atom(True, tokens, ntokens, ntf, ntl)
        symbols.symbol_value[iref] = symbols.symbol_save

    elif s == "unsetimage":

        # Get the command and break into tokens
        tokens = fstr(symbols.symbol_value[iref], 256)
        ntokens = parse_tokens()
        if first_token() == "y":
            add_image_atom(False, tokens, ntokens, ntf, ntl)
        symbols.symbol_value[iref] = symbols.symbol_save

    elif s == "molscalefc":

        # Get the command and break into tokens
        tokens = fstr(symbols.symbol_value[iref], 256)
        ntokens = parse_tokens()
        if first_token() == "y":
            set_molecule_scale_factors(tokens, ntokens, ntf, ntl)
        symbols.symbol_value[iref] = symbols.symbol_save

    elif s == "atlabel":

        ctemp, ierr = _read_word(symbols.symbol_value[iref], ctemp, 3)
        if ierr == 0:
            #
            # step through all the atoms and make sure any with the old
            # type are changed to the new
            #
            for ia in range(1, molecom.natom + 1):
                if molecom.itype[ia] == ifilecount:
                    molecom.alabel[ia] = ctemp
            #
            # finally change the type itself
            #
            molecom.type[ifilecount] = ctemp

    elif s == "epsilon":
        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom.epscor[ifilecount] = rtemp

    elif s == "sigma":
        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom.sigcor[ifilecount] = rtemp

    elif s == "atmass":
        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom.amassa[ifilecount] = rtemp

    elif s == "charge":
        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom.charge[ifilecount] = rtemp

    elif s == "qradius":
        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom.qradius[ifilecount] = rtemp

    elif s == "atsymbol":
        ctemp, ierr = _read_word(symbols.symbol_value[iref], ctemp, 3)
        if ierr == 0:
            #
            # need to check that this chemical symbol exists in the Periodic Table
            #
            #
            # test to see if this chemical symbol exists in the table
            #
            massno = 0
            ifind = 0
            btest, ifind = ato_manipulation_rev4.scatt_len(ctemp, massno, ifind)  # noqa: F841
            if ifind == 0:
                #
                # reference integer for atom label
                #
                irefl = file_and_directory_commands.ireff(symbols.nsymboluptodata + 1, ifilecount)
                units.write(6, fmt(_CC_F209, Ch(ctemp, 3), Ch(symbols.symbol_value[irefl], 256)))
                ierr = 1
            else:
                molecom.atsymbol[ifilecount] = ctemp

    elif s == "iexchange":
        itemp, ierr = _read_int(symbols.symbol_value[iref], itemp)
        if ierr == 0:
            molecom.iexchange[ifilecount] = itemp

    elif s == "natomsph":
        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom_mod.natomsphere[ifilecount] = rtemp

    elif s == "radminsph":
        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom_mod.radminsphere[ifilecount] = rtemp

    elif s == "radmaxsph":
        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0:
            molecom_mod.radmaxsphere[ifilecount] = rtemp

    elif s == "cgmolfac":
        rtemp, ierr = _read_real(symbols.symbol_value[iref], rtemp)
        if ierr == 0 and rtemp >= 1.0:
            molecom_mod.pseudofacmol[ifilecount] = rtemp
        else:
            molecom_mod.pseudofacmol[ifilecount] = 1.0

    #
    # reset the value of the symbol if an error occurred on trying to set it
    #
    if ierr != 0:
        symbols.symbol_value[iref] = symbols.symbol_save
        symbols.symbol_save = ""
        ierr = 0

    return ierr


_CB_F31 = "(/1x,'Give 2 atom type numbers, 2 atom numbers, and new bond: ')"


@routine
def change_bond():
    """process_symbols.f:change_bond -- changes the length of specified bonds (answers
    read from unit 5 until the first atom type number is 0)."""
    # NB: the READ targets are uninitialised before the first read; 0 is used for an
    # item left unchanged by a null value
    ncha1 = 0
    ncha2 = 0
    nat1 = 0
    nat2 = 0
    bonnew = f32(0.0)
    units.write(6, fmt(_CB_F31), advance=False)
    r = units.list_read(5)
    ncha1 = r.int(ncha1)
    ncha2 = r.int(ncha2)
    nat1 = r.int(nat1)
    nat2 = r.int(nat2)
    bonnew = r.real(bonnew)
    while ncha1 != 0:
        for ia in range(1, molecom.natom + 1):
            mmm = molecom.mref[ia]
            na1 = ia + 1 - clustercom.naf[mmm]
            if molecom.nbond[ia] > 0 and (na1 == nat1 or na1 == nat2):
                for ib in range(1, molecom.nbond[ia] + 1):
                    na2 = molecom.iabond[ia, ib] + 1 - clustercom.naf[mmm]
                    if (molecom.itype[ia] == ncha1
                            and molecom.itype[molecom.iabond[ia, ib]] == ncha2
                            and na1 == nat1 and na2 == nat2):
                        molecom.tbond[ia, ib] = bonnew
                    if (molecom.itype[ia] == ncha2
                            and molecom.itype[molecom.iabond[ia, ib]] == ncha1
                            and na1 == nat2 and na2 == nat1):
                        molecom.tbond[ia, ib] = bonnew
        units.write(6, fmt(_CB_F31))
        r = units.list_read(5)
        ncha1 = r.int(ncha1)
        ncha2 = r.int(ncha2)
        nat1 = r.int(nat1)
        nat2 = r.int(nat2)
        bonnew = r.real(bonnew)
    return


_CL_F35 = "(/'change_label> The following atom types are currently defined:-')"
_CL_F36 = ("(/'change_label> Give atom type number'"
           ",' and molecule type number to change (0 0 to exit): ')")
_CL_F1784 = "(/'change_label> Warning! No atoms of type ',a3,' at atom number ',i4)"
_CL_F1788 = "(/'change_label> There are ',i7,' atoms of type ',a3,' in molecules of type ',a3)"
_CL_F37 = ("(/'change_label> Current label for this type is ',a3"
           "//,' Give new label (a3): ')")
_CL_F39 = "(/'change_label> Give fraction of these atoms to change (1.0 for all): ')"
_CL_F1811 = "(/'change_label> Warning: No valid atoms to change labels found!!'/)"
_CL_F1814 = "(/'change_label> Changing labels on ',i7,' atoms'/)"


@routine
def change_label():
    """process_symbols.f:change_label -- changes the label of a random fraction of the
    atoms of one type in one molecule type."""
    a1 = ""  # character*3
    ats = ""  # character*3
    findtype = int_array(epsr_dimension.mcomp)  # noqa: F841 (declared, unused)
    changeselect = None
    changethis = None
    units.write(6, fmt(_CL_F35))
    for it in range(1, molecom.ntype + 1):
        units.write(6, lw(it, Ch(molecom.type[it], 3)))
    units.write(6, fmt("(/a)", "change_label> The following molecule types are currently defined:-"))
    for it in range(1, clustercom.nmoltype + 1):
        units.write(6, fmt("(i5,1x,a3)", it, Ch(clustercom.moltype[it], 3)))
    units.write(6, fmt(_CL_F36), advance=False)
    # NB: it and iamol are uninitialised before the read; 0 is kept for a null value
    it = 0
    iamol = 0
    r = units.list_read(5)
    it = r.int(it)
    iamol = r.int(iamol)
    if it < 1 or it > molecom.ntype or iamol < 1 or iamol > clustercom.nmoltype:
        return
    # First find out how many atoms of this type exist in the file
    nthistype = 0
    for im in range(1, clustercom.nmol + 1):
        if clustercom.mtype[im] == iamol:
            for ia in range(clustercom.naf[im], clustercom.nal[im] + 1):
                if molecom.itype[ia] == it:
                    nthistype = nthistype + 1
    if nthistype <= 0:
        units.write(6, fmt(_CL_F1784, Ch(molecom.type[it], 3), iamol))
        return
    else:
        units.write(6, fmt(_CL_F1788, nthistype, Ch(molecom.type[it], 3), Ch(clustercom.moltype[iamol], 3)))
    # Grab the list of atoms to (possibly) change
    if changeselect is not None:
        changeselect = None
    changeselect = int_array(nthistype)
    nthistype = 0
    for im in range(1, clustercom.nmol + 1):
        if clustercom.mtype[im] == iamol:
            for ia in range(clustercom.naf[im], clustercom.nal[im] + 1):
                if molecom.itype[ia] == it:
                    nthistype = nthistype + 1
                    changeselect[nthistype] = ia
    units.write(6, fmt(_CL_F37, Ch(molecom.type[it], 3)), advance=False)
    a1 = fstr(units.read_line(5), 3)  # 38 format(a)
    units.write(6, fmt(_CL_F39), advance=False)
    fracchange = units.list_read(5).real(f32(0.0))
    notochange = min(nthistype, ifix(fracchange * f32(nthistype)))
    if changethis is not None:
        changethis = None
    changethis = int_array(notochange)
    # Finally pick out the atoms to change at random from this list, remembering to reduce the list once it is changed
    # Get the list of atoms to change
    ichange = 0
    icount = 0
    while ichange < notochange:
        icount = icount + 1
        if icount > nthistype:
            icount = 1
        # Pick atom at random
        if (auxiliary_routines.ran1() < fracchange or fracchange >= 1.0) and changeselect[icount] > 0:
            ichange = ichange + 1
            changethis[ichange] = changeselect[icount]
            changeselect[icount] = 0  # Set this to zero to indicate this atom has already been selected
    nchange = ichange
    if nchange <= 0:
        units.write(6, fmt(_CL_F1811))
        return
    else:
        units.write(6, fmt(_CL_F1814, nchange))
    ichange = 0
    while ichange < nchange:
        ichange = ichange + 1
        ia = int(changethis[ichange])
        im = int(molecom.mref[ia])
        iaref = ia - clustercom.naf[im] + 1  # noqa: F841
        epsit = molecom.epscor[it]
        sigit = molecom.sigcor[it]
        amass = molecom.amassa[it]
        qt = molecom.charge[it]
        rqt = molecom.qradius[it]
        ats = molecom.atsymbol[it]
        iexch = int(molecom.iexchange[it])
        molecom.alabel[ia] = a1
        natomsph = molecom_mod.natomsphere[it]
        rhosph = molecom_mod.rhosphere[it]
        radminsph = molecom_mod.radminsphere[it]
        radmaxsph = molecom_mod.radmaxsphere[it]
        pseudofacsph = molecom_mod.pseudofacmol[it]
        ato_manipulation_rev4.ftype(ia)
        ato_manipulation_rev4.ctype(a1, ats, iexch, epsit, sigit, amass, qt, rqt, natomsph, rhosph, radminsph,
                                    radmaxsph, pseudofacsph)
    # Clean up the type lists
    cleanup_types()
    #
    # now check that all symbol_values have been redefined correctly for any type changes
    # that have been made
    #
    write_symbols_changeato(symbols.nsymboluptodata, symbols.nsymbol)

    return


_CT_F100 = "(/1x,'molecule type number',2x,'label',2x,'tethered status (T or F)',2x,'tethered atom'/)"
_CT_F101 = "(10x,i4,10x,a3,5x,l1,5x,i5)"  # (unused)
_CT_F102 = ("(/'change_tether> Give type number, new status (T or F), and tethering atom (0 for C.M.)"
            " ( to quit).')")


@routine
def change_tether(tokens, ntokens, ntf, ntl):
    """process_symbols.f:change_tether -- tethers specified molecule types."""
    moltetheredatomread = 0

    # Make a list of all current molecule types and their tethered status

    units.write(6, fmt(_CT_F100))
    for it in range(1, clustercom.nmoltype + 1):
        units.write(6, lw(it, Ch(clustercom.moltype[it], 3), bool(clustercom.moltethered[it]),
                          clustercom.moltetheredatom[it]))
    if ntokens > 4:
        # NB: no parse here: the words (ncf, ncl) used below are those of the previous parse
        # of sometext (Fortran behaviour, reproduced)
        commands.sometext = fstr(substr(tokens, ntf[2], ntl[4]), 3020)
    else:
        units.write(6, fmt(_CT_F102))
        commands.sometext = fstr(units.read_line(5), 3020)
        _parse_sometext(5)
    # NB: moltyperead is uninitialised if the read fails; it is only used when ierr = 0
    moltyperead, ierr = _read_int(_word(1), 0)
    if ierr == 0 and commands.nwords > 1:
        moltetheredread = _word(2) == "T" or _word(2) == "t"
        if moltetheredread and commands.nwords > 2:
            moltetheredatomread = read_internal(_word(3)).int(moltetheredatomread)
        else:
            moltetheredatomread = 0
        # Apply the new values
        for im in range(1, clustercom.nmol + 1):
            if clustercom.mtype[im] == moltyperead:
                # If this molecule is not tethered and we are making it tethered we need to set the tethering coordinates
                if not clustercom.tethered[im] and moltetheredread:
                    clustercom.phix[im] = clustercom.xadd[im]
                    clustercom.phiy[im] = clustercom.yadd[im]
                    clustercom.phiz[im] = clustercom.zadd[im]
                clustercom.tethered[im] = moltetheredread
                clustercom.tetheredatom[im] = moltetheredatomread
                # If not tethered then the thethering coordinates should be set to zero
                if not clustercom.tethered[im]:
                    clustercom.phix[im] = 0.0
                    clustercom.phiy[im] = 0.0
                    clustercom.phiz[im] = 0.0
                    clustercom.tetheredatom[im] = 0
        clustercom.moltethered[moltyperead] = moltetheredread
        clustercom.moltetheredatom[moltyperead] = moltetheredatomread
    #
    # now check that all symbol_values have been redefined correctly for any type changes
    # that have been made
    #
    write_symbols_changeato(symbols.nsymboluptodata, symbols.nsymbol)
    return


_SM_F100 = "(/1x,'molecule type number',2x,'label'/)"
_SM_F101 = "(10x,i4,10x,a3)"
_SM_F102 = ("(/'shift_molecules> Give type number, axis to shift along (x, y, or z),'/"
            ",' and amount to shift. (q to quit): ')")


@routine
def shift_molecules(tokens, ntokens, ntf, ntl):
    """process_symbols.f:shift_molecules -- applies a constant shift along a specified
    axis to all molecules of a specified type."""
    axisread = ""  # character*1
    axisshift = f32(0.0)

    # Make a list of all current molecule types and their tethered status

    units.write(6, fmt(_SM_F100))
    for it in range(1, clustercom.nmoltype + 1):
        units.write(6, fmt(_SM_F101, it, Ch(clustercom.moltype[it], 3)))
    if ntokens > 3:
        commands.sometext = fstr(substr(tokens, ntf[2], ntl[4]), 3020)
    else:
        units.write(6, fmt(_SM_F102), advance=False)
        commands.sometext = fstr(units.read_line(5), 3020)
    _parse_sometext(5)
    moltyperead, ierr = _read_int(_word(1), 0)
    if ierr == 0 and commands.nwords > 2:
        axisread = fstr(_word(2), 1)
        axisshift = read_internal(_word(3)).real(axisshift)
        # Apply the new values
        for im in range(1, clustercom.nmol + 1):
            if clustercom.mtype[im] == moltyperead:
                if axisread == "x":
                    clustercom.xadd[im] = clustercom.xadd[im] + axisshift
                if axisread == "y":
                    clustercom.yadd[im] = clustercom.yadd[im] + axisshift
                if axisread == "z":
                    clustercom.zadd[im] = clustercom.zadd[im] + axisshift
                if clustercom.tethered[im]:
                    if axisread == "x":
                        clustercom.phix[im] = clustercom.phix[im] + axisshift
                    if axisread == "y":
                        clustercom.phiy[im] = clustercom.phiy[im] + axisshift
                    if axisread == "z":
                        clustercom.phiz[im] = clustercom.phiz[im] + axisshift
    #
    # now check that all symbol_values have been redefined correctly for any type changes
    # that have been made
    #
    write_symbols_changeato(symbols.nsymboluptodata, symbols.nsymbol)
    return


_RM_F100 = "(/'rotate_molecules> Number of molecules: ',i0/,'molecule types:- '/)"
_RM_F101 = "(10x,i4,10x,a3)"
_RM_F102 = ("(/'rotate_molecules> Give molecule number or type, axis to rotate about (x, y, or z),'/"
            ",' and amount to rotate [deg.]. (q to quit): ')")


@routine
def rotate_molecules(tokens, ntokens, ntf, ntl):
    """process_symbols.f:rotate_molecules -- applies a rotation of a molecule about a
    specified axis to all molecules of a specified type (or one molecule number)."""
    moltyperead = ""  # character*3
    axisread = ""  # character*1
    axisshift = f32(0.0)

    pi = PI4ATAN  # 4.0*atan(1.0)
    degconv = pi / 180.0

    # Make a list of all current molecule types and their tethered status

    units.write(6, fmt(_RM_F100, clustercom.nmol))
    for it in range(1, clustercom.nmoltype + 1):
        units.write(6, fmt(_RM_F101, it, Ch(clustercom.moltype[it], 3)))
    if ntokens > 3:
        commands.sometext = fstr(substr(tokens, ntf[2], ntl[4]), 3020)
    else:
        units.write(6, fmt(_RM_F102), advance=False)
        commands.sometext = fstr(units.read_line(5), 3020)
    _parse_sometext(5)
    molnumberread, ierr = _read_int(_word(1), 0)
    if ierr != 0:  # Assume this is a molecule type, not a number
        moltyperead = fstr(_word(1), 3)
        molnumberread = 0
    else:
        moltyperead = ""
    if ierr == 0 and commands.nwords > 2:
        axisread = fstr(_word(2), 1)
        axisshift = read_internal(_word(3)).real(axisshift)
        phi = degconv * axisshift
        cphi = cosf(phi)
        sphi = sinf(phi)
        xa = molecom.xa
        ya = molecom.ya
        za = molecom.za
        # Apply the new values
        for im in range(1, clustercom.nmol + 1):
            if clustercom.moltype[clustercom.mtype[im]] == moltyperead or im == molnumberread:
                if axisread == "z":
                    # rotate about z axis
                    for ia in range(clustercom.naf[im], clustercom.nal[im] + 1):
                        xasave = xa[ia] * cphi - ya[ia] * sphi
                        ya[ia] = xa[ia] * sphi + ya[ia] * cphi
                        xa[ia] = xasave
                elif axisread == "y":
                    # rotate about y-axis
                    for ia in range(clustercom.naf[im], clustercom.nal[im] + 1):
                        zasave = za[ia] * cphi - xa[ia] * sphi
                        xa[ia] = za[ia] * sphi + xa[ia] * cphi
                        za[ia] = zasave
                elif axisread == "x":
                    # rotate about x-axis
                    for ia in range(clustercom.naf[im], clustercom.nal[im] + 1):
                        yasave = ya[ia] * cphi - za[ia] * sphi
                        za[ia] = ya[ia] * sphi + za[ia] * cphi
                        ya[ia] = yasave
    #
    # now check that all symbol_values have been redefined correctly for any type changes
    # that have been made
    #
    write_symbols_changeato(symbols.nsymboluptodata, symbols.nsymbol)
    return


_AM_F100 = "(/1x,'molecule type number',2x,'label'/)"
_AM_F101 = "(10x,i4,10x,a3)"
_AM_F102 = ("(/'align_molecules> Give type number, axis to align to (x,y,z,a,b,c),'/"
            ",' and 2 atom numbers to define molecular axis (q to quit): ')")


@routine
def align_molecules(tokens, ntokens, ntf, ntl):
    """process_symbols.f:align_molecules -- rotates all molecules of a specified type to
    align with a specified axis."""
    axisread = ""  # character*1
    axisatom1 = 0
    axisatom2 = 0
    # NB: labaxis is uninitialised in the Fortran if the axis letter is not one of
    # x,y,z,a,b,c; zeros are used here
    labaxis = real_array(3)

    # Make a list of all current molecule types
    units.write(6, fmt(_AM_F100))
    for it in range(1, clustercom.nmoltype + 1):
        units.write(6, fmt(_AM_F101, it, Ch(clustercom.moltype[it], 3)))
    if ntokens > 4:
        commands.sometext = fstr(substr(tokens, ntf[2], ntl[5]), 3020)
    else:
        units.write(6, fmt(_AM_F102), advance=False)
        commands.sometext = fstr(units.read_line(5), 3020)
    _parse_sometext(5)
    moltyperead, ierr = _read_int(_word(1), 0)
    if ierr == 0 and commands.nwords > 3:
        axisread = fstr(_word(2), 1)
        r = read_internal(substr(commands.sometext, commands.ncf[3], commands.ncl[4]))
        axisatom1 = r.int(axisatom1)
        axisatom2 = r.int(axisatom2)
        # Set up lab axis vector
        if axisread == "x" or axisread == "a":
            labaxis[1] = 1.0
            labaxis[2] = 0.0
            labaxis[3] = 0.0
        elif axisread == "y":
            labaxis[1] = 0.0
            labaxis[2] = 1.0
            labaxis[3] = 0.0
        elif axisread == "z":
            labaxis[1] = 0.0
            labaxis[2] = 0.0
            labaxis[3] = 1.0
        elif axisread == "b":
            labaxis[1] = lattice.blatx
            labaxis[2] = lattice.blaty
            labaxis[3] = 0.0
        elif axisread == "c":
            labaxis[1] = lattice.clatx
            labaxis[2] = lattice.claty
            labaxis[3] = lattice.clatz
        # Ensure this vector is a unit vector
        auxiliary_routines.normvec(labaxis)
        # Rotate specified molecules
        for im in range(1, clustercom.nmol + 1):
            if clustercom.mtype[im] == moltyperead:
                rotate_mol(im, labaxis, axisatom1, axisatom2)
    #
    # now check that all symbol_values have been redefined correctly for any type changes
    # that have been made
    #
    write_symbols_changeato(symbols.nsymboluptodata, symbols.nsymbol)
    return


@routine
def rotate_mol(im, labaxis, axisatom1, axisatom2):
    """process_symbols.f:rotate_mol -- rotates a molecule so the axis going from
    axisatom1 to axisatom2 coincides with labaxis.

    Labaxis is converted to a unit vector if it is not so already.
    """
    molaxis = real_array(3)

    # pvector is the vector about which the rotation occurs

    pvector = real_array(3)

    # xvector is the vector from centre of rotation to an individual atom
    # xpvector is the component of this vector along the axis of rotation
    # xrvector is the component of this vector perpendicular to the axis of rotation
    # xrpvector is the vector perpendicular to both axis of rotation and xrvector

    xvector = real_array(3)
    xpvector = real_array(3)
    xrvector = real_array(3)
    xrpvector = real_array(3)

    pi = r4(3.141592653)  # noqa: F841
    auxiliary_routines.normvec(labaxis)
    # First check specified axis atoms exist in this molecule.
    # If not nothing is done.

    #      write(6,*) im,axisatom1,axisatom2
    if axisatom1 < 1 or axisatom1 > clustercom.nato[im]:
        return
    if axisatom2 < 1 or axisatom2 > clustercom.nato[im]:
        return
    if axisatom1 == axisatom2:
        return
    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za
    nax1 = axisatom1 + int(clustercom.naf[im]) - 1
    nax2 = axisatom2 + int(clustercom.naf[im]) - 1
    # Axis defined by axisatom1 and axisatom2
    molaxis[1] = xa[nax2] - xa[nax1]
    molaxis[2] = ya[nax2] - ya[nax1]
    molaxis[3] = za[nax2] - za[nax1]
    auxiliary_routines.normvec(molaxis)
    #      write(6,*) im,(molaxis(j),j=1,3),(labaxis(j),j=1,3)
    # Cosine of required rotation
    cphi = auxiliary_routines.scapro(labaxis, molaxis)
    sphi = sqrtf(r4(1.0) - cphi * cphi)
    units.write(6, lw(im, cphi, sphi))
    # Axis of rotation is at right angles to the plane defined by molaxis and lab axis
    auxiliary_routines.vecpro(molaxis, labaxis, pvector)
    rab2 = auxiliary_routines.scapro(pvector, pvector)
    if rab2 > 0:
        rab = sqrtf(rab2)
        # Unit vector along axis of rotation
        for i in range(1, 4):
            pvector[i] = pvector[i] / rab
        # step through atoms in group
        for ia in range(1, clustercom.nato[im] + 1):
            #            if(ia.ne.axisatom1.and.ia.ne.axisatom2) then
            nat = ia + int(clustercom.naf[im]) - 1
            # first calculate vector from centre of rotation to each atom in group
            xvector[1] = xa[nat] - xa[nax1]
            xvector[2] = ya[nat] - ya[nax1]
            xvector[3] = za[nat] - za[nax1]
            # Vector parallel to axis of rotation
            frac = auxiliary_routines.scapro(pvector, xvector)
            for i in range(1, 4):
                xpvector[i] = frac * pvector[i]
            # Vector perpendicular to axis of rotation
            for i in range(1, 4):
                xrvector[i] = xvector[i] - xpvector[i]
            # Vector perpendicular to plane defined by axis of rotation and xrvector
            auxiliary_routines.vecpro(pvector, xrvector, xrpvector)
            # Generate rotated coordinates. The vector component parallel to the axis of rotation
            # is not rotated
            xa[nat] = xa[nax1] + xpvector[1]
            ya[nat] = ya[nax1] + xpvector[2]
            za[nat] = za[nax1] + xpvector[3]
            xa[nat] = xa[nat] + xrvector[1] * cphi + xrpvector[1] * sphi
            ya[nat] = ya[nat] + xrvector[2] * cphi + xrpvector[2] * sphi
            za[nat] = za[nat] + xrvector[3] * cphi + xrpvector[3] * sphi
            #            end if
    return


_AB_F41 = "(/1x,'Give 2 atom type numbers, 2 atom numbers,',' and bond length: ')"


@routine
def add_bond():
    """process_symbols.f:add_bond -- adds bonds between specified atoms (answers read
    from unit 5 until the first atom type number is not positive)."""
    # NB: the READ targets are uninitialised before the first read; 0 is kept for a null value
    ncha1 = 0
    ncha2 = 0
    nat1 = 0
    nat2 = 0
    bonnew = f32(0.0)
    units.write(6, fmt(_AB_F41), advance=False)
    r = units.list_read(5)
    ncha1 = r.int(ncha1)
    ncha2 = r.int(ncha2)
    nat1 = r.int(nat1)
    nat2 = r.int(nat2)
    bonnew = r.real(bonnew)
    while ncha1 > 0:
        for ia in range(1, molecom.natom + 1):
            mmm = molecom.mref[ia]
            na1 = ia + 1 - clustercom.naf[mmm]
            if molecom.itype[ia] == ncha1 and na1 == nat1:
                molecom.nbond[ia] = molecom.nbond[ia] + 1
                ib = molecom.nbond[ia]
                iref = clustercom.naf[mmm] + nat2 - 1
                molecom.iabond[ia, ib] = iref
                molecom.tbond[ia, ib] = bonnew
            if molecom.itype[ia] == ncha2 and na1 == nat2:
                molecom.nbond[ia] = molecom.nbond[ia] + 1
                ib = molecom.nbond[ia]
                iref = clustercom.naf[mmm] + nat1 - 1
                molecom.iabond[ia, ib] = iref
                molecom.tbond[ia, ib] = bonnew
        units.write(6, fmt(_AB_F41))
        r = units.list_read(5)
        ncha1 = r.int(ncha1)
        ncha2 = r.int(ncha2)
        nat1 = r.int(nat1)
        nat2 = r.int(nat2)
        bonnew = r.real(bonnew)
    return


_RB_F42 = "(1x,'Give 2 atom type numbers and 2 atom numbers',' of bond to delete')"


@routine
def remove_bond():
    """process_symbols.f:remove_bond -- removes the bonds between specified atoms
    (answers read from unit 5 until the first atom type number is not positive)."""
    # NB: the READ targets are uninitialised before the first read; 0 is kept for a null value
    ncha1 = 0
    ncha2 = 0
    nat1 = 0
    nat2 = 0
    units.write(6, fmt(_RB_F42))
    r = units.list_read(5)
    ncha1 = r.int(ncha1)
    ncha2 = r.int(ncha2)
    nat1 = r.int(nat1)
    nat2 = r.int(nat2)
    nbond = molecom.nbond
    iabond = molecom.iabond
    tbond = molecom.tbond
    while ncha1 > 0:
        for ia in range(1, molecom.natom + 1):
            mmm = molecom.mref[ia]
            na1 = ia - clustercom.naf[mmm] + 1
            if molecom.itype[ia] == ncha1 and na1 == nat1 and nbond[ia] > 0:
                found = False
                for ib in range(1, nbond[ia] + 1):
                    ibsav = ib
                    ibref = iabond[ia, ib]
                    na2 = ibref - clustercom.naf[mmm] + 1
                    if molecom.itype[ibref] == ncha2 and na2 == nat2:
                        found = True  # go to 43
                        break
                if found:  # (otherwise go to 44)
                    # 43 continue
                    nbond[ia] = nbond[ia] - 1
                    if nbond[ia] >= ibsav:
                        for ib in range(ibsav, nbond[ia] + 1):
                            ibref = ib + 1
                            iabond[ia, ib] = iabond[ia, ibref]
                            tbond[ia, ib] = tbond[ia, ibref]
                # 44 continue
            elif molecom.itype[ia] == ncha2 and na1 == nat2 and nbond[ia] > 0:
                found = False
                for ib in range(1, nbond[ia] + 1):
                    ibsav = ib
                    ibref = iabond[ia, ib]
                    na2 = ibref - clustercom.naf[mmm] + 1
                    if molecom.itype[ibref] == ncha1 and na2 == nat1:
                        found = True  # go to 45
                        break
                if found:  # (otherwise go to 46)
                    # 45 continue
                    nbond[ia] = nbond[ia] - 1
                    if nbond[ia] >= ibsav:
                        for ib in range(ibsav, nbond[ia] + 1):
                            ibref = ib + 1
                            iabond[ia, ib] = iabond[ia, ibref]
                            tbond[ia, ib] = tbond[ia, ibref]
                # 46 continue
        units.write(6, fmt(_RB_F42))
        r = units.list_read(5)
        ncha1 = r.int(ncha1)
        ncha2 = r.int(ncha2)
        nat1 = r.int(nat1)
        nat2 = r.int(nat2)
    return


_RT_F100 = "(/1x,'molecule type number',2x,'label',2x,'number',2x,'tethered status (T or F)'/)"
_RT_F101 = "(10x,i4,10x,a3,5x,i7,1x,l1)"
_RT_F102 = ("(/'remove_molecule_types> Give type number of molecule(s) to remove (q to quit).'"
            "/' and number to remove'//' If no number is given, all will be removed')")


@routine
def remove_molecule_types():
    """process_symbols.f:remove_molecule_types -- removes a specified type of molecule
    from the ato file."""
    findtype = int_array(epsr_dimension.mcomp)
    moltyperead = ""  # character(len=3)
    ierr = 0

    # Determine the number of molecules of each type

    for it in range(1, clustercom.nmoltype + 1):
        findtype[it] = 0
    for im in range(1, clustercom.nmol + 1):
        imoltype = 0
        ifound = 0
        while imoltype < clustercom.nmoltype and ifound == 0:
            imoltype = imoltype + 1
            if molecom.alabel[clustercom.naf[im]] == clustercom.moltype[imoltype]:
                ifound = imoltype
        if ifound > 0:
            findtype[imoltype] = findtype[imoltype] + 1

    # Make a list of all current molecule types and their tethered status

    # NB: at end of file on unit 5 the Fortran loops forever (readagain stays .true.);
    # reproduced
    readagain = True
    while readagain:
        units.write(6, fmt(_RT_F100))
        for it in range(1, clustercom.nmoltype + 1):
            units.write(6, fmt(_RT_F101, it, Ch(clustercom.moltype[it], 3), findtype[it],
                               bool(clustercom.moltethered[it])))
        units.write(6, fmt(_RT_F102))
        commands.sometext, ierr = _read_line5(commands.sometext, 3020)
        if ierr == 0 and len(commands.sometext) > 0:
            _parse_sometext(5)
            readagain = _word(1) != "q"
            if readagain:
                mtyperead = read_internal(_word(1)).int(0)
                ntoremove = -1
                if commands.nwords > 1:
                    ntoremove, ierr = _read_int(_word(2), ntoremove)
                if ntoremove < 0:
                    ntoremove = int(findtype[mtyperead])
                moltyperead = clustercom.moltype[mtyperead]
                im = 1
                nremove = 0
                while im <= clustercom.nmol:
                    if molecom.alabel[clustercom.naf[im]] == moltyperead and nremove < ntoremove:
                        units.write(6, lw(im, Ch(molecom.alabel[clustercom.naf[im]], 3), Ch(moltyperead, 3)))
                        remove_molecule(im)
                        nremove = nremove + 1
                        findtype[mtyperead] = findtype[mtyperead] - 1
                    else:
                        im = im + 1
        #            else
        #               readagain=.false.
    #
    # now check that all symbol_values have been redefined correctly for any type changes
    # that have been made
    #
    write_symbols_changeato(symbols.nsymboluptodata, symbols.nsymbol)
    return


@routine
def remove_molecule(iremove):
    """process_symbols.f:remove_molecule -- removes the specified molecule."""
    cc = clustercom
    mc = molecom
    mm = molecom_mod

    # Obviously can't remove a non-existent molecule. Also we will always leave at least one molecule
    # in file, otherwise the .ato file becomes unreadable.
    if iremove > 0 and iremove <= cc.nmol and cc.nmol > 1:
        # Get the rotational group number up to this point
        ig = 0
        im = 1
        while im < iremove:
            ig = ig + int(cc.ngrot[im])
            im = im + 1
        mc.natom = mc.natom - int(cc.nato[im])
        ia = int(cc.naf[im]) - 1
        while im < cc.nmol:
            cc.naf[im] = ia + 1
            cc.nato[im] = cc.nato[im + 1]
            cc.xadd[im] = cc.xadd[im + 1]
            cc.yadd[im] = cc.yadd[im + 1]
            cc.zadd[im] = cc.zadd[im + 1]
            cc.phix[im] = cc.phix[im + 1]
            cc.phiy[im] = cc.phiy[im + 1]
            cc.phiz[im] = cc.phiz[im + 1]
            cc.tethered[im] = cc.tethered[im + 1]
            ic = 0
            iaref = int(cc.naf[im + 1]) - 1
            while ic < cc.nato[im]:
                ic = ic + 1
                ia = ia + 1
                iaref = iaref + 1
                mc.alabel[ia] = mc.alabel[iaref]
                mc.xa[ia] = mc.xa[iaref]
                mc.ya[ia] = mc.ya[iaref]
                mc.za[ia] = mc.za[iaref]
                mc.nbond[ia] = mc.nbond[iaref]
                for ib in range(1, mc.nbond[ia] + 1):
                    mc.iabond[ia, ib] = mc.iabond[iaref, ib]
                    mc.tbond[ia, ib] = mc.tbond[iaref, ib]
            cc.nal[im] = ia
            # Now transfer the rotational groups
            igref = ig + int(cc.ngrot[im])
            cc.ngrot[im] = cc.ngrot[im + 1]
            ib = 0
            while ib < cc.ngrot[im]:
                ib = ib + 1
                ig = ig + 1
                igref = igref + 1
                mm.rotgrp[ig] = mm.rotgrp[igref]
                mm.labaxi[ig, 1] = mm.labaxi[igref, 1]
                mm.labaxi[ig, 2] = mm.labaxi[igref, 2]
                mm.nrot[ig] = mm.nrot[igref]
                mm.rotf[ig] = mm.rotf[igref]
                mm.rotl[ig] = mm.rotl[igref]
            im = im + 1
        cc.nmol = cc.nmol - 1
    # Finally clean up the type lists
    cleanup_types()
    return


@routine
def cleanup_types():
    """process_symbols.f:cleanup_types -- renews atom and molecule type lists, removing
    anything that is no longer needed."""
    cc = clustercom
    mc = molecom
    mm = molecom_mod
    findtype = int_array(epsr_dimension.mcomp)
    #
    # Check if any of the previous types no longer exist - if so remove them from the list
    #
    for it in range(1, mc.ntype + 1):
        findtype[it] = 0
    for im in range(1, cc.nmol + 1):
        for ia in range(cc.naf[im], cc.nal[im] + 1):
            ifound = 0
            it = 0
            while it < mc.ntype and ifound == 0:
                it = it + 1
                if mc.alabel[ia] == mc.type[it]:
                    ifound = it
            if ifound > 0:
                findtype[ifound] = findtype[ifound] + 1
    #
    # remove any types that are no longer needed. Start from the end to make sure we catch them all
    #
    it = mc.ntype
    while it > 0:
        if findtype[it] <= 0:
            #
            # remove this type
            #
            jt = it
            while jt < mc.ntype:
                mc.type[jt] = mc.type[jt + 1]
                mc.epscor[jt] = mc.epscor[jt + 1]
                mc.sigcor[jt] = mc.sigcor[jt + 1]
                mc.amassa[jt] = mc.amassa[jt + 1]
                mc.charge[jt] = mc.charge[jt + 1]
                mc.qradius[jt] = mc.qradius[jt + 1]
                mc.atsymbol[jt] = mc.atsymbol[jt + 1]
                mc.iexchange[jt] = mc.iexchange[jt + 1]
                mm.rhosphere[jt] = mm.rhosphere[jt + 1]
                mm.natomsphere[jt] = mm.natomsphere[jt + 1]
                mm.radminsphere[jt] = mm.radminsphere[jt + 1]
                mm.radmaxsphere[jt] = mm.radmaxsphere[jt + 1]
                jt = jt + 1
            mc.ntype = mc.ntype - 1
        it = it - 1
    # Finally reset the type numbers for each atom
    for ia in range(1, mc.natom + 1):
        it = 0
        ifound = 0
        while it < mc.ntype and ifound == 0:
            it = it + 1
            if mc.alabel[ia] == mc.type[it]:
                ifound = it
        mc.itype[ia] = ifound

    # Repeat this for the list of molecule types

    for it in range(1, cc.nmoltype + 1):
        findtype[it] = 0
    for im in range(1, cc.nmol + 1):
        imoltype = 0
        ifound = 0
        while imoltype < cc.nmoltype and ifound == 0:
            imoltype = imoltype + 1
            if mc.alabel[cc.naf[im]] == cc.moltype[imoltype]:
                ifound = imoltype
        if ifound > 0:
            findtype[imoltype] = findtype[imoltype] + 1
    for imoltype in range(1, cc.nmoltype + 1):
        units.write(6, fmt("(i4,1x,a3,1x,i7)", imoltype, Ch(cc.moltype[imoltype], 3), findtype[imoltype]))
    # remove any types no longer needed. Start from end to make sure we catch all of them.
    imoltype = cc.nmoltype
    while imoltype > 0:
        if findtype[imoltype] <= 0:
            jt = imoltype
            while jt < cc.nmoltype:
                cc.moltype[jt] = cc.moltype[jt + 1]
                cc.moltethered[jt] = cc.moltethered[jt + 1]
                cc.moltemplate[jt] = cc.moltemplate[jt + 1]
                cc.wholerot[jt] = cc.wholerot[jt + 1]
                cc.wholetrans[jt] = cc.wholetrans[jt + 1]
                jt = jt + 1
            cc.nmoltype = cc.nmoltype - 1
        imoltype = imoltype - 1
    # Finally reset the type numbers for each molecule
    for im in range(1, cc.nmol + 1):
        imoltype = 0
        ifound = 0
        while imoltype < cc.nmoltype and ifound == 0:
            imoltype = imoltype + 1
            if mc.alabel[cc.naf[im]] == cc.moltype[imoltype]:
                ifound = imoltype
        cc.mtype[im] = ifound

    return


_AW_F101 = ("(/'changeato> Type atom label of first atom in molecule (to identify molecule),'"
            "/'<Enter> to leave unchanged: ')")
_AW_F102 = "(/'changeato> Type two atom numbers (within molecule) to form rotational group: ')"


@routine
def add_whole_rot(tokens, ntokens, ntf, ntl):
    """process_symbols.f:add_whole_rot -- adds a rotational group which will rotate the
    whole molecule about a specified pair of atoms."""
    cc = clustercom
    mc = molecom
    mm = molecom_mod
    ctemp = ""  # character(len=3)
    ierr = 0

    # This routine will add a rotational group which will rotate the whole molecule about a specified pair of atoms
    if ntokens > 1:
        commands.sometext = fstr(substr(tokens, ntf[2], ntl[2]), 3020)
    else:
        units.write(6, fmt(_AW_F101), advance=False)
        commands.sometext = ""
        commands.sometext, ierr = _read_line5(commands.sometext, 3020)
    _parse_sometext(11)
    if commands.nwords > 0:
        ctemp = fstr(read_internal(_word(1)).word(ctemp), 3)
        # Find the molecule type this corresponds to
        imoltype = 0
        ifound = 0
        while imoltype < cc.nmoltype and ifound == 0:
            imoltype = imoltype + 1
            if ctemp == cc.moltype[imoltype]:
                ifound = imoltype
        if ifound > 0:
            rotatom1 = 0
            rotatom2 = 0
            if ntokens > 3:
                ierr = 0
                try:
                    r = read_internal(substr(tokens, ntf[3], ntl[4]))
                    rotatom1 = r.int(rotatom1)
                    rotatom2 = r.int(rotatom2)
                except FortranEOF:
                    ierr = -1
                except FortranIOError as e:
                    ierr = e.iostat
            else:
                units.write(6, fmt(_AW_F102), advance=False)
                ierr = 0
                try:
                    r = units.list_read(5)
                    rotatom1 = r.int(rotatom1)
                    rotatom2 = r.int(rotatom2)
                except FortranEOF:
                    ierr = -1
                except FortranIOError as e:
                    ierr = e.iostat  # noqa: F841 (iostat result unused, as in the Fortran)
            # Only proceed if both atom numbers are .gt. 0 and not the same value
            if rotatom1 > 0 and rotatom2 > 0 and rotatom1 != rotatom2:
                # Step through each molecule and stop if the molecule is the same type as that specified
                ig = 0
                for im in range(1, cc.nmol + 1):
                    ig = ig + int(cc.ngrot[im])
                    if cc.mtype[im] == ifound and rotatom1 <= cc.nato[im] and rotatom2 <= cc.nato[im]:
                        # Set the rotational axis atom to the actual atom number for this molecule
                        rotatom1this = rotatom1 + int(cc.naf[im]) - 1
                        rotatom2this = rotatom2 + int(cc.naf[im]) - 1
                        # Create a new rotational group. First need to shift all the existing groups beyond ig up by 1
                        ig = ig + 1
                        ih = mm.ngrrot
                        ih1 = mm.ngrrot + 1
                        mm.ngrrot = ih1
                        lib_molecom_mod.testsize_rotgrps()
                        while ih >= ig:
                            mm.rotgrp[ih1] = mm.rotgrp[ih]
                            mm.labaxi[ih1, 1] = mm.labaxi[ih, 1]
                            mm.labaxi[ih1, 2] = mm.labaxi[ih, 2]
                            mm.nrot[ih1] = mm.nrot[ih]
                            mm.rotf[ih1] = mm.rotf[ih]
                            mm.rotl[ih1] = mm.rotl[ih]
                            # Do any dihedral groups.
                            if mm.rotgrp[ih] == "DIH":
                                mc.ndihedralang[ih1] = mc.ndihedralang[ih]
                                for k in range(1, mc.ndihedralang[ih1] + 1):
                                    mc.dihedralang[k, ih1] = mc.dihedralang[k, ih]
                                    for j in range(1, 3):
                                        mc.labrotdih[j, k, ih1] = mc.labrotdih[j, k, ih]
                            ih1 = ih
                            ih = ih - 1
                        # The rotational group will be all the atoms in the molecule except rotatom1 and rotatom2
                        cc.ngrot[im] = cc.ngrot[im] + 1
                        #                     call testsize_rotgrps()
                        mm.rotgrp[ig] = "ROT"
                        mm.labaxi[ig, 1] = rotatom1this
                        mm.labaxi[ig, 2] = rotatom2this
                        mm.rotf[ig] = mm.nrottot + 1
                        # Number of atoms in rotational group is equal to number of atoms in molecule except the two forming the axis
                        mm.nrot[ig] = cc.nato[im] - 2
                        mm.nrottot = mm.nrottot + int(mm.nrot[ig])
                        mm.rotl[ig] = mm.nrottot
                        lib_molecom_mod.testsize_labrot()
                        labref = int(mm.rotf[ig]) - 1
                        units.write(6, lw(ig, mm.rotf[ig], mm.rotl[ig], rotatom1this, rotatom2this, cc.naf[im],
                                          cc.nal[im]))
                        for ia in range(cc.naf[im], cc.nal[im] + 1):
                            if ia != rotatom1this and ia != rotatom2this:
                                labref = labref + 1
                                mm.labrot[labref] = ia
                                units.write(6, lw(labref, ia, mm.labrot[labref]))
    return


_AI_F101 = ("(/'changeato> Type atom label of first atom in molecule (to identify molecule),'"
            "/'<Enter> to leave unchanged: ')")
_AI_F102 = ("(/'changeato> Type atom number and image atom number to set (within same molecule)"
            " (<Enter> to quit): ')")
_AI_F103 = "('changeato> Remove connection between atom ',i5,' and its image ',i5,' (y or n)? ')"  # (unused)


@routine
def add_image_atom(ido, tokens, ntokens, ntf, ntl):
    """process_symbols.f:add_image_atom -- adds (ido) or removes (.not.ido) an image
    atom, joining one end of a long molecule to its tail on the opposite side of the
    simulation box."""
    cc = clustercom
    mc = molecom
    ctemp = ""  # character(len=3)
    yorno = ""  # character(len=1)
    ierr = 0  # noqa: F841

    # This routine will add an image atom, joining one end of a long molecule to its tail on the opposite side of the simulation box.
    if ntokens > 1:
        commands.sometext = fstr(substr(tokens, ntf[2], ntl[2]), 3020)
    else:
        units.write(6, fmt(_AI_F101), advance=False)
        commands.sometext = ""
        commands.sometext, ierr = _read_line5(commands.sometext, 3020)
    _parse_sometext(11)
    if commands.nwords > 0:
        ctemp = fstr(read_internal(_word(1)).word(ctemp), 3)
        # Find the molecule type this corresponds to
        imoltype = 0
        ifound = 0
        while imoltype < cc.nmoltype and ifound == 0:
            imoltype = imoltype + 1
            if ctemp == cc.moltype[imoltype]:
                ifound = imoltype
        if ifound > 0 and ido:
            rotatom1s = -1
            rotatom2s = -1
            # NB: when the atom numbers come from tokens (ntokens > 3) the Fortran loop never
            # ends (the same tokens are re-read each pass); reproduced
            while rotatom1s != 0 and rotatom2s != 0:
                if ntokens > 3:
                    commands.sometext = fstr(substr(tokens, ntf[3], ntl[4]), 3020)
                else:
                    units.write(6, fmt(_AI_F102), advance=False)
                    commands.sometext, ierr = _read_line5(commands.sometext, 3020)
                _parse_sometext(11)
                if commands.nwords > 1:
                    r = read_internal(commands.sometext)
                    rotatom1s = r.int(rotatom1s)
                    rotatom2s = r.int(rotatom2s)
                else:
                    rotatom1s = 0
                    rotatom2s = 0
                # Only proceed if both atom numbers are .gt. 0 and not the same value
                if rotatom1s > 0 and rotatom2s > 0 and rotatom1s != rotatom2s:
                    # Step through each molecule and stop if the molecule is the same type as that specified
                    ig = 0  # noqa: F841
                    for im in range(1, cc.nmol + 1):
                        if cc.mtype[im] == ifound and rotatom1s <= cc.nato[im] and rotatom2s <= cc.nato[im]:
                            # Set the rotational axis atom to the actual atom number for this molecule
                            rotatom1 = rotatom1s + int(cc.naf[im]) - 1
                            rotatom2 = rotatom2s + int(cc.naf[im]) - 1
                            units.write(6, lw(im, cc.mtype[im], ifound, rotatom1s, rotatom2s, rotatom1, rotatom2,
                                              cc.naf[im]))
                            mc.imageatom[rotatom1] = rotatom2
                            mc.imageatom[rotatom2] = rotatom1
        elif ifound > 0 and not ido:
            # Make a list of the image atoms for this molecular type
            im = 0
            found = False
            while im < cc.nmol and not found:
                im = im + 1
                found = cc.mtype[im] == ifound
            if found:
                nadd = 1 - int(cc.naf[im])
                for ia in range(cc.naf[im], cc.nal[im] + 1):
                    if mc.imageatom[ia] > 0:
                        iaremove = ia + nadd
                        iaremoveimage = int(mc.imageatom[ia]) + nadd
                        #                     write(6,103,advance='no') iaremove,iaremoveimage
                        #                     read(5,'(a)',iostat=ierr) yorno
                        yorno = "y"
                        if yorno == "y" or yorno == "Y":
                            # Do this for all the other molecules of this type. im should be the first of its type
                            # found in the .ato file
                            im1 = im - 1
                            while im1 < cc.nmol:
                                im1 = im1 + 1
                                # NB: the Fortran tests mtype(im), not mtype(im1); reproduced
                                if cc.mtype[im] == ifound:
                                    nadd1 = int(cc.naf[im1]) - 1
                                    mc.imageatom[iaremove + nadd1] = 0
                                    mc.imageatom[iaremoveimage + nadd1] = 0
    return


_MM_F600 = "(1x,'potential ',a3,1x,4(f10.5,1x),a3)"


@routine
def make_mol_files():
    """process_symbols.f:make_mol_files -- writes a .mol template file for each molecule
    type that does not have one yet."""
    cc = clustercom
    mc = molecom
    int_as_str1 = ""  # character(len=3)
    int_as_str2 = ""
    int_as_str3 = ""
    # Variables for checking
    moltemplatecheck = ""  # character(len=256)
    moltemplateshort = ""  # character(len=256)
    lenmoltemplatecheck = 0
    a = real_array(3)
    a1 = real_array(3)
    a2 = real_array(3)
    b = real_array(3)
    c1perp = real_array(3)
    c2perp = real_array(3)
    csign = real_array(3)
    xa = mc.xa
    ya = mc.ya
    za = mc.za

    pi = PI4ATAN  # 4.0*atan(1.0)
    piconv = f32(180) / pi
    nbondsfound = 0
    ndihedralsfound = 0
    for imoltype in range(1, cc.nmoltype + 1):
        lenmoltemplateshort = auxiliary_routines.indexlastdot(cc.moltemplate[imoltype]) - 1
        if lenmoltemplateshort <= 0:
            lenmoltemplateshort = len(cc.moltemplate[imoltype])
        moltemplateshort = fstr(substr(cc.moltemplate[imoltype], 1, lenmoltemplateshort) + ".mol", 256)
        lenmoltemplateshort = lenmoltemplateshort + 4
        moltemplatecheck, lenmoltemplatecheck = file_and_directory_commands.get_full_filename(
            moltemplateshort, lenmoltemplateshort, moltemplatecheck, lenmoltemplatecheck)
        units.write(6, fmt("(/1x,a)", substr(moltemplatecheck, 1, lenmoltemplatecheck)))
        # Try to open this file. If it exists, we won't do anything.
        try:
            units.open(10, moltemplatecheck, status="old")
            ierr = 0
        except FortranIOError as e:
            ierr = e.iostat
        if ierr != 0:
            # Create a new template file based on .ato filename plus label of first atom in molecule
            cc.moltemplate[imoltype] = fstr(
                substr(symbols.fname_inp, 1, auxiliary_routines.indexlastdot(symbols.fname_inp) - 1)
                + substr(cc.moltype[imoltype], 1, len(cc.moltype[imoltype])), 256)
            lenmoltemplateshort = auxiliary_routines.indexlastdot(symbols.fname_inp) - 1 + len(cc.moltype[imoltype])
            units.write(6, fmt("(/1x,a,a)", "make_mol_files> molecule template name: ",
                               substr(cc.moltemplate[imoltype], 1, lenmoltemplateshort)))
            moltemplatecheck = fstr(
                substr(symbols.fnameinp, 1, auxiliary_routines.indexlastdot(symbols.fnameinp) - 1)
                + substr(cc.moltype[imoltype], 1, len(cc.moltype[imoltype])) + ".mol", 256)
            lenmoltemplatecheck = (auxiliary_routines.indexlastdot(symbols.fnameinp) - 1
                                   + len(cc.moltype[imoltype]) + 4)
            units.write(6, fmt('(/1x,"make_mol_files> creating .mol file: ",a)',
                               substr(moltemplatecheck, 1, lenmoltemplatecheck)))
            im = 0
            found = False
            while im < cc.nmol and not found:
                im = im + 1
                found = mc.alabel[cc.naf[im]] == cc.moltype[imoltype]
            if found:  # Use this molecule to form the template
                # The next bit is a replica of ghemical_routines_busch.f (readmole) with appropriate changes related to atom labelling.
                try:
                    units.open(10, moltemplatecheck, status="unknown")
                    ierr = 0
                except FortranIOError as e:
                    ierr = e.iostat
                # Signal that this is a .gmol file and that the atom numbers listed are the EPSR atom numbers
                if ierr == 0:
                    units.write(10, lw(".gmol ", 0))
                    iasub = 1 - int(cc.naf[im])
                    # Write atoms and their bonding arrangements.
                    for ia1 in range(cc.naf[im], cc.nal[im] + 1):
                        ia = ia1 + iasub
                        int_as_str1 = fstr(fmt("(i3.3)", ia), 3)
                        units.write(10, lw(" atom ", ia, " " + pad(int_as_str1, 3) + " ", xa[ia1], ya[ia1], za[ia1],
                                           mc.nbond[ia1],
                                           *[int(mc.iabond[ia1, j]) + iasub for j in range(1, mc.nbond[ia1] + 1)]))
                    # Write the bond list.
                    for ia1 in range(cc.naf[im], cc.nal[im] + 1):
                        int_as_str1 = fstr(fmt("(i3.3)", ia1 + iasub), 3)
                        ib1 = 0
                        while ib1 < mc.nbond[ia1]:
                            ib1 = ib1 + 1
                            ia2 = int(mc.iabond[ia1, ib1])
                            if ia2 > ia1:
                                nbondsfound = nbondsfound + 1
                                int_as_str2 = fstr(fmt("(i3.3)", ia2 + iasub), 3)
                                units.write(10, fmt('(" bond ",a3,1x,a3,1x,f10.5)', Ch(int_as_str1, 3),
                                                    Ch(int_as_str2, 3), mc.tbond[ia1, ib1]))
                    # Write the bond angle list
                    for ia1 in range(cc.naf[im], cc.nal[im] + 1):
                        int_as_str1 = fstr(fmt("(i3.3)", ia1 + iasub), 3)
                        for ib1 in range(1, mc.nbond[ia1] + 1):
                            ia2 = int(mc.iabond[ia1, ib1])
                            int_as_str2 = fstr(fmt("(i3.3)", ia2 + iasub), 3)
                            a[1] = xa[ia2] - xa[ia1]
                            a[2] = ya[ia2] - ya[ia1]
                            a[3] = za[ia2] - za[ia1]
                            ia2ref = ia2 + iasub  # noqa: F841
                            for ib2 in range(1, mc.nbond[ia1] + 1):
                                ia3 = int(mc.iabond[ia1, ib2])
                                if ia3 > ia2:
                                    int_as_str3 = fstr(fmt("(i3.3)", ia3 + iasub), 3)
                                    b[1] = xa[ia3] - xa[ia1]
                                    b[2] = ya[ia3] - ya[ia1]
                                    b[3] = za[ia3] - za[ia1]
                                    angle = piconv * acosf(auxiliary_routines.scapro(a, b)
                                                           / sqrtf(auxiliary_routines.scapro(a, a)
                                                                   * auxiliary_routines.scapro(b, b)))
                                    units.write(10, fmt("(a,a,1x,f10.5)", " angle ",
                                                        pad(int_as_str2, 3) + " " + pad(int_as_str1, 3) + " "
                                                        + pad(int_as_str3, 3), angle))
                    # Set up all possible dihedral angles
                    for ia1 in range(cc.naf[im], cc.nal[im] + 1):
                        ib1 = 0
                        while ib1 < mc.nbond[ia1]:
                            ib1 = ib1 + 1
                            ia2 = int(mc.iabond[ia1, ib1])
                            # Define 2 axis atoms
                            if ia2 > ia1:
                                b[1] = xa[ia2] - xa[ia1]
                                b[2] = ya[ia2] - ya[ia1]
                                b[3] = za[ia2] - za[ia1]
                                # Get neighbours of ia1
                                ib2 = 0
                                while ib2 < mc.nbond[ia1]:
                                    ib2 = ib2 + 1
                                    ia1ref = int(mc.iabond[ia1, ib2])
                                    if ia1ref != ia2:
                                        a1[1] = xa[ia1ref] - xa[ia1]
                                        a1[2] = ya[ia1ref] - ya[ia1]
                                        a1[3] = za[ia1ref] - za[ia1]
                                        # Calculate the vector perpendicular to the plane formed by these two vectors
                                        auxiliary_routines.vecpro(a1, b, c1perp)
                                        # Normalise this vector
                                        amp = sqrtf(auxiliary_routines.scapro(c1perp, c1perp))
                                        for i in range(1, 4):
                                            c1perp[i] = c1perp[i] / amp
                                        # Now search the bonds of ia2 for the 4th atom
                                        ib3 = 0
                                        while ib3 < mc.nbond[ia2]:
                                            ib3 = ib3 + 1
                                            ia2ref = int(mc.iabond[ia2, ib3])
                                            if ia2ref != ia1 and ia2ref != ia1ref:
                                                a2[1] = xa[ia2ref] - xa[ia2]
                                                a2[2] = ya[ia2ref] - ya[ia2]
                                                a2[3] = za[ia2ref] - za[ia2]
                                                # Calculate the vector perpendicular to the plane formed by these two vectors
                                                auxiliary_routines.vecpro(a2, b, c2perp)
                                                # Normalise this vector
                                                amp = sqrtf(auxiliary_routines.scapro(c2perp, c2perp))
                                                for i in range(1, 4):
                                                    c2perp[i] = c2perp[i] / amp
                                                # Cosine of the angle between them
                                                cosphi = auxiliary_routines.scapro(c1perp, c2perp)
                                                if cosphi > 1.0:
                                                    cosphi = f32(1.0)
                                                if cosphi < -1.0:
                                                    cosphi = f32(-1.0)
                                                # Sine of the angle between them - see torangles for this
                                                auxiliary_routines.vecpro(c1perp, c2perp, csign)
                                                sinphi = (auxiliary_routines.scapro(csign, b)
                                                          / sqrtf(auxiliary_routines.scapro(b, b)))
                                                angle = piconv * atan2f(sinphi, cosphi)
                                                units.write(10, lw("dihedral ", ia1ref, ia1, ia2, ia2ref, angle))
                                                ndihedralsfound = ndihedralsfound + 1
                    # Write the changelabel commands
                    for ia1 in range(cc.naf[im], cc.nal[im] + 1):
                        int_as_str1 = fstr(fmt("(i3.3)", ia1 + iasub), 3)
                        units.write(10, lw("changelabel " + pad(int_as_str1, 3) + " ", Ch(mc.alabel[ia1], 3)))
                    # Write potential parameters
                    for it in range(1, mc.ntype + 1):
                        # Check that this type exists in the molecule
                        found = False
                        ia = int(cc.naf[im]) - 1
                        while ia < cc.nal[im] and not found:
                            ia = ia + 1
                            found = mc.alabel[ia] == mc.type[it]
                        if found:
                            units.write(10, fmt(_MM_F600, Ch(mc.type[it], 3), mc.epscor[it], mc.sigcor[it],
                                                mc.amassa[it], mc.charge[it], Ch(mc.atsymbol[it], 3)))
                    units.write(10, lw("temperature ", 300))
                    units.write(10, lw("vibtemp ", f32(65.0)))
                    units.write(10, lw("angtemp ", f32(1.0)))
                    units.write(10, lw("dihtemp ", f32(10.0)))
                    units.write(10, lw("ecoredcore ", f32(0.0), f32(0.0)))
                    units.write(10, lw("density ", epsr16.rho))
                    # Done.
                    units.close(10)
                    units.write(6, fmt("(/1x,a)", "make_mol_files> *** Remember to run fmole on the new .ato file. ***"))
    if nbondsfound > 0:
        mc.stepmi = r4(0.1)
    else:
        mc.stepmi = f32(0.0)
    if ndihedralsfound > 0:
        mc.stepri = r4(0.1)
    else:
        mc.stepri = f32(0.0)
    # Ensure the bond angle and dihedral weights have been set
    mc.tempangreq = f32(1.0)
    mc.tempdihreq = f32(10.0)
    mc.tempangwt = f32(1.0)
    mc.tempdihwt = f32(10.0)
    return


_SF_F100 = "(/1x,'molecule type number',2x,'label',2x,'number',2x,'relative width'/)"
_SF_F101 = "(10x,i4,10x,a3,5x,i0,1x,f10.5)"
_SF_F102 = ("(/'set_molecule_scale_factors> Give type number of molecule(s) to set (q to quit).'"
            "/' and relative width of Poisson function (0 - 1)')")


@routine
def set_molecule_scale_factors(tokens, ntokens, ntf, ntl):
    """process_symbols.f:set_molecule_scale_factors -- sets the molecular scale factor
    for a particular molecule type."""
    cc = clustercom
    mm = molecom_mod
    findtype = int_array(epsr_dimension.mcomp)
    moltyperead = ""  # character(len=3)
    # NB: rfract is uninitialised in the Fortran when its read fails; 0.0 is used here
    rfract = f32(0.0)  # Relative width of Poisson function. 0.0 means all molecules have the same size, scale factor = 1.0
    ierr = 0

    thrd = r4(1.0) / r4(3.0)
    # Determine the number of molecules of each type
    for it in range(1, cc.nmoltype + 1):
        findtype[it] = 0
    for im in range(1, cc.nmol + 1):
        imoltype = 0
        ifound = 0
        while imoltype < cc.nmoltype and ifound == 0:
            imoltype = imoltype + 1
            if molecom.alabel[cc.naf[im]] == cc.moltype[imoltype]:
                ifound = imoltype
        if ifound > 0:
            findtype[imoltype] = findtype[imoltype] + 1

    # Make a list of all current molecule types and their current scale factor
    units.write(6, fmt(_SF_F100))
    for it in range(1, cc.nmoltype + 1):
        units.write(6, fmt(_SF_F101, it, Ch(cc.moltype[it], 3), findtype[it], mm.rfractmol[it]))
    if ntokens > 2:
        commands.sometext = fstr(substr(tokens, ntf[2], ntl[3]), 3020)
        ierr = 0
    else:
        units.write(6, fmt(_SF_F102))
        commands.sometext, ierr = _read_line5(commands.sometext, 3020)
    if ierr == 0 and len(commands.sometext) > 0:
        _parse_sometext(5)
        # NB: mtyperead is uninitialised if this read fails; it is then not used
        mtyperead, ierr = _read_int(_word(1), 0)
        if commands.nwords > 1 and ierr == 0:
            rfract, ierr = _read_real(_word(2), rfract)
            if rfract > 1.0:
                rfract = f32(1.0)
            if rfract > 0.0:
                pfroutines_mod.getpfprobability(rfract)
            mm.rfractmol[mtyperead] = rfract
            ntoset = int(findtype[mtyperead])
            moltyperead = cc.moltype[mtyperead]
            im = 0
            nset = 0
            while im < cc.nmol:
                im = im + 1
                if molecom.alabel[cc.naf[im]] == moltyperead and nset < ntoset:
                    if rfract > 0.0:
                        mm.scalefacmol[im] = powf(pfroutines_mod.getrratio(), thrd)  # Assume the MASS of the particle is distributed normally
                    else:
                        mm.scalefacmol[im] = 1.0
                    units.write(6, lw(im, Ch(molecom.alabel[cc.naf[im]], 3), " ", Ch(moltyperead, 3),
                                      mm.scalefacmol[im]))
                    nset = nset + 1
    return
