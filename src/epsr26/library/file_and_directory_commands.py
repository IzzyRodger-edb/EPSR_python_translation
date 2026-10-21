"""Translation of fortran_src/file_and_directory_commands.f.

Command-line arguments, the current/home folders, directory listing, file search and
the input-file and symbol-table helpers used by every EPSR program.

FORMATs ending in ``$`` (no advance) are written without the ``$`` and with
``advance=False`` (fformat does not support ``$``).

Character comparisons that involve a substring (blank-padded) use ``_feq``/``_fne``,
which ignore trailing blanks as Fortran's ``.eq.``/``.ne.`` do.
"""

import os

from ..common import commands, directory_name, symbols, system_commands, version
from ..fortran import (
    Ch,
    FortranEOF,
    FortranIOError,
    findex,
    findex_padded,
    fmt,
    fstr,
    idiv,
    lw,
    runtime,
    substr,
    units,
)
from . import parse, stringutility
from . import system_commands as system_commands_routines
from ..trace import routine


def _feq(a, b):
    """Fortran ``a .eq. b`` for character values (the shorter is blank-padded)."""
    return a.rstrip(" ") == b.rstrip(" ")


def _fne(a, b):
    """Fortran ``a .ne. b`` for character values."""
    return a.rstrip(" ") != b.rstrip(" ")


def _getcwd():
    """``ierr = getcwd(name)`` (gfortran GETCWD function): (status, name)."""
    try:
        return 0, os.getcwd()
    except OSError as e:
        return e.errno or 1, ""


def _read_a(unit, old, n):
    """``read(unit,'(a)',iostat=ierr) var`` for a ``character*n`` var: (var, ierr).

    On end of file or error the variable keeps its old value.
    """
    try:
        return fstr(units.read_line(unit), n), 0
    except FortranEOF:
        return old, -1
    except FortranIOError as e:
        return old, e.iostat


_F600 = "(/'Binaries folder: ',a)"
_F100_get_command_from_args = "(/'Current folder for input files is: ',a)"
_F102_get_command_from_args = "(a)"


@routine
def get_command_from_args():
    """file_and_directory_commands.f:get_command_from_args.

    Get the command into sometext from the command line arguments
    (argument 1: the current folder; argument 2: symbol_setup; the rest: the command).
    """
    # use String_Utility
    # character*256 text
    # integer*4 ilensometext

    # Get the current version

    version.version = fstr(version.date, 20)
    version.ilenversion = version.ilendate
    #
    # get the system commands for the particular operating system in use
    #
    system_commands_routines.get_system_commands()

    # Add on the path separator to the system_binaries.

    sc = system_commands
    sc.system_binaries = fstr(substr(sc.system_binaries, 1, sc.ilensystem_binaries)
                              + substr(sc.system_termination, 1, sc.ilensystem_termination), 256)
    sc.ilensystem_binaries = sc.ilensystem_binaries + sc.ilensystem_termination

    units.write(6, fmt(_F600, substr(sc.system_binaries, 1, sc.ilensystem_binaries)))
    #
    # set initial directory name - this is called the home directory
    #
    get_directory_name()
    directory_name.ncharhomedir = directory_name.nchardir
    directory_name.home_directory = directory_name.directory_name

    # Get input from any supplied arguments

    narg = runtime.iargc()
    iarg = 0
    commands.sometext = ""
    ilensometext = 0
    while iarg < narg:

        iarg = iarg + 1
        text = fstr(runtime.getarg(iarg), 256)

        # Parse the text
        # NB: text is character*256 but parse is told nstsize; parse stops after 5 blanks,
        # so only an argument of more than 250 characters would read past text.

        commands.nwords = parse.parse(text, commands.nwords, commands.ncf, commands.ncl, 5,
                                      commands.nstsize, commands.nworddim)
        nwords = commands.nwords
        ncf = commands.ncf
        ncl = commands.ncl

        # Set the current directory - should be the first argument

        if iarg == 1:

            directory_name.directory_name = fstr(substr(text, ncf[1], ncl[nwords]), 256)
            directory_name.nchardir = int(ncl[nwords]) - int(ncf[1]) + 1
            units.write(6, fmt(_F100_get_command_from_args,
                               substr(directory_name.directory_name, 1, directory_name.nchardir)))

        elif iarg == 2:

            # The second argument will be normally be symbol_setup: these are always lower case.

            commands.sometext = fstr(substr(text, ncf[1], ncl[nwords]), commands.nstsize)
            ilensometext = int(ncl[nwords]) - int(ncf[1]) + 1
            symbols.symbol_setup = fstr(stringutility.strlowcase(substr(commands.sometext, 1, ilensometext)), 20)
            symbols.ilensetup = ilensometext

        else:

            # Otherwise accumulate the arguments into a single text line separated by spaces

            commands.sometext = fstr(substr(commands.sometext, 1, ilensometext) + " "
                                     + substr(text, ncf[1], ncl[nwords]), commands.nstsize)
            ilensometext = ilensometext + int(ncl[nwords]) - int(ncf[1]) + 2

    units.write(6, fmt(_F102_get_command_from_args, substr(commands.sometext, 1, ilensometext)))

    # Read the command

    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    units.write(6, lw(commands.nwords))

    return


_F100_get_directory_name = "(/'Home folder is: ',a)"


@routine
def get_directory_name():
    """file_and_directory_commands.f:get_directory_name.

    sets the directory where input will be read from.
    """
    ierr, cwd = _getcwd()  # noqa: F841 (ierr is not used in the Fortran either)
    commands.sometext = fstr(cwd, commands.nstsize)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)

    # concatenate the words into the directory name

    nlast = commands.nwords
    _, nlast, directory_name.directory_name, directory_name.nchardir = parse.concatenate(
        1, nlast, directory_name.directory_name, directory_name.nchardir,
        commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)

    # To cope with spaces in folder names we need to insert the directory name between quotes.
    # These cannot include the first (if present) and last system termination character.

    # Check whether the first character is a system termination character.

    #      if(directory_name(1:ilensystem_termination).eq.system_termination(1:ilensystem_termination)) then
    #         directory_name=system_termination(1:ilensystem_termination)//'"'
    #     *//directory_name(ilensystem_termination+1:nchardir)
    #      else
    #         directory_name='"'//directory_name(ilensystem_termination+1:nchardir)
    #      endif
    #      nchardir=nchardir+1

    # check that the last characters are the correct termination sequence. If not add it.

    dn = directory_name
    sc = system_commands
    ntest = dn.nchardir - sc.ilensystem_termination + 1
    if _fne(substr(dn.directory_name, ntest, dn.nchardir),
            substr(sc.system_termination, 1, sc.ilensystem_termination)):
        dn.directory_name = fstr(substr(dn.directory_name, 1, dn.nchardir)
                                 #    *   //'"'
                                 + substr(sc.system_termination, 1, sc.ilensystem_termination), 256)
        dn.nchardir = dn.nchardir + sc.ilensystem_termination
    else:
        dn.directory_name = fstr(substr(dn.directory_name, 1, ntest - 1)
                                 #     *   //'"'
                                 + substr(sc.system_termination, 1, sc.ilensystem_termination), 256)
    #      nchardir=nchardir+1

    units.write(6, fmt(_F100_get_directory_name, substr(dn.directory_name, 1, dn.nchardir)))

    return


_F101_set_directory_name = "(/'No such folder! : ',a)"
_F100_set_directory_name = "(/'Current folder for input files: ',a)"


def _back_one_directory(dn, sc):
    """The '..' block of set_directory_name (it occurs twice, identically):

    Working from the end of the directory name search for the first occurrence of
    a system termination character (but ignoring the very last one). If none are found,
    do nothing.
    """
    ntest = dn.nchardir - sc.ilensystem_termination
    ipos = ntest + sc.ilensystem_termination - 1
    while (_fne(substr(dn.directory_name, ntest, ipos), substr(sc.system_termination, 1, sc.ilensystem_termination))
           and ntest > 0):
        ntest = ntest - 1
        ipos = ipos - 1
    if ntest > 0:
        dn.nchardir = ipos
        dn.directory_name = fstr(substr(dn.directory_name, 1, dn.nchardir), 256)


@routine
def set_directory_name():
    """file_and_directory_commands.f:set_directory_name.

    sets the current directory

    if no directory is given the current directory is shown
    """
    # Temporary character store
    # character*256 temp_name,stringcheck,dir_temp,dir_test
    # character*20 testseq
    # integer*4 nchartemp,ilentestseq,ipos,ilenstringcheck,nchardirtemp,nchardirtest
    # logical fileexists
    dn = directory_name
    sc = system_commands
    sometext = commands.sometext
    ncf = commands.ncf
    ncl = commands.ncl
    nwords = commands.nwords
    temp_name = ""
    nchartemp = 0
    #
    # Store the current directory in a temporary variable

    dir_temp = fstr(substr(dn.directory_name, 1, dn.nchardir), 256)
    nchardirtemp = dn.nchardir
    #
    # check that a name was specified - if not do nothing
    #
    if nwords > 1:
        if (substr(sometext, ncf[2], ncl[2]) == "home"
                or substr(sometext, ncf[2], ncl[2]) == "."):

            # If the second word is 'home'or '.'
            # the directory name is reset to the home directory
            # (which is the directory where EPSRshell was first started)

            dn.directory_name = dn.home_directory
            dn.nchardir = dn.ncharhomedir

        elif substr(sometext, ncf[2], ncl[2]) == "..":

            # If the second word is '..' the directory name goes one level back up the tree

            # Working from the end of the directory name search for the first occurrence of
            # a system termination character (but ignoring the very last one). If none are found,
            # do nothing.

            _back_one_directory(dn, sc)

        else:

            nlast = nwords

            # Place the remaining words into a store

            _, nlast, temp_name, nchartemp = parse.concatenate(
                2, nlast, temp_name, nchartemp, nwords, ncf, ncl, sometext, 256)

            # Search for system_termination characters

            icount = 0
            testseq = fstr(substr(sc.system_termination, 1, sc.ilensystem_termination), 20)
            ilentestseq = sc.ilensystem_termination
            ifirst = 1
            iadd = findex(substr(temp_name, ifirst, nchartemp), substr(testseq, 1, ilentestseq))

            while iadd > 0:

                ilast = ifirst + iadd - 2

                # Check what is the string before this. If ilast.lt.ifirst there is nothing before this
                # so start a new directory name.

                # Otherwise check whether the characters before the system_termination are '.' or '..'

                if ilast >= ifirst:

                    stringcheck = fstr(substr(temp_name, ifirst, ilast), 256)
                    ilenstringcheck = ilast - ifirst + 1

                    if _feq(substr(temp_name, ifirst, ilast), "."):

                        # Single '.' means backup to home directory

                        dn.directory_name = dn.home_directory
                        dn.nchardir = dn.ncharhomedir

                    elif _feq(substr(temp_name, ifirst, ilast), ".."):

                        # Double '..' means backup one directory

                        # Working from the end of the directory name search for the first occurrence of
                        # a system termination character (but ignoring the very last one). If none are found,
                        # do nothing.

                        _back_one_directory(dn, sc)

                    else:

                        # Otherwise add the specified string to the existing directory name

                        dn.directory_name = fstr(substr(dn.directory_name, 1, dn.nchardir)
                                                 + substr(stringcheck, 1, ilenstringcheck)
                                                 + substr(sc.system_termination, 1, sc.ilensystem_termination), 256)
                        dn.nchardir = dn.nchardir + ilenstringcheck + sc.ilensystem_termination

                elif ifirst == 1:

                    # In the event that we encounter a system_termination string at the beginnning,
                    # we start a new directory name

                    dn.directory_name = fstr(substr(sc.system_termination, 1, sc.ilensystem_termination), 256)
                    dn.nchardir = sc.ilensystem_termination

                else:

                    # Do nothing
                    pass

                icount = icount + 1
                ifirst = ilast + sc.ilensystem_termination + 1
                #     write(6,*) icount,ifirst,ilast

                if ifirst <= nchartemp:

                    iadd = findex(substr(temp_name, ifirst, nchartemp), substr(testseq, 1, ilentestseq))

                    # If there is still text in the word, we need one more pass to process it

                    if iadd == 0:
                        iadd = nchartemp - ifirst + 2

                else:

                    iadd = 0

            # In case there is still text in the word, we add it to the current directory name

            if ifirst <= nchartemp:

                dn.directory_name = fstr(substr(dn.directory_name, 1, dn.nchardir)
                                         + substr(temp_name, ifirst, nchartemp)
                                         + substr(sc.system_termination, 1, sc.ilensystem_termination), 256)
                dn.nchardir = dn.nchardir + nchartemp - ifirst + 1 + sc.ilensystem_termination

        # check that the last characters are the correct termination sequence. If not add it.

        ntest = dn.nchardir - sc.ilensystem_termination + 1
        if _fne(substr(dn.directory_name, ntest, dn.nchardir),
                substr(sc.system_termination, 1, sc.ilensystem_termination)):
            dn.directory_name = fstr(substr(dn.directory_name, 1, dn.nchardir)
                                     + substr(sc.system_termination, 1, sc.ilensystem_termination), 256)
            dn.nchardir = dn.nchardir + sc.ilensystem_termination

    # Now check that this folder exists

    dir_test = fstr(substr(dn.directory_name, 1, dn.nchardir) + ".", 256)
    nchardirtest = dn.nchardir + 1
    # Check whether this file exists
    fileexists = units.inquire_exists(substr(dir_test, 1, nchardirtest))
    #      write(6,*) dir_test(1:nchardirtest),' ',fileexists
    if not fileexists:
        units.write(6, fmt(_F101_set_directory_name, substr(dn.directory_name, 1, dn.nchardir)))
        dn.directory_name = fstr(substr(dir_temp, 1, nchardirtemp), 256)
        dn.nchardir = nchardirtemp

    units.write(6, fmt(_F100_set_directory_name, substr(dn.directory_name, 1, dn.nchardir)))

    return


@routine
def list_directory():
    """file_and_directory_commands.f:list_directory.

    shows the files in the current input directory
    """
    # character*256 command,text
    # integer*4 nchartext
    dn = directory_name
    sc = system_commands
    text = ""
    nchartext = 0
    if commands.nwords > 1:
        #
        # concatenate the remaining words into the filename list
        #
        nlast = commands.nwords
        _, nlast, text, nchartext = parse.concatenate(
            2, nlast, text, nchartext, commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
        command = fstr(substr(sc.system_ls, 1, sc.ilensystem_ls)
                       + " " + substr(sc.system_quotes, 1, sc.ilensystem_quotes)
                       + substr(dn.directory_name, 1, dn.nchardir) + substr(text, 1, nchartext)
                       + substr(sc.system_quotes, 1, sc.ilensystem_quotes), 256)

    else:

        command = fstr(substr(sc.system_ls, 1, sc.ilensystem_ls)
                       + " " + substr(sc.system_quotes, 1, sc.ilensystem_quotes)
                       + substr(dn.directory_name, 1, dn.nchardir)
                       + substr(sc.system_quotes, 1, sc.ilensystem_quotes), 256)

    ierr = runtime.system(command)  # noqa: F841
    return


# 100   format(/'make_directory> Type the directory name to create: ',$)
_F100_make_directory = "(/'make_directory> Type the directory name to create: ')"


@routine
def make_directory():
    """file_and_directory_commands.f:make_directory.

    shows the files in the current input directory (sic: creates a folder)
    """
    # character*256 command,text,fnamefull
    # integer*4 nchartext,ilenfnamefull
    sc = system_commands
    text = ""
    nchartext = 0
    fnamefull = ""
    ilenfnamefull = 0
    if commands.nwords > 1:
        #
        # concatenate the remaining words into the filename list
        #
        nlast = commands.nwords
        _, nlast, text, nchartext = parse.concatenate(
            2, nlast, text, nchartext, commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
    else:
        units.write(6, fmt(_F100_make_directory), advance=False)
        commands.sometext = fstr(units.read_line(5), commands.nstsize)
        commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                      commands.nstsize, commands.nworddim)
        _, commands.nwords, text, nchartext = parse.concatenate(
            1, commands.nwords, text, nchartext, commands.nwords, commands.ncf, commands.ncl,
            commands.sometext, 256)

    fnamefull, ilenfnamefull = get_full_filename(text, nchartext, fnamefull, ilenfnamefull)
    command = fstr(substr(sc.system_md, 1, sc.ilensystem_md)
                   + " "
                   + substr(sc.system_quotes, 1, sc.ilensystem_quotes)
                   + substr(fnamefull, 1, ilenfnamefull)
                   + substr(sc.system_quotes, 1, sc.ilensystem_quotes), 256)

    ierr = runtime.system(command)  # noqa: F841

    return


_F100_edit_file = "(/'edit_file> ',a)"


@routine
def edit_file():
    """file_and_directory_commands.f:edit_file.

    shows the files in the current input directory (sic: opens a file in the editor)
    """
    # character*256 command,text,fnamefull
    # integer*4 nchartext,ilenfnamefull
    sc = system_commands
    text = ""
    nchartext = 0
    fnamefull = ""
    ilenfnamefull = 0
    if commands.nwords > 1:
        #
        # concatenate the remaining words into the filename list
        #
        nlast = commands.nwords
        _, nlast, text, nchartext = parse.concatenate(
            2, nlast, text, nchartext, commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
        fnamefull, ilenfnamefull = get_full_filename(text, nchartext, fnamefull, ilenfnamefull)
    else:
        fnamefull = ""
        ilenfnamefull = 0

    units.write(6, fmt(_F100_edit_file, substr(fnamefull, 1, ilenfnamefull)))

    if ilenfnamefull > 0:
        command = fstr(substr(sc.system_edit, 1, sc.ilensystem_edit)
                       + " "
                       + substr(sc.system_quotes, 1, sc.ilensystem_quotes)
                       + substr(fnamefull, 1, ilenfnamefull)
                       + substr(sc.system_quotes, 1, sc.ilensystem_quotes), 256)

        ilencommand = sc.ilensystem_edit + ilenfnamefull + 1 + 2 * sc.ilensystem_quotes

    else:

        command = fstr(substr(sc.system_edit, 1, sc.ilensystem_edit), 256)
        ilencommand = sc.ilensystem_edit

    units.write(6, fmt(_F100_edit_file, substr(command, 1, ilencommand)))

    ierr = runtime.system(command)  # noqa: F841

    return


_F100_delete_file = "(/'delete_file> No files specified to delete. ')"


@routine
def delete_file():
    """file_and_directory_commands.f:delete_file.

    shows the files in the current input directory (sic: deletes a file)
    """
    # character*256 command,text,fnamefull
    # integer*4 nchartext,ilenfnamefull
    sc = system_commands
    text = ""
    nchartext = 0
    fnamefull = ""
    ilenfnamefull = 0
    if commands.nwords > 1:
        #
        # concatenate the remaining words into the filename list
        #
        nlast = commands.nwords
        _, nlast, text, nchartext = parse.concatenate(
            2, nlast, text, nchartext, commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
        fnamefull, ilenfnamefull = get_full_filename(text, nchartext, fnamefull, ilenfnamefull)
        command = fstr(substr(sc.system_del, 1, sc.ilensystem_del)
                       + " "
                       + substr(sc.system_quotes, 1, sc.ilensystem_quotes)
                       + substr(fnamefull, 1, ilenfnamefull)
                       + substr(sc.system_quotes, 1, sc.ilensystem_quotes), 256)

        ierr = runtime.system(command)  # noqa: F841

    else:

        units.write(6, fmt(_F100_delete_file))

    return


@routine
def get_full_filename(fnameshort, ilenshort, fnamelong, ilenlong):
    """file_and_directory_commands.f:get_full_filename.

    checks to see whether a directory path is present in fnameshort, and if not
    adds it. Output in fnamelong.
    Returns (fnamelong, ilenlong).
    """
    # character*256 fnameshort,fnamelong
    # integer*4 ilenshort,ilenlong,ilentermination
    sc = system_commands
    dn = directory_name
    fnameshort = fstr(fnameshort, 256)
    #
    # search for directory termination characters in fnameshort
    #
    ilentermination = findex_padded(fnameshort, substr(sc.system_termination, 1, sc.ilensystem_termination), 256)

    if dn.nchardir > 0 and ilentermination == 0:
        fnamelong = fstr(substr(dn.directory_name, 1, dn.nchardir) + substr(fnameshort, 1, ilenshort), 256)
        ilenlong = dn.nchardir + ilenshort
    else:
        fnamelong = fstr(substr(fnameshort, 1, ilenshort), 256)
        ilenlong = ilenshort

    return fnamelong, ilenlong


@routine
def get_nterms(nin, fname, extin, nterms, ierr):
    """file_and_directory_commands.f:get_nterms.

    returns the number of columns in a file. Lines containing a # are ignored.
    It is assumed the maximum number of columns in a file is 201 (EPSR or Gudrun
    files). If this number is reached the nextfile in the sequence is opened if
    possible
    Returns (fname, nterms, ierr).
    """
    from . import auxiliary_routines  # imported here: auxiliary_routines.f belongs to another package

    # character*256 fname
    # character*2 ext
    # character*10 extin
    # integer*4 ilenextin,nterms
    fname = fstr(fname, 256)
    extin = fstr(extin, 10)

    ilenextin = findex_padded(extin, " ", 10) - 1
    units.write(6, lw(Ch(fname, 256), " ", substr(extin, 1, ilenextin)))
    ntest = 100
    nterms = 0
    ifilecount = 0
    ierr = 0
    while ntest >= 100 and ierr == 0:
        ifilecount = ifilecount + 1
        #
        # generate the version number
        #
        ext = fstr(fmt("(i2.2)", ifilecount), 2)
        #
        # strip any existing extension off the input filename
        #
        ilenfname = auxiliary_routines.indexlastdot(fname) - 1
        if ilenfname <= 0:
            ilenfname = findex_padded(fname, " ", 256) - 1
        #
        # add the expected extension back on
        #
        fname = fstr(substr(fname, 1, ilenfname) + substr(extin, 1, ilenextin), 256)
        ilenfname = ilenfname + ilenextin
        #
        # strip the last two characters off the filename *assume it is of the form
        # xxxxx.y01, etc.) and create the sequence of versions in turn
        #
        ilenfname = ilenfname - 2
        fname = fstr(substr(fname, 1, ilenfname) + substr(ext, 1, 2), 256)
        #
        units.write(6, lw(substr(fname, 1, ilenfname + 2)))
        #
        # try to open this file
        #
        try:
            units.open(nin, fname, status="old")
            ierr = 0
        except FortranIOError as e:
            ierr = e.iostat
        if ierr == 0:
            #
            # read lines of text until one occurs without a # in it
            #
            nhash = 1
            ierr = 0
            while nhash != 0 and ierr == 0:
                commands.sometext, ierr = _read_a(nin, commands.sometext, commands.nstsize)
                nhash = findex(commands.sometext, "#")
            #
            # find out how many words this line contains
            #
            if nhash == 0:
                commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                              commands.nstsize, commands.nworddim)
                #
                # number of terms = (ncolumns-1)/2
                #
                ntest = idiv(commands.nwords - 1, 2)
                nterms = nterms + ntest
            units.close(nin)
    return fname, nterms, ierr


_F101_setup_input_file = "(/'setup_input_file> File class: \"',a,'\"; file extension: \"',a,'\"')"
_F100_setup_input_file = "(/'setup input file> No input file specified - routine will exit'/)"
_F105_setup_input_file = "(/'setup_input_file> Full filename = ',a)"
_F102_setup_input_file = "(/'setup_input_file> Reading input file: \"',a,'\"')"
_F104_setup_input_file = "(/'setup_input_file> Run name in input file is different from ','filename specified')"
_F103_setup_input_file = "(/'setup_input_file> Problems with specified input file: ',a,//'   - will use default values')"


@routine
def setup_input_file(nin, nout, ierrout):
    """file_and_directory_commands.f:setup_input_file.

    sets up the EPSR input files: find the input file name in sometext (after the
    symbol_setup word) or ask for it, open it on unit nin and read the run name.
    Returns (nout, ierrout).
    """
    from . import auxiliary_routines  # imported here: auxiliary_routines.f belongs to another package

    # input and output stream number
    # integer*4 nin,nout
    # integer*4 ierrout,ilenext,ilenuptodot
    # character*1 iaccept
    # character*20 ext
    # character*256 text
    sy = symbols
    iaccept = ""
    #
    # generate the extension from the last dot onwards
    #
    text = fstr(sy.extinp, 256)
    sy.ilenextinpdot = auxiliary_routines.indexlastdot(text)
    sy.extdot = fstr(substr(sy.extinp, sy.ilenextinpdot, sy.ilenextinp), 20)
    sy.ilenextdot = sy.ilenextinp - sy.ilenextinpdot + 1
    #     write(6,*) extdot(1:ilenextdot),ilenextinp,ilenextinpdot,ilenextdot
    #
    # write out the value of the input class and symbol_setup to ensure we are
    # indeed accessing the class of file requested
    #
    units.write(nout, fmt(_F101_setup_input_file, substr(sy.symbol_setup, 1, sy.ilensetup),
                          substr(sy.extinp, 1, sy.ilenextinp)))
    #
    # parse the supplied text
    #
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    #
    # generate the name of the input file
    #
    iw = 1
    while iw <= commands.nwords and _fne(substr(commands.sometext, commands.ncf[iw], commands.ncl[iw]),
                                         substr(sy.symbol_setup, 1, sy.ilensetup)):
        iw = iw + 1
    if iw >= commands.nwords:
        #
        # no filename found so prompt for it from the user
        #
        ext = fstr(sy.extinp, 20)
        ilenext = sy.ilenextinp
        #     write(6,*) ext(1:ilenext)
        commands.nwords = 0
        commands.sometext = ""
        while commands.nwords <= 0:
            iaccept = search_for_files(ext, ilenext, iaccept, len_iaccept=1)
            commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                          commands.nstsize, commands.nworddim)
            if commands.nwords <= 0:
                units.write(nout, fmt(_F100_setup_input_file))
                ierrout = 1
                return nout, ierrout
        iw = 0
    iw = iw + 1
    iw, commands.nwords, sy.fname_inp, sy.ilenfname_inp = parse.concatenate(
        iw, commands.nwords, sy.fname_inp, sy.ilenfname_inp,
        commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
    #
    # if this filename has a .ext we do nothing and will attempt to read it as is,
    # warts and all!
    #
    text = fstr(sy.fname_inp, 256)
    ilenuptodot = findex(sy.fname_inp, ".") - 1
    #     write(6,*) fname_inp(1:ilenfname_inp),ilenfname_inp,ilenuptodot
    if ilenuptodot < 0:
        #
        # .ext not supplied - therefore we need to add on the correct one
        #
        ilenuptodot = sy.ilenfname_inp
        sy.fname_inp = fstr(substr(sy.fname_inp, 1, ilenuptodot) + substr(sy.extinp, 1, sy.ilenextinp), 256)
        sy.ilenfname_inp = ilenuptodot + sy.ilenextinp
    #
    # create the full name of the input file by adding on the directory name if not present
    #
    sy.fnameinp, sy.ilenfnameinp = get_full_filename(sy.fname_inp, sy.ilenfname_inp, sy.fnameinp, sy.ilenfnameinp)
    units.write(6, fmt(_F105_setup_input_file, substr(sy.fnameinp, 1, sy.ilenfnameinp)))
    #
    # open this file if possible
    #
    try:
        units.open(nin, sy.fnameinp, status="old")
        ierrout = 0
    except FortranIOError as e:
        ierrout = e.iostat
    if ierrout == 0:
        #
        #  If O.K. read in the data from the specified input file
        #
        units.write(nout, fmt(_F102_setup_input_file, substr(sy.fname_inp, 1, sy.ilenfname_inp)))
        #
        # find the number of characters in this filename up to the last dot
        #
        text = fstr(sy.fname_inp, 256)
        ilenuptodot = auxiliary_routines.indexlastdot(text) - 1
        #
        # read the name of the run - this should correspond to fname_inp but without the
        # last .extension. It will be used to construct the names of input
        # and output files (except .ato and .pcof files which may be named independently
        # from the name of the run)
        #
        try:
            r = units.list_read(nin)
            sy.fname_run = fstr(r.word(sy.fname_run), 256)
            ierrout = 0
        except FortranEOF:
            ierrout = -1
        except FortranIOError as e:
            ierrout = e.iostat
        sy.ilenfname_run = len(sy.fname_run)
        #
        # we don't do the following two checks for plot input files
        #
        #
        # remove the .extension if supplied
        #
        if substr(sy.symbol_setup, 1, 4) != "plot":
            #
            # check that this name corresponds to filename (without the .extension) specified.
            # If not signal the error but carry on, using the correct name when writing the file
            #
            if _fne(substr(sy.fname_run, 1, sy.ilenfname_run), substr(sy.fname_inp, 1, ilenuptodot)):
                units.write(nout, fmt(_F104_setup_input_file))
                sy.fname_run = fstr(substr(sy.fname_inp, 1, ilenuptodot), 256)
                sy.ilenfname_run = ilenuptodot
        #
        # generate the long form of fname_run
        #
        sy.fnamerun, sy.ilenfnamerun = get_full_filename(sy.fname_run, sy.ilenfname_run,
                                                         sy.fnamerun, sy.ilenfnamerun)
        #
        # signal that if an error occurred the attempt to open the input file
        # was unsuccessful, but we can still read commands from the terminal
        #
        if ierrout != 0:
            ierrout = -1
        return nout, ierrout
    else:
        sy.fname_run = fstr(substr(sy.fname_inp, 1, ilenuptodot), 256)
        sy.ilenfname_run = ilenuptodot
        units.write(nout, fmt(_F103_setup_input_file, substr(sy.fname_inp, 1, sy.ilenfname_inp)))
        #
        # generate the long form of fname_run
        #
        sy.fnamerun, sy.ilenfnamerun = get_full_filename(sy.fname_run, sy.ilenfname_run,
                                                         sy.fnamerun, sy.ilenfnamerun)
        #
        # signify that the input file could not be opened or read, but a run
        # name has been successfully defined with the default values
        #
        ierrout = -1
        return nout, ierrout


# 123   format(/'Filename: ',a,'? (Type y to accept, u to go back, e to exit) ',$)
_F123_search_for_files = "(/'Filename: ',a,'? (Type y to accept, u to go back, e to exit) ')"
_F125_search_for_files = "(/1x,'No files of extension \"',a,'\" found in directory \"',a,'\"'/)"
# 124   format(/1x,'No files selected...'//1x,'Type the required filename (without extension): ',$)
_F124_search_for_files = "(/1x,'No files selected...'//1x,'Type the required filename (without extension): ')"


@routine
def search_for_files(ext, ilenext, iaccept, len_iaccept=1):
    """file_and_directory_commands.f:search_for_files.

    list any files with the required extension in the current directory and let the
    user pick one (into sometext).  iaccept is 'y' if a file was selected or typed.
    ``len_iaccept`` is the declared length of the caller's iaccept (``character(len=*)``).
    Returns iaccept.
    """
    # character*256 fname
    # character(len=*) ext
    # character(len=*) iaccept
    # integer*4 ilenext,ilenstar,ilenextin
    # integer*4 ifirstchar,idirpath
    sc = system_commands
    dn = directory_name
    nstsize = commands.nstsize
    fname = ""
    ilenfname = 0
    #
    # first set up the command to do a directory listing into a file
    #
    commands.sometext = ""
    #      write(6,*) sometext(1:len_trim(sometext))
    commands.sometext = fstr(substr(sc.system_ls, 1, sc.ilensystem_ls) + " "
                             + substr(sc.system_quotes, 1, sc.ilensystem_quotes)
                             + substr(dn.directory_name, 1, dn.nchardir) + "*" + substr(ext, 1, ilenext)
                             + substr(sc.system_quotes, 1, sc.ilensystem_quotes)
                             + " > dirlist.txt", nstsize)
    #      write(6,*) sometext(1:len_trim(sometext))
    #
    # get the directory listing into a file, then asks the user to select a file
    # if found then iaccept is set to 'y'. If none found or none specified then iaccept is set to 'n'
    #
    ierr = runtime.system(commands.sometext)
    #
    # open this file and list it line by line
    #
    if ierr == 0:
        try:
            units.open(10, "dirlist.txt", status="old")
            ierr = 0
        except FortranIOError as e:
            ierr = e.iostat
    #
    # strip out the * in ext if present
    #
    ilenstar = findex(ext, "*")
    if ilenstar > 0:
        ilenextin = ilenstar - 1  # noqa: F841
    else:
        ilenextin = ilenext  # noqa: F841
    iaccept = ""
    ifound = 0
    #
    # number of lines read
    #
    nline = 0
    while ierr == 0 and iaccept != "y" and iaccept != "e":
        commands.sometext = ""
        commands.nwords = 0
        commands.sometext, ierr = _read_a(10, commands.sometext, nstsize)
        nline = nline + 1
        commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                      nstsize, commands.nworddim)
        if commands.nwords > 0:
            ifound = ifound + 1
            _, commands.nwords, fname, ilenfname = parse.concatenate(
                1, commands.nwords, fname, ilenfname, commands.nwords, commands.ncf, commands.ncl,
                commands.sometext, 256)
            #           write(6,*) ifound,nwords

            # Strip out everything up to the last directory path character sequence

            idirpath = 1
            ifirstchar = 1
            while idirpath > 0 and ifirstchar <= ilenfname:
                idirpath = findex(substr(fname, ifirstchar, ilenfname),
                                  substr(sc.system_termination, 1, sc.ilensystem_termination))
                if idirpath > 0:
                    ifirstchar = ifirstchar + idirpath + sc.ilensystem_termination - 1
            fname = fstr(substr(fname, ifirstchar, ilenfname), 256)
            ilenfname = ilenfname - ifirstchar + 1
            units.write(6, fmt(_F123_search_for_files, substr(fname, 1, ilenfname)), advance=False)
            #
            # save the filename in case it is accepted
            #
            commands.sometext = fstr(substr(fname, 1, ilenfname), nstsize)
            iaccept = fstr(units.read_line(5), len_iaccept)
            #      write(6,*) 'search_for_files> ',iaccept(1:len_trim(iaccept))
            if iaccept == "u":
                #     write(6,*) nline
                #
                # request to backspace through file
                #
                units.backspace(10)
                units.backspace(10)
                nline = nline - 2
                #
                # If nline is .lt. 0 it means we have attempted to skip past the beginning of file.
                # Therefore read to end of file, counting the number of lines and reposition before last
                # line
                #
                if nline < 0:
                    #
                    # ensure file is positioned at beginning
                    #
                    units.rewind(10)
                    nline = 0
                    ierr = 0
                    while ierr == 0:
                        commands.sometext = ""
                        commands.nwords = 0
                        commands.sometext, ierr = _read_a(10, commands.sometext, nstsize)
                        nline = nline + 1
                        commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf,
                                                      commands.ncl, 5, nstsize, commands.nworddim)
                    #     write(6,*) ierr,nline,nwords,sometext(ncf(1):ncl(1))
                    units.backspace(10)
                    nline = nline - 1
                    #
                    # If the line that caused ierr to be non-zero contains no words, we need to backspace
                    # again.
                    #
                    if commands.nwords <= 0:
                        units.backspace(10)
                        nline = nline - 1
                    ierr = 0
        if ifound > 0 and ierr != 0:
            units.rewind(10)
            ifound = 0
            ierr = 0
            nline = 0
    units.close(10)
    if ifound == 0:
        units.write(6, fmt(_F125_search_for_files, substr(ext, 1, ilenext), substr(dn.directory_name, 1, dn.nchardir)))
        iaccept = fstr("n", len_iaccept)
    if iaccept != "y":
        units.write(6, fmt(_F124_search_for_files), advance=False)
        commands.sometext = ""
        commands.sometext, ierr = _read_a(5, commands.sometext, nstsize)
        #
        # strip out any leading spaces
        #
        commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                      nstsize, commands.nworddim)
        if commands.nwords > 0:
            iaccept = fstr("y", len_iaccept)
            _, commands.nwords, fname, ilenfname = parse.concatenate(
                1, commands.nwords, fname, ilenfname, commands.nwords, commands.ncf, commands.ncl,
                commands.sometext, 256)
            commands.sometext = fstr(substr(fname, 1, ilenfname), nstsize)
        else:
            iaccept = fstr("n", len_iaccept)
    return iaccept


@routine
def get_a_line(nin, nspcword, nchartext, nfirst, nlast, ifilecount):
    """file_and_directory_commands.f:get_a_line.

    gets a line of text from the input file and parses it into words
    ncf and ncl contain the first and last character numbers of each word
    nchartext is the maximum number of characters in a sometext,
    nworddim is dimension of the
    word array, nspcword is the number of spaces between words needed to decide that the
    end of text has been reached.

    The words are read into the symbols if they exist, otherwise the defaults are adopted
    Returns (nfirst, nlast).
    """
    # integer*4 nspcword,nchartext
    # integer*4 nfirst,nlast,ifilecount
    #
    # make sure allowed limits are not exceeded
    #
    if nfirst < 1:
        nfirst = 1
    if nlast > symbols.nsymbol:
        nlast = symbols.nsymbol
    #
    # read a line of text
    #
    commands.sometext, ierr = _read_a(nin, commands.sometext, commands.nstsize)
    #
    # parse it into words
    #
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, nspcword,
                                  nchartext, commands.nworddim)
    i = 0
    icount = nfirst - 1
    while icount < nlast:
        i = i + 1
        icount = icount + 1
        #
        # set up the reference value for this symbol
        #
        iref = ireff(icount, ifilecount)
        if i <= commands.nwords:
            symbols.symbol_value[iref] = fstr(substr(commands.sometext, commands.ncf[i], commands.ncl[i]), 256)
        else:
            symbols.symbol_value[iref] = symbols.symbol_default[icount]
    return nfirst, nlast


@routine
def ireff(icount, ifilecount):
    """file_and_directory_commands.f:ireff.

    calculates the reference integer for the symbol table
    Returns ireff.
    """
    #
    # set up the reference value for this symbol
    #
    if icount > symbols.nsymboldata:
        ireff_ = icount + (symbols.nfiles - 1) * symbols.nsymbolplus
        #         ifilecount=nfiles
    elif icount > symbols.nsymboluptodata:
        ireff_ = icount + (ifilecount - 1) * symbols.nsymbolplus
    else:
        ireff_ = icount
        #         ifilecount=0
    return ireff_


# 101   format(/'search_for_symbol> searching for symbol ',a)
# 100   format(/'search_for_symbol> text ',a,' not found in symbol table')
# 102   format(/'search_for_symbol> found symbol ',a,' number ',i3)


@routine
def search_for_symbol(text, ifound):
    """file_and_directory_commands.f:search_for_symbol.

    returns the symbol number corresponding to a specified text string
    Returns ifound.
    """
    # integer*4 icount,ifound
    # character*20 text
    text = fstr(text, 20)

    #
    # first check that the current value of ifound does not already correspond
    # to the specified symbol
    #
    #     write(6,101) text

    if ifound > 0 and ifound <= symbols.nsymbol:
        if symbols.symbol[ifound] == text:
            return ifound
    ifound = 0
    icount = 0
    while icount < symbols.nsymbol and ifound == 0:
        icount = icount + 1
        if symbols.symbol[icount] == text:
            ifound = icount
    if ifound == 0:
        #        write(6,100) text
        pass
    else:
        #        write(6,102) symbol(ifound),ifound
        pass

    return ifound
