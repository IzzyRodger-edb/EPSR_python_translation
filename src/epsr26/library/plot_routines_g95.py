"""Translation of ``subroutine read_plot3djmol`` from fortran_src/plot_routines_g95.f.

plot3djmol.f (``run_plot3djmol``) and setup_generic.f (``setup``) call
``read_plot3djmol``, but no active source file defines it: the active
plot_routines.f has no such routine, and plot_routines_g95.f (an older variant of
plot_routines.f, otherwise not used) is the only file in fortran_src/ that does.  Only
this routine is translated from it; everything it calls is taken from the active files
(``plot_routines.find_plot_par``, ``parse.concatenate``, ...), which is how the
comparison binaries in work/N2/fbuild were linked (see docs/translation_notes/N2.md).

It reads a plot3djmol input file in the old format (values one per line, in symbol
order, no keywords).
"""

from ..common import commands, symbols
from ..fortran import FortranEOF, FortranIOError, FortranStop, f32, fstr, iostat_of, lw, read_internal, substr, units
from ..trace import routine
from . import parse


def _const_nlast(nlast, nwords):
    """concatenate assigns its nlast dummy (``if(nlast.gt.nwords) nlast=nwords``).  When the
    actual argument is a literal constant, gfortran keeps it in read-only memory and the
    assignment kills the program with SIGSEGV (seen with the reference build: an old-format
    file whose chsymbol/rthetaphi line has fewer than 4 words).  Reproduced as an error
    stop with the status of a SIGSEGV death (139).  (The 1,1 / 3,3 / 4,4 calls cannot
    reach it: they are guarded by nwords.)"""
    if nlast > nwords:
        raise FortranStop(139, "Program received signal SIGSEGV: Segmentation fault - invalid memory reference.")


@routine
def read_plot3djmol(nin, nout, ierrout):
    """plot_routines_g95.f:read_plot3djmol -- sets up the plot3d input files.

    Returns (nout, ierrout).
    """
    # use molecule_coordinates_mod   (none of its variables is used here)
    from . import file_and_directory_commands, plot_routines

    # integer nin,nout,ierr,ierrout
    # real realtest
    cm = commands
    sy = symbols
    # NB: ilen/ilen1/ierr2/nvaryt are set before use (uninitialised in the Fortran)
    nvaryt = 0
    #
    # return if an error occurred on trying to open the input file
    #
    if ierrout != 0:
        units.close(nin)
        return nout, ierrout
    #
    #  If O.K. read in the data from the specified input file
    #
    sy.symbol_value[1] = fstr(substr(sy.fname_run, 1, sy.ilenfname_run), 256)
    #
    # commands to read the rest of the input file
    #
    icount = 0
    ifilecount = 1
    ierr = 0
    while icount < sy.nsymbol and ierr == 0:
        icount = icount + 1
        iref = file_and_directory_commands.ireff(icount, ifilecount)
        sy.symbol_save = sy.symbol_value[iref]
        if (sy.symbol[icount] != 'mirrorsym' and sy.symbol[icount] != 'tlucency'
                and sy.symbol[icount] != 'scalefact' and sy.symbol[icount] != 'rotatemol'):
            try:
                cm.sometext = fstr(units.read_line(nin), cm.nstsize)
                ierr = 0
            except (FortranEOF, FortranIOError) as e:
                ierr = iostat_of(e)
            units.write(6, lw(icount, iref, ifilecount, substr(cm.sometext, 1, len(cm.sometext))))
            if ierr == 0:
                cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
                if cm.nwords > 0:
                    if icount == sy.nsymboluptodata:
                        _, cm.nwords, sy.symbol_value[iref], ilen = parse.concatenate(
                            1, cm.nwords, sy.symbol_value[iref], 0, cm.nwords, cm.ncf, cm.ncl, cm.sometext)
                        try:
                            sy.nfiles = read_internal(sy.symbol_value[iref]).int(sy.nfiles)
                            ierr = 0
                        except (FortranEOF, FortranIOError) as e:
                            ierr = iostat_of(e)
                    elif icount == sy.nsymboluptodata + 1:
                        _, _, sy.symbol_value[iref], ilen = parse.concatenate(
                            1, 1, sy.symbol_value[iref], 0, cm.nwords, cm.ncf, cm.ncl, cm.sometext)
                        icount = icount + 1
                        iref = file_and_directory_commands.ireff(icount, ifilecount)
                        _const_nlast(4, cm.nwords)
                        _, _, sy.symbol_value[iref], ilen1 = parse.concatenate(
                            2, 4, sy.symbol_value[iref], 0, cm.nwords, cm.ncf, cm.ncl, cm.sometext)
                    elif sy.symbol[icount] == 'rotelev':
                        _const_nlast(2, cm.nwords)
                        _, _, sy.symbol_value[iref], ilen = parse.concatenate(
                            1, 2, sy.symbol_value[iref], 0, cm.nwords, cm.ncf, cm.ncl, cm.sometext)
                        icount = icount + 1
                        iref = file_and_directory_commands.ireff(icount, ifilecount)
                        if cm.nwords > 2:
                            _, _, sy.symbol_value[iref], ilen = parse.concatenate(
                                3, 3, sy.symbol_value[iref], 0, cm.nwords, cm.ncf, cm.ncl, cm.sometext)
                        # Convert this value (scale factor) to a real number to decide whether
                        # to rotate the molecule if appropriate
                        realtest = f32(0.0)
                        try:
                            realtest = read_internal(sy.symbol_value[iref]).real(realtest)
                            ierr = 0
                        except (FortranEOF, FortranIOError) as e:
                            ierr = iostat_of(e)
                        if ierr != 0:
                            realtest = f32(1.0)
                        # Make sure a scale factor .gt. 0.0 will be used
                        if realtest <= 0.0:
                            sy.symbol_value[iref] = '1.0'
                        icount = icount + 1
                        iref = file_and_directory_commands.ireff(icount, ifilecount)
                        if cm.nwords > 3:
                            _, _, sy.symbol_value[iref], ilen = parse.concatenate(
                                4, 4, sy.symbol_value[iref], 0, cm.nwords, cm.ncf, cm.ncl, cm.sometext)
                        icount = icount + 1
                        iref = file_and_directory_commands.ireff(icount, ifilecount)
                        # Symbol: rotatemol
                        if realtest > 0.0:
                            sy.symbol_value[iref] = 'y'
                        else:
                            sy.symbol_value[iref] = 'no'
                    else:
                        _, cm.nwords, sy.symbol_value[iref], ilen = parse.concatenate(
                            1, cm.nwords, sy.symbol_value[iref], 0, cm.nwords, cm.ncf, cm.ncl, cm.sometext)
                    # if symbol_value .eq. nvary has changed need to change symbol_text for the
                    # next line.
                    if sy.symbol[icount] == 'nvary':
                        try:
                            nvaryt = read_internal(sy.symbol_value[iref]).int(nvaryt)
                        except (FortranEOF, FortranIOError):
                            pass  # ierr2
                        if nvaryt == 1:
                            sy.symbol_text[iref + 1] = 'phim, thetam, chim (degrees - separated by spaces)'
                        elif nvaryt == 2:
                            sy.symbol_text[iref + 1] = 'thetal, phil, chim (degrees - separated by spaces)'
                        elif nvaryt == 3:
                            sy.symbol_text[iref + 1] = 'thetal, phil, phim (degrees - separated by spaces)'
            if icount >= sy.nsymboldata and ifilecount < sy.nfiles:
                ifilecount = ifilecount + 1
                icount = sy.nsymboluptodata
    units.close(nin)
    #
    # determine the number of coefficients in the coefficients file and whether the
    # molecules are identical or distinct
    #
    nout, ierrout = plot_routines.find_plot_par(nin, nout, ierrout)
    return nout, ierrout
