"""Translation of fortran_src/alignmole.f: ``program standalone``, ``subroutine
alignmole`` and ``function displacement``.

Aligns every molecule of a .ato file along the x, y or z axis by Monte Carlo rotations
(``rotatm``) that minimise the mean square distance of its atoms (all, or a list of atom
numbers) from the axis, and rewrites the .ato file.  See docs/translation_notes/N2.md.
"""

from ..common import clustercom, commands, molecom, symbols
from ..fortran import (Ch, FortranEOF, FortranIOError, f32, findex, fmt, fstr, idiv, int_array, lw,
                       read_internal, real_array, substr, units)
from ..library import ato_read_and_save, file_and_directory_commands, molfit, parse
from ..trace import routine

_F499 = "(/'alignmole> No .ato file specified! - routine will exit'/)"
_F100 = "(/'alignmole> Type axis to align to (x, y, or z): ')"
_F101 = "(/'alignmole> Give number of iterations: ')"
_F104 = "(/'alignmole> Give list of atom numbers to define molecular axis' /,' (Type \"all\" to use all atoms in each molecule) : ')"
_F599 = "(/'Problems with .ato file ',a/'Routine will exit now'/)"
_F103 = "(/'alignmole> Molecule: ',i4,' Initial energy: ',e12.5)"
_F102 = "(/'alignmole> Iteration ',i7,' Energy: ',e12.5,' Reject rate: ',f10.5)"


def main():
    """alignmole.f:standalone -- 2009/3/30 Attempt to make this a standalone program."""
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    alignmole()


@routine
def alignmole():
    """alignmole.f:alignmole -- Align a specified molecule along a particular axis."""
    # character(len=256)                  :: readfile
    # character(len=4)                    :: ext
    # character(len=1)                    :: iaccept,alignaxis
    # real                                :: displacement
    # logical                             :: useall
    # integer                             :: niter,iter,naxisatom,j
    # integer, dimension(:), allocatable  :: axisatom   !atom numbers used to define molecular axis
    # real, dimension(:), allocatable     :: xr,yr,zr
    cm = commands
    cc = clustercom
    axisatom = None
    xr = yr = zr = None
    iaccept = ''
    niter = 0
    #
    # generate the name of the .ato file
    #
    ext = '.ato'
    ilenext = 4
    if cm.nwords < 2:
        cm.nwords = 0
        while cm.nwords <= 0:
            iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)
            cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
            if cm.nwords <= 0:
                units.write(6, fmt(_F499))
                return
        readfile = fstr(substr(cm.sometext, cm.ncf[1], cm.ncl[1]), 256)
        ilenreadfile = int(cm.ncl[1]) - int(cm.ncf[1]) + 1
    else:
        readfile = fstr(substr(cm.sometext, cm.ncf[2], cm.ncl[2]), 256)
        ilenreadfile = int(cm.ncl[2]) - int(cm.ncf[2]) + 1
    if cm.nwords > 2:
        alignaxis = fstr(substr(cm.sometext, cm.ncf[3], cm.ncf[3]), 1)
    else:
        units.write(6, fmt(_F100), advance=False)
        alignaxis = fstr(units.read_line(5), 1)
    if cm.nwords > 3:
        niter = read_internal(substr(cm.sometext, cm.ncf[4], cm.ncl[4])).int(niter)
    else:
        units.write(6, fmt(_F101), advance=False)
        niter = units.list_read(5).int(niter)

    # Get the atom numbers to be used to align the molecules. These numbers have to be common
    # to all molecules in the box. Hence best to use this routine for boxes containing only
    # one type of molecule.

    naxisatom = 0
    nf = 5
    # See if the word ' all ' appears on its own in sometext - if so we will use all the atoms
    # in each molecule to align it
    useall = True
    if cm.nwords >= nf:
        useall = findex(substr(cm.sometext, cm.ncf[nf], cm.ncl[nf]), 'all') > 0
    if axisatom is not None:
        axisatom = None
    while axisatom is None and not useall:
        naxisatom = cm.nwords - nf + 1
        if naxisatom >= 2:  # Must have at least 2 axis atoms
            axisatom = int_array(naxisatom)
            ierr = 0
            try:
                r = read_internal(substr(cm.sometext, cm.ncf[nf], cm.ncl[cm.nwords]))
                for j in range(1, naxisatom + 1):
                    axisatom[j] = r.int(axisatom[j])
            except FortranEOF:
                ierr = -1
            except FortranIOError as e:
                ierr = e.iostat
            if ierr != 0:  # If an error occurred on read, then will have input values again
                axisatom = None
        if axisatom is None:
            units.write(6, fmt(_F104), advance=False)
            cm.sometext = fstr(units.read_line(5), cm.nstsize)
            useall = findex(cm.sometext, 'all') > 0
            cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 10, cm.nstsize, cm.nworddim)
            nf = 1
    #
    # strip out any extension that was supplied - extension must be .ato, otherwise
    # nothing is done.
    #
    ilendot = findex(readfile, '.') - 1
    if ilendot > 0:
        readfile = fstr(substr(readfile, 1, ilendot), 256)
        ilenreadfile = ilendot
    #
    # create the full name of the input file by adding on the directory name if not present
    #
    symbols.fnameato, symbols.ilenfnameato = file_and_directory_commands.get_full_filename(
        readfile, ilenreadfile, symbols.fnameato, symbols.ilenfnameato)
    symbols.fnameato = fstr(substr(symbols.fnameato, 1, symbols.ilenfnameato) + '.ato', 256)
    symbols.ilenfnameato = symbols.ilenfnameato + 4
    #
    # initialise number of atom and molecule types to zero
    #
    molecom.ntype = 0
    cc.nmoltype = 0
    #
    # read .ato file
    #
    molecom.sizefactor = f32(1.0)
    # rmaxcell, cellmax, upol, rho, amubar, amubar2 (real) and ndist (integer) are locals
    rmaxcell, cellmax, upol, rho, amubar, amubar2, ndist, ierr = ato_read_and_save.minitc(
        symbols.fnameato, f32(0.0), f32(0.0), f32(0.0), f32(0.0), f32(0.0), f32(0.0), 0, 0)
    if rho == 0.0 or ierr != 0:
        units.write(6, fmt(_F599, Ch(symbols.fnameato, 256)))
        return
    units.write(6, lw(cc.dmax, molecom.dmaxmol))

    # Now check that the specified axis atoms actually exist for all the molecules.
    # If not we use all atoms for all molecules

    if axisatom is not None:
        im = 0
        while im < cc.nmol and ierr == 0:
            im = im + 1
            ia = 0
            while ia < naxisatom and ierr == 0:
                ia = ia + 1
                if axisatom[ia] < 1 or axisatom[ia] > cc.nato[im]:
                    ierr = 1
        if ierr != 0:
            axisatom = None
            naxisatom = 0
            useall = True
    if naxisatom > 0:
        units.write(6, lw(naxisatom, *[axisatom[j] for j in range(1, naxisatom + 1)]))

    # Save the rotational step size

    steprs = molecom.stepr
    if molecom.stepr <= 0.0:
        molecom.stepr = f32(0.1)

    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za
    for im in range(1, cc.nmol + 1):
        # NB: xr, yr, zr are allocated but never assigned in the Fortran, so a rejected
        # move "restores" uninitialised heap memory; zero-filled arrays are used here
        # (as in the Rust port, alignmole.rs).  The test cases (ex47, ex48) do not
        # depend on this.
        if xr is None:
            xr = real_array(int(cc.nato[im]))
            yr = real_array(int(cc.nato[im]))
            zr = real_array(int(cc.nato[im]))
        elif xr.shape[0] - 1 < cc.nato[im]:
            xr = yr = zr = None
            xr = real_array(int(cc.nato[im]))
            yr = real_array(int(cc.nato[im]))
            zr = real_array(int(cc.nato[im]))
        ntrials = 0
        nreject = 0
        naccum = idiv(niter, 10)
        iter_ = 0

        # Calculate the initial displacement from the specified axis

        energold = displacement(im, alignaxis, useall, axisatom, naxisatom)
        units.write(6, fmt(_F103, im, energold))
        while iter_ < niter:
            iter_ = iter_ + 1
            ntrials = ntrials + 1
            # Random rotation about any axis
            molfit.rotatm(im, molecom.stepr)
            energnew = displacement(im, alignaxis, useall, axisatom, naxisatom)
            if energnew > energold:
                # Reject the move
                nreject = nreject + 1
                # Restore the coordinates
                nat = int(cc.naf[im]) - 1
                for ia in range(1, int(cc.nato[im]) + 1):
                    nat = nat + 1
                    xa[nat] = xr[ia]
                    ya[nat] = yr[ia]
                    za[nat] = zr[ia]
            else:
                # If accepted, save the energy
                energold = energnew
            # Write out a message and reappraise step size when appropriate
            # NB: niter < 10 gives naccum = 0, an integer division by zero (SIGFPE) in the
            # Fortran; idiv raises ZeroDivisionError here.
            if naccum * idiv(iter_, naccum) == iter_:
                rat = f32(nreject) / f32(ntrials)
                units.write(6, fmt(_F102, iter_, energold, rat))
                rat = rat / f32(0.75)
                if rat > f32(1.1):
                    rat = f32(1.1)
                if rat < f32(0.9):
                    rat = f32(0.9)
                molecom.stepr = molecom.stepr / rat
    #
    # set the rotational step sizes back to what they were on input
    #
    molecom.stepr = steprs
    #
    # write out the new molecule coordinates
    #
    ato_read_and_save.moloutc(symbols.fnameato)
    return


@routine
def displacement(im, alignaxis, useall, axisatom, naxisatom):
    """alignmole.f:displacement -- mean square distance of the atoms of molecule im (all,
    or the naxisatom atoms listed in axisatom) from the axis alignaxis."""
    # real                    :: displacement,sum,diff
    # character(len=1)        :: alignaxis
    # integer                 :: im,ia,ir,naxisatom
    # integer, dimension(*)   :: axisatom
    # logical                 :: useall
    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za
    naf = clustercom.naf
    nal = clustercom.nal
    sum_ = f32(0.0)
    displacement_ = f32(0.0)
    #      write(6,*) useall,naxisatom,im,nato(im)
    if useall or naxisatom <= 0:
        if alignaxis == 'x':
            for ia in range(int(naf[im]), int(nal[im]) + 1):
                diff = ya[ia]
                sum_ = sum_ + diff * diff
                diff = za[ia]
                sum_ = sum_ + diff * diff
        elif alignaxis == 'y':
            for ia in range(int(naf[im]), int(nal[im]) + 1):
                diff = xa[ia]
                sum_ = sum_ + diff * diff
                diff = za[ia]
                sum_ = sum_ + diff * diff
        else:
            for ia in range(int(naf[im]), int(nal[im]) + 1):
                diff = xa[ia]
                sum_ = sum_ + diff * diff
                diff = ya[ia]
                sum_ = sum_ + diff * diff
        displacement_ = sum_ / f32(clustercom.nato[im])
    else:
        if alignaxis == 'x':
            for ir in range(1, naxisatom + 1):
                ia = int(axisatom[ir]) + int(naf[im]) - 1
                diff = ya[ia]
                sum_ = sum_ + diff * diff
                diff = za[ia]
                sum_ = sum_ + diff * diff
        elif alignaxis == 'y':
            for ir in range(1, naxisatom + 1):
                ia = int(axisatom[ir]) + int(naf[im]) - 1
                diff = xa[ia]
                sum_ = sum_ + diff * diff
                diff = za[ia]
                sum_ = sum_ + diff * diff
        else:
            for ir in range(1, naxisatom + 1):
                ia = int(axisatom[ir]) + int(naf[im]) - 1
                diff = xa[ia]
                sum_ = sum_ + diff * diff
                diff = ya[ia]
                sum_ = sum_ + diff * diff
        displacement_ = sum_ / f32(naxisatom)
    #      write(6,*) im,displacement
    return displacement_
