"""Translation of fortran_src/sdfcube.f: ``program standalone`` (sdfcube),
``subroutine run_sdfcube``, ``initcube``, ``outcube`` and ``sdfcube`` -- the spatial density
function of second molecules around a central molecule, accumulated directly on a cubic
grid in the frame of the central molecule and written as a Gaussian .cube file.

The data of ``module sdfcube_mod`` (lines 1-13 of the Fortran file) is in
common/sdfcube_mod.py; ``common/sdfcube.py`` is the include file sdfcube.inc (imported
here as ``sdfcube_inc``, since ``sdfcube`` is the name of a routine of this file).
form_eulerangles is in library/sharm_routines_rev1.py.  ``setup_symbols_sdfcube`` is in
setup_symbols_generic.f (library/setup_symbols_generic.py, imported inside the routine:
another work package).

Quirks reproduced (docs/translation_notes/Q.md): in outcube ``fc1=1.0-fc1`` reads the
uninitialised fc1, which gfortran -O2 folds to ``fc1=1.0`` (Rust sdfcube.rs); for identical
central and second molecules ``averagemol`` is read from the ``axiss1`` symbol (the Fortran
has no ``iref=iref+2`` in that branch).
"""

from ..common import clustercom, commands, lattice, molecom, sdfcube_mod, sharm_sdf, symbols
from ..common import sdfcube as sdfcube_inc  # sdfcube.inc
from ..fortran import (
    Ch,
    FortranIOError,
    acosf,
    char_array,
    f32,
    findex,
    fmt,
    fstr,
    int_array,
    iostat_of,
    lw,
    nint,
    r4,
    read_internal,
    real_array,
    substr,
    units,
)
from ..library import (
    auxiliary_routines,
    file_and_directory_commands,
    minimagec,
    parse,
    read_and_save_files,
    reallocation_routines,
    sharm_routines_rev1,
)
from ..library import ato_read_and_save
from ..trace import routine

_F201 = "(/'sharm> Atom type ',a,' does not exist in .ato file ',a)"
_F333 = "(1x,'Number of type 1 molecules = ',i10/1x, 'Number of type 2 molecules = ',i10)"
_F502 = "(2(1x,i5))"
_F193 = "(1x,'rho = ',e13.6)"
_FXYZ = "(a3,1x,f10.5,1x,f10.5,1x,f10.5,1x,a3)"

# real, parameter :: pi=3.141592653,bohr=0.52917721092 (initcube, outcube)
_BOHR = r4(0.52917721092)


def main():
    """sdfcube.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_sdfcube()


def _parse(iref, nspc):
    """call parse(symbol_value(iref),nwords,ncf,ncl,nspc,256,nworddim)"""
    commands.nwords = parse.parse(symbols.symbol_value[iref], commands.nwords, commands.ncf, commands.ncl, nspc,
                                  256, commands.nworddim)


def _word(iref, ic):
    """symbol_value(iref)(ncf(ic):ncl(ic))"""
    return substr(symbols.symbol_value[iref], commands.ncf[ic], commands.ncl[ic])


@routine
def run_sdfcube():
    """sdfcube.f:run_sdfcube -- performs calculation of the spatial density function on a
    cube from .ato file.

    uses arithmetic centre of molecule rather than a particular atom to define
    molecular centre.
    """
    from ..library import setup_symbols_generic

    D = sdfcube_mod
    H = sdfcube_inc
    sv = symbols.symbol_value
    #
    # setup the SHARM symbols
    #
    setup_symbols_generic.setup_symbols_sdfcube()
    #
    # now setup the input file
    #
    nin = 11
    nout = 6
    # NB: ierrout is not initialised in the Fortran; setup_input_file sets it
    ierrout = 0
    nout, ierrout = file_and_directory_commands.setup_input_file(nin, nout, ierrout)
    #
    # if an error occurred on setup, exit now
    #
    if ierrout != 0:
        return
    #
    # read the input file
    #
    nout, ierrout = read_and_save_files.read_input_file(nin, nout, ierrout)
    #
    # if an error occurred on read, exit now
    #
    if ierrout != 0:
        return
    # Set up the name of the file where the CUBE data will be stored
    D.lenfname_cube = findex(symbols.fnameinp, '.SDFCUBE.dat')
    if D.lenfname_cube <= 0:
        return  # Apparently the wrong input file specified!
    D.fname_cube = fstr(substr(symbols.fnameinp, 1, D.lenfname_cube - 1) + '.SDFCUBE.cube', 256)
    D.lenfname_cube = len(D.fname_cube)
    iref = 1
    #
    # set up number of pixels
    #
    iref = iref + 1
    H.nz = read_internal(sv[iref]).int(H.nz)
    #
    # set up rmax of calculation
    #
    iref = iref + 1
    H.zmax = read_internal(sv[iref]).real(H.zmax)
    #
    # set up number of configurations already stored
    #
    iref = iref + 1
    H.nsumt = read_internal(sv[iref]).int(H.nsumt)
    #
    # for molecule 1 read number and type of atoms used to define the molecular
    # centre
    #
    iref = iref + 1
    _parse(iref, 5)
    H.ncent1 = 0
    #
    # skip through the atom types to find out the type number of the specified atom
    #
    for ic in range(1, commands.nwords + 1):
        it = 0
        ifound = 0
        while it < molecom.ntype and ifound == 0:
            it = it + 1
            if molecom.type[it] == _word(iref, ic):
                ifound = it
        if ifound <= 0:
            units.write(6, fmt(_F201, _word(iref, ic), substr(symbols.fname_ato, 1, symbols.ilenfname_ato)))
            return
        else:
            if H.ncent1 < H.mcentre:
                H.ncent1 = H.ncent1 + 1
                H.lctyp1[H.ncent1] = ifound
    #
    # for molecule 1 read primary axis type (1=X, 2=Y, 3=Z), and
    # the number and atom numbers of atoms used to define the molecule's nzt1
    # axis
    #
    iref = iref + 1
    _parse(iref, 5)
    H.nzt1 = fstr(_word(iref, 1), 2)
    H.ncz1 = commands.nwords - 1
    if H.ncz1 > H.mcentre:
        H.ncz1 = H.mcentre
    for i in range(1, H.ncz1 + 1):
        iread = i + 1
        H.lztyp1[i] = read_internal(_word(iref, iread)).int(H.lztyp1[i])
    #
    # for molecule 1 read secondary axis type (1=X, 2=Y, 3=Z), and
    # the number and atom numbers of atoms used to define the
    # molecule's nyt1 axis
    #
    iref = iref + 1
    _parse(iref, 5)
    H.nyt1 = fstr(_word(iref, 1), 2)
    H.ncy1 = commands.nwords - 1
    if H.ncy1 > H.mcentre:
        H.ncy1 = H.mcentre
    for i in range(1, H.ncy1 + 1):
        iread = i + 1
        H.lytyp1[i] = read_internal(_word(iref, iread)).int(H.lytyp1[i])

    # for molecule 2 read number and type of atoms used to define the molecular
    # centre. Also check whether this molecule is identical to the first

    iref = iref + 1
    _parse(iref, 5)
    H.ncent2 = 0
    #
    # skip through the atom types to find out the type number of the specified atom
    #
    for ic in range(1, commands.nwords + 1):
        it = 0
        ifound = 0
        while it < molecom.ntype and ifound == 0:
            it = it + 1
            if molecom.type[it] == _word(iref, ic):
                ifound = it
        if ifound == 0:
            units.write(6, fmt(_F201, _word(iref, ic), substr(symbols.fname_ato, 1, symbols.ilenfname_ato)))
            return
        else:
            if H.ncent2 < H.mcentre:
                H.ncent2 = H.ncent2 + 1
                H.lctyp2[H.ncent2] = ifound

    # test for identity of central and second molecules

    H.ldistinct = 0
    if H.ncent2 == H.ncent1:
        H.ldistinct = 0
        for it in range(1, H.ncent1 + 1):
            if H.lctyp2[it] != H.lctyp1[it]:
                H.ldistinct = 1
    else:
        H.ldistinct = 1

    units.write(6, lw(H.ncent1, H.ncent2, H.ldistinct))

    if H.ldistinct != 0:

        # if molecules are different, then repeat this for the second molecule

        # for molecule 2 read primary axis type (1=X, 2=Y, 3=Z), and
        # the number and atom numbers of atoms used to define the molecule's nzt1
        # axis

        iref = iref + 1
        _parse(iref, 5)
        H.nzt2 = fstr(_word(iref, 1), 2)
        H.ncz2 = commands.nwords - 1
        if H.ncz2 > H.mcentre:
            H.ncz2 = H.mcentre
        for i in range(1, H.ncz2 + 1):
            iread = i + 1
            H.lztyp2[i] = read_internal(_word(iref, iread)).int(H.lztyp2[i])

        # for molecule 2 read secondary axis type (1=X, 2=Y, 3=Z), and
        # the number and atom numbers of atoms used to define the
        # molecule's nyt1 axis

        iref = iref + 1
        _parse(iref, 5)
        H.nyt2 = fstr(_word(iref, 1), 2)
        H.ncy2 = commands.nwords - 1
        if H.ncy2 > H.mcentre:
            H.ncy2 = H.mcentre
        for i in range(1, H.ncy2 + 1):
            iread = i + 1
            H.lytyp2[i] = read_internal(_word(iref, iread)).int(H.lytyp2[i])
    else:

        # otherwise set all the values to the same as the central molecule

        H.ncent2 = H.ncent1
        H.nzt2 = H.nzt1
        H.nyt2 = H.nyt1
        if H.ncent2 > 0:
            for ic in range(1, H.ncent1 + 1):
                H.lctyp2[ic] = H.lctyp1[ic]
        H.ncz2 = H.ncz1
        for i in range(1, H.ncz1 + 1):
            H.lztyp2[i] = H.lztyp1[i]
        H.ncy2 = H.ncy1
        for i in range(1, abs(H.ncy1) + 1):
            H.lytyp2[i] = H.lytyp1[i]
        # NB: no iref=iref+2 here in the Fortran, so averagemol is read from axiss1 below
    D.performtest = False
    iref = iref + 1
    D.outputaveragemol = sv[iref] == 'y'
    # c     write(6,*) ncent1,ncent2,ncz1,ncy1,ncz2,ncy2
    #
    # only perform the calculation when nsum is zero or a multiple of nmult
    #
    nmult = 5
    nsum = H.nsumt // nmult
    nt = nsum * nmult
    if nt == H.nsumt:
        #
        # normalise atomic coordinates to box size
        #
        ato_read_and_save.normalise_coordinates()
        #
        # set up step size and z values
        #
        H.zstep = H.zmax / f32(H.nz)
        initcube(nin)
        sdfcube(nin, symbols.fnameinp)
        outcube(nin)
    H.nsumt = H.nsumt + 1
    #
    # write input file
    #
    #
    # make sure all current values are stored
    #
    sv[4] = fstr(lw(H.nsumt), 256)
    read_and_save_files.save_files_others(14)
    units.close(14)
    return


@routine
def initcube(nin):
    """sdfcube.f:initcube -- set up the cube grid (reading the step size from an existing
    cube file when nsumt > 0) and zero the density array."""
    D = sdfcube_mod
    H = sdfcube_inc
    # If nsumt.gt.0 we need to read in the current values from the specified cube file
    if H.nsumt > 0:  # Make sure our step sizes are the same as in the existing file
        try:
            units.open(nin, D.fname_cube, status='old')
            ierr = 0
        except FortranIOError as e:
            ierr = iostat_of(e)
        if ierr == 0:
            # Read two blank lines
            units.skip_record(nin)
            units.skip_record(nin)
            # Number of atoms in central molecule
            r = units.list_read(nin)
            D.natosave = r.int(D.natosave)
            r.end()
            # Number of pixels along each of x,y,z - we will assume this is the same number for
            # each axis
            r = units.list_read(nin)
            H.nz0 = r.int(H.nz0)
            H.zstep = r.real(H.zstep)  # zstep in Bohr units
            r.end()
            H.zstep = H.zstep * _BOHR  # Convert to Angstroms
            units.close(nin)
            H.nz = (H.nz0 - 1) // 2
            H.nx = H.nz
            H.ny = H.nz
            H.nz0 = 2 * H.nz + 1  # Always an odd number of pixels along each axis (nz0=2*nz+1)
            H.nx0 = H.nz0
            H.ny0 = H.nz0
            H.xstep = H.zstep
            H.ystep = H.zstep
        else:  # If file could not be opened, reset the distributions
            H.nsumt = 0
            H.nx = H.nz
            H.ny = H.nz
            H.nz0 = 2 * H.nz + 1  # Always an odd number of pixels along each axis (nz0=2*nz+1)
            H.nx0 = H.nz0
            H.ny0 = H.nz0
            H.xstep = H.zstep
            H.ystep = H.zstep
    else:
        H.nsumt = 0
        H.nx = H.nz
        H.ny = H.nz
        H.nz0 = 2 * H.nz + 1  # Always an odd number of pixels along each axis (nz0=2*nz+1)
        H.nx0 = H.nz0
        H.ny0 = H.nz0
        H.xstep = H.zstep
        H.ystep = H.zstep
    # Allocate the density array
    if D.rhocube is not None:
        D.rhocube = None
    D.rhocube = real_array(H.nz0, H.nz0, H.nz0)
    # Set all values to zero
    rhocube = D.rhocube
    for ix in range(1, H.nx0 + 1):
        for iy in range(1, H.ny0 + 1):
            for iz in range(1, H.nz0 + 1):
                rhocube[ix, iy, iz] = 0.0
    return


@routine
def outcube(nout):
    """sdfcube.f:outcube -- average the new density with the one in the existing cube
    file (nsumt > 0) and write the cube file.  Returns nout."""
    D = sdfcube_mod
    H = sdfcube_inc
    rhocube = D.rhocube
    rhocubet = real_array(10)

    # If nsumt.gt.0 we need to read in the current values from the specified cube file
    if H.nsumt > 0:  # Make sure our step sizes are the same as in the existing file
        try:
            units.open(nout, D.fname_cube, status='old')
            ierr = 0
        except FortranIOError as e:
            ierr = iostat_of(e)
        if ierr == 0:
            # Read two blank lines
            units.skip_record(nout)
            units.skip_record(nout)
            # Number of atoms in central molecule
            natosavet = 0
            r = units.list_read(nout)
            natosavet = r.int(natosavet)
            r.end()
            # Number of pixels along each of x,y,z - we will assume this is the same number for
            # each axis
            units.skip_record(nout)
            units.skip_record(nout)
            units.skip_record(nout)
            # Coordinates of central molecule - these will be rewritten so ignore current values
            for i in range(1, natosavet + 1):
                units.skip_record(nout)
            # Get the existing density distribution and add in the new values
            fc = r4(1.0) / f32(H.nsumt)
            # NB: fc1=1.0-fc1 reads the uninitialised fc1; gfortran -O2 folds it to 1.0
            # (checked in the reference object code, Rust sdfcube.rs)
            fc1 = r4(1.0)
            for ix in range(1, H.nx0 + 1):
                for iy in range(1, H.nz0 + 1):
                    for iz in range(1, H.nz0 + 1, 10):
                        itmax = min(9, H.nz0 - iz)
                        r = units.list_read(nout)
                        for it in range(0, itmax + 1):
                            rhocubet[1 + it] = r.real(rhocubet[1 + it])
                        r.end()
                        for it in range(0, itmax + 1):
                            izref = iz + it
                            rhocube[ix, iy, izref] = fc1 * rhocubet[1 + it] + fc * rhocube[ix, iy, izref]
            units.close(nout)
        else:  # If file could not be opened, reset the distributions
            H.nsumt = 0
    # Write out the revised values
    try:
        units.open(nout, D.fname_cube, status='unknown')
        ierr = 0
    except FortranIOError as e:
        ierr = iostat_of(e)
    if ierr == 0:
        units.write(nout, lw('# Cube file output by sdfcube'))
        units.write(nout, lw('Outer loop x, middle loop y, inner loop z, filname: ' + D.fname_cube))
        units.write(nout, lw(D.natosave, -(f32(H.nx) * H.xstep / _BOHR), -(f32(H.ny) * H.ystep / _BOHR),
                             -(f32(H.nz) * H.zstep / _BOHR)))
        units.write(nout, lw(H.nx0, H.xstep / _BOHR, f32(0.0), f32(0.0)))
        units.write(nout, lw(H.ny0, f32(0.0), H.ystep / _BOHR, f32(0.0)))
        units.write(nout, lw(H.nz0, f32(0.0), f32(0.0), H.zstep / _BOHR))
        for is_ in range(1, D.natosave + 1):
            iatomicno = auxiliary_routines.get_atomicnumber(D.atsymbolsave[is_])
            units.write(nout, lw(iatomicno, f32(iatomicno), D.xasave[is_] / _BOHR, D.yasave[is_] / _BOHR,
                                 D.zasave[is_] / _BOHR))
        for ix in range(1, H.nx0 + 1):
            for iy in range(1, H.ny0 + 1):
                for iz in range(1, H.nz0 + 1, 10):
                    itmax = min(9, H.nz0 - iz)
                    units.write(nout, fmt('(10(1x,e11.4))', *[rhocube[ix, iy, iz + it] for it in range(0, itmax + 1)]))
        units.close(nout)
    return nout


@routine
def sdfcube(nout, fname_inp):
    """sdfcube.f:sdfcube -- molecular centres and Euler angles, the averaged central and
    second molecules (xyz files), and the histogram of second-molecule centres in the frame
    of each central molecule (rhocube), normalised to the density.

    Returns nout.
    """
    D = sdfcube_mod
    H = sdfcube_inc
    M = molecom
    C = clustercom
    # real amass,amasssum
    # logical :: centralmolecule,found,foundatom
    # character(len=256) :: fname_xyz,textline
    # (xstepd,ystepd,zstepd are declared real here but are the sdfcube.inc COMMON variables)
    xacentral = yacentral = zacentral = None  # real, dimension(:), allocatable
    atsymbolcentral = attypecentral = None  # character(len=3), dimension(:), allocatable
    centraldefinition = None  # logical, dimension(:), allocatable
    xat = yat = zat = None  # real, dimension(:), allocatable
    natocentral = 0

    sharm_sdf.pi = r4(3.141592653)
    # Arrays to store Euler angles
    if D.nm is not None:
        D.nm = D.phim = D.cthetm = D.chim = D.xcent = D.ycent = D.zcent = None
    nmol = C.nmol
    D.nm = int_array(nmol)
    D.phim = real_array(nmol)
    D.cthetm = real_array(nmol)
    D.chim = real_array(nmol)
    D.xcent = real_array(nmol)
    D.ycent = real_array(nmol)
    D.zcent = real_array(nmol)
    #
    # set up coordinates of Laboratory X, Y and Z axes
    #
    alab = sharm_sdf.alab
    blab = sharm_sdf.blab
    clab = sharm_sdf.clab
    alab[1] = 1.0
    alab[2] = 0.0
    alab[3] = 0.0
    blab[1] = 0.0
    blab[2] = 1.0
    blab[3] = 0.0
    clab[1] = 0.0
    clab[2] = 0.0
    clab[3] = 1.0
    # Renormalised step sizes
    dmax = C.dmax
    H.xstepd = H.xstep / dmax
    H.ystepd = H.ystep / dmax
    H.zstepd = H.zstep / dmax
    #
    # set up coordinates of molecular centres and Euler angles for each molecule
    # in laboratory coordinate frame
    #
    nsmmo1 = 0
    nsmmo2 = 0
    for imol1 in range(1, nmol + 1):
        # define origin of this molecule
        xo1 = f32(0.0)
        yo1 = f32(0.0)
        zo1 = f32(0.0)
        amasssum = f32(0.0)
        ncref = 0  # noqa: F841
        D.nm[imol1] = 0
        #
        # test to see if this is molecule 1 or molecule 2. Need to check that all the specified
        # types occur - not just one
        #
        found = True
        ic = 0
        while ic < H.ncent1 and found:
            ic = ic + 1
            foundatom = False
            ia1 = C.naf[imol1] - 1
            while ia1 < C.nal[imol1] and not foundatom:
                ia1 = ia1 + 1
                foundatom = M.itype[ia1] == H.lctyp1[ic]
            if ic == 1:
                found = foundatom
            else:
                found = found and foundatom
        if found:
            D.nm[imol1] = 1
        else:  # molecule 1 not found, try molecule 2
            found = True
            ic = 0
            while ic < H.ncent2 and found:
                ic = ic + 1
                foundatom = False
                ia1 = C.naf[imol1] - 1
                while ia1 < C.nal[imol1] and not foundatom:
                    ia1 = ia1 + 1
                    foundatom = M.itype[ia1] == H.lctyp2[ic]
                if ic == 1:
                    found = foundatom
                else:
                    found = found and foundatom
            if found:
                D.nm[imol1] = 2
        if D.nm[imol1] == 1:
            nsmmo1 = nsmmo1 + 1
        if D.nm[imol1] == 2:
            nsmmo2 = nsmmo2 + 1
        if D.nm[imol1] == 1 and H.ncent1 > 0:
            for ia1 in range(C.naf[imol1], C.nal[imol1] + 1):
                for ic in range(1, H.ncent1 + 1):
                    if M.itype[ia1] == H.lctyp1[ic]:
                        amass = M.amassa[M.itype[ia1]]
                        xo1 = xo1 + M.xa[ia1] * amass
                        yo1 = yo1 + M.ya[ia1] * amass
                        zo1 = zo1 + M.za[ia1] * amass
                        amasssum = amasssum + amass
        elif D.nm[imol1] == 2 and H.ncent2 > 0:
            for ia1 in range(C.naf[imol1], C.nal[imol1] + 1):
                for ic in range(1, H.ncent2 + 1):
                    if M.itype[ia1] == H.lctyp2[ic]:
                        amass = M.amassa[M.itype[ia1]]
                        xo1 = xo1 + M.xa[ia1] * amass
                        yo1 = yo1 + M.ya[ia1] * amass
                        zo1 = zo1 + M.za[ia1] * amass
                        amasssum = amasssum + amass
        else:
            for ia1 in range(C.naf[imol1], C.nal[imol1] + 1):
                amass = M.amassa[M.itype[ia1]]
                xo1 = xo1 + M.xa[ia1] * amass
                yo1 = yo1 + M.ya[ia1] * amass
                zo1 = zo1 + M.za[ia1] * amass
                amasssum = amasssum + amass
        if amasssum > 0.0 and D.nm[imol1] != 0:
            xo1 = xo1 / amasssum
            yo1 = yo1 / amasssum
            zo1 = zo1 / amasssum
            # c     write(6,991) xo1,yo1,zo1
            D.xcent[imol1] = C.xadd[imol1] + xo1
            D.ycent[imol1] = C.yadd[imol1] + yo1
            D.zcent[imol1] = C.zadd[imol1] + zo1
            # Form atomic coordinates relative to this centre
            xat = reallocation_routines.reallocate1d_r(xat, C.nato[imol1])
            yat = reallocation_routines.reallocate1d_r(yat, C.nato[imol1])
            zat = reallocation_routines.reallocate1d_r(zat, C.nato[imol1])
            for ia in range(1, C.nato[imol1] + 1):
                iaref = ia + C.naf[imol1] - 1
                xat[ia] = M.xa[iaref] - xo1
                yat[ia] = M.ya[iaref] - yo1
                zat[ia] = M.za[iaref] - zo1
            if D.nm[imol1] == 1:
                D.phim[imol1], D.cthetm[imol1], D.chim[imol1] = sharm_routines_rev1.form_eulerangles(
                    imol1, xat, yat, zat, H.nzt1, H.nyt1, H.ncz1, H.ncy1, H.lztyp1, H.lytyp1, alab, blab, clab,
                    D.phim[imol1], D.cthetm[imol1], D.chim[imol1])
            elif D.nm[imol1] == 2:
                D.phim[imol1], D.cthetm[imol1], D.chim[imol1] = sharm_routines_rev1.form_eulerangles(
                    imol1, xat, yat, zat, H.nzt2, H.nyt2, H.ncz2, H.ncy2, H.lztyp2, H.lytyp2, alab, blab, clab,
                    D.phim[imol1], D.cthetm[imol1], D.chim[imol1])
        else:
            D.nm[imol1] = 0
        # 499   format(1x,i5,4(1x,e12.5))
    units.write(6, fmt(_F333, nsmmo1, nsmmo2))
    # Write out coordinates of centralmolecule. We will use the average coordinates of all
    # these molecules.
    im = 0
    ncentral = 0
    centralmolecule = False
    while im < nmol:
        im = im + 1
        centralmolecule = D.nm[im] == 1
        if centralmolecule:  # Write out coordinates of this molecule
            if ncentral == 0 or D.outputaveragemol:
                ncentral = ncentral + 1
                if ncentral == 1:  # First time allocate the arrays to store summed coordinates
                    natocentral = int(C.nato[im])
                    xacentral = real_array(natocentral)
                    yacentral = real_array(natocentral)
                    zacentral = real_array(natocentral)
                    atsymbolcentral = char_array(natocentral)
                    attypecentral = char_array(natocentral)
                    centraldefinition = [False] * (natocentral + 1)
                # Shift coordinates to centre used by SHARM/SDF
                xshift = C.xadd[im] - D.xcent[im]
                yshift = C.yadd[im] - D.ycent[im]
                zshift = C.zadd[im] - D.zcent[im]
                for ia in range(C.naf[im], C.nal[im] + 1):
                    # Temporary coordinates
                    xt = M.xa[ia] + xshift
                    yt = M.ya[ia] + yshift
                    zt = M.za[ia] + zshift
                    # Need to rotate these to align with the laboratory axes - performed in reverse
                    # order about the laboratory axes
                    xr, yr, zr = auxiliary_routines.rotateeuler(xt, yt, zt, -D.chim[im], -acosf(D.cthetm[im]),
                                                                -D.phim[im], None, None, None)
                    iaref = ia - C.naf[im] + 1
                    if ncentral == 1:
                        xacentral[iaref] = xr
                        yacentral[iaref] = yr
                        zacentral[iaref] = zr
                        atsymbolcentral[iaref] = M.atsymbol[M.itype[ia]]
                        attypecentral[iaref] = M.type[M.itype[ia]]
                        # Check whether this atom is part of the central molecule definition
                        centraldefinition[iaref] = _central_definition(ia, iaref)
                    else:
                        xacentral[iaref] = xacentral[iaref] + xr
                        yacentral[iaref] = yacentral[iaref] + yr
                        zacentral[iaref] = zacentral[iaref] + zr
    # Form average coordinates for central molecule
    if ncentral > 0:
        for ia in range(1, natocentral + 1):
            xacentral[ia] = xacentral[ia] * dmax / f32(ncentral)
            yacentral[ia] = yacentral[ia] * dmax / f32(ncentral)
            zacentral[ia] = zacentral[ia] * dmax / f32(ncentral)
        # Save these coordinates to output to the cube file
        if D.xasave is not None:
            D.xasave = D.yasave = D.zasave = D.atsymbolsave = None
        D.natosave = natocentral
        D.xasave = real_array(D.natosave)
        D.yasave = real_array(D.natosave)
        D.zasave = real_array(D.natosave)
        D.atsymbolsave = char_array(D.natosave)
        for ia in range(1, natocentral + 1):
            D.xasave[ia] = xacentral[ia]
            D.yasave[ia] = yacentral[ia]
            D.zasave[ia] = zacentral[ia]
            D.atsymbolsave[ia] = atsymbolcentral[ia]
        fname_xyz = fstr(substr(fname_inp, 1, auxiliary_routines.indexlastdot(fname_inp) - 1) + '-first.xyz', 256)
        _write_xyz(nout, fname_xyz, natocentral, atsymbolcentral, xacentral, yacentral, zacentral, attypecentral,
                   centraldefinition)
        if H.ldistinct == 0:  # Also write this out as the second molecule
            fname_xyz = fstr(substr(fname_inp, 1, auxiliary_routines.indexlastdot(fname_inp) - 1) + '-second.xyz',
                             256)
            _write_xyz(nout, fname_xyz, natocentral, atsymbolcentral, xacentral, yacentral, zacentral, attypecentral,
                       centraldefinition)
    # Now repeat this for second molecules if they are different from first molecules
    # Write out coordinates of centralmolecule. We will use the average coordinates of all
    # these molecules.
    im = 0
    ncentral = 0
    centralmolecule = False
    while im < nmol:
        im = im + 1
        centralmolecule = D.nm[im] == 2
        if centralmolecule:  # Write out coordinates of this molecule
            if ncentral == 0 or D.outputaveragemol:
                ncentral = ncentral + 1
                if ncentral == 1:  # First time allocate the arrays to store summed coordinates
                    natocentral = int(C.nato[im])
                    xacentral = real_array(natocentral)
                    yacentral = real_array(natocentral)
                    zacentral = real_array(natocentral)
                    atsymbolcentral = char_array(natocentral)
                    attypecentral = char_array(natocentral)
                    centraldefinition = [False] * (natocentral + 1)
                # Shift coordinates to centre used by SHARM/SDF
                xshift = C.xadd[im] - D.xcent[im]
                yshift = C.yadd[im] - D.ycent[im]
                zshift = C.zadd[im] - D.zcent[im]
                for ia in range(C.naf[im], C.nal[im] + 1):
                    # Temporary coordinates
                    xt = M.xa[ia] + xshift
                    yt = M.ya[ia] + yshift
                    zt = M.za[ia] + zshift
                    # Need to rotate these to align with the laboratory axes - performed in reverse
                    # order about the laboratory axes
                    xr, yr, zr = auxiliary_routines.rotateeuler(xt, yt, zt, -D.chim[im], -acosf(D.cthetm[im]),
                                                                -D.phim[im], None, None, None)
                    iaref = ia - C.naf[im] + 1
                    if ncentral == 1:
                        xacentral[iaref] = xr
                        yacentral[iaref] = yr
                        zacentral[iaref] = zr
                        atsymbolcentral[iaref] = M.atsymbol[M.itype[ia]]
                        attypecentral[iaref] = M.type[M.itype[ia]]
                        # Check whether this atom is part of the central molecule definition
                        # (the Fortran tests the first-molecule definition here too)
                        centraldefinition[iaref] = _central_definition(ia, iaref)
                    else:
                        xacentral[iaref] = xacentral[iaref] + xr
                        yacentral[iaref] = yacentral[iaref] + yr
                        zacentral[iaref] = zacentral[iaref] + zr
    # Form average coordinates for central molecule
    if ncentral > 0:
        for ia in range(1, natocentral + 1):
            xacentral[ia] = xacentral[ia] * dmax / f32(ncentral)
            yacentral[ia] = yacentral[ia] * dmax / f32(ncentral)
            zacentral[ia] = zacentral[ia] * dmax / f32(ncentral)
        fname_xyz = fstr(substr(fname_inp, 1, auxiliary_routines.indexlastdot(fname_inp) - 1) + '-second.xyz', 256)
        _write_xyz(nout, fname_xyz, natocentral, atsymbolcentral, xacentral, yacentral, zacentral, attypecentral,
                   centraldefinition)
    #
    # initialize counters for number of molecules of each kind in the system
    #
    if H.ldistinct == 0:
        if D.nm[nmol] == 1 or D.nm[nmol] == 2:
            lsum1 = 1
            lsum2 = 1
        else:
            lsum1 = 0
            lsum2 = 0
    else:
        lsum1 = 0
        lsum2 = 0
        if D.nm[nmol] == 1:
            lsum1 = 1
        if D.nm[nmol] == 2:
            lsum2 = 1
    units.write(6, lw(nmol, D.nm[nmol], lsum1, lsum2))
    nm = D.nm
    xcent = D.xcent
    ycent = D.ycent
    zcent = D.zcent
    rhocube = D.rhocube
    xstepd = H.xstepd
    ystepd = H.ystepd
    zstepd = H.zstepd
    nx, ny, nz = H.nx, H.ny, H.nz
    nx0, ny0, nz0 = H.nx0, H.ny0, H.nz0
    for imol1 in range(1, nmol - 1 + 1):
        if nm[imol1] != 0:
            units.write(6, lw(imol1))
            # c     call system_clock(time3,rate,max)
            #
            # increment counters for number of molecules of types 1 and 2 in the system
            #
            if H.ldistinct == 0:
                if nm[imol1] == 1 or nm[imol1] == 2:
                    lsum1 = lsum1 + 1
                    lsum2 = lsum2 + 1
            else:
                if nm[imol1] == 1:
                    lsum1 = lsum1 + 1
                if nm[imol1] == 2:
                    lsum2 = lsum2 + 1
            xc1 = xcent[imol1]
            yc1 = ycent[imol1]
            zc1 = zcent[imol1]
            #
            # calc now look at second molecules
            #
            nm2f = imol1 + 1
            for imol2 in range(nm2f, nmol + 1):
                if imol1 != imol2 and nm[imol2] != 0:
                    iflag1 = 0
                    iflag2 = 0
                    nm1 = nm[imol1]
                    nm2 = nm[imol2]
                    if H.ldistinct == 0 and nm1 == 1 and nm2 == 1:
                        iflag1 = 1
                        iflag2 = 1
                    if H.ldistinct == 1 and nm1 == 1 and nm2 == 2:
                        iflag1 = 1
                        iflag2 = 0
                    if H.ldistinct == 1 and nm1 == 2 and nm2 == 1:
                        iflag1 = 0
                        iflag2 = 1
                    xc2 = xcent[imol2]
                    yc2 = ycent[imol2]
                    zc2 = zcent[imol2]
                    xdif = xc2 - xc1
                    ydif = yc2 - yc1
                    zdif = zc2 - zc1
                    # Apply periodic boundary conditions
                    xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                    if iflag1 == 1:
                        # Need to rotate coordinates of this molecule onto the frame of the first
                        # molecule
                        xdifr, ydifr, zdifr = auxiliary_routines.rotateeuler(
                            xdif, ydif, zdif, -D.chim[imol1], -acosf(D.cthetm[imol1]), -D.phim[imol1], None, None, None)
                        xdifr, ydifr, zdifr = minimagec.minimagec(xdifr, ydifr, zdifr)
                        #                     xdifr=xdif
                        ix = nx + nint(xdifr / xstepd) + 1
                        if ix > 0 and ix <= nx0:
                            iy = ny + nint(ydifr / ystepd) + 1
                            if iy > 0 and iy <= ny0:
                                iz = nz + nint(zdifr / zstepd) + 1
                                if iz > 0 and iz <= nz0:
                                    #  write(6,*) xstepd,ystepd,zstepd
                                    #  write(6,*) imol1,imol2,xdifr,ydifr,zdifr,ix,iy,iz
                                    rhocube[ix, iy, iz] = rhocube[ix, iy, iz] + r4(1.0)
                    # We add on the case when molecule2 is type 1 and at the origin
                    if iflag2 == 1:
                        xdifr, ydifr, zdifr = auxiliary_routines.rotateeuler(
                            -xdif, -ydif, -zdif, -D.chim[imol2], -acosf(D.cthetm[imol2]), -D.phim[imol2], None, None,
                            None)
                        xdifr, ydifr, zdifr = minimagec.minimagec(xdifr, ydifr, zdifr)
                        #                     xdifr=-xdif
                        ix = nx + nint(xdifr / xstepd) + 1
                        if ix > 0 and ix <= nx0:
                            iy = ny + nint(ydifr / ystepd) + 1
                            if iy > 0 and iy <= ny0:
                                iz = nz + nint(zdifr / zstepd) + 1
                                if iz > 0 and iz <= nz0:
                                    #  write(6,*) xstepd,ystepd,zstepd
                                    #  write(6,*) imol2,imol1,xdifr,ydifr,zdifr,ix,iy,iz
                                    rhocube[ix, iy, iz] = rhocube[ix, iy, iz] + r4(1.0)
                        #                     write(6,*) imol2,imol1,xdifr,ydifr,zdifr,ix,iy,iz
            # c     call system_clock(time4,rate,max)
            # c     write(6,*) rate,max,time4-time3,nsec
    units.write(6, fmt(_F502, lsum1, lsum2))
    # renormalize g(r)'s to density of second atoms in sphere
    vol = lattice.alatx * (lattice.blaty * lattice.clatz)  # vol: a local implicit real here
    rho = f32(lsum1) * f32(lsum2) / vol
    units.write(6, fmt(_F193, rho))
    # volume of element
    gfac = xstepd * ystepd * zstepd
    gfac1 = r4(1.0) / (gfac * rho)
    for ix in range(1, nx0 + 1):
        for iy in range(1, ny0 + 1):
            for iz in range(1, nz0 + 1):
                rhocube[ix, iy, iz] = rhocube[ix, iy, iz] * gfac1
    return nout


def _write_xyz(nout, fname_xyz, natocentral, atsym, xs, ys, zs, attype, centraldef):
    """The coordinate-file block of sdfcube (repeated three times in the Fortran):
    open(nout,...,iostat=ierr); if(ierr.eq.0) then ... close(nout)."""
    try:
        units.open(nout, fname_xyz, status='unknown')
        ierr = 0
    except FortranIOError as e:
        ierr = iostat_of(e)
    if ierr == 0:
        # Write number of atoms in molecule
        units.write(nout, fmt('(i6)', natocentral))
        # Write information line - not used for plotting
        units.write(nout, fmt('(a)', 'Central molecule coordinate file: ' + fname_xyz))
        for ia in range(1, natocentral + 1):
            textline = fstr(fmt(_FXYZ, Ch(atsym[ia], 3), xs[ia], ys[ia], zs[ia], Ch(attype[ia], 3)), 256)
            # Append a comment if this atom is part of the central molecule definition
            if centraldef[ia]:
                textline = fstr(textline + '  #CENTRALMOLECULE', 256)
            units.write(nout, fmt('(a)', textline))
        units.close(nout)
    return ierr


def _central_definition(ia, iaref):
    """The test 'Check whether this atom is part of the central molecule definition' of
    sdfcube (always with the first molecule's ncent1, lctyp1, ncz1, lztyp1, ncy1, lytyp1).
    Returns found."""
    H = sdfcube_inc
    found = False
    ic = 0
    while not found and ic < H.ncent1:
        ic = ic + 1
        found = molecom.itype[ia] == H.lctyp1[ic]
    iz = 0
    while not found and iz < H.ncz1:
        iz = iz + 1
        found = iaref == H.lztyp1[iz]
    iy = 0
    while not found and iy < H.ncy1:
        iy = iy + 1
        found = iaref == H.lytyp1[iy]
    return bool(found)
