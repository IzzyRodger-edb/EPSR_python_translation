"""Translation of fortran_src/readpdb.f: ``program standalone``, ``subroutine run_readpdb``,
``subroutine readpdb`` and ``function upcase`` -- convert a .pdb file into a .ato file.

The reference binary is the *patched* readpdb (tools/build_patched_readxyz_readpdb.sh):
the original calls ``ctype`` with 11 of its 13 arguments and no explicit interface and
crashes at that call; the patch adds an INTERFACE with the last five arguments
``optional``.  The Python reproduces the patched program: the 11 positional arguments
fill a1 ... radminsph, radmaxsph and pseudofacsph are absent.

Quirks reproduced (docs/rust_port/TOOLS.md): ``dihedral(im, ...)`` is called with
absolute atom numbers (correct only for molecule 1); a molecule of one atom never has
``nal`` set (it keeps the value already in the array, normally 0).
"""

from ..common import chemical_symbols, clustercom, commands, epsr_dimension, lattice, molecom, molecom_mod, symbols
from ..fortran import (
    PI4ATAN,
    Ch,
    FortranEOF,
    FortranIOError,
    char_array,
    f32,
    findex,
    fmt,
    fstr,
    int_array,
    iostat_of,
    lw,
    r4,
    read_internal,
    real_array,
    setsub,
    sqrtf,
    substr,
    units,
)
from ..fortran.listio import read_formatted
from ..library import ato_manipulation_rev4
from ..library import ato_read_and_save
from ..library import auxiliary_routines
from ..library import file_and_directory_commands
from ..library import molecom_mod as molecom_mod_lib
from ..library import parse
from ..trace import routine

# FORMAT statements of run_readpdb ('$' dropped: written with advance=False)
_F499 = "(/'readpdb> No .pdb file specified! - routine will exit'/)"
_F599 = "(/'readpdb> Problems with .pdb file ',a/'Routine will exit now'/)"
_F115 = "(/'readpdb> Specify the chemical symbol for atom type ',a3,': ')"
# FORMAT statements of readpdb
_F201 = "(/' readpdb> Following molecule types defined initially:-'/)"
_F202 = "(1x,i3,1x,a3,1x,a10,2(1x,e12.5))"
_F48 = ("(/' readpdb> Number of van der Waals radii read in: ',i3"
        "//' Give van der Waals factor (suggest 0.7). This is used to test whether two atoms are bonded or not: ')")
_F49 = "(/' readpdb> Create dihedral angles (y or no): ')"
_FROT = "(i5/(10(1x,i6)))"


def main():
    """readpdb.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_readpdb()


@routine
def run_readpdb():
    """readpdb.f:run_readpdb -- read a .pdb file and write the corresponding .ato file.

    2006-11-15 Output options are (1) GNUplot, (2) PGPlot, (3) ATOMS.  Also atoms to be
    plotted will be defined by their types, not atomic symbol.
    """
    ncf = commands.ncf
    ncl = commands.ncl
    # character*256 readfile ; character*4 ext ; character*1 iaccept
    iaccept = ""  # NB: not initialised in the Fortran; search_for_files sets it
    # character*3020 sometextsave ; character*4 yesorno
    # integer iread,nwordssave,ncfsave(nworddim),nclsave(nworddim)
    rho = f32(0.0)  # implicitly typed real local
    ierr = 0

    # First identify the number of arguments read in

    commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 5, commands.nstsize,
                                  commands.nworddim)

    # generate the name of the .pdb file

    ext = ".pdb"
    ilenext = 4
    if commands.nwords < 2:
        commands.nwords = 0
        while commands.nwords <= 0:
            iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)
            commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 5, commands.nstsize,
                                          commands.nworddim)
            if commands.nwords <= 0:
                units.write(6, fmt(_F499))
                return
        readfile = fstr(substr(commands.sometext, ncf[1], ncl[1]), 256)
        ilenreadfile = ncl[1] - ncf[1] + 1
    else:
        readfile = fstr(substr(commands.sometext, ncf[2], ncl[2]), 256)
        ilenreadfile = ncl[2] - ncf[2] + 1
    #
    # strip out any extension that was supplied - extension must be .xyz, otherwise
    # nothing is done.
    #
    ilendot = auxiliary_routines.indexlastdot(readfile) - 1
    if ilendot > 0:
        readfile = fstr(substr(readfile, 1, ilendot), 256)
        ilenreadfile = ilendot
    #
    # create the full name of the input file by adding on the directory name if not present
    #
    symbols.fnameato, symbols.ilenfnameato = file_and_directory_commands.get_full_filename(
        readfile, ilenreadfile, symbols.fnameato, symbols.ilenfnameato)
    symbols.fnameato = fstr(substr(symbols.fnameato, 1, symbols.ilenfnameato) + ".pdb", 256)
    symbols.ilenfnameato = symbols.ilenfnameato + 4
    #
    # initialise number of atom and molecule types to zero
    #
    molecom.ntype = 0
    clustercom.nmoltype = 0
    #
    # read .pdb file
    #
    molecom.sizefactor = f32(1.0)
    rho, ierr = readpdb(symbols.fnameato, rho, ierr)
    if rho == 0.0 or ierr != 0:
        units.write(6, fmt(_F599, Ch(symbols.fnameato, 256)))
        return
    #     call normalise_coordinates
    #
    # step through the atom types and make sure the chemical symbol has been specified
    #
    for it in range(1, molecom.ntype + 1):
        if molecom.atsymbol[it] == "XXX" or molecom.atsymbol[it] == "":
            units.write(6, fmt(_F115, Ch(molecom.type[it], 3)), advance=False)
            molecom.atsymbol[it] = fstr(units.list_read(5).word(molecom.atsymbol[it]), 3)

    # Ensure the coordinates of each atom in each molecule are relative to the centre of mass

    for im in range(1, clustercom.nmol + 1):

        ato_manipulation_rev4.centreofmass(im)

    # Replace the .pdb extension with the .ato extension

    ilendot = auxiliary_routines.indexlastdot(symbols.fnameato) - 1
    symbols.fnameato = fstr(substr(symbols.fnameato, 1, ilendot) + ".ato", 256)
    ato_read_and_save.moloutc(symbols.fnameato)

    return


def _read_word(text, old, n):
    """``read(text,*,iostat=ierr) c`` into a character*n variable: (value, iostat)."""
    try:
        return fstr(read_internal(text).word(old), n), 0
    except (FortranEOF, FortranIOError) as e:
        return old, iostat_of(e)


def _read_int(text, old):
    """``read(text,*,iostat=ierr) i``: (value, iostat)."""
    try:
        return read_internal(text).int(old), 0
    except (FortranEOF, FortranIOError) as e:
        return old, iostat_of(e)


def _read_reals_fmt(fmtstr, text, olds):
    """``read(text,'(nfw.d)',iostat=ierr) a,b,c``: (values, iostat)."""
    try:
        vals = read_formatted(fmtstr, text)
        return [f32(v) for v in vals[:len(olds)]], 0
    except (FortranEOF, FortranIOError, ValueError) as e:
        return list(olds), (iostat_of(e) if not isinstance(e, ValueError) else 5010)


def _do_end(first, last):
    """The value of a DO variable after ``do v=first,last`` completes (Fortran)."""
    return max(first, last + 1)


@routine
def readpdb(readfile, rho, ierr):
    """readpdb.f:readpdb -- read the CRYST1, ATOM and CONECT records of a .pdb file into
    the molecule arrays, then build bonds, bond angles and (optionally) dihedral groups.

    Returns (rho, ierr).
    """
    # use molecom_mod
    # INTERFACE get_rotgrp / get_whole_rotgrp (allocatable intent(out) labrott): the
    #   Python routines return (labrott, nrott) (docs/SIGNATURE_EXCEPTIONS.md)
    # INTERFACE upcase
    matom = epsr_dimension.matom
    mmol = epsr_dimension.mmol
    near = epsr_dimension.near
    melements = chemical_symbols.melements
    chsymbol = chemical_symbols.chsymbol
    chatwt = chemical_symbols.chatwt
    ncf = commands.ncf
    ncl = commands.ncl
    mol = molecom
    cc = clustercom
    lat = lattice
    mm = molecom_mod
    # character*256 readfile
    # integer ilenreadfile,ierr,iline,natomexpect
    # real vanderWaalsFactor,rho,xdif,ydif,zdif,pi,piconv,alphaaxes,betaaxes,gammaaxes,dtest,dvanderwaals
    vanderWaalsFactor = f32(0.0)  # NB: uninitialised in the Fortran; reset to 0.7 if the read fails
    # real, dimension(melements) :: vanderWaalsRadius
    vanderWaalsRadius = real_array(melements)
    # integer, allocatable :: iatomicno,ialabel
    iatomicno = None
    ialabel = None
    # character*2 chsymbolsave ; character(len=2), allocatable :: chsymboltemp
    chsymbolsave = ""  # NB: uninitialised in the Fortran until read
    chsymboltemp = None
    # character*3 a1,ats,yesorno ; character*4 alabelsave ; character*10 junk
    alabelsave = ""
    # integer iexch,j,k,ir,ia,ja,ifirst ; integer, dimension(4) :: isecond
    isecond = int_array(4)
    # external cosd,sind  -> EPSR's own (auxiliary_routines)
    cosd = auxiliary_routines.cosd
    sind = auxiliary_routines.sind
    # logical found,dodihedrals
    # temporary rotational groups
    nrott1 = nrott2 = 0
    labrott1 = labrott2 = None
    # !Temporary arrays for moving atoms around
    nbondmove = imageatommove = ialabelmove = None
    xamove = yamove = zamove = None
    chsymbolmove = alabelmove = None
    # implicitly typed real locals
    xat = yat = zat = f32(0.0)
    alpha = beta = gam = f32(0.0)
    cell = rmax = f32(0.0)  # noqa: F841

    pi = PI4ATAN
    piconv = r4(180.0) / pi
    ilenreadfile = len(readfile)  # noqa: F841
    try:
        units.open(10, readfile, status="old")
        ierr = 0
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    # !      write(6,*) ' readpdb> Open file ',ierr
    if ierr != 0:
        return rho, ierr  # Any errors opening file, return
    iline = 0
    natomexpect = 0  # noqa: F841
    mol.sizefactor = f32(1.0)
    lat.cubic = False
    cc.nmol = 0
    mol.simtk = f32(300.0)
    mol.tempb = f32(0.65e2)
    mol.tempang = f32(1.0)
    mol.tempdih = f32(10.0)
    lat.phib = f32(90.0)
    lat.thetac = f32(0.0)
    lat.phic = f32(0.0)
    mol.tol = f32(0.0)
    mol.tolout = f32(0.0)
    # !Ensure no atom moves are attempted
    mol.stepmi = f32(0.0)
    mol.stepri = f32(0.0)
    mol.stepr = f32(0.0)
    mol.stepm = f32(0.0)
    mol.tempangreq = f32(1.0)
    mol.tempdihreq = f32(10.0)
    ia = 0
    im = 0
    ig = 0  # noqa: F841
    mol.natom = 0
    mol.natomr = 0
    mm.ngrrot = 0
    mm.nrottot = 0
    mol.ntype = 0
    # ! Reset array size counters
    mm.ngrrotmax = 0
    mm.nrottotmax = 0
    imsave = 0
    # if(allocated(chsymboltemp)) deallocate(chsymboltemp,ialabel)
    chsymboltemp = char_array(matom)
    ialabel = int_array(matom)
    nbond = mol.nbond
    imageatom = mol.imageatom
    mref = mol.mref
    xa = mol.xa
    ya = mol.ya
    za = mol.za
    naf = cc.naf
    nal = cc.nal
    nato = cc.nato
    while ierr == 0 and ia < matom and im < mmol:
        try:
            commands.sometext = fstr(units.read_line(10), commands.nstsize)
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = iostat_of(e)
        # !         write(6,*) ' readpdb> Read line: ',iline,ierr
        if ierr == 0:
            sometext = commands.sometext
            # !            call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
            iline = iline + 1
            if findex(substr(sometext, 1, 6), "CRYST1") > 0:
                (lat.alat, lat.blat, lat.clat), ierr = _read_reals_fmt(
                    "(3f9.3)", substr(sometext, 7, 33), (lat.alat, lat.blat, lat.clat))
                (alpha, beta, gam), ierr = _read_reals_fmt("(3f7.2)", substr(sometext, 34, 54), (alpha, beta, gam))
                while ierr != 0:
                    units.write(6, lw("readpdb> Need a,b,c axes and 3 angles (alpha,beta,gamma):"))
                    try:
                        r = units.list_read(5)
                        lat.alat = r.real(lat.alat)
                        lat.blat = r.real(lat.blat)
                        lat.clat = r.real(lat.clat)
                        alpha = r.real(alpha)
                        beta = r.real(beta)
                        gam = r.real(gam)
                        ierr = 0
                    except (FortranEOF, FortranIOError) as e:
                        ierr = iostat_of(e)
                # ! Convert angles to internal spherical polar coordinates
                lat.phib, lat.thetac, lat.phic = auxiliary_routines.convertalphabetagamma(
                    alpha, beta, gam, lat.phib, lat.thetac, lat.phic)
                # set up lattice parameters for applying the minimum image convention
                lat.alatx = lat.alat
                lat.blatx = lat.blat * cosd(lat.phib)
                lat.blaty = lat.blat * sind(lat.phib)
                lat.clatx = lat.clat * sind(lat.thetac) * cosd(lat.phic)
                lat.claty = lat.clat * sind(lat.thetac) * sind(lat.phic)
                lat.clatz = lat.clat * cosd(lat.thetac)
                # define the Cartesian coordinates of the lattice vectors
                lat.alatc[1] = lat.alatx
                lat.alatc[2] = 0.0
                lat.alatc[3] = 0.0
                lat.blatc[1] = lat.blatx
                lat.blatc[2] = lat.blaty
                lat.blatc[3] = 0.0
                lat.clatc[1] = lat.clatx
                lat.clatc[2] = lat.claty
                lat.clatc[3] = lat.clatz
                lat.bxax = lat.blatx / lat.alatx
                lat.cxax = lat.clatx / lat.alatx
                lat.cyby = lat.claty / lat.blaty
                # calculate the perpendicular distance between the faces of the lattice
                perpab = auxiliary_routines.perpd(lat.alatc, lat.blatc, lat.clatc)
                perpbc = auxiliary_routines.perpd(lat.blatc, lat.clatc, lat.alatc)
                perpca = auxiliary_routines.perpd(lat.clatc, lat.alatc, lat.blatc)
                units.write(6, lw(perpab, perpbc, perpca))
                # set up maximum radius to be used for Fourier transforms and to truncate
                # the potential. This is defined as the minimum perpendicular distance between
                # planes.
                #
                # G(r) and potentials will be calculated to rmax, or maximum allowed distance
                # based on number of points and step size, whichever is smaller.
                cc.dmax = min(perpab, perpbc, perpca)
                units.write(6, lw(cc.dmax))
                cell = cc.dmax
                #
                # setup molecular size factor
                #
                mol.dmaxmol = cc.dmax * mol.sizefactor
                #
                # determine maximum distance across box
                #
                # G(r) and potentials will be calculated to rmax, or maximum allowed distance
                # based on number of points and step size, whichever is smaller.
                #
                xmax = abs(lat.alatx) + abs(lat.blatx) + abs(lat.clatx)
                ymax = abs(lat.blaty) + abs(lat.claty)
                zmax = abs(lat.clatz)
                rmax = f32(0.5) * sqrtf(xmax * xmax + ymax * ymax + zmax * zmax)
                units.write(6, lw(rmax, cell))
                #
                # convert temperature in K to kJ/mole
                #
                mol.simt = f32(0.008314) * mol.simtk
                #
                # calculate volume of box
                #
                cc.volume = lat.alatx * lat.blaty * lat.clatz
            elif findex(substr(sometext, 1, 6), "ATOM") > 0:
                # !Otherwise we read the list of atomic coordinates
                ia, ierr = _read_int(substr(sometext, 7, 11), ia)
                chsymbolsave, ierr = _read_word(substr(sometext, 13, 14), chsymbolsave, 2)
                alabelsave, ierr = _read_word(substr(sometext, 13, 16), alabelsave, 4)
                im, ierr = _read_int(substr(sometext, 23, 26), im)
                (xat, yat, zat), ierr = _read_reals_fmt("(3f8.3)", substr(sometext, 31, 54), (xat, yat, zat))
                units.write(6, fmt("(i5,a2,i5,a4)", ia, Ch(chsymbolsave, 2), im, Ch(alabelsave, 4)))
                if ierr == 0:
                    mol.natom = mol.natom + 1
                    natom = mol.natom
                    nbond[natom] = 0
                    imageatom[natom] = 0
                    ialabel[natom] = ia  # The value of ia read in is kept as a (temporary) label to identify this atom uniquely
                    chsymboltemp[natom] = chsymbolsave
                    mol.alabel[natom] = fstr(alabelsave, 3)
                    xa[natom] = xat
                    ya[natom] = yat
                    za[natom] = zat
                    #
                    # do not count atoms which begin with 'Q' or 'q' - these are assumed
                    # invisible to the scattering radiation and experiment (although they
                    # can be used to calculate other structural quantities if necessary
                    #
                    if im != imsave:
                        cc.nmol = cc.nmol + 1
                        nmol = cc.nmol
                        naf[nmol] = natom
                        cc.xadd[nmol] = 0.0
                        cc.yadd[nmol] = 0.0
                        cc.zadd[nmol] = 0.0
                        nato[nmol] = 1
                        cc.ngrot[nmol] = 0
                        mm.ngrrot = mm.ngrrot + int(cc.ngrot[nmol])
                        cc.tethered[nmol] = False
                        cc.phix[nmol] = 0.0
                        cc.phiy[nmol] = 0.0
                        cc.phiz[nmol] = 0.0
                    else:
                        nal[cc.nmol] = natom
                        nato[cc.nmol] = nato[cc.nmol] + 1
                    mref[natom] = cc.nmol
                    imsave = im
            elif findex(substr(sometext, 1, 6), "CONECT") > 0 and len(sometext) > 15:
                # It will be assumed any atom numbers read in here have already been defined in a
                # !corresponding ATOM line. Otherwise this line will be ignored.
                # !Determine the number of numbers on this line
                nsecond = int((len(sometext) - 6) // 5) - 1
                ifirst = read_internal(substr(sometext, 7, 11)).int(0)
                # !Get the internal atom number associated with this label
                found = False
                ia = 0
                while not found and ia < mol.natom:
                    ia = ia + 1
                    found = ialabel[ia] == ifirst
                if found:
                    ifirst = ia
                ir1 = 12
                for i in range(1, nsecond + 1):
                    ir2 = ir1 + 4
                    isecond[i] = read_internal(substr(sometext, ir1, ir2)).int(isecond[i])
                    # !Get the internal atom number associated with this label
                    found = False
                    ia = 0
                    while not found and ia < mol.natom:
                        ia = ia + 1
                        found = ialabel[ia] == isecond[i]
                    if found:
                        isecond[i] = ia
                    ir1 = ir2 + 1
                items = [" readpdb> CONECT: ", ifirst, mref[ifirst]]
                for i in range(1, nsecond + 1):
                    items += [isecond[i], mref[isecond[i]]]
                units.write(6, lw(*items))
                i = 1
                while i <= nsecond and ifirst <= mol.natom and isecond[i] <= mol.natom:
                    if mref[ifirst] != mref[isecond[i]]:
                        # !Need to merge molecules
                        immin = int(min(mref[ifirst], mref[isecond[i]]))
                        immax = int(max(mref[ifirst], mref[isecond[i]]))
                        # !Put all the atoms in immax into immin
                        # !Check whether atom numbers in immax are contiguous with those in immin
                        if naf[immax] != (nal[immin] + 1):
                            # !If not contiguous, need to reorder atoms in list between immin and immax
                            # !First save the data in immax
                            n = int(nato[immax])
                            nbondmove = int_array(n)
                            imageatommove = int_array(n)
                            ialabelmove = int_array(n)
                            chsymbolmove = char_array(n)
                            alabelmove = char_array(n)
                            xamove = real_array(n)
                            yamove = real_array(n)
                            zamove = real_array(n)
                            for ia in range(1, n + 1):
                                iaref = ia + naf[immax] - 1
                                units.write(6, lw(immax, ia, iaref, ialabel[iaref]))
                                nbondmove[ia] = nbond[iaref]
                                imageatommove[ia] = imageatom[iaref]
                                ialabelmove[ia] = ialabel[iaref]
                                chsymbolmove[ia] = chsymboltemp[iaref]
                                alabelmove[ia] = mol.alabel[iaref]
                                xamove[ia] = xa[iaref]
                                yamove[ia] = ya[iaref]
                                zamove[ia] = za[iaref]
                            # !Now move atoms, starting from the one above naf(immax) and going backwards
                            ia = int(naf[immax]) - 1
                            mrefsave = int(mref[ia])
                            nal[mrefsave] = ia + nato[immax]
                            while ia > nal[immin]:
                                iaref = ia + int(nato[immax])
                                units.write(6, lw(ia, iaref, ialabel[ia]))
                                nbond[iaref] = nbond[ia]
                                imageatom[iaref] = imageatom[ia]
                                ialabel[iaref] = ialabel[ia]
                                chsymboltemp[iaref] = chsymboltemp[ia]
                                mol.alabel[iaref] = mol.alabel[ia]
                                xa[iaref] = xa[ia]
                                ya[iaref] = ya[ia]
                                za[iaref] = za[ia]
                                if mref[ia] != mrefsave:
                                    # !Remember to reset the integers defining the atom numbers for each molecule.
                                    mrefsave = int(mref[ia])
                                    nal[mrefsave] = iaref
                                else:
                                    naf[mrefsave] = iaref
                                mref[iaref] = mrefsave  # The molecule numbers for this region do not change
                                ia = ia - 1
                            # !Finally insert the molecule in the space created
                            for ia in range(1, n + 1):
                                iaref = ia + nal[immin]
                                units.write(6, lw(ia, iaref, ialabelmove[ia]))
                                nbond[iaref] = nbondmove[ia]
                                imageatom[iaref] = imageatommove[ia]
                                ialabel[iaref] = ialabelmove[ia]
                                chsymboltemp[iaref] = chsymbolmove[ia]
                                mol.alabel[iaref] = alabelmove[ia]
                                xa[iaref] = xamove[ia]
                                ya[iaref] = yamove[ia]
                                za[iaref] = zamove[ia]
                                mref[iaref] = immin
                            ia = _do_end(1, n)
                            nbondmove = imageatommove = ialabelmove = chsymbolmove = alabelmove = None
                            xamove = yamove = zamove = None
                        nato[immin] = nato[immin] + nato[immax]
                        nal[immin] = nal[immin] + nato[immax]
                        # !Since we haven't defined any rotational or dihedral groups at this stage we don't need to worry about them
                        for ia in range(naf[immin], nal[immin] + 1):
                            mref[ia] = immin
                        ia = _do_end(int(naf[immin]), int(nal[immin]))
                        # !Need to remove 1 molecule from system if immax is less than nmol
                        im = immax
                        while im < cc.nmol:
                            im1 = im + 1
                            nato[im] = nato[im1]
                            naf[im] = naf[im1]
                            nal[im] = nal[im1]
                            for ia in range(naf[im], nal[im] + 1):
                                mref[ia] = im
                            ia = _do_end(int(naf[im]), int(nal[im]))
                            im = im + 1
                            # !No need to worry about molecular coordinates, etc., at this stage as none have been defined.
                        cc.nmol = cc.nmol - 1
                        for im in range(1, immax + 1):
                            units.write(6, lw(im, naf[im], nal[im], mref[naf[im]], mref[nal[im]]))
                        im = _do_end(1, immax)
                    i = i + 1
    units.close(10)
    if mol.natom > 0:
        ierr = 0  # Reset the error indicator so that calling program does not think read failed.
    for im in range(1, cc.nmol + 1):
        xsum = f32(0.0)
        ysum = f32(0.0)
        zsum = f32(0.0)
        for ib in range(naf[im], nal[im] + 1):
            xsum = xsum + xa[ib]
            ysum = ysum + ya[ib]
            zsum = zsum + za[ib]
        cc.xadd[im] = xsum / f32(nato[im])
        cc.yadd[im] = ysum / f32(nato[im])
        cc.zadd[im] = zsum / f32(nato[im])
        for ib in range(naf[im], nal[im] + 1):
            xa[ib] = xa[ib] - cc.xadd[im]
            ya[ib] = ya[ib] - cc.yadd[im]
            za[ib] = za[ib] - cc.zadd[im]
    # !Set up the atomic numbers of all the atoms
    iatomicno = int_array(mol.natom)
    mol.natomr = 0
    for ia in range(1, mol.natom + 1):
        ir = 0
        found = False
        while not found and ir < melements:
            ir = ir + 1
            # !It is assumed the atom labels are also the chemical symbols
            found = chsymboltemp[ia] == upcase(chsymbol[ir], len_string=2)  # For PDB files the first 1 or 2 characters of the atom name are the atomic symbol in upper case
            if found:
                iatomicno[ia] = ir
        if not found:
            # ! If no valid atomic symbol found. This implies this atom is a dummy, non-scattering atom, so we convert it a Q atom of type oxygen
            iatomicno[ia] = 8
            mol.alabel[ia] = setsub(mol.alabel[ia], 1, 1, "Q", 3)
        # set up 'type' of atom - this is used in assigning neutron weights
        ato_manipulation_rev4.ftype(ia)
        if mol.qtype[ia] != 1:
            mol.natomr = mol.natomr + 1
    rho = f32(mol.natomr) / cc.volume
    units.write(6, lw(" readpdb> rho = ", rho, mol.natomr, mol.natom, cc.nmol))
    for it in range(1, mol.ntype + 1):
        units.write(6, lw(it, Ch(mol.type[it], 3)))
    chsymboltemp = None

    #
    # make a list of the atom types read in from the current file
    #
    mol.ntypefile = 0
    ia = 0
    while ia < mol.natom:
        ia = ia + 1
        itypefile = 0
        ifound = 0
        while itypefile < mol.ntypefile and ifound == 0:
            itypefile = itypefile + 1
            if mol.alabel[ia] == mol.typefile[itypefile]:
                ifound = itypefile
        if ifound == 0:
            mol.ntypefile = mol.ntypefile + 1
            mol.typefile[mol.ntypefile] = mol.alabel[ia]
            mol.iatomicnotype[mol.ntypefile] = iatomicno[ia]
            units.write(6, lw(mol.ntypefile, Ch(mol.typefile[mol.ntypefile], 3),
                              mol.iatomicnotype[mol.ntypefile]))
    # !Use these values to set the potential parameters
    for it in range(1, mol.ntypefile + 1):
        iat = mol.iatomicnotype[it]
        if chsymbol[iat] == "H":
            epsilon = f32(0)
            sigma = f32(0)
            amass = f32(2)
        elif chsymbol[iat] == "O":
            epsilon = f32(0.65)
            sigma = f32(3.1)
            amass = chatwt[iat]
        elif chsymbol[iat] == "N":
            epsilon = f32(0.7)
            sigma = f32(3.2)
            amass = chatwt[iat]
        elif chsymbol[iat] == "C":
            epsilon = f32(0.8)
            sigma = f32(3.7)
            amass = chatwt[iat]
        else:
            epsilon = f32(0.8)
            sigma = f32(3.2)
            amass = chatwt[iat]
        # patched reference: 11 positional arguments, radmaxsph and pseudofacsph absent
        ato_manipulation_rev4.ctype(mol.typefile[it], chsymbol[iat], 0, epsilon, sigma, amass, f32(1.0), f32(0.0),
                                    natomsph=f32(0.5) * sigma, rhosph=f32(1.0), radminsph=f32(0.0))
    #
    # make a list of the molecule types read in from the current file
    #
    im = 0
    cc.nmoltype = 0
    while im < cc.nmol:
        im = im + 1
        imoltype = 0
        ifound = 0
        while imoltype < cc.nmoltype and ifound == 0:
            imoltype = imoltype + 1
            if mol.alabel[naf[im]] == cc.moltype[imoltype]:
                ifound = imoltype
                cc.mtype[im] = ifound
        if ifound == 0:
            cc.nmoltype = cc.nmoltype + 1
            nmt = cc.nmoltype
            cc.moltype[nmt] = mol.alabel[naf[im]]
            cc.moltethered[nmt] = cc.tethered[im]
            #
            # initialise the template filename
            #
            cc.moltemplate[nmt] = fstr(fmt('("moltype",i2.2)', nmt), 256)
            cc.mtype[im] = nmt

            # Initialise the rotational and translational step sizes

            cc.interntrans[nmt] = mol.stepmi
            cc.internrot[nmt] = mol.stepri
            cc.wholerot[nmt] = mol.stepr
            cc.wholetrans[nmt] = mol.stepm

    # Write these out as a check

    units.write(6, fmt(_F201))

    for imoltype in range(1, cc.nmoltype + 1):

        units.write(6, fmt(_F202, imoltype, Ch(cc.moltype[imoltype], 3), substr(cc.moltemplate[imoltype], 1, 10),
                           cc.wholerot[imoltype], cc.wholetrans[imoltype]))

    # ! Now setup the bonding parameters within each molecule
    # ! Get the list of van der Waals radii
    try:
        units.open(10, "vanderWaalsRadii.txt", status="old")
        ierr = 0
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if ierr != 0:
        units.write(6, lw(" readpdb> Unable to open vanderWaalsRadii.txt"))
        return rho, ierr
    try:
        commands.sometext = fstr(units.read_line(10), commands.nstsize)
        ierr = 0
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    icount = 0
    while ierr == 0:
        commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 11, commands.nstsize,
                                      commands.nworddim)
        if commands.nwords > 2:  # Ignore lines with insufficient data
            # !The first word is the atomic number.
            # !The second is the atomic symbol - ignored in this case.
            # !The third is the van der Waals radius.
            ir = read_internal(substr(commands.sometext, ncf[1], ncl[1])).int(0)
            if ir > 0 and ir <= melements:
                vanderWaalsRadius[ir] = read_internal(substr(commands.sometext, ncf[3], ncl[3])).real(
                    vanderWaalsRadius[ir])
                if ir == 1:
                    vanderWaalsRadius[ir] = 0.0  # The van der Waals radius for H is set to zero.
                icount = icount + 1
        try:
            commands.sometext = fstr(units.read_line(10), commands.nstsize)
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = iostat_of(e)
    units.close(10)
    if icount > 0:
        ierr = 0
    units.write(6, fmt(_F48, icount), advance=False)
    try:
        vanderWaalsFactor = units.list_read(5).real(vanderWaalsFactor)
        ierr = 0
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if ierr != 0 or vanderWaalsFactor <= 0.0:
        vanderWaalsFactor = f32(0.7)
    units.write(6, fmt(_F49), advance=False)
    yesorno = fstr(units.read_line(5), 3)
    dodihedrals = findex(yesorno, "y") > 0
    ato_read_and_save.normalise_coordinates()
    iabond = mol.iabond
    tbond = mol.tbond
    nbondp = mol.nbondp
    itype = mol.itype
    iatomicnotype = mol.iatomicnotype
    dmax = cc.dmax
    for im in range(1, cc.nmol + 1):
        for ia in range(naf[im], nal[im] - 1 + 1):
            for ja in range(ia + 1, nal[im] + 1):
                xdif = xa[ia] - xa[ja]
                ydif = ya[ia] - ya[ja]
                zdif = za[ia] - za[ja]
                dtest = dmax * sqrtf(xdif * xdif + ydif * ydif + zdif * zdif)
                dvanderwaals = vanderWaalsFactor * (vanderWaalsRadius[iatomicnotype[itype[ia]]]
                                                    + vanderWaalsRadius[iatomicnotype[itype[ja]]])
                found = dtest < dvanderwaals
                # !               write(6,*) im,ia,ja,dtest,dvanderwaals,found
                if found:
                    if nbond[ia] < near and nbond[ja] < near:
                        nbond[ia] = nbond[ia] + 1
                        iabond[ia, nbond[ia]] = ja
                        tbond[ia, nbond[ia]] = dtest / dmax
                        nbond[ja] = nbond[ja] + 1
                        iabond[ja, nbond[ja]] = ia
                        tbond[ja, nbond[ja]] = dtest / dmax
    # !Now check for bond angles within each molecule
    for ia in range(1, mol.natom + 1):
        nbondp[ia] = nbond[ia]
    ir = 0
    for im in range(1, cc.nmol + 1):
        for ia in range(naf[im], nal[im] + 1):
            ib1 = 0
            while ib1 < nbondp[ia] - 1:  # Don't want to include any new bond angles that are already defined
                ib1 = ib1 + 1
                iaref1 = int(iabond[ia, ib1])
                ib2 = ib1
                while ib2 < nbondp[ia]:
                    ib2 = ib2 + 1
                    iaref2 = int(iabond[ia, ib2])
                    # !Need to check that iabond(ia,ib2) does not already exist in iaref1's neighbour list
                    ic = 0
                    found = False
                    while not found and ic < nbond[iaref1]:
                        ic = ic + 1
                        found = iaref2 == iabond[iaref1, ic]
                    if not found:
                        xdif = xa[iaref1] - xa[iaref2]
                        ydif = ya[iaref1] - ya[iaref2]
                        zdif = za[iaref1] - za[iaref2]
                        dtest = sqrtf(xdif * xdif + ydif * ydif + zdif * zdif)
                        if nbond[iaref1] < near and nbond[iaref2] < near:
                            nbond[iaref1] = nbond[iaref1] + 1
                            nbond[iaref2] = nbond[iaref2] + 1
                            iabond[iaref1, nbond[iaref1]] = iaref2
                            iabond[iaref2, nbond[iaref2]] = iaref1
                            tbond[iaref1, nbond[iaref1]] = -dtest
                            tbond[iaref2, nbond[iaref2]] = -dtest
        # !Dihedral angles do...
        if dodihedrals and nato[im] > 3:  # Can only do this if molecule contains at least 4 atoms
            # !Also dihedrals can only be created between atoms that do not have bonds
            for iaref1 in range(naf[im], nal[im] + 1):  # ia1ref will be the first axis atom
                ib1 = 0
                while ib1 < nbondp[iaref1]:  # Look at the first neighbours of iaref1 to determine the second atom to form the rotation axis (iaref2)
                    ib1 = ib1 + 1
                    iaref2 = int(iabond[iaref1, ib1])
                    # ! Need to check whether a dihedral group for this pair of axis atoms has already been defined.
                    ir1 = iaref1
                    ir2 = iaref2
                    iref = 0
                    ifounddih = 0
                    while iref < ir and ifounddih == 0:
                        iref = iref + 1
                        if (mm.rotgrp[iref] == "DIH"
                                and ((ir1 == mm.labaxi[iref, 1] and ir2 == mm.labaxi[iref, 2])
                                     or (ir2 == mm.labaxi[iref, 1] and ir1 == mm.labaxi[iref, 2]))):
                            ifounddih = iref
                    ib2 = 0
                    while ib2 < nbondp[iaref1]:  # Look at other first neighbours of iaref1 to determine the second axis atom (iaref2)
                        ib2 = ib2 + 1
                        ibref1 = int(iabond[iaref1, ib2])  # (ibref1 must be a direct bond to iaref1, not via a bond angle, hence us nbondp and not nbond)
                        ib3 = 0
                        while ib3 < nbondp[iaref2] and ibref1 != iaref2:  # Look at the first neighbours of iaref2 to determine second atom in the dihedral (ibref2)
                            ib3 = ib3 + 1
                            ibref2 = int(iabond[iaref2, ib3])
                            ic = 0
                            found = ibref1 == ibref2 or ibref2 == iaref1
                            while ic < nbond[ibref2] and not found:  # Need to check that ibref2 has no bonds with ibref1
                                ic = ic + 1
                                found = ibref1 == iabond[ibref2, ic]
                            if ifounddih > 0 and not found:
                                # !Need to check whether a dihedral group with these second atoms already exists
                                iref = 0
                                while iref < mol.ndihedralang[ifounddih] and not found:
                                    iref = iref + 1
                                    found = ((ibref1 == mol.labrotdih[1, iref, ifounddih]
                                              and ibref2 == mol.labrotdih[2, iref, ifounddih])
                                             or (ibref2 == mol.labrotdih[1, iref, ifounddih]
                                                 and ibref1 == mol.labrotdih[2, iref, ifounddih]))
                            if not found:
                                dread = ato_manipulation_rev4.dihedral(im, iaref1, iaref2, ibref1, ibref2)
                                units.write(6, lw(ifounddih, iaref1, iaref2, ibref1, ibref2, dread))
                                dread = dread / piconv  # Need to convert angles to radians as these are normalised coordinates
                                if ifounddih == 0:
                                    # ! Start a new rotational group:
                                    # ! Get the two possible rotational groups into temporary arrays
                                    labrott1, nrott1 = ato_manipulation_rev4.get_rotgrp(ir1, ir2, labrott1, nrott1)
                                    labrott2, nrott2 = ato_manipulation_rev4.get_rotgrp(ir2, ir1, labrott2, nrott2)
                                    units.write(6, fmt(_FROT, nrott1, *[labrott1[j] for j in range(1, nrott1 + 1)]))
                                    units.write(6, fmt(_FROT, nrott2, *[labrott2[j] for j in range(1, nrott2 + 1)]))
                                    ir = ir + 1
                                    mol.ndihedralang[ir] = 1
                                    mol.dihedralang[mol.ndihedralang[ir], ir] = dread
                                    molecom_mod_lib.testsize_rotgrps(ir)
                                    mm.rotgrp[ir] = "DIH"
                                    ifounddih = ir  # Need to record the fact that this is the current found dihedral group
                                    # !Use the smallest of the two groups
                                    if nrott1 < nrott2:
                                        mm.labaxi[ir, 1] = ir1
                                        mm.labaxi[ir, 2] = ir2
                                        mm.nrot[ir] = nrott1
                                        mm.rotf[ir] = mm.nrottot + 1
                                        mm.nrottot = mm.nrottot + int(mm.nrot[ir])
                                        mm.rotl[ir] = mm.nrottot
                                        molecom_mod_lib.testsize_labrot()
                                        for i in range(1, mm.nrot[ir] + 1):
                                            labref = mm.rotf[ir] + i - 1
                                            mm.labrot[labref] = labrott1[i]
                                        # !Remember to preserve the linking to the axis atoms!!
                                        mol.labrotdih[1, mol.ndihedralang[ir], ir] = ibref1
                                        mol.labrotdih[2, mol.ndihedralang[ir], ir] = ibref2
                                    else:
                                        mm.labaxi[ir, 1] = ir2
                                        mm.labaxi[ir, 2] = ir1
                                        mm.nrot[ir] = nrott2
                                        mm.rotf[ir] = mm.nrottot + 1
                                        mm.nrottot = mm.nrottot + int(mm.nrot[ir])
                                        mm.rotl[ir] = mm.nrottot
                                        molecom_mod_lib.testsize_labrot()
                                        for i in range(1, mm.nrot[ir] + 1):
                                            labref = mm.rotf[ir] + i - 1
                                            mm.labrot[labref] = labrott2[i]
                                        # !Remember to preserve the linking to the axis atoms!!
                                        mol.labrotdih[1, mol.ndihedralang[ir], ir] = ibref2
                                        mol.labrotdih[2, mol.ndihedralang[ir], ir] = ibref1
                                else:
                                    # ! Add the new dihedral angle to the existing list
                                    if mol.ndihedralang[ifounddih] < near:
                                        mol.ndihedralang[ifounddih] = mol.ndihedralang[ifounddih] + 1
                                        nd = mol.ndihedralang[ifounddih]
                                        mol.dihedralang[nd, ifounddih] = dread
                                        if ir1 == mm.labaxi[ifounddih, 1] and ir2 == mm.labaxi[ifounddih, 2]:
                                            mol.labrotdih[1, nd, ifounddih] = ibref1
                                            mol.labrotdih[2, nd, ifounddih] = ibref2
                                        elif ir2 == mm.labaxi[ifounddih, 1] and ir1 == mm.labaxi[ifounddih, 2]:
                                            mol.labrotdih[1, nd, ifounddih] = ibref2
                                            mol.labrotdih[2, nd, ifounddih] = ibref1
        cc.ngrot[im] = ir - mm.ngrrot
        mm.ngrrot = ir
    # !      do ia=1,natom
    # !         nbond(ia)=nbondp(ia)
    # !      end do
    ato_read_and_save.unnormalise_coordinates()

    return rho, ierr


@routine
def upcase(string, len_string=None):
    """readpdb.f:upcase -- the string in upper case (``character(len=len(string))``).

    ``len_string`` is the declared length of the actual argument (the result keeps it;
    trailing blanks are dropped as for every Python character value).
    """
    n = len(string) if len_string is None else len_string
    s = string.ljust(n)[:n]
    out = []
    for i in range(1, n + 1):
        j = ord(s[i - 1])
        if j >= ord("a") and j <= ord("z"):
            out.append(chr(ord(s[i - 1]) - 32))
        else:
            out.append(s[i - 1])
    return fstr("".join(out), n)
