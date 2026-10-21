"""Translation of fortran_src/get_molecules.f."""

from ..common import chemical_symbols, clustercom, commands, epsr_dimension, molecom
from ..fortran import (
    Ch,
    FortranEOF,
    FortranIOError,
    f32,
    fmt,
    fstr,
    int_array,
    logical_array,
    lw,
    read_internal,
    real_array,
    sqrtf,
    substr,
    units,
)
from . import minimagec, parse, reallocation_routines
from ..trace import routine

_F48 = ("(/1x,'get_molecules> Number of van der Waals radii read in: ',i3"
        "//' Give van der Waals factor (suggested value 0.8, but check): ')")
_F24 = "(/1x,'get_molecules> There are ',i7,' molecules in this file'/)"
_F27 = "(1x,'Molecule type ',i2,' has label ',a3)"
_F32 = "(/1x,'get_molecules> There are ',i7,' atoms and ',i2,' types of atom in this file'/)"
_F33 = "(1x,'Atom type ',i2,' has label ',a3)"
_F381 = "(/1x,'get_molecules> Give pairs of atom types needed to define bonds (e.g. 3 4 4 2), 0 to exit): ')"
_F57 = "(' get_molecules> For atom type ',i2,1x,i5,' atoms found')"
_F168 = "(/1x,'get_molecules> Initial molecule no,, number of atoms, atom list'/)"
_F237 = "(/1x,'get_molecules> Final molecule no,, number of atoms, atom list'/)"
_F389 = "(/1x,'get_molecules> There are now ',i7,' molecules and ',i2,' types of molecule in this file'/)"
_F392 = "(1x,'Molecule type ',i2,' has label ',a3,1x,a,4(1x,f10.5))"
_F413 = "(/1x,'get_molecules> *** Use the label command to give atoms distinct names if needed ***')"


@routine
def get_molecules():
    """get_molecules.f:get_molecules -- rebuild molecules from van der Waals contacts.

    resets the value of a specified parameter in the .ato file if requested
    by calling program.  Atoms of the requested type pairs closer than the scaled sum of
    their van der Waals radii are bonded and their molecules merged; the .ato arrays are
    then rewritten molecule by molecule and the molecule types regenerated.
    """
    # Imported here: these files belong to other packages.
    from . import ato_manipulation_rev4, ato_read_and_save

    # use reallocation_routines
    near = epsr_dimension.near
    melements = chemical_symbols.melements
    nstsize = commands.nstsize
    nworddim = commands.nworddim
    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za
    naf = clustercom.naf
    nal = clustercom.nal
    nato = clustercom.nato
    xadd = clustercom.xadd
    yadd = clustercom.yadd
    zadd = clustercom.zadd
    mref = molecom.mref
    itype = molecom.itype
    nbond = molecom.nbond
    iabond = molecom.iabond
    tbond = molecom.tbond

    atomtypelist = None  # integer, dimension(:), allocatable
    natomslist = None
    atomslist = None  # integer, dimension(:,:), allocatable
    vanderWaalsFactor = f32(0.0)
    # NB: vanderWaalsRadius is an uninitialised local in the Fortran; elements not read
    # from vanderWaalsRadii.txt are taken as 0.
    vanderWaalsRadius = real_array(melements)

    # Arrays for new molecule definitions (allocatable)
    natonew = nafnew = nalnew = nbondnew = mrefnew = imageatomnew = itypenew = None
    iabondnew = None
    tbondnew = None
    xanew = yanew = zanew = None
    alabelnew = None
    tetherednew = None
    tetheredatomnew = None
    # Arrays to store old molecule types
    itold = irold = wrold = wtold = None  # Step sizes for current molecules
    moltypeold = None
    moltemplateold = None
    # NB: ir is read from the file before it is set; 0 is its value if the read is null.
    ir = 0

    ato_read_and_save.normalise_coordinates()

    # Get the list of van der Waals radii
    try:
        units.open(10, "vanderWaalsRadii.txt", status="old")
        ierr = 0
    except FortranIOError as e:
        ierr = e.iostat
    if ierr != 0:
        units.write(6, lw("get_molecules> Unable to open vanderWaalsRadii.txt"))
        return
    try:
        commands.sometext = fstr(units.read_line(10), nstsize)
        ierr = 0
    except FortranEOF:
        ierr = -1
    except FortranIOError as e:
        ierr = e.iostat
    icount = 0
    while ierr == 0:
        commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 11,
                                      nstsize, nworddim)
        if commands.nwords > 2:  # Ignore lines with insufficient data
            # The first word is the atomic number.
            # The second is the atomic symbol - ignored in this case.
            # The third is the van der Waals radius.
            ir = read_internal(substr(commands.sometext, commands.ncf[1], commands.ncl[1])).int(ir)
            if ir > 0 and ir <= melements:
                vanderWaalsRadius[ir] = read_internal(
                    substr(commands.sometext, commands.ncf[3], commands.ncl[3])).real(vanderWaalsRadius[ir])
                icount = icount + 1
        try:
            commands.sometext = fstr(units.read_line(10), nstsize)
            ierr = 0
        except FortranEOF:
            ierr = -1
        except FortranIOError as e:
            ierr = e.iostat
    units.close(10)
    units.write(6, fmt(_F48, icount), advance=False)
    r = units.list_read(5)
    vanderWaalsFactor = r.real(vanderWaalsFactor)
    ierr = 0
    nmoltype = clustercom.nmoltype
    # Allocate space to save molecule template information
    itold = reallocation_routines.reallocate1d_r(itold, nmoltype)
    irold = reallocation_routines.reallocate1d_r(irold, nmoltype)
    wrold = reallocation_routines.reallocate1d_r(wrold, nmoltype)
    wtold = reallocation_routines.reallocate1d_r(wtold, nmoltype)
    moltypeold = reallocation_routines.reallocate1d_c(moltypeold, 3, nmoltype)
    moltemplateold = reallocation_routines.reallocate1d_c(moltemplateold, 256, nmoltype)
    # And save it
    nmoltypeold = nmoltype
    for it in range(1, nmoltype + 1):
        itold[it] = clustercom.interntrans[it]
        irold[it] = clustercom.internrot[it]
        wrold[it] = clustercom.wholerot[it]
        wtold[it] = clustercom.wholetrans[it]
        moltypeold[it] = clustercom.moltype[it]
        moltemplateold[it] = clustercom.moltemplate[it]
    # For this exercise we need to set the coordinates of each atom to their absolute values
    # (because we are going to alter the molecule reference numbers)
    for im in range(1, clustercom.nmol + 1):
        for ia in range(naf[im], nal[im] + 1):
            xa[ia] = xa[ia] + xadd[im]
            ya[ia] = ya[ia] + yadd[im]
            za[ia] = za[ia] + zadd[im]
        xadd[im] = 0.0
        yadd[im] = 0.0
        zadd[im] = 0.0
    units.write(6, fmt(_F24, clustercom.nmol))
    for it in range(1, nmoltypeold + 1):
        units.write(6, fmt(_F27, it, Ch(moltypeold[it], 3)))
    units.write(6, fmt(_F32, molecom.natom, molecom.ntype))
    for it in range(1, molecom.ntype + 1):
        units.write(6, fmt(_F33, it, Ch(molecom.type[it], 3)))
    units.write(6, fmt(_F381), advance=False)
    commands.sometext = fstr(units.read_line(5), nstsize)
    units.write(6, lw())
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  nstsize, nworddim)
    if atomtypelist is not None:
        atomtypelist = None
    if commands.nwords < 2:  # Need a minimum of two numbers to generate bonds
        return
    else:
        # Only use the numbers in pairs so ignore any extra words at end
        commands.nwords = 2 * (commands.nwords // 2)
        lenatomtypelist = commands.nwords
        atomtypelist = int_array(commands.nwords)
        rd = read_internal(substr(commands.sometext, commands.ncf[1], commands.ncl[commands.nwords]))
        for it in range(1, commands.nwords + 1):
            atomtypelist[it] = rd.int(atomtypelist[it])
    natomslist = reallocation_routines.reallocate1d_i(natomslist, lenatomtypelist)
    for it in range(1, lenatomtypelist + 1):
        # Make a list of atoms of each specified type. First count the number, then save the atom numbers
        natomslist[it] = 0
    for ia in range(1, molecom.natom + 1):
        for it in range(1, lenatomtypelist + 1):
            if itype[ia] == atomtypelist[it]:
                natomslist[it] = natomslist[it] + 1
    # Get the largest of these lists
    maxatoms = 0
    for it in range(1, lenatomtypelist + 1):
        maxatoms = max(maxatoms, int(natomslist[it]))
        natomslist[it] = 0
    atomslist = reallocation_routines.reallocate2d_i(atomslist, maxatoms, lenatomtypelist)
    # Get the values
    for ia in range(1, molecom.natom + 1):
        for it in range(1, lenatomtypelist + 1):
            if itype[ia] == atomtypelist[it]:
                natomslist[it] = natomslist[it] + 1
                atomslist[natomslist[it], it] = ia
    for it in range(1, lenatomtypelist + 1):
        units.write(6, fmt(_F57, atomtypelist[it], natomslist[it]))
        units.write(6, lw(*[atomslist[j, it] for j in range(1, natomslist[it] + 1)]))
    # Create an array to store the sequence of atom numbers for each new molecule
    # Make sure this is large enough to cope with the existing molecules
    maxatoms = 0
    for im in range(1, clustercom.nmol + 1):
        maxatoms = max(maxatoms, int(nato[im]))
    newmolatomlist = None
    newmolatomlist = reallocation_routines.reallocate2d_i(newmolatomlist, maxatoms, clustercom.nmol)
    # Populate this list
    units.write(6, fmt(_F168))
    for im in range(1, clustercom.nmol + 1):
        for ia in range(1, nato[im] + 1):
            iaref = ia + int(naf[im]) - 1
            newmolatomlist[ia, im] = iaref
    #         write(6,*) 'get_molecules> ',im,nato(im),(newmolatomlist(ia,im),ia=1,nato(im))
    # Look for bonded atom pairs between the specified atom types, using atom type pairs in list
    it = 1
    while it < lenatomtypelist:
        jt = it + 1
        for il in range(1, natomslist[it] + 1):
            ia = int(atomslist[il, it])
            xi = xa[ia]  # Coordinates of this atom
            yi = ya[ia]
            zi = za[ia]
            for jl in range(1, natomslist[jt] + 1):
                ja = int(atomslist[jl, jt])
                xj = xa[ja]  # Coordinates of this atom
                yj = ya[ja]
                zj = za[ja]
                if ia != ja and nato[mref[ja]] > 0:  # Don't include atoms from molecules which no longer have any atoms
                    xdif = (xj - xi)
                    ydif = (yj - yi)
                    zdif = (zj - zi)
                    xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                    dtest = clustercom.dmax * sqrtf(xdif * xdif + ydif * ydif + zdif * zdif)
                    im = int(mref[ia])
                    jm = int(mref[ja])
                    dvanderwaals = vanderWaalsFactor * (
                        vanderWaalsRadius[molecom.iatomicnotype[itype[ia]]]
                        + vanderWaalsRadius[molecom.iatomicnotype[itype[ja]]])
                    found = dtest < dvanderwaals  # Need to check these atoms EVEN if they are already part of the same molecule
                    #                   write(6,*) ia,ja,dtest,dvanderwaals
                    if found:
                        units.write(6, lw(ia, ja, im, jm, nato[im], nato[jm], dtest, dvanderwaals))
                        if nbond[ia] < near and nbond[ja] < near:
                            # need to check that ja is not already a neighbour of ia
                            found = False
                            ib = 0
                            while ib < nbond[ia] and not found:
                                ib = ib + 1
                                found = iabond[ia, ib] == ja
                            if not found:
                                nbond[ia] = nbond[ia] + 1
                                iabond[ia, nbond[ia]] = ja
                                tbond[ia, nbond[ia]] = dtest / clustercom.dmax
                                nbond[ja] = nbond[ja] + 1
                                iabond[ja, nbond[ja]] = ia
                                tbond[ja, nbond[ja]] = dtest / clustercom.dmax
                                # Only do the next bit if im.ne.jm
                                if im != jm:
                                    maxatoms = max(maxatoms, int(nato[im] + nato[jm]))
                                    newmolatomlist = reallocation_routines.reallocate2d_i(
                                        newmolatomlist, maxatoms, clustercom.nmol)
                                    # Need to give all the atoms contained in molecule mref(ja) the new molecule number
                                    ireflist = int(nato[im])
                                    for jb in range(1, nato[jm] + 1):  # No. of atoms in each molecule is kept up to date
                                        jat = int(newmolatomlist[jb, jm])  # Use newmolatomlist to get the atoms in this molecule because naf and nal are not reliable here
                                        ireflist = ireflist + 1
                                        mref[jat] = im
                                        newmolatomlist[ireflist, im] = jat
                                    nato[im] = nato[im] + nato[jm]
                                    nato[jm] = 0
                    #                      write(6,*) ia,ja,im,jm,nato(im),nato(jm)
        it = jt + 1  # Get next pair
    # Determine the maximum number of atoms in any molecule, and the number of molecules with non-zero
    # numbers of atoms
    maxatoms = 0
    nmolnew = 0
    units.write(6, fmt(_F237))
    for im in range(1, clustercom.nmol + 1):
        maxatoms = max(maxatoms, int(nato[im]))
        #         write(6,*) ' get_molecules> ',im,nato(im),(newmolatomlist(ia,im),ia=1,nato(im))
        if nato[im] > 0:
            nmolnew = nmolnew + 1
            units.write(6, lw("get_molecules> ", im, nato[im],
                              *[newmolatomlist[ia, im] for ia in range(1, nato[im] + 1)]))
    # Create arrays needed for temporary storage during resorting of atoms
    natom = molecom.natom
    natonew = reallocation_routines.reallocate1d_i(natonew, nmolnew)
    nafnew = reallocation_routines.reallocate1d_i(nafnew, nmolnew)
    nalnew = reallocation_routines.reallocate1d_i(nalnew, nmolnew)
    nbondnew = reallocation_routines.reallocate1d_i(nbondnew, natom)
    itypenew = reallocation_routines.reallocate1d_i(itypenew, natom)
    iabondnew = reallocation_routines.reallocate2d_i(iabondnew, natom, near)
    tbondnew = reallocation_routines.reallocate2d_r(tbondnew, natom, near)
    mrefnew = reallocation_routines.reallocate1d_i(mrefnew, natom)
    xanew = reallocation_routines.reallocate1d_r(xanew, natom)
    yanew = reallocation_routines.reallocate1d_r(yanew, natom)
    zanew = reallocation_routines.reallocate1d_r(zanew, natom)
    alabelnew = reallocation_routines.reallocate1d_c(alabelnew, 3, natom)
    imageatomnew = reallocation_routines.reallocate1d_i(imageatomnew, natom)
    if tetherednew is not None:
        tetherednew = None
        tetheredatomnew = None
    tetherednew = logical_array(nmolnew)
    tetheredatomnew = int_array(nmolnew)
    # Create the new molecules
    imnew = 0
    ianew = 0
    for im in range(1, clustercom.nmol + 1):
        if nato[im] > 0:
            imnew = imnew + 1
            natonew[imnew] = nato[im]
            tetherednew[imnew] = clustercom.tethered[im]
            tetheredatomnew[imnew] = clustercom.tetheredatom[im]
            nafnew[imnew] = ianew + 1
            for iaref in range(1, nato[im] + 1):
                ia = int(newmolatomlist[iaref, im])
                ianew = ianew + 1
                alabelnew[ianew] = molecom.alabel[ia]
                xanew[ianew] = xa[ia]
                yanew[ianew] = ya[ia]
                zanew[ianew] = za[ia]
                imageatomnew[ianew] = molecom.imageatom[ia]
                nbondnew[ianew] = nbond[ia]
                itypenew[ianew] = itype[ia]
                for ib in range(1, nbondnew[ianew] + 1):
                    # Search for the correct bond reference atom
                    found = False
                    is_ = 0  # Fortran variable 'is'
                    while is_ < nato[im] and not found:
                        is_ = is_ + 1
                        found = iabond[ia, ib] == newmolatomlist[is_, im]
                    iabondnew[ianew, ib] = is_ + int(nafnew[imnew]) - 1
                    tbondnew[ianew, ib] = tbond[ia, ib]
                mrefnew[ianew] = imnew
            nalnew[imnew] = ianew
    nmolnew = imnew
    # Finally overwrite the existing .ato file arrays
    clustercom.nmol = nmolnew
    for im in range(1, clustercom.nmol + 1):
        nato[im] = natonew[im]
        naf[im] = nafnew[im]
        nal[im] = nalnew[im]
        clustercom.tethered[im] = tetherednew[im]
        clustercom.tetheredatom[im] = tetheredatomnew[im]
        xadd[im] = 0.0
        yadd[im] = 0.0
        zadd[im] = 0.0
        clustercom.phix[im] = 0.0
        clustercom.phiy[im] = 0.0
        clustercom.phiz[im] = 0.0
        for ia in range(naf[im], nal[im] + 1):
            molecom.alabel[ia] = alabelnew[ia]
            xa[ia] = xanew[ia]
            ya[ia] = yanew[ia]
            za[ia] = zanew[ia]
            molecom.imageatom[ia] = imageatomnew[ia]  # Need to set this value once new box is complete
            nbond[ia] = nbondnew[ia]
            itype[ia] = itypenew[ia]
            for ib in range(1, nbond[ia] + 1):
                iabond[ia, ib] = iabondnew[ia, ib]
                tbond[ia, ib] = tbondnew[ia, ib]
            mref[ia] = mrefnew[ia]
        # For each atom in a molecule the coordinates should be relative to the centre of mass - not the periodic image
        for ia in range(naf[im], nal[im] + 1):
            x1 = xa[ia] - xa[naf[im]]
            y1 = ya[ia] - ya[naf[im]]
            z1 = za[ia] - za[naf[im]]
            x2 = x1
            y2 = y1
            z2 = z1
            x2, y2, z2 = minimagec.minimagec(x2, y2, z2)
            # write(6,*) ia,dmax*sqrt(x1*x1+y1*y1+z1*z1),dmax*sqrt(x2*x2+y2*y2+z2*z2)
            if abs(x2) < abs(x1):
                xa[ia] = x2 + xa[naf[im]]
            if abs(y2) < abs(y1):
                ya[ia] = y2 + ya[naf[im]]
            if abs(z2) < abs(z1):
                za[ia] = z2 + za[naf[im]]
        ato_manipulation_rev4.centreofmass(im)
        xadd[im], yadd[im], zadd[im] = minimagec.minimagec(xadd[im], yadd[im], zadd[im])  # Place molecule at its minimum image position
        if clustercom.tethered[im]:  # Set the tethering points if required.
            clustercom.phix[im] = xadd[im]
            clustercom.phiy[im] = yadd[im]
            clustercom.phiz[im] = zadd[im]
    ato_read_and_save.unnormalise_coordinates()
    #
    # get_molecules will normally be called either after makelattice - generate a crystalline lattice from a set of crystal
    # coordinates - or after readxyz - generate a .ato file from a .xyz file. In the former case the atoms may have sensible
    # names, but in the latter case they certainly won't since the xyz file consists of the chemical symbol of each atom
    # and its (absolute) xyz coordinates. In that case the chemical symbol becomes the atom label.
    # Hence, in either case, for molecular mixtures, there may be at this point more than one molecule type that begins with
    # the same atom label, causing confusion to EPSR, which requires the first atom of each molecule type to have
    # a distinct label. Equally there may be several atoms with the same label, e.g. C, but which have different bonding
    # patterns.
    # Make a new list of moltypes. We need to check the old list to be sure that existing moltypes
    # are not overwritten, unless they do not exist.
    #
    # make a list of the molecule types read in from the current file
    #
    clustercom.nmoltype = 0
    im = 0
    while im < clustercom.nmol:
        im = im + 1
        # First check whether this molecule type already occurs in the new list
        imoltype = 0
        ifoundnew = 0
        while imoltype < clustercom.nmoltype and ifoundnew == 0:
            imoltype = imoltype + 1
            if molecom.alabel[naf[im]] == clustercom.moltype[imoltype]:
                ifoundnew = imoltype
        if ifoundnew == 0:  # Does not already occur, so now check the old list
            imoltype = 0
            ifoundold = 0
            while imoltype < nmoltypeold and ifoundold == 0:
                imoltype = imoltype + 1
                if molecom.alabel[naf[im]] == moltypeold[imoltype]:
                    ifoundold = imoltype
            if ifoundold == 0:  # Not found in old list either so create a new mol type
                clustercom.nmoltype = clustercom.nmoltype + 1
                nmt = clustercom.nmoltype
                clustercom.moltype[nmt] = molecom.alabel[naf[im]]
                # initialise the template filename
                clustercom.moltemplate[nmt] = fstr(fmt('("moltype",i2.2)', nmt), 256)
                # Initialise the rotational and translational step sizes
                clustercom.interntrans[nmt] = 0.0
                clustercom.internrot[nmt] = 0.0
                clustercom.wholerot[nmt] = molecom.stepr
                clustercom.wholetrans[nmt] = molecom.stepm
            else:  # Found in old list
                # getting values from old list
                clustercom.nmoltype = clustercom.nmoltype + 1
                nmt = clustercom.nmoltype
                clustercom.moltype[nmt] = molecom.alabel[naf[im]]
                # Copy the template filename
                clustercom.moltemplate[nmt] = moltemplateold[ifoundold]
                # Copy the rotational and translational step sizes
                clustercom.interntrans[nmt] = itold[ifoundold]
                clustercom.internrot[nmt] = irold[ifoundold]
                clustercom.wholerot[nmt] = wrold[ifoundold]
                clustercom.wholetrans[nmt] = wtold[ifoundold]
    units.write(6, fmt(_F389, clustercom.nmol, clustercom.nmoltype))
    for it in range(1, clustercom.nmoltype + 1):
        units.write(6, fmt(_F392, it, Ch(clustercom.moltype[it], 3), substr(clustercom.moltemplate[it], 1, 10),
                           clustercom.interntrans[it], clustercom.internrot[it], clustercom.wholerot[it],
                           clustercom.wholetrans[it]))
    units.write(6, fmt(_F413))
    return
