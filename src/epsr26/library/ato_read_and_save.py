"""Translation of fortran_src/ato_read_and_save.f: read and write .ato configuration files.

* ``minitc`` reads a .ato file (molecules, atoms, bonds, rotational groups, atom-type
  parameters, the saved random-number-generator state and the molecule templates) and
  sets up the box, densities, pair tables and dipole moments;
* ``normalise_coordinates`` / ``unnormalise_coordinates`` convert lengths to and from
  box-fraction units;
* ``moloutc`` writes the .ato file.

Imports of other library modules (and of the Fortran-module data owned by other work
packages) are made inside the routines, which avoids import cycles
(ato_manipulation_rev4 <-> ato_read_and_save) and lets this module import before those
modules exist.
"""

from ..common import clustercom, commands, lattice, molecom
from ..fortran import (
    PI4ATAN,
    Ch,
    FortranEOF,
    FortranIOError,
    cosd_intrinsic,
    f32,
    findex,
    findex_padded,
    fmt,
    fstr,
    int_array,
    lw,
    r4,
    read_internal,
    sind_intrinsic,
    sqrtf,
    substr,
    units,
)
from ..trace import routine
from . import parse

_F599_minitc = "(/'Problems with file ',a/)"
_F201_minitc = "(/'minitc> Following molecule types defined initially:-'/)"
_F202_minitc = "(1x,i3,1x,a3,1x,a10,2(1x,e12.5))"
_F203_minitc = "(1x,i3,1x,a3,1x,a10,5(1x,e12.5))"
_F917 = "(1x,'no. of molecules to read = ',4(1x,i0),1x,e13.6,1x,i0)"
_F916_minitc = "(1x,'Atomic fraction (actual and coarse grained)',i4,1x,a,' = ',e12.5,1x,e12.5)"
_F997 = "(1x,'Molecule ',i5,' dipole moment',e13.6)"


def _read_a(unit, old, n):
    """``read(unit,'(a)',iostat=ierr) var`` for a character*n var -> (var, ierr).

    Private helper (not a Fortran routine): on an error or end of file the variable
    keeps its old value and ierr is the gfortran iostat.
    """
    try:
        return fstr(units.read_line(unit), n), 0
    except FortranEOF:
        return old, -1
    except FortranIOError as e:
        return old, e.iostat


@routine
def minitc(fname1, rmax, cell, upol, rho, amubar, amubar2, ndist, ierr):
    """ato_read_and_save.f:minitc -- read a .ato file and set up the box.

    this is a modified version of minitx from epsr15x to include a purely cubic box
    Returns (rmax, cell, upol, rho, amubar, amubar2, ndist, ierr).
    """
    from ..common import molecom_mod
    from . import auxiliary_routines, bondnumbers, reallocation_routines
    from . import ato_manipulation_rev4
    from . import molecom_mod as molecom_mod_routines

    # character*3 a1,ats
    # integer*4 iexch,j,k
    # external cosd,sind  -> EPSR's own cosd/sind (auxiliary_routines)
    # real pi,piconv,natomsph,rhosph,radminsph,radmaxsph,volsph,pseudofacsph
    # NB: implicitly typed real locals that are read with iostat= (and so could be left
    # unset by a failed read) start at 0 here; in the Fortran they are uninitialised.
    epsit = f32(0.0)
    sigt = f32(0.0)
    amass = f32(0.0)
    qt = f32(0.0)
    rqt = f32(0.0)
    natomsph = f32(0.0)
    radminsph = f32(0.0)
    radmaxsph = f32(0.0)
    pseudofacsph = f32(0.0)
    imolread = 0

    pi = PI4ATAN
    piconv = pi / r4(180.0)  # noqa: F841

    # Create the lock file and write set to true (= locked)

    #      write(6,'(a)') fname1(1:len_trim(fname1))//'.lck'
    #      open(10,file=fname1(1:len_trim(fname1))//'.lck',status='unknown',iostat=ierr)
    #      write(10,'(a)') 'T'
    #      close (10)

    try:
        units.open(10, fname1, status="old")
        ierr = 0
    except FortranIOError as e:
        ierr = e.iostat
    #
    # if there is a problem with this file, print out an error message and exit
    #
    if ierr != 0:
        ilenfname = findex(fname1, ".")
        ilenfname = ilenfname + findex(substr(fname1, ilenfname, 256), " ") - 1
        units.write(6, fmt(_F599_minitc, substr(fname1, 1, ilenfname)))
        return rmax, cell, upol, rho, amubar, amubar2, ndist, ierr

    # Decide whether this is a cubic or a non-cubic simulation box

    lattice.cubic = True
    commands.sometext, ierr = _read_a(10, commands.sometext, commands.nstsize)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    if commands.nwords != 3:

        lattice.cubic = False
        # read no. of atoms in molecule, temperature in Kelvin, and lattice vectors
        ierr = 0
        try:
            r = read_internal(commands.sometext)
            clustercom.nmol = r.int(clustercom.nmol)
            molecom.simtk = r.real(molecom.simtk)
        except FortranEOF:
            ierr = -1
        except FortranIOError as e:
            ierr = e.iostat
        ierr = 0
        try:
            r = units.list_read(10)
            lattice.alat = r.real(lattice.alat)
            lattice.blat = r.real(lattice.blat)
            lattice.clat = r.real(lattice.clat)
        except FortranEOF:
            ierr = -1
        except FortranIOError as e:
            ierr = e.iostat
        ierr = 0
        try:
            r = units.list_read(10)
            lattice.phib = r.real(lattice.phib)
            lattice.thetac = r.real(lattice.thetac)
            lattice.phic = r.real(lattice.phic)
        except FortranEOF:
            ierr = -1
        except FortranIOError as e:
            ierr = e.iostat

    else:

        lattice.cubic = True
        ierr = 0
        try:
            r = read_internal(commands.sometext)
            clustercom.nmol = r.int(clustercom.nmol)
            clustercom.dmax = r.real(clustercom.dmax)
            molecom.simtk = r.real(molecom.simtk)
        except FortranEOF:
            ierr = -1
        except FortranIOError as e:
            ierr = e.iostat
        lattice.alat = clustercom.dmax
        lattice.blat = clustercom.dmax
        lattice.clat = clustercom.dmax
        lattice.phib = r4(90.0)
        lattice.thetac = r4(0.0)
        lattice.phic = r4(0.0)

    units.write(6, fmt("('minitc> cubic: ',L10)", lattice.cubic))
    # Assign initial scale factors for each molecule if not defined on input
    molecom_mod.scalefacmol = reallocation_routines.reallocate1d_r(molecom_mod.scalefacmol, clustercom.nmol)
    molecom_mod.scalefacmol[:] = 1.0
    #
    # read fitting tolerence, step sizes, and temperatures for harmonic terms
    #
    commands.sometext, ierr = _read_a(10, commands.sometext, commands.nstsize)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    ierr = 0
    try:
        r = read_internal(commands.sometext)
        molecom.tol = r.real(molecom.tol)
        molecom.stepmi = r.real(molecom.stepmi)
        molecom.stepri = r.real(molecom.stepri)
        molecom.stepr = r.real(molecom.stepr)
        molecom.stepm = r.real(molecom.stepm)
        molecom.tempb = r.real(molecom.tempb)
    except FortranEOF:
        ierr = -1
    except FortranIOError as e:
        ierr = e.iostat
    molecom.tempangreq = r4(1.0)
    molecom.tempdihreq = r4(10.0)
    molecom.tempangwt = r4(1.0)
    molecom.tempdihwt = r4(1.0)
    ncf = commands.ncf
    ncl = commands.ncl
    if commands.nwords == 7:
        ierr = 0
        try:
            r = read_internal(substr(commands.sometext, ncf[7], ncl[7]))
            molecom.tempdihreq = r.real(molecom.tempdihreq)
        except FortranEOF:
            ierr = -1
        except FortranIOError as e:
            ierr = e.iostat
    if commands.nwords > 7:
        ierr = 0
        try:
            r = read_internal(substr(commands.sometext, ncf[7], ncl[8]))
            molecom.tempangreq = r.real(molecom.tempangreq)
            molecom.tempdihreq = r.real(molecom.tempdihreq)
        except FortranEOF:
            ierr = -1
        except FortranIOError as e:
            ierr = e.iostat
    if commands.nwords > 9:
        ierr = 0
        try:
            r = read_internal(substr(commands.sometext, ncf[9], ncl[10]))
            molecom.tempangwt = r.real(molecom.tempangwt)
            molecom.tempdihwt = r.real(molecom.tempdihwt)
        except FortranEOF:
            ierr = -1
        except FortranIOError as e:
            ierr = e.iostat
    molecom.tempang = molecom.tempangwt
    molecom.tempdih = molecom.tempdihwt * molecom.tempb
    molecom.tolout = molecom.tol
    ia = 0
    ig = 0
    molecom.natom = 0
    molecom.natomr = 0
    molecom_mod.ngrrot = 0
    molecom_mod.nrottot = 0
    # Reset array size counters
    molecom_mod.ngrrotmax = 0
    molecom_mod.nrottotmax = 0
    for im in range(1, clustercom.nmol + 1):
        commands.sometext, ierr = _read_a(10, commands.sometext, commands.nstsize)
        commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 10,
                                      commands.nstsize, commands.nworddim)
        ierr = 0
        try:
            r = read_internal(commands.sometext)
            clustercom.nato[im] = r.int(clustercom.nato[im])
            clustercom.xadd[im] = r.real(clustercom.xadd[im])
            clustercom.yadd[im] = r.real(clustercom.yadd[im])
            clustercom.zadd[im] = r.real(clustercom.zadd[im])
            clustercom.phix[im] = r.real(clustercom.phix[im])
            clustercom.phiy[im] = r.real(clustercom.phiy[im])
            clustercom.phiz[im] = r.real(clustercom.phiz[im])
        except FortranEOF:
            ierr = -1
        except FortranIOError as e:
            ierr = e.iostat
        clustercom.tetheredatom[im] = 0
        if commands.nwords > 7:
            clustercom.tethered[im] = substr(commands.sometext, ncf[8], ncf[8]) == "T"
            if ncl[8] > ncf[8]:
                ierr = 0
                try:
                    r = read_internal(substr(commands.sometext, ncf[8] + 1, ncl[8]))
                    clustercom.tetheredatom[im] = r.int(clustercom.tetheredatom[im])
                except FortranEOF:
                    ierr = -1
                except FortranIOError as e:
                    ierr = e.iostat
                if ierr != 0:
                    clustercom.tetheredatom[im] = 0
        else:
            clustercom.tethered[im] = False
        if commands.nwords > 8:
            clustercom.nstepreject[im] = read_internal(substr(commands.sometext, ncf[9], ncl[9])).int(
                clustercom.nstepreject[im])
        else:
            clustercom.nstepreject[im] = 0
        if commands.nwords > 9:
            molecom_mod.scalefacmol[im] = read_internal(substr(commands.sometext, ncf[10], ncl[10])).real(
                molecom_mod.scalefacmol[im])
        else:
            molecom_mod.scalefacmol[im] = 1.0
        #         write(6,*) im,nwords,nstepreject(im),scalefacmol(im)
        molecom.natom = molecom.natom + int(clustercom.nato[im])
        clustercom.naf[im] = ia + 1
        # set a flag to determine if this molecule has any Q-type atoms
        iflag = 0  # noqa: F841
        for ib in range(1, int(clustercom.nato[im]) + 1):
            ia = ia + 1
            molecom.mref[ia] = im
            # for each atom read its label
            commands.sometext, ierr = _read_a(10, commands.sometext, commands.nstsize)
            commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 10,
                                          commands.nstsize, commands.nworddim)
            # read(sometext(ncf(1):ncl(1)),'(a3)',iostat=ierr) alabel(ia): the first 3
            # characters of the word (blank padded if shorter)
            molecom.alabel[ia] = fstr(substr(commands.sometext, ncf[1], ncl[1]), 3)
            ierr = 0
            # Read the image atom if present.
            if commands.nwords > 2:
                ierr = 0
                try:
                    r = read_internal(substr(commands.sometext, ncf[3], ncl[3]))
                    molecom.imageatom[ia] = r.int(molecom.imageatom[ia])
                except FortranEOF:
                    ierr = -1
                except FortranIOError as e:
                    ierr = e.iostat
                # As a precaution, don't allow unphysical values
                if molecom.imageatom[ia] < 0:
                    molecom.imageatom[ia] = 0
                if molecom.imageatom[ia] > clustercom.nato[im]:
                    molecom.imageatom[ia] = 0
            else:
                molecom.imageatom[ia] = 0
            # set up 'type' of atom - this is used in assigning neutron weights
            ato_manipulation_rev4.ftype(ia)
            #
            # do not count atoms which begin with 'Q' or 'q' - these are assumed
            # invisible to the scattering radiation and experiment (although they
            # can be used to calculate other structural quantities if necessary
            #
            if molecom.qtype[ia] != 1:
                molecom.natomr = molecom.natomr + 1
            else:
                iflag = 1  # noqa: F841
            # read its starting coordinates and sigmas
            ierr = 0
            try:
                r = units.list_read(10)
                molecom.xa[ia] = r.real(molecom.xa[ia])
                molecom.ya[ia] = r.real(molecom.ya[ia])
                molecom.za[ia] = r.real(molecom.za[ia])
            except FortranEOF:
                ierr = -1
            except FortranIOError as e:
                ierr = e.iostat
            # read the number of atoms to which it must be bound, the atom numbers of
            # those neighbours and the anticipated bond length
            ierr = 0
            try:
                r = units.list_read(10)
                molecom.nbond[ia] = r.int(molecom.nbond[ia])
                for j in range(1, int(molecom.nbond[ia]) + 1):
                    molecom.iabond[ia, j] = r.int(molecom.iabond[ia, j])
                    molecom.tbond[ia, j] = r.real(molecom.tbond[ia, j])
            except FortranEOF:
                ierr = -1
            except FortranIOError as e:
                ierr = e.iostat
        clustercom.nal[im] = ia
        #
        # reset positions of Q-type atoms
        #
        #     if(iflag.eq.1) call resetq(im)
        # read in number of internal rotation groups
        ierr = 0
        try:
            r = units.list_read(10)
            clustercom.ngrot[im] = r.int(clustercom.ngrot[im])
        except FortranEOF:
            ierr = -1
        except FortranIOError as e:
            ierr = e.iostat
        clustercom.nrf[im] = ig + 1
        if clustercom.ngrot[im] > 0:
            molecom_mod.ngrrot = molecom_mod.ngrrot + int(clustercom.ngrot[im])
            molecom_mod_routines.testsize_rotgrps()
            for ih in range(1, int(clustercom.ngrot[im]) + 1):
                ig = ig + 1
                # save molecule in which this rotational group occurs
                molecom_mod.mrot[ig] = im
                #                  write(6,*) 'minitc> ',ig,ngrrot,ngrrotmax,size(rotgrp)
                ierr = 0
                try:
                    molecom_mod.rotgrp[ig] = fstr(units.read_fmt(10, "(1x,a3)")[0], 3)
                except FortranEOF:
                    ierr = -1
                except FortranIOError as e:
                    ierr = e.iostat
                # for rotational groups read atom numbers of two atoms to form axis of rotation
                ierr = 0
                try:
                    r = units.list_read(10)
                    molecom_mod.labaxi[ig, 1] = r.int(molecom_mod.labaxi[ig, 1])
                    molecom_mod.labaxi[ig, 2] = r.int(molecom_mod.labaxi[ig, 2])
                except FortranEOF:
                    ierr = -1
                except FortranIOError as e:
                    ierr = e.iostat
                # For DIH groups read number of dihedral angles defined for a given axis, and
                # two atom numbers to calculate the dihedral angle and the required dihedral angle
                if molecom_mod.rotgrp[ig] == "DIH":
                    ierr = 0
                    try:
                        r = units.list_read(10)
                        molecom.ndihedralang[ig] = r.int(molecom.ndihedralang[ig])
                    except FortranEOF:
                        ierr = -1
                    except FortranIOError as e:
                        ierr = e.iostat
                    for k in range(1, int(molecom.ndihedralang[ig]) + 1):
                        ierr = 0
                        try:
                            r = units.list_read(10)
                            for j in range(1, 3):
                                molecom.labrotdih[j, k, ig] = r.int(molecom.labrotdih[j, k, ig])
                            molecom.dihedralang[k, ig] = r.real(molecom.dihedralang[k, ig])
                        except FortranEOF:
                            ierr = -1
                        except FortranIOError as e:
                            ierr = e.iostat
                else:
                    molecom.ndihedralang[ig] = 0
                # Read number of atoms in group and their numbers
                # Read the first line
                commands.sometext, ierr = _read_a(10, commands.sometext, commands.nstsize)
                commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                              commands.nstsize, commands.nworddim)
                molecom_mod.nrot[ig] = read_internal(substr(commands.sometext, ncf[1], ncl[1])).int(
                    molecom_mod.nrot[ig])
                molecom_mod.rotf[ig] = molecom_mod.nrottot + 1
                molecom_mod.nrottot = molecom_mod.nrottot + int(molecom_mod.nrot[ig])
                molecom_mod.rotl[ig] = molecom_mod.nrottot
                molecom_mod_routines.testsize_labrot()
                # Read the remaining words on this line. The first has already been read so we start from 2
                iw = 1
                labref = int(molecom_mod.rotf[ig])
                while iw < commands.nwords:
                    iw = iw + 1
                    labref = int(molecom_mod.rotf[ig]) + iw - 2
                    molecom_mod.labrot[labref] = read_internal(substr(commands.sometext, ncf[iw], ncl[iw])).int(
                        molecom_mod.labrot[labref])
                # Read the remaining groups directly if needed
                if labref < molecom_mod.rotl[ig]:
                    labref = labref + 1
                    ierr = 0
                    try:
                        r = units.list_read(10)
                        for j in range(labref, int(molecom_mod.rotl[ig]) + 1):
                            molecom_mod.labrot[j] = r.int(molecom_mod.labrot[j])
                    except FortranEOF:
                        ierr = -1
                    except FortranIOError as e:
                        ierr = e.iostat
        clustercom.nrl[im] = ig
    #
    # check bond atom numbers and rotational group atom numbers are correct
    # if not convert to correct values
    #
    bondnumbers.bondnumbers()
    #
    # make a list of the atom types read in from the current file
    #
    molecom.ntypefile = 0
    ia = 0
    while ia < molecom.natom:
        ia = ia + 1
        itypefile = 0
        ifound = 0
        while itypefile < molecom.ntypefile and ifound == 0:
            itypefile = itypefile + 1
            if molecom.alabel[ia] == molecom.typefile[itypefile]:
                ifound = itypefile
        if ifound == 0:
            molecom.ntypefile = molecom.ntypefile + 1
            molecom.typefile[molecom.ntypefile] = molecom.alabel[ia]
    #
    # make a list of the molecule types read in from the current file
    #
    nmoltypefirst = clustercom.nmoltype
    im = 0
    while im < clustercom.nmol:
        im = im + 1
        imoltype = 0
        ifound = 0
        while imoltype < clustercom.nmoltype and ifound == 0:
            imoltype = imoltype + 1
            if molecom.alabel[clustercom.naf[im]] == clustercom.moltype[imoltype]:
                ifound = imoltype
                clustercom.mtype[im] = ifound
        if ifound == 0:
            clustercom.nmoltype = clustercom.nmoltype + 1
            clustercom.moltype[clustercom.nmoltype] = molecom.alabel[clustercom.naf[im]]
            clustercom.moltethered[clustercom.nmoltype] = clustercom.tethered[im]
            clustercom.moltetheredatom[clustercom.nmoltype] = clustercom.tetheredatom[im]
            #
            # initialise the template filename
            #
            clustercom.moltemplate[clustercom.nmoltype] = fstr(fmt('("moltype",i2.2)', clustercom.nmoltype), 256)
            clustercom.mtype[im] = clustercom.nmoltype

            # Initialise the rotational and translational step sizes

            clustercom.wholerot[clustercom.nmoltype] = molecom.stepr
            clustercom.wholetrans[clustercom.nmoltype] = molecom.stepm

    # Initialise the relative broadening of each molecule type if coarse graining
    molecom_mod.rfractmol = reallocation_routines.reallocate1d_r(molecom_mod.rfractmol, clustercom.nmoltype)
    molecom_mod.rfractmol[:] = 0.0

    # Write these out as a check

    units.write(6, fmt(_F201_minitc))

    for imoltype in range(1, clustercom.nmoltype + 1):
        units.write(6, fmt(_F202_minitc, imoltype, Ch(clustercom.moltype[imoltype], 3),
                           substr(clustercom.moltemplate[imoltype], 1, 10), clustercom.wholerot[imoltype],
                           clustercom.wholetrans[imoltype]))

    #      do it=1,ntypefile
    #            write(6,'(a,i5,1x,a)') 'minitc> ',it,typefile(it)
    #      end do

    #
    # read Lennard-Jones, mass and charge parameters for each atom type.
    # rhosphere,radminsphere and radmaxsphere (used in coarse graining) are allocated in ctype.
    #
    for it in range(1, molecom.ntypefile + 1):
        commands.sometext, ierr = _read_a(10, commands.sometext, commands.nstsize)
        #
        # parse this text into words
        #
        commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 8,
                                      commands.nstsize, commands.nworddim)
        a1 = fstr(substr(commands.sometext, ncf[1], ncl[1]), 3)
        if commands.nwords > 1:
            ats = fstr(substr(commands.sometext, ncf[2], ncl[2]), 3)
        else:
            ats = "XXX"
        if commands.nwords > 2:
            iexch = read_internal(substr(commands.sometext, ncf[3], ncl[3])).int(0)
        else:
            iexch = -1
        commands.sometext, ierr = _read_a(10, commands.sometext, commands.nstsize)
        commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 11,
                                      commands.nstsize, commands.nworddim)
        ierr = 0
        try:
            r = read_internal(substr(commands.sometext, ncf[1], ncl[5]))
            epsit = r.real(epsit)
            sigt = r.real(sigt)
            amass = r.real(amass)
            qt = r.real(qt)
            rqt = r.real(rqt)
        except FortranEOF:
            ierr = -1
        except FortranIOError as e:
            ierr = e.iostat
        # Get the coarse graining parameters if present
        if commands.nwords > 5:
            ierr = 0
            try:
                # If this is negative it represents a density, not the number of atoms in the sphere.
                natomsph = read_internal(substr(commands.sometext, ncf[6], ncl[6])).real(natomsph)
            except FortranEOF:
                ierr = -1
            except FortranIOError as e:
                ierr = e.iostat
        else:
            natomsph = r4(1.0)
        if commands.nwords > 7:
            ierr = 0
            try:
                r = read_internal(substr(commands.sometext, ncf[7], ncl[8]))
                radminsph = r.real(radminsph)
                radmaxsph = r.real(radmaxsph)
            except FortranEOF:
                ierr = -1
            except FortranIOError as e:
                ierr = e.iostat
            if radmaxsph <= 0.0:
                radmaxsph = r4(0.5) * sigt  # In case default value is read in, set this to radius of atom
        else:
            radminsph = r4(0.0)
            radmaxsph = r4(0.5) * sigt
        if commands.nwords > 8:
            ierr = 0
            try:
                pseudofacsph = read_internal(substr(commands.sometext, ncf[9], ncl[9])).real(pseudofacsph)
            except FortranEOF:
                ierr = -1
            except FortranIOError as e:
                ierr = e.iostat
            if pseudofacsph < 1.0:
                pseudofacsph = r4(1.0)
        else:
            pseudofacsph = r4(1.0)
        # Get the density of this atom type
        volsph = r4(4.0) * pi * (radmaxsph * radmaxsph * radmaxsph - radminsph * radminsph * radminsph) / r4(3.0)
        if abs(natomsph) > 0.0:
            if natomsph < 0.0:
                rhosph = abs(natomsph)
                natomsph = rhosph * volsph / pseudofacsph
            elif volsph > 0.0:
                rhosph = (natomsph * pseudofacsph) / volsph
            else:
                rhosph = r4(0.0)
        else:
            rhosph = r4(0.0)
        # assign type number
        ato_manipulation_rev4.ctype(a1, ats, iexch, epsit, sigt, amass, qt, rqt, natomsph, rhosph, radminsph,
                                    radmaxsph, pseudofacsph)
    #
    # read core energy and diameter - used in fmole for non-bonded atoms
    #
    ierr = 0
    try:
        r = units.list_read(10)
        molecom.ecore = r.real(molecom.ecore)
        molecom.dcore = r.real(molecom.dcore)
    except FortranEOF:
        ierr = -1
    except FortranIOError as e:
        ierr = e.iostat
    #
    # initialise random number generator
    #
    if ierr == 0:
        auxiliary_routines.initran1(10)
        #
        # read the names of the template files if present. Normally these should be
        # specified in the order in which the different molecule types occur in the
        # .ato file.
        #
        icount = 0
        while icount < (clustercom.nmoltype - nmoltypefirst) and ierr == 0:
            icount = icount + 1
            commands.sometext, ierr = _read_a(10, commands.sometext, commands.nstsize)
            #
            # parse this text into words
            #
            if ierr == 0:
                commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 11,
                                              commands.nstsize, commands.nworddim)
                if commands.nwords > 1:
                    #
                    # get the molecular type this line belongs to
                    #
                    ierr = 0
                    try:
                        imolread = read_internal(substr(commands.sometext, ncf[1], ncl[1])).int(imolread)
                    except FortranEOF:
                        ierr = -1
                    except FortranIOError as e:
                        ierr = e.iostat
                    imoltype = nmoltypefirst + imolread
                    clustercom.moltemplate[imoltype] = fstr(substr(commands.sometext, ncf[2], ncl[2]), 256)
                    clustercom.interntrans[imoltype] = 0.0
                    clustercom.internrot[imoltype] = 0.0
                    clustercom.wholerot[imoltype] = 0.0
                    clustercom.wholetrans[imoltype] = 0.0
                    # Read the rotational and translational step sizes for each molecule type if specified
                    if commands.nwords > 2:
                        ierr = 0
                        try:
                            clustercom.interntrans[imoltype] = read_internal(
                                substr(commands.sometext, ncf[3], ncl[3])).real(clustercom.interntrans[imoltype])
                        except FortranEOF:
                            ierr = -1
                        except FortranIOError as e:
                            ierr = e.iostat
                    if commands.nwords > 3:
                        ierr = 0
                        try:
                            clustercom.internrot[imoltype] = read_internal(
                                substr(commands.sometext, ncf[4], ncl[4])).real(clustercom.internrot[imoltype])
                        except FortranEOF:
                            ierr = -1
                        except FortranIOError as e:
                            ierr = e.iostat
                    if commands.nwords > 4:
                        ierr = 0
                        try:
                            clustercom.wholerot[imoltype] = read_internal(
                                substr(commands.sometext, ncf[5], ncl[5])).real(clustercom.wholerot[imoltype])
                        except FortranEOF:
                            ierr = -1
                        except FortranIOError as e:
                            ierr = e.iostat
                    if commands.nwords > 5:
                        ierr = 0
                        try:
                            clustercom.wholetrans[imoltype] = read_internal(
                                substr(commands.sometext, ncf[6], ncl[6])).real(clustercom.wholetrans[imoltype])
                        except FortranEOF:
                            ierr = -1
                        except FortranIOError as e:
                            ierr = e.iostat
                    if commands.nwords > 6:
                        ierr = 0
                        try:
                            molecom_mod.rfractmol[imoltype] = read_internal(
                                substr(commands.sometext, ncf[7], ncl[7])).real(molecom_mod.rfractmol[imoltype])
                        except FortranEOF:
                            ierr = -1
                        except FortranIOError as e:
                            ierr = e.iostat

        # set the error code back to 0 so that if the template file has not been specified
        # routines which do not need this template can still run

        ierr = 0
        units.close(10)

    else:

        units.close(10)
        return rmax, cell, upol, rho, amubar, amubar2, ndist, ierr

    # Write moltypes out as a check

    units.write(6, fmt(_F201_minitc))

    for imoltype in range(1, clustercom.nmoltype + 1):
        units.write(6, fmt(_F203_minitc, imoltype, Ch(clustercom.moltype[imoltype], 3),
                           substr(clustercom.moltemplate[imoltype], 1, 10),
                           clustercom.interntrans[imoltype], clustercom.internrot[imoltype],
                           clustercom.wholerot[imoltype], clustercom.wholetrans[imoltype],
                           molecom_mod.rfractmol[imoltype]))

    # set up lattice parameters for applying the minimum image convention

    lattice.alatx = lattice.alat
    lattice.blatx = lattice.blat * auxiliary_routines.cosd(lattice.phib)
    lattice.blaty = lattice.blat * auxiliary_routines.sind(lattice.phib)
    lattice.clatx = lattice.clat * auxiliary_routines.sind(lattice.thetac) * auxiliary_routines.cosd(lattice.phic)
    lattice.claty = lattice.clat * auxiliary_routines.sind(lattice.thetac) * auxiliary_routines.sind(lattice.phic)
    lattice.clatz = lattice.clat * auxiliary_routines.cosd(lattice.thetac)
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
    cell = clustercom.dmax
    #
    # setup molecular size factor
    #
    molecom.dmaxmol = clustercom.dmax * molecom.sizefactor
    #
    # determine maximum distance across box
    #
    # G(r) and potentials will be calculated to rmax, or maximum allowed distance
    # based on number of points and step size, whichever is smaller.
    #
    xmax = abs(lattice.alatx) + abs(lattice.blatx) + abs(lattice.clatx)
    ymax = abs(lattice.blaty) + abs(lattice.claty)
    zmax = abs(lattice.clatz)
    rmax = r4(0.5) * sqrtf(xmax * xmax + ymax * ymax + zmax * zmax)
    units.write(6, lw(rmax, cell))
    #
    # convert temperature in K to kJ/mole
    #
    molecom.simt = r4(0.008314) * molecom.simtk
    units.write(6, fmt(_F917, clustercom.nmol, molecom.natom, molecom.natomr, molecom.ntype, clustercom.dmax))
    #
    # calculate atomic fractions
    #
    molecom_mod.natomcg = r4(0.0)
    molecom_mod.natomrcg = r4(0.0)
    for it in range(1, molecom.ntype + 1):
        molecom.cfrac[it] = 0.0
        molecom.cfraccg[it] = 0.0
    for ia in range(1, molecom.natom + 1):
        it = int(molecom.itype[ia])
        molecom_mod.natomcg = molecom_mod.natomcg + molecom_mod.natomsphere[it] * molecom_mod.pseudofacmol[it]
        if molecom.qtype[ia] != 1:
            molecom_mod.natomrcg = molecom_mod.natomrcg + molecom_mod.natomsphere[it] * molecom_mod.pseudofacmol[it]
        # Beware! This gets confusing because the natomr, cfrac, rho
        # represents the fraction of spheres, whereas natomrcg, cfraccg, rhocg represents the fraction of atoms.
        molecom.cfraccg[it] = molecom.cfraccg[it] + molecom_mod.natomsphere[it] * molecom_mod.pseudofacmol[it]
        molecom.cfrac[it] = molecom.cfrac[it] + r4(1.0)
    #
    # calculate volume of box and atomic number density
    #
    clustercom.volume = lattice.alatx * lattice.blaty * lattice.clatz
    rho = f32(molecom.natomr) / clustercom.volume
    molecom_mod.rhocg = molecom_mod.natomrcg / clustercom.volume
    units.write(6, lw(clustercom.volume, molecom_mod.natomcg, molecom_mod.natomrcg, rho, molecom_mod.rhocg))
    if molecom_mod.natomrcg > 0.0:
        for it in range(1, molecom.ntype + 1):
            molecom.cfraccg[it] = molecom.cfraccg[it] / molecom_mod.natomrcg
            molecom.cfrac[it] = molecom.cfrac[it] / f32(molecom.natomr)
    elif molecom_mod.natomcg > 0.0:
        for it in range(1, molecom.ntype + 1):
            molecom.cfraccg[it] = molecom.cfraccg[it] / molecom_mod.natomcg
            molecom.cfrac[it] = molecom.cfrac[it] / f32(molecom.natom)
    molecom.sumcharge = r4(0.0)
    for it in range(1, molecom.ntype + 1):
        units.write(6, fmt(_F916_minitc, it, Ch(molecom.type[it], 3), molecom.cfrac[it], molecom.cfraccg[it]))
        molecom.sumcharge = molecom.sumcharge + molecom.cfrac[it] * molecom.charge[it]
    #
    # calculate mass of each molecule
    #
    for im in range(1, clustercom.nmol + 1):
        clustercom.amasst[im] = 0.0
        for ia in range(int(clustercom.naf[im]), int(clustercom.nal[im]) + 1):
            if molecom.qtype[ia] == 0:
                amass = molecom.amassa[molecom.itype[ia]]
                clustercom.amasst[im] = clustercom.amasst[im] + amass
        clustercom.amasstsqrt[im] = sqrtf(clustercom.amasst[im])
        #
        # set maximum distance of each atom from centre of mass
        #
        clustercom.moldiam[im] = 0.0
        for ia in range(int(clustercom.naf[im]), int(clustercom.nal[im]) + 1):
            sum_ = molecom.xa[ia] * molecom.xa[ia] + molecom.ya[ia] * molecom.ya[ia] + molecom.za[ia] * molecom.za[ia]
            if sum_ > clustercom.moldiam[im]:
                clustercom.moldiam[im] = sum_
        clustercom.moldiam[im] = sqrtf(clustercom.moldiam[im])
        #           write(6,*) im,moldiam(im)
    ip = 0
    for it in range(1, molecom.ntype + 1):
        for jt in range(it, molecom.ntype + 1):
            ip = ip + 1
            molecom.n1extr[ip] = it
            molecom.n2extr[ip] = jt
            molecom.npairref[it, jt] = ip
            molecom.npairref[jt, it] = ip
            ect = sqrtf(molecom.epscor[it] * molecom.epscor[jt])
            molecom.ecor[ip] = ect
            molecom.qcor[ip] = molecom.charge[it] * molecom.charge[jt] * f32(1390)
            dct = r4(0.5) * (molecom.sigcor[it] + molecom.sigcor[jt])
            molecom.dcor[ip] = dct
            molecom.dscor[ip] = dct * dct
            am1 = molecom.amassa[it]
            am2 = molecom.amassa[jt]
            #                  amassr(ip)=tempb*simt*sqrt(am1*am2/(am1+am2))
            molecom.amassr[ip] = sqrtf(am1 * am2 / (am1 + am2))
            #                 write(6,*) it,jt,sigcor(it),sigcor(jt),dcor(ip)
    molecom.npairs = ip
    ndist = molecom.npairs
    units.write(6, fmt("(5(1x,i0))", clustercom.nmol, molecom.natomr, molecom.ntype, molecom.npairs,
                       molecom_mod.rhosphere.shape[0] - 1))
    #
    #     calculate dipole moments of molecules. Currently it is assumed the
    # the default dipole moment is that of an isolated water molecule
    #
    amu0 = r4(1.85)
    amu0 = amu0 / r4(4.80273)
    amusum = r4(0.0)
    amubar = r4(0.0)
    molecom.nbar = 0
    for im in range(1, clustercom.nmol + 1):
        if clustercom.nato[im] > 1:
            qx = r4(0.0)
            qy = r4(0.0)
            qz = r4(0.0)
            for ia in range(int(clustercom.naf[im]), int(clustercom.nal[im]) + 1):
                qt = molecom.charge[molecom.itype[ia]]
                qx = qx + qt * molecom.xa[ia]
                qy = qy + qt * molecom.ya[ia]
                qz = qz + qt * molecom.za[ia]
            amu2 = qx * qx + qy * qy + qz * qz
            amu = sqrtf(amu2)
            amud = amu * r4(4.80273)
            molecom.nbar = molecom.nbar + 1
            amubar = amubar + amud
            amubar2 = amubar2 + amud * amud
            #                 write(6,997) im,amu*4.80273
            amudif = amu - amu0
            amusum = amusum + amudif * amudif
    if molecom.nbar > 0:
        amubar = amubar / f32(molecom.nbar)
        amubar2 = amubar2 / f32(molecom.nbar) - amubar * amubar
        if amubar2 > 0.0:
            amubar2 = sqrtf(amubar2)
        else:
            amubar2 = r4(0.0)
        amusum = amusum / f32(clustercom.nmol)
        #
        # calculate static energy correction due to induced dipoles (c.f. Berendsen)
        #
        upol = r4(0.5) * amusum * f32(1390) / r4(1.445)
    else:
        amubar = r4(0.0)
        upol = r4(0.0)
    return rmax, cell, upol, rho, amubar, amubar2, ndist, ierr


@routine
def normalise_coordinates():
    """ato_read_and_save.f:normalise_coordinates -- lengths -> box-fraction units.

    sets up atom and molecule coordinates normalised to box dimension
    Returns nothing.
    """
    from ..common import coordconstraints_mod, molecom_mod, threebody_mod
    from . import ato_manipulation_rev4

    fractcoord = ato_manipulation_rev4.fractcoord
    dmax = clustercom.dmax
    #
    # fractional step sizes
    #
    molecom.stepmi = fractcoord(molecom.stepmi, molecom.dmaxmol)
    molecom.stepm = fractcoord(molecom.stepm, dmax)
    lattice.alatx = fractcoord(lattice.alatx, dmax)
    lattice.blatx = fractcoord(lattice.blatx, dmax)
    lattice.clatx = fractcoord(lattice.clatx, dmax)
    lattice.blaty = fractcoord(lattice.blaty, dmax)
    lattice.claty = fractcoord(lattice.claty, dmax)
    lattice.clatz = fractcoord(lattice.clatz, dmax)
    lattice.alat = fractcoord(lattice.alat, dmax)
    lattice.blat = fractcoord(lattice.blat, dmax)
    lattice.clat = fractcoord(lattice.clat, dmax)
    for i in range(1, 4):
        lattice.alatc[i] = fractcoord(lattice.alatc[i], dmax)
        lattice.blatc[i] = fractcoord(lattice.blatc[i], dmax)
        lattice.clatc[i] = fractcoord(lattice.clatc[i], dmax)
    if molecom.stepm > 0.25:
        molecom.stepm = r4(0.25)
    for im in range(1, clustercom.nmol + 1):
        #
        # set up fractional coordinates
        #
        clustercom.xadd[im] = fractcoord(clustercom.xadd[im], dmax)
        clustercom.yadd[im] = fractcoord(clustercom.yadd[im], dmax)
        clustercom.zadd[im] = fractcoord(clustercom.zadd[im], dmax)
        if clustercom.tethered[im]:
            clustercom.phix[im] = fractcoord(clustercom.phix[im], dmax)
            clustercom.phiy[im] = fractcoord(clustercom.phiy[im], dmax)
            clustercom.phiz[im] = fractcoord(clustercom.phiz[im], dmax)
        #
        # set maximum distance of each atom from centre of mass
        #
        clustercom.moldiam[im] = fractcoord(clustercom.moldiam[im], molecom.dmaxmol)
        for ia in range(int(clustercom.naf[im]), int(clustercom.nal[im]) + 1):
            molecom.xa[ia] = fractcoord(molecom.xa[ia], molecom.dmaxmol)
            molecom.ya[ia] = fractcoord(molecom.ya[ia], molecom.dmaxmol)
            molecom.za[ia] = fractcoord(molecom.za[ia], molecom.dmaxmol)
            if molecom.nbond[ia] > 0:
                for ib in range(1, int(molecom.nbond[ia]) + 1):
                    molecom.tbond[ia, ib] = fractcoord(molecom.tbond[ia, ib], molecom.dmaxmol)
    for ip in range(1, molecom.npairs + 1):
        molecom.amassr[ip] = molecom.amassr[ip] * molecom.dmaxmol
    # If any dihedral angles are defined, convert them to radians
    pi = PI4ATAN
    piconv = pi / r4(180.0)
    ir = 0
    while ir < molecom_mod.ngrrot:
        ir = ir + 1
        if molecom_mod.rotgrp[ir] == "DIH":
            for k in range(1, int(molecom.ndihedralang[ir]) + 1):
                molecom.dihedralang[k, ir] = molecom.dihedralang[k, ir] * piconv
    # If any three-body potentials are defined, normalise the maximum distances
    if threebody_mod.nthreeb > 0:
        for ir in range(1, threebody_mod.nthreeb + 1):
            for j in range(1, 3):
                threebody_mod.threebrmax[j, ir] = fractcoord(threebody_mod.threebrmax[j, ir], dmax)
    # If any coordination constraints are defined, normalise specified bond lengths
    if coordconstraints_mod.ncoordconstraints > 0:
        for ir in range(1, coordconstraints_mod.ncoordconstraints + 1):
            cc = coordconstraints_mod.coordconstraints[ir]
            for j in range(1, cc.ns + 1):
                units.write(6, lw(ir, j, cc.ds[j], dmax))
                cc.ds[j] = fractcoord(cc.ds[j], dmax)
                cc.dsmax[j] = fractcoord(cc.dsmax[j], dmax)
            #            call fractcoord(coordconstraints(ir)%ccdistort,dmax)
            #            if(coordconstraints(ir)%hasthreebody) then
            #               call fractcoord(coordconstraints(ir)%threebodydistort,dmax)
            #            endif
    return


@routine
def unnormalise_coordinates():
    """ato_read_and_save.f:unnormalise_coordinates -- box-fraction units -> lengths.

    sets up atom and molecule coordinates in unnormalised units
    Returns nothing.
    """
    from ..common import coordconstraints_mod, molecom_mod, threebody_mod
    from . import ato_manipulation_rev4

    wholecoord = ato_manipulation_rev4.wholecoord
    dmax = clustercom.dmax
    #
    # fractional step sizes
    #
    molecom.stepmi = wholecoord(molecom.stepmi, molecom.dmaxmol)
    if molecom.stepm > 0.25:
        molecom.stepm = r4(0.25)
    molecom.stepm = wholecoord(molecom.stepm, dmax)
    lattice.alatx = wholecoord(lattice.alatx, dmax)
    lattice.blatx = wholecoord(lattice.blatx, dmax)
    lattice.clatx = wholecoord(lattice.clatx, dmax)
    lattice.blaty = wholecoord(lattice.blaty, dmax)
    lattice.claty = wholecoord(lattice.claty, dmax)
    lattice.clatz = wholecoord(lattice.clatz, dmax)
    lattice.alat = wholecoord(lattice.alat, dmax)
    lattice.blat = wholecoord(lattice.blat, dmax)
    lattice.clat = wholecoord(lattice.clat, dmax)
    for i in range(1, 4):
        lattice.alatc[i] = wholecoord(lattice.alatc[i], dmax)
        lattice.blatc[i] = wholecoord(lattice.blatc[i], dmax)
        lattice.clatc[i] = wholecoord(lattice.clatc[i], dmax)
    for im in range(1, clustercom.nmol + 1):
        #
        # set up fractional coordinates
        #
        clustercom.xadd[im] = wholecoord(clustercom.xadd[im], dmax)
        clustercom.yadd[im] = wholecoord(clustercom.yadd[im], dmax)
        clustercom.zadd[im] = wholecoord(clustercom.zadd[im], dmax)
        if clustercom.tethered[im]:
            clustercom.phix[im] = wholecoord(clustercom.phix[im], dmax)
            clustercom.phiy[im] = wholecoord(clustercom.phiy[im], dmax)
            clustercom.phiz[im] = wholecoord(clustercom.phiz[im], dmax)
        #
        # set maximum distance of each atom from centre of mass
        #
        clustercom.moldiam[im] = wholecoord(clustercom.moldiam[im], molecom.dmaxmol)
        for ia in range(int(clustercom.naf[im]), int(clustercom.nal[im]) + 1):
            molecom.xa[ia] = wholecoord(molecom.xa[ia], molecom.dmaxmol)
            molecom.ya[ia] = wholecoord(molecom.ya[ia], molecom.dmaxmol)
            molecom.za[ia] = wholecoord(molecom.za[ia], molecom.dmaxmol)
            if molecom.nbond[ia] > 0:
                for ib in range(1, int(molecom.nbond[ia]) + 1):
                    molecom.tbond[ia, ib] = wholecoord(molecom.tbond[ia, ib], molecom.dmaxmol)
    for ip in range(1, molecom.npairs + 1):
        molecom.amassr[ip] = molecom.amassr[ip] / molecom.dmaxmol
    # If any dihedral angles are defined, convert them to degrees
    pi = PI4ATAN
    piconv = r4(180.0) / pi
    ir = 0
    while ir < molecom_mod.ngrrot:
        ir = ir + 1
        if molecom_mod.rotgrp[ir] == "DIH":
            for k in range(1, int(molecom.ndihedralang[ir]) + 1):
                molecom.dihedralang[k, ir] = molecom.dihedralang[k, ir] * piconv
    # If any three-body potentials are defined, unnormalise the maximum distances
    if threebody_mod.nthreeb > 0:
        for ir in range(1, threebody_mod.nthreeb + 1):
            for j in range(1, 3):
                threebody_mod.threebrmax[j, ir] = wholecoord(threebody_mod.threebrmax[j, ir], dmax)
    # If any coordination constraints are defined, unnormalise specified bond lengths
    if coordconstraints_mod.ncoordconstraints > 0:
        for ir in range(1, coordconstraints_mod.ncoordconstraints + 1):
            cc = coordconstraints_mod.coordconstraints[ir]
            for j in range(1, cc.ns + 1):
                cc.ds[j] = wholecoord(cc.ds[j], dmax)
                cc.dsmax[j] = wholecoord(cc.dsmax[j], dmax)
            #            call wholecoord(coordconstraints(ir)%ccdistort,dmax)
            #            if(coordconstraints(ir)%hasthreebody) then
            #               call wholecoord(coordconstraints(ir)%threebodydistort,dmax)
            #            endif
    return


_F207 = "(1x,i7,7(1x,e15.8))"
_F208 = "(1x,6(1x,e15.8))"
_F206 = "(1x,10(1x,e12.5))"
_F201 = "(1x,a3,2(2x,i5))"
_F202 = "(8(1x,f14.8))"
_F203 = "(1x,i4,(5(1x,i4,1x,e9.3)))"
_F204 = "(1x,14(i4,1x))"
_F209 = "(1x,2(i4,1x),f10.3)"
_F205 = "(1x,i32)"
_F211 = "(1x,a3,1x,a3,1x,i2)"
_F212 = "(1x,i3,1x,a,4(1x,e13.6),1x,f10.6)"
_F200 = "(1x,i0,6(1x,e15.8),2x,a6,1x,i0,1x,f10.5,1x,i5)"
_F599_moloutc = "(/'moloutc> Problems with file ',a/)"
_F598 = "(/'moloutc> WARNING! ROT groups defined BUT step size is ZERO!')"
_F915 = "(' Atomic number density (coarse grain and actual) = ',2(1x,e14.7))"
_F916_moloutc = "(1x,a3,2x,i1,2x,e12.5,10x,e12.5,2x,a3)"


@routine
def moloutc(fname1):
    """ato_read_and_save.f:moloutc -- write the configuration to a .ato file.

    Also writes the random number generator state (saveran1), the molecule templates
    and the atomic composition, and unlocks fname1.lck if it exists.  cosd/sind are
    not declared external here, so they are gfortran's intrinsics (unlike minitc).
    Returns nothing.
    """
    from ..common import molecom_mod
    from . import auxiliary_routines, reallocation_routines

    # integer, dimension(:), allocatable :: moltypeswritten
    # character*256 fname1
    # character*6 tetheredtext
    # logical found
    # If molecular scale factors have not been assigned, assign them
    #      write(6,*) 'moloutc> ',nmol,size(scalefacmol),allocated(scalefacmol)
    #      write(6,*) 'moloutc> ',ntype,size(rhosphere),allocated(rhosphere)
    if molecom_mod.scalefacmol is None:
        molecom_mod.scalefacmol = reallocation_routines.reallocate1d_r(molecom_mod.scalefacmol, clustercom.nmol)
        molecom_mod.scalefacmol[:] = 1.0
    elif molecom_mod.scalefacmol.shape[0] - 1 < clustercom.nmol:
        molecom_mod.scalefacmol = reallocation_routines.reallocate1d_r(molecom_mod.scalefacmol, clustercom.nmol)
        molecom_mod.scalefacmol[:] = 1.0
    try:
        units.open(10, fname1, status="unknown")
        ierr = 0
    except FortranIOError as e:
        ierr = e.iostat
    #
    # Rewind in case a previous routine has left it open
    #
    units.rewind(10)
    #
    # if there is a problem with this file, print out an error message and exit
    #
    if ierr != 0:
        ilenfname = findex(fname1, ".")
        ilenfname = ilenfname + findex(substr(fname1, ilenfname, 256), " ") - 1
        units.write(6, fmt(_F599_moloutc, substr(fname1, 1, ilenfname)))
        clustercom.nmol = 0
        return

    #     Write a warning message if ROT groups are defined, but internal rotations are
    #     switched off.

    if molecom_mod.ngrrot > 0 and molecom.stepri <= 0.0:
        units.write(6, fmt(_F598))

    ig = 0
    #
    # write no. of molecules and tolerance on bond distances
    #
    if lattice.cubic:

        units.write(10, fmt(_F207, clustercom.nmol, clustercom.dmax, molecom.simtk))

    else:

        units.write(10, fmt(_F207, clustercom.nmol, molecom.simtk))
        units.write(10, fmt(_F208, lattice.alat, lattice.blat, lattice.clat))
        units.write(10, fmt(_F208, lattice.phib, lattice.thetac, lattice.phic))

    # set up lattice parameters for calculating volume and density
    lattice.alatx = lattice.alat
    lattice.blatx = lattice.blat * cosd_intrinsic(lattice.phib)
    lattice.blaty = lattice.blat * sind_intrinsic(lattice.phib)
    lattice.clatx = lattice.clat * sind_intrinsic(lattice.thetac) * cosd_intrinsic(lattice.phic)
    lattice.claty = lattice.clat * sind_intrinsic(lattice.thetac) * sind_intrinsic(lattice.phic)
    lattice.clatz = lattice.clat * cosd_intrinsic(lattice.thetac)
    units.write(10, fmt(_F206, molecom.tolout, molecom.stepmi, molecom.stepri, molecom.stepr, molecom.stepm,
                        molecom.tempb, molecom.tempangreq, molecom.tempdihreq, molecom.tempangwt, molecom.tempdihwt))
    molecom.natomr = 0
    molecom.natom = 0
    molecom_mod.natomcg = f32(0)
    molecom_mod.natomrcg = f32(0)
    for it in range(1, molecom.ntype + 1):
        molecom.cfrac[it] = 0.0
        molecom.cfraccg[it] = 0.0
    for im in range(1, clustercom.nmol + 1):
        if clustercom.tetheredatom[im] < 1 or clustercom.tetheredatom[im] > clustercom.nato[im]:
            clustercom.tetheredatom[im] = 0
        if clustercom.tethered[im]:
            tetheredtext = fstr(fmt("(l1,i5.5)", bool(clustercom.tethered[im]), clustercom.tetheredatom[im]), 6)
        else:
            tetheredtext = fstr(fmt("(l1)", bool(clustercom.tethered[im])), 6)
        units.write(10, fmt(_F200, clustercom.nato[im], clustercom.xadd[im], clustercom.yadd[im],
                            clustercom.zadd[im], clustercom.phix[im], clustercom.phiy[im], clustercom.phiz[im],
                            Ch(tetheredtext, 6), clustercom.nstepreject[im], molecom_mod.scalefacmol[im], im))
        #
        # integer to convert from absolute atom numbers to relative atom numbers
        #
        molecom.natom = molecom.natom + int(clustercom.nato[im])
        nadd = 1 - int(clustercom.naf[im])

        for ia in range(int(clustercom.naf[im]), int(clustercom.nal[im]) + 1):
            # for each atom write its label
            if molecom.imageatom[ia] > 0:
                units.write(10, fmt(_F201, Ch(molecom.alabel[ia], 3), (ia + nadd), (int(molecom.imageatom[ia]) + nadd)))
            else:
                units.write(10, fmt(_F201, Ch(molecom.alabel[ia], 3), (ia + nadd), 0))

            # write its starting coordinates and sigmas
            units.write(10, fmt(_F202, molecom.xa[ia], molecom.ya[ia], molecom.za[ia]))
            #
            # write the number of atoms to which it must be bound, the atom numbers of
            # those neighbours and the anticipated bond length
            #
            vals = [molecom.nbond[ia]]
            for j in range(1, int(molecom.nbond[ia]) + 1):
                vals.append(int(molecom.iabond[ia, j]) + nadd)
                vals.append(molecom.tbond[ia, j])
            units.write(10, fmt(_F203, *vals))
            it = int(molecom.itype[ia])
            molecom_mod.natomcg = molecom_mod.natomcg + molecom_mod.natomsphere[it] * molecom_mod.pseudofacmol[it]
            if molecom.qtype[ia] == 0:
                molecom.natomr = molecom.natomr + 1
                molecom_mod.natomrcg = (molecom_mod.natomrcg
                                        + molecom_mod.natomsphere[it] * molecom_mod.pseudofacmol[it])
            molecom.cfrac[it] = molecom.cfrac[it] + r4(1.0)
            molecom.cfraccg[it] = molecom.cfraccg[it] + molecom_mod.natomsphere[it] * molecom_mod.pseudofacmol[it]
        # write number of internal rotation groups
        units.write(10, fmt(_F204, clustercom.ngrot[im]))
        if clustercom.ngrot[im] > 0:
            for ih in range(1, int(clustercom.ngrot[im]) + 1):
                ig = ig + 1
                units.write(10, fmt(_F201, Ch(molecom_mod.rotgrp[ig], 3)))
                # for rotational groups read atom numbers of two atoms to form axis of rotation
                units.write(10, fmt(_F204, (int(molecom_mod.labaxi[ig, 1]) + nadd),
                                    (int(molecom_mod.labaxi[ig, 2]) + nadd)))
                # For DIH groups read two atom numbers to calculate the dihedral angle and the required dihedral angle
                if molecom_mod.rotgrp[ig] == "DIH":
                    units.write(10, fmt(_F204, molecom.ndihedralang[ig]))
                    for k in range(1, int(molecom.ndihedralang[ig]) + 1):
                        units.write(10, fmt(_F209, *[(int(molecom.labrotdih[j, k, ig]) + nadd) for j in range(1, 3)],
                                            molecom.dihedralang[k, ig]))
                # write number of atoms in group and their numbers
                units.write(10, fmt(_F204, molecom_mod.nrot[ig],
                                    *[(int(molecom_mod.labrot[j]) + nadd)
                                      for j in range(int(molecom_mod.rotf[ig]), int(molecom_mod.rotl[ig]) + 1)]))
    if molecom_mod.natomsphere is None:
        molecom_mod.natomsphere = reallocation_routines.reallocate1d_r(molecom_mod.natomsphere, molecom.ntype)
        molecom_mod.rhosphere = reallocation_routines.reallocate1d_r(molecom_mod.rhosphere, molecom.ntype)
        molecom_mod.radminsphere = reallocation_routines.reallocate1d_r(molecom_mod.radminsphere, molecom.ntype)
        molecom_mod.radmaxsphere = reallocation_routines.reallocate1d_r(molecom_mod.radmaxsphere, molecom.ntype)
        molecom_mod.pseudofacmol = reallocation_routines.reallocate1d_r(molecom_mod.pseudofacmol, molecom.ntype)
        molecom_mod.natomsphere[:] = 1.0
        molecom_mod.rhosphere[:] = 0.0
        molecom_mod.radminsphere[:] = 0.0
        for it in range(1, molecom.ntype + 1):
            molecom_mod.radmaxsphere[it] = r4(0.5) * molecom.sigcor[it]
        molecom_mod.pseudofacmol[:] = 1.0
    elif molecom_mod.natomsphere.shape[0] - 1 < molecom.ntype:
        molecom_mod.natomsphere = reallocation_routines.reallocate1d_r(molecom_mod.natomsphere, molecom.ntype)
        molecom_mod.rhosphere = reallocation_routines.reallocate1d_r(molecom_mod.rhosphere, molecom.ntype)
        molecom_mod.radminsphere = reallocation_routines.reallocate1d_r(molecom_mod.radminsphere, molecom.ntype)
        molecom_mod.radmaxsphere = reallocation_routines.reallocate1d_r(molecom_mod.radmaxsphere, molecom.ntype)
        molecom_mod.pseudofacmol = reallocation_routines.reallocate1d_r(molecom_mod.pseudofacmol, molecom.ntype)
        molecom_mod.natomsphere[:] = 1.0
        molecom_mod.rhosphere[:] = 0.0
        molecom_mod.radminsphere[:] = 0.0
        for it in range(1, molecom.ntype + 1):
            molecom_mod.radmaxsphere[it] = r4(0.5) * molecom.sigcor[it]
        molecom_mod.pseudofacmol[:] = 1.0
    for it in range(1, molecom.ntype + 1):
        units.write(10, fmt(_F211, Ch(molecom.type[it], 3), Ch(molecom.atsymbol[it], 3), molecom.iexchange[it]))
        units.write(10, fmt("(5(1x,f10.6),1x,e13.6,2(1x,f10.4),1x,e13.6)", molecom.epscor[it], molecom.sigcor[it],
                            molecom.amassa[it], molecom.charge[it], molecom.qradius[it],
                            molecom_mod.natomsphere[it], molecom_mod.radminsphere[it],
                            molecom_mod.radmaxsphere[it], molecom_mod.pseudofacmol[it]))
    # write core energy and diameter which represents the closest overlap of ANY
    # non-bonded atom
    units.write(10, fmt(_F202, molecom.ecore, molecom.dcore))
    #
    # save random number generator parameters
    #
    auxiliary_routines.saveran1(10)
    #
    # write out the molecular template files if defined
    # In the order in which they occur in the .ato file!
    #
    moltypeswritten = int_array(clustercom.nmoltype)
    #      do imoltype=1,nmoltype
    #         moltypeswritten(imoltype)=0
    #      end do
    # In case the relative broadening for each molecule has not been defined, define them
    if molecom_mod.rfractmol is None:
        molecom_mod.rfractmol = reallocation_routines.reallocate1d_r(molecom_mod.rfractmol, clustercom.nmoltype)
        molecom_mod.rfractmol[:] = 0.0
    elif molecom_mod.rfractmol.shape[0] - 1 < clustercom.nmoltype:
        molecom_mod.rfractmol = reallocation_routines.reallocate1d_r(molecom_mod.rfractmol, clustercom.nmoltype)
        molecom_mod.rfractmol[:] = 0.0
    nmoltypeswritten = 0
    for im in range(1, clustercom.nmol + 1):
        # Search for a moltype corresponding to this molecule, but if already written, ignore.
        found = False
        imoltype = 0
        # First find the molecule type of this molecule
        while imoltype < clustercom.nmoltype and not found:
            imoltype = imoltype + 1
            found = molecom.alabel[clustercom.naf[im]] == clustercom.moltype[imoltype]
        # Now search the list of written molecules
        found = False
        imoltypewritten = 0
        while imoltypewritten < nmoltypeswritten and not found:
            imoltypewritten = imoltypewritten + 1
            found = imoltype == moltypeswritten[imoltypewritten]
        if not found:
            nmoltypeswritten = nmoltypeswritten + 1
            ilen = findex_padded(clustercom.moltemplate[imoltype], " ", 256) - 1
            units.write(10, fmt(_F212, nmoltypeswritten, substr(clustercom.moltemplate[imoltype], 1, ilen),
                                clustercom.interntrans[imoltype], clustercom.internrot[imoltype],
                                clustercom.wholerot[imoltype], clustercom.wholetrans[imoltype],
                                molecom_mod.rfractmol[imoltype]))
            moltypeswritten[nmoltypeswritten] = imoltype
    #
    # calculate volume of box and atomic number density
    #
    clustercom.volume = lattice.alatx * lattice.blaty * lattice.clatz
    rho = f32(molecom.natomr) / clustercom.volume
    molecom_mod.rhocg = molecom_mod.natomrcg / clustercom.volume
    if molecom.natomr > 0:
        for it in range(1, molecom.ntype + 1):
            molecom.cfrac[it] = molecom.cfrac[it] / f32(molecom.natomr)
            molecom.cfraccg[it] = molecom.cfraccg[it] / molecom_mod.natomrcg
    #      write(6,915) rho
    units.write(10, fmt(_F915, molecom_mod.rhocg, rho))
    # Write out the atomic composition in the form required by Gudrun
    for it in range(1, molecom.ntype + 1):
        #         if(qtypea(it).eq.0) write(6,916) atsymbol(it),0,cfrac(it),type(it)
        if molecom.qtypea[it] == 0:
            units.write(10, fmt(_F916_moloutc, Ch(molecom.atsymbol[it], 3), 0, molecom.cfraccg[it],
                                molecom.cfrac[it], Ch(molecom.type[it], 3)))

    units.close(10)

    # Open the lock file if it exists and set to false (= unlocked)

    try:
        units.open(10, fname1 + ".lck", status="old")
        ierr = 0
    except FortranIOError as e:
        ierr = e.iostat
    if ierr == 0:
        units.write(10, fmt("(a)", "F"))
        units.close(10)

    return
