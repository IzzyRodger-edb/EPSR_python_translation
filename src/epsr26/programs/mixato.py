"""Translation of fortran_src/mixato.f: ``program standalone`` and ``subroutine mixato``
-- mixes up to mcomp .ato files (each containing one molecule type).
"""

from ..common import clustercom, commands, directory_name, epsr_dimension, lattice, molecom, molecom_mod, symbols
from ..fortran import (
    Ch,
    FortranEOF,
    FortranIOError,
    cosd_intrinsic,
    f32,
    findex,
    fmt,
    fstr,
    iostat_of,
    lw,
    powf,
    r4,
    read_internal,
    sind_intrinsic,
    substr,
    units,
)
from ..library import ato_manipulation_rev4
from ..library import ato_read_and_save
from ..library import file_and_directory_commands
from ..library import molecom_mod as molecom_mod_lib
from ..library import parse
from ..library import reallocation_routines
from ..trace import routine

# FORMAT statements ('$' dropped: written with advance=False)
_F580 = "(/'mixato> How many .ATO files do you want to mix? ')"
_F699 = "(//'mixato> Search for .ato file ',i2)"
_F700 = "(/'mixato> No .ato files found in the current directory - try another directory.'/)"
_F499 = "(/'mixato> No .ato file specified - try again!'/)"
_F604 = "(/'mixato> ',a,'.ato contains more than 1 type of molecule'//'mixato> This is not allowed - start again')"
_F36 = "(/1x,i3,' new atom types in file ',a/)"
_F33 = "(1x,'Atom type ',i2,' has label ',a3)"
_F581 = "(/'mixato> How many of these molecules do you want in the mixture? ')"
_F946 = "(1x,'No. of molecules requested .gt. than allowed dimension:',i6)"
_F599 = "(/'Problems with .ato file ',a,' Try again.'/)"
_F601 = "(/'mixato> Tried to define more molecules than allowed by dimension'/)"
_F602 = "(/'mixato> Tried to define more atoms than allowed by dimension'/)"
_F34 = "(/'mixato> Give atomic number density (per A**3) of mixture: ')"
_F35 = "(/'mixato> Type name of file to put mixture in: ')"
_F101 = "(a)"


def main():
    """mixato.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    mixato()


@routine
def mixato():
    """mixato.f:mixato -- mixes up to mcomp .ato files. If any of these contain more than
    one molecule, mixato will fail."""
    # use reallocation_routines
    # use molecom_mod
    mcomp = epsr_dimension.mcomp
    near = epsr_dimension.near
    mmol = epsr_dimension.mmol
    matom = epsr_dimension.matom
    mmix = mcomp  # noqa: F841  parameter(mmix=mcomp,ma=near)
    ncf = commands.ncf
    ncl = commands.ncl
    # character*256 readfile ; character*4 ext ; character*1 iaccept
    iaccept = ""  # NB: not initialised in the Fortran; search_for_files sets it
    # integer, allocatable :: nmincl,natos,ngrots,nmoltypes(mmix)
    nmincl = natos = ngrots = nmoltypes = None
    # integer, allocatable :: nbondss,nrots,ndihedralangs (near,mmix)
    nbondss = nrots = ndihedralangs = None
    # integer, allocatable :: iabonds,labaxis,labrots
    iabonds = labaxis = labrots = None
    # integer, allocatable :: labrotdihs(2,near,near,mmix)
    labrotdihs = None
    # Save coordinates and bonds
    xadds = yadds = zadds = phixs = phiys = phizs = None
    xas = yas = zas = None
    tbonds = dihedralangs = None
    # Save step sizes and vibrational temperature
    # integer natomax,ngrotmax,nrotmax,ndihmax
    # real stepmis,stepris,steprs,stepms,tempbs,natomrnew
    # character(len=3), allocatable :: alabels,rotgrps
    alabels = rotgrps = None
    # implicitly typed locals passed to minitc (epsr16.inc is not included)
    rmaxcell = cellmax = upol = rho = amubar = amubar2 = f32(0.0)
    ndist = 0
    nmix = 0  # NB: uninitialised in the Fortran until read

    ierr = 1
    if commands.nwords > 1:
        try:
            nmix = read_internal(substr(commands.sometext, ncf[2], ncl[2])).int(nmix)
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = iostat_of(e)
    else:
        ierr = 1
    while ierr != 0:
        units.write(6, fmt(_F580), advance=False)
        try:
            nmix = units.list_read(5).int(nmix)
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = iostat_of(e)
    nmincl = reallocation_routines.reallocate1d_i(nmincl, nmix)
    natos = reallocation_routines.reallocate1d_i(natos, nmix)
    ngrots = reallocation_routines.reallocate1d_i(ngrots, nmix)
    nmoltypes = reallocation_routines.reallocate1d_i(nmoltypes, nmix)
    xadds = reallocation_routines.reallocate1d_r(xadds, nmix)
    yadds = reallocation_routines.reallocate1d_r(yadds, nmix)
    zadds = reallocation_routines.reallocate1d_r(zadds, nmix)
    phixs = reallocation_routines.reallocate1d_r(phixs, nmix)
    phiys = reallocation_routines.reallocate1d_r(phiys, nmix)
    phizs = reallocation_routines.reallocate1d_r(phizs, nmix)

    nm1 = 0
    natomax = 0
    ngrotmax = 0
    nrotmax = 0
    ndihmax = 0
    #
    # set the number of atom and molecule types to zero at the outset
    #
    molecom.ntype = 0
    clustercom.nmoltype = 0
    ntypefirst = molecom.ntype
    nmoltypefirst = clustercom.nmoltype

    # Set the maximum step sizes to zero

    stepmis = f32(0)
    stepris = f32(0)
    steprs = f32(0)
    stepms = f32(0)
    tempbs = f32(0)
    tempdihs = f32(0)

    imix = 0
    while imix < nmix:
        imix = imix + 1
        units.write(6, fmt(_F699, imix))
        #
        # generate the name of each .ato file in turn
        #
        ext = ".ato"
        ilenext = 4
        commands.nwords = 0
        while commands.nwords <= 0:
            iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)
            #
            # if no .ato files are available in this directory, abandon this attempt at mixato
            #
            if iaccept != "y":
                units.write(6, fmt(_F700))
                return
            commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 5, commands.nstsize,
                                          commands.nworddim)
            if commands.nwords <= 0:
                units.write(6, fmt(_F499))
        readfile = fstr(substr(commands.sometext, ncf[1], ncl[1]), 256)
        ilenreadfile = ncl[1] - ncf[1] + 1
        #
        # strip out any extension that was supplied - extension must be .ato, otherwise
        # nothing is done.
        #
        ilendot = findex(readfile, ".") - 1
        if ilendot > 0:
            readfile = fstr(substr(readfile, 1, ilendot), 256)
            ilenreadfile = ilendot
        #
        # create the full name of the input file by adding on the directory name if not present
        #
        symbols.fnameato, symbols.ilenfnameato = file_and_directory_commands.get_full_filename(
            readfile, ilenreadfile, symbols.fnameato, symbols.ilenfnameato)
        symbols.fnameato = fstr(substr(symbols.fnameato, 1, symbols.ilenfnameato) + ".ato", 256)
        symbols.ilenfnameato = symbols.ilenfnameato + 4
        #
        # read .ato file
        #
        molecom.sizefactor = f32(1.0)
        (rmaxcell, cellmax, upol, rho, amubar, amubar2, ndist,
         ierr) = ato_read_and_save.minitc(symbols.fnameato, rmaxcell, cellmax, upol, rho, amubar, amubar2,
                                          ndist, ierr)
        if ierr == 0:
            ntypefile = molecom.ntype - ntypefirst
            nmoltypefile = clustercom.nmoltype - nmoltypefirst
            #
            # because we have to read the entire file to get all the relevant information, a
            # .ato file which contains more than 1 molecule type cannot be mixed
            #
            if nmoltypefile > 1:
                units.write(6, fmt(_F604, substr(readfile, 1, ilenreadfile)))
                return
            ntypefirst = molecom.ntype
            nmoltypefirst = clustercom.nmoltype
            units.write(6, fmt(_F36, ntypefile, substr(symbols.fnameato, 1, symbols.ilenfnameato)))
            for it in range(molecom.itype[1], molecom.itype[molecom.natom] + 1):
                units.write(6, fmt(_F33, it, Ch(molecom.type[it], 3)))
            units.write(6, fmt(_F581), advance=False)
            nmincl[imix] = units.list_read(5).int(nmincl[imix])
            #
            # generate the extra molecules based on the last one read in
            #
            nm1 = nm1 + nmincl[imix]
            if nm1 > mmol:
                units.write(6, fmt(_F946, mmol))
                return
            #
            # save the parameters of the first molecule
            #
            # Use the largest of the step sizes

            if molecom.stepmi > stepmis:
                stepmis = molecom.stepmi
            if molecom.stepri > stepris:
                stepris = molecom.stepri
            if molecom.stepr > steprs:
                steprs = molecom.stepr
            if molecom.stepm > stepms:
                stepms = molecom.stepm
            if molecom.tempb > tempbs:
                tempbs = molecom.tempb
            if molecom.tempdih > tempdihs:
                tempdihs = molecom.tempdih

            natos[imix] = clustercom.nato[1]
            xadds[imix] = clustercom.xadd[1]
            yadds[imix] = clustercom.yadd[1]
            zadds[imix] = clustercom.zadd[1]
            phixs[imix] = clustercom.phix[1]
            phiys[imix] = clustercom.phiy[1]
            phizs[imix] = clustercom.phiz[1]
            if clustercom.nato[1] > natomax:
                natomax = int(clustercom.nato[1])
                alabels = reallocation_routines.reallocate2d_c(alabels, 3, natomax, nmix)
                xas = reallocation_routines.reallocate2d_r(xas, natomax, nmix)
                yas = reallocation_routines.reallocate2d_r(yas, natomax, nmix)
                zas = reallocation_routines.reallocate2d_r(zas, natomax, nmix)
                nbondss = reallocation_routines.reallocate2d_i(nbondss, natomax, nmix)
                iabonds = reallocation_routines.reallocate3d_i(iabonds, natomax, near, nmix)
                tbonds = reallocation_routines.reallocate3d_r(tbonds, natomax, near, nmix)
            for ia in range(clustercom.naf[1], clustercom.nal[1] + 1):
                alabels[ia, imix] = molecom.alabel[ia]
                xas[ia, imix] = molecom.xa[ia]
                yas[ia, imix] = molecom.ya[ia]
                zas[ia, imix] = molecom.za[ia]
                nbondss[ia, imix] = molecom.nbond[ia]
                if molecom.nbond[ia] > 0:
                    for ib in range(1, molecom.nbond[ia] + 1):
                        iabonds[ia, ib, imix] = molecom.iabond[ia, ib]
                        tbonds[ia, ib, imix] = molecom.tbond[ia, ib]
            #
            # rotational groups
            #
            ngrots[imix] = clustercom.ngrot[1]
            if clustercom.ngrot[1] > ngrotmax:
                ngrotmax = int(clustercom.ngrot[1])
                rotgrps = reallocation_routines.reallocate2d_c(rotgrps, 3, ngrotmax, nmix)
                nrots = reallocation_routines.reallocate2d_i(nrots, ngrotmax, nmix)
                labaxis = reallocation_routines.reallocate3d_i(labaxis, ngrotmax, 2, nmix)
                ndihedralangs = reallocation_routines.reallocate2d_i(ndihedralangs, ngrotmax, nmix)
                dihedralangs = reallocation_routines.reallocate3d_r(dihedralangs, ndihmax, ngrotmax, nmix)
                labrotdihs = reallocation_routines.reallocate4d_i(labrotdihs, 2, ndihmax, ngrotmax, nmix)
                labrots = reallocation_routines.reallocate3d_i(labrots, ngrotmax, nrotmax, nmix)
            if clustercom.ngrot[1] > 0:
                for ig in range(1, clustercom.ngrot[1] + 1):
                    rotgrps[ig, imix] = molecom_mod.rotgrp[ig]
                    labaxis[ig, 1, imix] = molecom_mod.labaxi[ig, 1]
                    labaxis[ig, 2, imix] = molecom_mod.labaxi[ig, 2]
                    nrots[ig, imix] = molecom_mod.nrot[ig]
                    if molecom_mod.nrot[ig] > nrotmax:
                        nrotmax = int(molecom_mod.nrot[ig])
                        labrots = reallocation_routines.reallocate3d_i(labrots, ngrotmax, nrotmax, nmix)
                    for j in range(1, molecom_mod.nrot[ig] + 1):
                        labref = molecom_mod.rotf[ig] + j - 1
                        labrots[ig, j, imix] = molecom_mod.labrot[labref]
                    # !Check if there are dihedral angles defined
                    if molecom_mod.rotgrp[ig] == "DIH":
                        ndihedralangs[ig, imix] = molecom.ndihedralang[ig]
                        if molecom.ndihedralang[ig] > ndihmax:
                            ndihmax = int(molecom.ndihedralang[ig])
                            dihedralangs = reallocation_routines.reallocate3d_r(dihedralangs, ndihmax, ngrotmax,
                                                                                nmix)
                            labrotdihs = reallocation_routines.reallocate4d_i(labrotdihs, 2, ndihmax, ngrotmax,
                                                                              nmix)
                        for k in range(1, molecom.ndihedralang[ig] + 1):
                            dihedralangs[k, ig, imix] = molecom.dihedralang[k, ig]
                            for j in range(1, 3):
                                labrotdihs[j, k, ig, imix] = molecom.labrotdih[j, k, ig]
            # !Find out if the supplied template file actually exists.
            nmt = clustercom.nmoltype
            try:
                units.open(35, clustercom.moltemplate[nmt] + ".ato", status="old")
                ierr = 0
            except (FortranEOF, FortranIOError) as e:
                ierr = iostat_of(e)
            units.close(35)
            # !write(6,*) 'mixato> Found or not? ',...
            if ierr != 0:
                units.write(6, lw("mixato> Did not find template file: ", clustercom.moltemplate[nmt] + ".ato"))
                # !Create a name based on the name of the current file being read from
                clustercom.moltemplate[nmt] = fstr(substr(readfile, 1, ilenreadfile), 256)
            units.write(6, lw())
            units.write(6, lw(" mixato> moltemplate file: ", clustercom.moltemplate[nmt] + ".ato"))

        else:
            molecom.ntype = ntypefirst
            clustercom.nmoltype = nmoltypefirst
            units.write(6, fmt(_F599, Ch(symbols.fnameato, 256)))
            imix = imix - 1
    #
    # O.K. so now generate the new molecular coordinates based on the values for the
    # individual molecules read in
    #
    ia = 0
    im = 0
    ig = 0
    molecom.natom = 0
    molecom.natomr = 0
    molecom_mod.natomcg = f32(0)
    molecom_mod.natomrcg = f32(0)
    molecom_mod.ngrrot = 0
    molecom_mod.nrottot = 0
    for imix in range(1, nmix + 1):
        for it in range(1, nmincl[imix] + 1):
            if im < mmol:
                im = im + 1
            else:
                units.write(6, fmt(_F601))
                return
            clustercom.xadd[im] = xadds[imix]
            clustercom.yadd[im] = yadds[imix]
            clustercom.zadd[im] = zadds[imix]
            clustercom.nato[im] = natos[imix]
            nadd = ia
            clustercom.naf[im] = ia + 1
            for is_ in range(1, clustercom.nato[im] + 1):
                if ia < matom:
                    ia = ia + 1
                else:
                    units.write(6, fmt(_F602))
                    return
                molecom.alabel[ia] = alabels[is_, imix]
                ato_manipulation_rev4.ftype(ia)
                itia = molecom.itype[ia]
                molecom_mod.natomcg = molecom_mod.natomcg + molecom_mod.natomsphere[itia] * molecom_mod.pseudofacmol[
                    itia]
                #
                # do not count atoms which begin with 'Q' or 'q' - these are assumed
                # invisible to the scattering radiation and experiment (although they
                # can be used to calculate other structural quantities if necessary
                #
                if molecom.qtype[ia] != 1:
                    molecom.natomr = molecom.natomr + 1
                    molecom_mod.natomrcg = molecom_mod.natomrcg + molecom_mod.natomsphere[itia] * \
                        molecom_mod.pseudofacmol[itia]
                molecom.xa[ia] = xas[is_, imix]
                molecom.ya[ia] = yas[is_, imix]
                molecom.za[ia] = zas[is_, imix]
                molecom.nbond[ia] = nbondss[is_, imix]
                if molecom.nbond[ia] > 0:
                    for ib in range(1, molecom.nbond[ia] + 1):
                        molecom.iabond[ia, ib] = iabonds[is_, ib, imix] + nadd
                        molecom.tbond[ia, ib] = tbonds[is_, ib, imix]
            clustercom.nal[im] = ia
            #
            # rotational groups
            #
            clustercom.ngrot[im] = ngrots[imix]
            if clustercom.ngrot[im] > 0:
                molecom_mod.ngrrot = molecom_mod.ngrrot + int(clustercom.ngrot[im])
                molecom_mod_lib.testsize_rotgrps()
                for ih in range(1, clustercom.ngrot[im] + 1):
                    ig = ig + 1
                    # save molecule in which this rotational group occurs
                    molecom_mod.mrot[ig] = im
                    molecom_mod.rotgrp[ig] = rotgrps[ih, imix]
                    # for rotational groups read atom numbers of two atoms to form axis of rotation
                    for j in range(1, 3):
                        molecom_mod.labaxi[ig, j] = labaxis[ih, j, imix] + nadd
                    # read number of atoms in group and their numbers
                    molecom_mod.nrot[ig] = nrots[ih, imix]
                    molecom_mod.rotf[ig] = molecom_mod.nrottot + 1
                    molecom_mod.nrottot = molecom_mod.nrottot + int(molecom_mod.nrot[ig])
                    molecom_mod.rotl[ig] = molecom_mod.nrottot
                    molecom_mod_lib.testsize_labrot()
                    for j in range(1, molecom_mod.nrot[ig] + 1):
                        labref = molecom_mod.rotf[ig] + j - 1
                        molecom_mod.labrot[labref] = labrots[ih, j, imix] + nadd
                    if molecom_mod.rotgrp[ig] == "DIH":
                        molecom.ndihedralang[ig] = ndihedralangs[ih, imix]
                        for k in range(1, molecom.ndihedralang[ig] + 1):
                            molecom.dihedralang[k, ig] = dihedralangs[k, ih, imix]
                            for j in range(1, 3):
                                molecom.labrotdih[j, k, ig] = labrotdihs[j, k, ih, imix] + nadd
    clustercom.nmol = im
    molecom.natom = ia
    molecom_mod.ngrrot = ig
    units.write(6, fmt(_F34), advance=False)
    rhoatom = units.list_read(5).real(f32(0.0))
    clustercom.dmax = molecom_mod.natomrcg / rhoatom
    clustercom.dmax = powf(clustercom.dmax, r4(1.0) / r4(3.0))
    #
    # normalise all the coordinates to this new value of dmax
    #
    molecom.dmaxmol = molecom.sizefactor * clustercom.dmax

    # ! Set up the lattice dimensions for calculating the density in moloutc.

    lattice.cubic = True
    lattice.alat = clustercom.dmax
    lattice.blat = clustercom.dmax
    lattice.clat = clustercom.dmax
    lattice.phib = f32(90.0)
    lattice.thetac = f32(0.0)
    lattice.phic = f32(0.0)
    lattice.alatx = lattice.alat
    lattice.blatx = lattice.blat * cosd_intrinsic(lattice.phib)
    lattice.blaty = lattice.blat * sind_intrinsic(lattice.phib)
    lattice.clatx = lattice.clat * sind_intrinsic(lattice.thetac) * cosd_intrinsic(lattice.phic)
    lattice.claty = lattice.clat * sind_intrinsic(lattice.thetac) * sind_intrinsic(lattice.phic)
    lattice.clatz = lattice.clat * cosd_intrinsic(lattice.thetac)

    # Set the step sizes to the largest found in all the ato files read in

    molecom.stepmi = stepmis
    molecom.stepri = f32(0.3)
    molecom.stepr = f32(0.3)
    molecom.stepm = f32(1.0)
    molecom.tempb = tempbs
    molecom.tempdih = tempdihs
    #
    # get the new filename to put the mixture in
    #
    commands.nwords = 0
    while commands.nwords <= 0:
        units.write(6, fmt(_F35), advance=False)
        try:
            commands.sometext = fstr(units.read_line(5), commands.nstsize)
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = iostat_of(e)
        commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 5, commands.nstsize,
                                      commands.nworddim)
        if commands.nwords <= 0:
            units.write(6, fmt(_F499))
    readfile = fstr(substr(commands.sometext, ncf[1], ncl[1]), 256)
    ilenreadfile = ncl[1] - ncf[1] + 1
    #
    # strip out any extension that was supplied - extension must be .ato, otherwise
    # nothing is done.
    #
    ilendot = findex(readfile, ".") - 1
    if ilendot > 0:
        readfile = fstr(substr(readfile, 1, ilendot), 256)
        ilenreadfile = ilendot
    symbols.fnameato = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                            + substr(readfile, 1, ilenreadfile) + ".ato", 256)
    symbols.ilenfnameato = directory_name.nchardir + ilenreadfile + 4
    #
    # write the new .ato file
    #
    ato_read_and_save.moloutc(symbols.fnameato)
    return
