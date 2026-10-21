"""Translation of fortran_src/addato.f: ``program standalone`` (``main``) and ``addato``.

adds up to 20 .ato files to an existing .ato file.  See docs/translation_notes/N1.md.
"""

from ..common import clustercom
from ..common import commands
from ..common import directory_name
from ..common import lattice
from ..common import molecom
from ..common import molecom_mod
from ..common import symbols
from ..common.epsr_dimension import matom, mmol, near
from ..fortran import Ch
from ..fortran import FortranEOF
from ..fortran import FortranIOError
from ..fortran import PI4ATAN
from ..fortran import f32
from ..fortran import findex
from ..fortran import fmt
from ..fortran import fstr
from ..fortran import int_array
from ..fortran import iostat_of
from ..fortran import lw
from ..fortran import pad
from ..fortran import r4
from ..fortran import read_internal
from ..fortran import substr
from ..fortran import units
from ..library import ato_manipulation_rev4
from ..library import ato_read_and_save
from ..library import auxiliary_routines
from ..library import file_and_directory_commands
from ..library import minimagec
from ..library import molfit
from ..library import parse
from ..library import reallocation_routines
from ..library import molecom_mod as molecom_mod_routines
from ..trace import routine

# 580 format(/'addato> How many .ATO files do you want to add? ',$)
_F580 = "(/'addato> How many .ATO files do you want to add? ')"
_F699 = "(//'addato> Search for .ato file to be ADDED ')"
_F700 = "(/'addato> No .ato files found in the current directory - try another directory.'/)"
_F499 = "(/'addato> No .ato file specified - try again!'/)"
_F604 = "(/'addato> ',a,'.ato contains more than 1 type of molecule'//'addato> This is not allowed - try again')"
_F36 = "(/1x,i3,' new atom types in file ',a/)"
_F33 = "(1x,'Atom type ',i2,' has label ',a3)"
_F581 = "(/'addato> How many of these molecules do you want in the mixture? ')"
_F946 = "(1x,'No. of molecules requested .gt. than allowed dimension:',i6)"
_F599 = "(/'addato> Problems with .ato file ',a,' Try again.'/)"
_F698 = "(//'addato> Search for .ato file to be ADDED to')"
_F601 = "(/'addato> Tried to define more molecules than allowed by dimension'/)"
_F697 = "(/'addato> Could find no insertion point after ',i4,' attempts.')"
_F602 = "(/'addato> Tried to define more atoms than allowed by dimension'/)"
_F603 = "(/'addato> Tried to define more rotational groups than allowed by dimension'/)"
# 35 format(/'addato> Type name of file to put mixture in: ',$)
_F35 = "(/'addato> Type name of file to put mixture in: ')"
_F101 = "(a)"


def main():
    """addato.f:standalone -- Program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    addato()


@routine
def addato():
    """addato.f:addato -- add molecules from up to 20 .ato files to an existing .ato file.

    If any of these contain more than one molecule, addato will use only the first molecule
    encountered. Only files containing one type of molecule can be added.
    The box size of the .ato file to which these molecules cannot be changed, so it is up to
    the user to add the correct number of molecules to get the correct density. The added
    molecules will be placed at random, except that they will not overlap with any tethered
    molecules - the degree of overlap is determined by dcore which is read in from the file
    being added to.
    """
    # use reallocation_routines
    # use molecom_mod
    # parameter(mmix=20,ma=near)
    m = molecom
    c = clustercom
    cm = commands
    mm = molecom_mod
    lat = lattice
    ra = reallocation_routines

    # character*256 readfile
    # character*4 ext
    # character*1 iaccept
    readfile = ""
    iaccept = ""
    # NB: ext and ilenext are only set inside the first loop; if it does not run (nmix <= 0)
    # the Fortran uses them uninitialised.  Blank / 0 here.
    ext = ""
    ilenext = 0
    nmix = 0  # NB: uninitialised in the Fortran if the first read fails; 0 here
    ilenreadfile = 0

    nmincl = natos = ngrots = None  # integer(:) allocatable (mmix)
    nbondss = nrots = itypes = ndihedralangs = None  # integer(:,:) allocatable
    iabonds = labaxis = labrots = None  # integer(:,:,:) allocatable
    labrotdihs = None  # integer(:,:,:,:) allocatable

    # Save coordinates and bonds

    xadds = yadds = zadds = phixs = phiys = phizs = None  # real(:) allocatable (mmix)
    xas = yas = zas = None  # real(:,:) allocatable (ma,mmix)
    tbonds = dihedralangs = None  # real(:,:,:) allocatable

    # Save step sizes and vibrational temperature

    # integer natomax,ngrotmax,nrotmax,ndihmax
    # real stepmis,stepris,steprs,stepms,tempbs,rminsqr,xdif,ydif,zdif,xdifc,ydifc,zdifc,pi,piconv

    alabels = rotgrps = None  # character(len=3)(:,:) allocatable

    # logical :: overlap
    tnumber = None  # integer(:) allocatable
    # integer :: it,itimes,ntimes,itn,ntn,im,jm,iat,jat
    tethereds = None  # logical(:) allocatable

    rmaxcell = f32(0.0)
    cellmax = f32(0.0)
    upol = f32(0.0)
    rho = f32(0.0)
    amubar = f32(0.0)
    amubar2 = f32(0.0)
    ndist = 0

    ierr = 1
    if cm.nwords > 1:
        try:
            nmix = read_internal(substr(cm.sometext, cm.ncf[2], cm.ncl[2])).int(nmix)
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
    nmincl = ra.reallocate1d_i(nmincl, nmix)
    natos = ra.reallocate1d_i(natos, nmix)
    ngrots = ra.reallocate1d_i(ngrots, nmix)
    xadds = ra.reallocate1d_r(xadds, nmix)
    yadds = ra.reallocate1d_r(yadds, nmix)
    zadds = ra.reallocate1d_r(zadds, nmix)
    phixs = ra.reallocate1d_r(phixs, nmix)
    phiys = ra.reallocate1d_r(phiys, nmix)
    phizs = ra.reallocate1d_r(phizs, nmix)
    tethereds = ra.reallocate1d_l(tethereds, nmix)

    nm1 = 0
    nmixmax = nmix
    natomax = 0
    ngrotmax = 0
    nrotmax = 0
    ndihmax = 0
    #
    # set the number of atom and molecule types to zero at the outset
    #
    m.ntype = 0
    c.nmoltype = 0
    ntypefirst = m.ntype
    nmoltypefirst = c.nmoltype

    # Set the maximum step sizes to zero

    stepmis = f32(0)
    stepris = f32(0)
    steprs = f32(0)
    stepms = f32(0)
    tempbs = f32(0)
    tempdihs = f32(0)

    imix = 0
    while imix < nmix:
        units.write(6, fmt(_F699))
        #
        # generate the name of each .ato file in turn
        #
        ext = fstr(".ato", 4)
        ilenext = 4
        cm.nwords = 0
        while cm.nwords <= 0:
            iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)
            #
            # if no .ato files are available in this directory, abandon this attempt at addato
            #
            if iaccept != "y":
                units.write(6, fmt(_F700))
                return
            cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
            if cm.nwords <= 0:
                units.write(6, fmt(_F499))
        readfile = fstr(substr(cm.sometext, cm.ncf[1], cm.ncl[1]), 256)
        ilenreadfile = cm.ncl[1] - cm.ncf[1] + 1
        #
        # strip out any extension that was supplied - extension must be .ato, otherwise
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
        symbols.fnameato = fstr(substr(symbols.fnameato, 1, symbols.ilenfnameato) + ".ato", 256)
        symbols.ilenfnameato = symbols.ilenfnameato + 4
        #
        # read .ato file
        #
        m.sizefactor = f32(1.0)
        (rmaxcell, cellmax, upol, rho, amubar, amubar2, ndist, ierr) = ato_read_and_save.minitc(
            symbols.fnameato, rmaxcell, cellmax, upol, rho, amubar, amubar2, ndist, ierr)
        if ierr == 0:
            ntypefile = m.ntype - ntypefirst
            nmoltypefile = c.nmoltype - nmoltypefirst
            #
            # because we have to read the entire file to get all the relevant information, a
            # .ato file which contains more than 1 molecule type cannot be mixed
            #
            if nmoltypefile > 1:
                units.write(6, fmt(_F604, substr(readfile, 1, ilenreadfile)))
            else:
                ntypefirst = m.ntype
                nmoltypefirst = c.nmoltype
                units.write(6, fmt(_F36, ntypefile, substr(symbols.fnameato, 1, symbols.ilenfnameato)))
                for it in range(m.itype[1], m.itype[m.natom] + 1):
                    units.write(6, fmt(_F33, it, Ch(m.type[it], 3)))

                # If a molecule type is tethered, we use all the molecules, otherwise we use only
                # the first molecule and specify the number of times we want this molecule

                if c.moltethered[c.nmoltype]:
                    nmix = nmix + c.nmol - 1
                    nmoladd = 1
                    immax = c.nmol
                else:
                    units.write(6, fmt(_F581), advance=False)
                    nmoladd = units.list_read(5).int(0)
                    immax = 1
                im = 0
                irgref = 0
                while im < immax:
                    im = im + 1
                    imix = imix + 1
                    if imix > nmixmax:
                        nmixmax = nmixmax + 10
                        nmincl = ra.reallocate1d_i(nmincl, nmixmax)
                        natos = ra.reallocate1d_i(natos, nmixmax)
                        ngrots = ra.reallocate1d_i(ngrots, nmixmax)
                        xadds = ra.reallocate1d_r(xadds, nmixmax)
                        yadds = ra.reallocate1d_r(yadds, nmixmax)
                        zadds = ra.reallocate1d_r(zadds, nmixmax)
                        phixs = ra.reallocate1d_r(phixs, nmixmax)
                        phiys = ra.reallocate1d_r(phiys, nmixmax)
                        phizs = ra.reallocate1d_r(phizs, nmixmax)
                        tethereds = ra.reallocate1d_l(tethereds, nmixmax)
                    nmincl[imix] = nmoladd
                    nm1 = nm1 + int(nmincl[imix])
                    natos[imix] = c.nato[im]
                    xadds[imix] = c.xadd[im]
                    yadds[imix] = c.yadd[im]
                    zadds[imix] = c.zadd[im]
                    phixs[imix] = c.phix[im]
                    phiys[imix] = c.phiy[im]
                    phizs[imix] = c.phiz[im]
                    tethereds[imix] = c.tethered[im]
                    iaadd = int(c.naf[im]) - 1
                    if c.nato[im] > natomax:
                        natomax = int(c.nato[im])
                        alabels = ra.reallocate2d_c(alabels, 3, natomax, nmixmax)
                        itypes = ra.reallocate2d_i(itypes, natomax, nmixmax)
                        xas = ra.reallocate2d_r(xas, natomax, nmixmax)
                        yas = ra.reallocate2d_r(yas, natomax, nmixmax)
                        zas = ra.reallocate2d_r(zas, natomax, nmixmax)
                        nbondss = ra.reallocate2d_i(nbondss, natomax, nmixmax)
                        iabonds = ra.reallocate3d_i(iabonds, natomax, near, nmixmax)
                        tbonds = ra.reallocate3d_r(tbonds, natomax, near, nmixmax)
                    for ir in range(1, natos[imix] + 1):
                        ia = iaadd + ir
                        alabels[ir, imix] = m.alabel[ia]
                        itypes[ir, imix] = m.itype[ia]
                        xas[ir, imix] = m.xa[ia]
                        yas[ir, imix] = m.ya[ia]
                        zas[ir, imix] = m.za[ia]
                        nbondss[ir, imix] = m.nbond[ia]
                        if m.nbond[ia] > 0:
                            for ib in range(1, m.nbond[ia] + 1):
                                iabonds[ir, ib, imix] = m.iabond[ia, ib]
                                tbonds[ir, ib, imix] = m.tbond[ia, ib]
                    #
                    # rotational groups
                    #
                    ngrots[imix] = c.ngrot[im]
                    if c.ngrot[im] > ngrotmax:
                        ngrotmax = int(c.ngrot[im])
                        rotgrps = ra.reallocate2d_c(rotgrps, 3, ngrotmax, nmixmax)
                        nrots = ra.reallocate2d_i(nrots, ngrotmax, nmixmax)
                        labaxis = ra.reallocate3d_i(labaxis, ngrotmax, 2, nmixmax)
                        ndihedralangs = ra.reallocate2d_i(ndihedralangs, ngrotmax, nmix)
                        dihedralangs = ra.reallocate3d_r(dihedralangs, ndihmax, ngrotmax, nmix)
                        labrotdihs = ra.reallocate4d_i(labrotdihs, 2, ndihmax, ngrotmax, nmix)
                        labrots = ra.reallocate3d_i(labrots, ngrotmax, nrotmax, nmix)
                    if c.ngrot[im] > 0:
                        for ig in range(1, c.ngrot[im] + 1):
                            irgref = irgref + 1
                            rotgrps[ig, imix] = mm.rotgrp[irgref]
                            labaxis[ig, 1, imix] = mm.labaxi[irgref, 1]
                            labaxis[ig, 2, imix] = mm.labaxi[irgref, 2]
                            nrots[ig, imix] = mm.nrot[irgref]
                            if mm.nrot[irgref] > nrotmax:
                                nrotmax = int(mm.nrot[irgref])
                                labrots = ra.reallocate3d_i(labrots, ngrotmax, nrotmax, nmix)
                            # NB: the Fortran indexes nrot, rotf and rotgrp with ig here (not irgref)
                            for j in range(1, mm.nrot[ig] + 1):
                                labref = int(mm.rotf[ig]) + j - 1
                                labrots[ig, j, imix] = mm.labrot[labref]
                            # Check if there are dihedral angles defined
                            if mm.rotgrp[ig] == "DIH":
                                ndihedralangs[ig, imix] = m.ndihedralang[ig]
                                if m.ndihedralang[ig] > ndihmax:
                                    ndihmax = int(m.ndihedralang[ig])
                                    dihedralangs = ra.reallocate3d_r(dihedralangs, ndihmax, ngrotmax, nmix)
                                    labrotdihs = ra.reallocate4d_i(labrotdihs, 2, ndihmax, ngrotmax, nmix)
                                for k in range(1, m.ndihedralang[ig] + 1):
                                    dihedralangs[k, ig, imix] = m.dihedralang[k, ig]
                                    for j in range(1, 3):
                                        labrotdihs[j, k, ig, imix] = m.labrotdih[j, k, ig]

                # Save the template file name for output later
                # First find out if the supplied template file actually exists.
                # (the file name is the blank-padded character*256 moltemplate followed by '.ato')
                try:
                    units.open(35, pad(c.moltemplate[c.nmoltype], 256) + ".ato", status="old")
                    ierr = 0
                except (FortranEOF, FortranIOError) as e:
                    ierr = iostat_of(e)
                units.close(35)
                if ierr != 0:
                    # Create a name based on the name of the current file being read from
                    c.moltemplate[c.nmoltype] = fstr(substr(readfile, 1, ilenreadfile), 256)
                units.write(6, lw())
                units.write(6, lw(" addato> moltemplate file: ", c.moltemplate[c.nmoltype] + ".ato"))

                # Use the largest of the step sizes

                if m.stepmi > stepmis:
                    stepmis = m.stepmi
                if m.stepri > stepris:
                    stepris = m.stepri
                if m.stepr > steprs:
                    steprs = m.stepr
                if m.stepm > stepms:
                    stepms = m.stepm
                if m.tempb > tempbs:
                    tempbs = m.tempb
                if m.tempdih > tempdihs:
                    tempdihs = m.tempdih

                if nm1 > mmol:
                    units.write(6, fmt(_F946, mmol))
                    return
        else:
            m.ntype = ntypefirst
            c.nmoltype = nmoltypefirst
            units.write(6, fmt(_F599, Ch(symbols.fnameato, 256)))

    # Now read the .ato file you want to ADD these molecules to

    while imix < nmix + 1:
        units.write(6, fmt(_F698))
        cm.nwords = 0
        while cm.nwords <= 0:
            iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)
            #
            # if no .ato files are available in this directory, abandon this attempt at addato
            #
            if iaccept != "y":
                units.write(6, fmt(_F700))
                return
            cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
            if cm.nwords <= 0:
                units.write(6, fmt(_F499))
        readfile = fstr(substr(cm.sometext, cm.ncf[1], cm.ncl[1]), 256)
        ilenreadfile = cm.ncl[1] - cm.ncf[1] + 1
        #
        # strip out any extension that was supplied - extension must be .ato, otherwise
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
        symbols.fnameato = fstr(substr(symbols.fnameato, 1, symbols.ilenfnameato) + ".ato", 256)
        symbols.ilenfnameato = symbols.ilenfnameato + 4
        #
        # read .ato file
        #
        m.sizefactor = f32(1.0)
        (rmaxcell, cellmax, upol, rho, amubar, amubar2, ndist, ierr) = ato_read_and_save.minitc(
            symbols.fnameato, rmaxcell, cellmax, upol, rho, amubar, amubar2, ndist, ierr)
        if ierr == 0:
            imix = imix + 1
            ntypefile = m.ntype - ntypefirst
            nmoltypefile = c.nmoltype - nmoltypefirst  # noqa: F841 (set but not used, as in the Fortran)
            units.write(6, fmt(_F36, ntypefile, substr(symbols.fnameato, 1, symbols.ilenfnameato)))
            for it in range(m.itype[1], m.itype[m.natom] + 1):
                units.write(6, fmt(_F33, it, Ch(m.type[it], 3)))

            # Use the largest of the step sizes

            if m.stepmi > stepmis:
                stepmis = m.stepmi
            if m.stepri > stepris:
                stepris = m.stepri
            if m.stepr > steprs:
                steprs = m.stepr
            if m.stepm > stepms:
                stepms = m.stepm
            if m.tempb > tempbs:
                tempbs = m.tempb
            if m.tempdih > tempdihs:
                tempdihs = m.tempdih

        else:
            m.ntype = ntypefirst
            c.nmoltype = nmoltypefirst
            units.write(6, fmt(_F599, Ch(symbols.fnameato, 256)))

    #
    # O.K. so now generate the new molecular coordinates based on the values for the
    # individual molecules read in
    #
    # First make a list of all tethered molecules: these will be used to test for overlaps

    ntn = 0
    for im in range(1, c.nmol + 1):
        if c.tethered[im]:
            ntn = ntn + 1
    # Gather the tethered molecule numbers
    # if(allocated(tnumber)) deallocate(tnumber)
    tnumber = int_array(ntn)  # allocate(tnumber(ntn))
    itn = 0
    for im in range(1, c.nmol + 1):
        if c.tethered[im]:
            itn = itn + 1
            tnumber[itn] = im
    units.write(6, lw("addato> Number of tethered molecules at start = ", ntn))

    # To be safe normalise the existing coordinates, add normalised coordinates, then
    # unnormalise them at the end.
    ato_read_and_save.normalise_coordinates()
    pi = PI4ATAN
    piconv = pi / f32(180.0)

    im = c.nmol
    ia = m.natom
    ig = mm.ngrrot
    # Box distance squared to test for overlaps
    rminsqr = f32(1.0) / (c.dmax * c.dmax)
    for imix in range(1, nmix + 1):
        for iloop in range(1, nmincl[imix] + 1):
            if im < mmol:
                im = im + 1
            else:
                units.write(6, fmt(_F601))
                return
            if not tethereds[imix]:
                # Set the tethered flag to false
                c.tethered[im] = False
                overlap = True
                itimes = 0
                # Allow 1000 attempts before giving up
                ntimes = 1000
                while overlap and itimes < ntimes:
                    itimes = itimes + 1
                    # Choose a position at random
                    fa = (auxiliary_routines.ran1() - f32(0.5))
                    fb = (auxiliary_routines.ran1() - f32(0.5))
                    fc = (auxiliary_routines.ran1() - f32(0.5))

                    # Convert to Cartesian coordinates. Note that by using the a,b,c axes
                    # formalism we avoid the need to use normalised coordinates for the
                    # mininum image convention.

                    c.zadd[im] = fc * lat.clatz
                    c.yadd[im] = fb * lat.blaty + fc * lat.claty
                    c.xadd[im] = fa * lat.alatx + fb * lat.blatx + fc * lat.clatx
                    c.phix[im] = f32(0.0)
                    c.phiy[im] = f32(0.0)
                    c.phiz[im] = f32(0.0)

                    # ensure all molecules are in the central box ( + or - 0.5)

                    c.xadd[im], c.yadd[im], c.zadd[im] = minimagec.minimagec(c.xadd[im], c.yadd[im], c.zadd[im])
                    # Step through the tethered atoms test for overlap
                    itn = 0
                    overlap = False
                    while not overlap and itn < ntn:
                        itn = itn + 1
                        jm = int(tnumber[itn])
                        xdifc = c.xadd[im] - c.xadd[jm]
                        ydifc = c.yadd[im] - c.yadd[jm]
                        zdifc = c.zadd[im] - c.zadd[jm]
                        # Compare the first atoms of each molecule for overlaps
                        iat = 0
                        while iat < natos[imix] and not overlap:
                            iat = iat + 1
                            it = int(itypes[iat, imix])
                            jat = int(c.naf[jm]) - 1
                            while jat < c.nal[jm] and not overlap:
                                jat = jat + 1
                                jt = int(m.itype[jat])
                                xdif = xdifc + xas[iat, imix] / m.dmaxmol - m.xa[jat]
                                ydif = ydifc + yas[iat, imix] / m.dmaxmol - m.ya[jat]
                                zdif = zdifc + zas[iat, imix] / m.dmaxmol - m.za[jat]
                                xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                                dtest = xdif * xdif + ydif * ydif + zdif * zdif
                                overlap = (m.ecor[m.npairref[it, jt]] > 0.0
                                           and dtest < rminsqr * m.dscor[m.npairref[it, jt]])
                                if overlap:
                                    units.write(6, lw(im, jm, iat, jat, it, jt, m.npairref[it, jt], dtest,
                                                      m.dscor[m.npairref[it, jt]], rminsqr))
            else:
                itimes = 1  # noqa: F841 (set but not used, as in the Fortran)
                ntimes = 1
                overlap = False
                c.tethered[im] = tethereds[imix]
                c.xadd[im] = xadds[imix] / c.dmax
                c.yadd[im] = yadds[imix] / c.dmax
                c.zadd[im] = zadds[imix] / c.dmax
                c.phix[im] = phixs[imix] / c.dmax
                c.phiy[im] = phiys[imix] / c.dmax
                c.phiz[im] = phizs[imix] / c.dmax

            if overlap:

                units.write(6, fmt(_F697, ntimes))
                im = im - 1

            else:

                c.nato[im] = natos[imix]
                nadd = ia
                c.naf[im] = ia + 1
                for is_ in range(1, c.nato[im] + 1):
                    if ia < matom:
                        ia = ia + 1
                    else:
                        units.write(6, fmt(_F602))
                        return
                    m.alabel[ia] = alabels[is_, imix]
                    units.write(6, lw(ia, is_, imix, Ch(alabels[is_, imix], 3)))
                    ato_manipulation_rev4.ftype(ia)
                    #
                    # do not count atoms which begin with 'Q' or 'q' - these are assumed
                    # invisible to the scattering radiation and experiment (although they
                    # can be used to calculate other structural quantities if necessary
                    #
                    if m.qtype[ia] != 1:
                        m.natomr = m.natomr + 1
                    m.xa[ia] = xas[is_, imix] / m.dmaxmol
                    m.ya[ia] = yas[is_, imix] / m.dmaxmol
                    m.za[ia] = zas[is_, imix] / m.dmaxmol
                    m.nbond[ia] = nbondss[is_, imix]
                    if m.nbond[ia] > 0:
                        for ib in range(1, m.nbond[ia] + 1):
                            m.iabond[ia, ib] = iabonds[is_, ib, imix] + nadd
                            m.tbond[ia, ib] = tbonds[is_, ib, imix] / m.dmaxmol
                c.nal[im] = ia
                #
                # rotational groups
                #
                c.ngrot[im] = ngrots[imix]
                if c.ngrot[im] > 0:
                    mm.ngrrot = mm.ngrrot + int(c.ngrot[im])
                    molecom_mod_routines.testsize_rotgrps()
                    for ih in range(1, c.ngrot[im] + 1):
                        if ig < matom:
                            ig = ig + 1
                        else:
                            units.write(6, fmt(_F603))
                            return
                        # save molecule in which this rotational group occurs
                        mm.mrot[ig] = im
                        mm.rotgrp[ig] = rotgrps[ih, imix]
                        # for rotational groups read atom numbers of two atoms to form axis of rotation
                        for j in range(1, 3):
                            mm.labaxi[ig, j] = labaxis[ih, j, imix] + nadd
                        # read number of atoms in group and their numbers
                        mm.nrot[ig] = nrots[ih, imix]
                        mm.rotf[ig] = mm.nrottot + 1
                        mm.nrottot = mm.nrottot + int(mm.nrot[ig])
                        mm.rotl[ig] = mm.nrottot
                        molecom_mod_routines.testsize_labrot()
                        for j in range(1, mm.nrot[ig] + 1):
                            labref = int(mm.rotf[ig]) + j - 1
                            mm.labrot[labref] = labrots[ih, j, imix] + nadd
                        if mm.rotgrp[ig] == "DIH":
                            m.ndihedralang[ig] = ndihedralangs[ih, imix]
                            for k in range(1, m.ndihedralang[ig] + 1):
                                m.dihedralang[k, ig] = dihedralangs[k, ih, imix] * piconv
                                for j in range(1, 3):
                                    m.labrotdih[j, k, ig] = labrotdihs[j, k, ih, imix] + nadd
                # perform random whole molecule rotation
                if not c.tethered[im]:
                    steprs = m.stepr
                    if steprs > 0.0:
                        m.stepr = f32(1.0)
                    molfit.rotatm(im, m.stepr)
                    m.stepr = steprs
    c.nmol = im
    m.natom = ia

    # Unnormalise the coordinates

    ato_read_and_save.unnormalise_coordinates()

    # Leave the step sizes to those found in last file

    m.stepmi = stepmis
    m.stepri = r4(0.3)
    m.stepr = r4(0.3)
    m.stepm = f32(1.0)
    m.tempb = tempbs
    m.tempdih = tempdihs
    #
    # get the new filename to put the mixture in
    #
    cm.nwords = 0
    while cm.nwords <= 0:
        units.write(6, fmt(_F35), advance=False)
        try:
            cm.sometext = fstr(units.read_line(5), cm.nstsize)
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = iostat_of(e)
        cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
        if cm.nwords <= 0:
            units.write(6, fmt(_F499))
    readfile = fstr(substr(cm.sometext, cm.ncf[1], cm.ncl[1]), 256)
    ilenreadfile = cm.ncl[1] - cm.ncf[1] + 1
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
