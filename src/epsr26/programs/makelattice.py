"""Translation of fortran_src/makelattice.f: ``program standalone`` and ``subroutine
makelattice``.

Makes an EPSR .ato file from a file of unit cell coordinates and potential parameters:
every atom becomes a one-atom molecule tethered to its lattice site, in a box of
na x nb x nc unit cells.  See docs/translation_notes/N2.md.
"""

import numpy as np

from ..common import clustercom, commands, directory_name, lattice, molecom, symbols
from ..fortran import (Ch, FortranEOF, FortranIOError, PI4ATAN, char_array, cosd_intrinsic, f32, findex, fmt,
                       fstr, lw, r4, read_internal, real_array, sind_intrinsic, substr, units)
from ..library import ato_manipulation_rev4, ato_read_and_save, auxiliary_routines, file_and_directory_commands, parse
from ..trace import routine

_F100 = "(/'makelattice> Type name of file containing unit cell coordinates (incl. extension): ')"
_F101 = "(a)"
_F499 = "(/'makelatice> No file specified! - try again!'/)"
_F606 = "('makelattice> Give number of unit cells along c-axis: ')"
_F607 = "('makelattice> Give number of unit cells along b and c-axes: ')"
_F608 = "('makelattice> Give number of unit cells along a, b and c-axes: ')"
_F602 = "('makelattice> Problems with file: ',a)"
# 603 format('makelattice> Lattice constants, angles and derived angles:'
#     13(/3(1x,f10.5)))       (the 1 is the continuation mark)
_F603 = "('makelattice> Lattice constants, angles and derived angles:',3(/3(1x,f10.5)))"
_F604 = "('makelattice> Insufficient data supplied in input file.')"
_F605 = "('makelattice> Number of atom types found: ',i3)"
# 35 format(/'makelattice> Type name of file to put lattice in: ',$)
_F35 = "(/'makelattice> Type name of file to put lattice in: ')"
_F609 = "('makelattice> Output being written to file: ',a)"
_F610 = "('makelattice> ',1x,i5,' molecules and ',i5,' atoms created')"
_F611 = "('makelattice> Centre of box at: ',3(1x,e12.5))"
_F201 = "(/'makelattice> Following molecule types defined initially:-'/)"
_F202 = "(1x,i3,1x,a3,1x,a10,2(1x,e12.5))"


def main():
    """makelattice.f:standalone -- 2009/3/30 Attempt to make this a standalone program."""
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    makelattice()


def _read_text(unit, old, n):
    """``read(unit,'(a)',iostat=ierr) text``: returns (text, ierr); text unchanged on error."""
    try:
        return fstr(units.read_line(unit), n), 0
    except FortranEOF:
        return old, -1
    except FortranIOError as e:
        return old, e.iostat


@routine
def makelattice():
    """makelattice.f:makelattice -- Makes an EPSR18 ato file from a set of supplied lattice
    coordinates.

    The unit cell coordinates are supplied in a file alongside other potential
    parameters - this filename will be given to the created .ato file.
    """
    # character(len=256)  :: readfile
    # integer             :: ilenreadfile,ierr,i,iexch,nwordstest
    # real                :: pi,radconv,alpha,beta,gamma,alphar,betar,sinphib,cosphib,cosalpha,cosbeta,top,bottom,sum
    # real, dimension(:), allocatable :: xunit,yunit,zunit !x,y,z coordinates in unit cell
    # integer             :: nunitcell,na,nb,nc,ka,kb,kc !Number of atoms in unit cell and number of steps along a, b and c axes
    # character(len=3), dimension(:), allocatable  :: atsymbolunit !Atom symbols in unit cell
    # character(len=3)    :: a1,ats
    # real                :: epsit,sigt,amass,qt,rqt,acent,bcent,ccent,fa,fb,fc,xav,yav,zav,xbox,ybox,zbox
    cm = commands
    cc = clustercom
    lt = lattice
    dn = directory_name
    na = nb = nc = 0
    epsit = sigt = amass = qt = rqt = f32(0.0)
    alpha = beta = gam = f32(0.0)
    iexch = 0

    pi = PI4ATAN
    radconv = pi / r4(180.0)  # noqa: F841
    degconv = r4(180.0) / pi  # noqa: F841
    #
    # generate the name of the file containing the unit cell coordinates
    #
    nwordstest = 2
    if cm.nwords < 2:
        cm.nwords = 0
        nwordstest = 1
        while cm.nwords <= 0:
            units.write(6, fmt(_F100), advance=False)
            cm.sometext = fstr(units.read_line(5), cm.nstsize)
            cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
            if cm.nwords <= 0:
                units.write(6, fmt(_F499))
        readfile = fstr(substr(cm.sometext, cm.ncf[1], cm.ncl[1]), 256)
        ilenreadfile = int(cm.ncl[cm.nwords]) - int(cm.ncf[1]) + 1
    else:
        readfile = fstr(substr(cm.sometext, cm.ncf[2], cm.ncl[2]), 256)
        ilenreadfile = int(cm.ncl[2]) - int(cm.ncf[2]) + 1
    units.write(6, lw(cm.nwords, nwordstest))
    if cm.nwords > nwordstest:
        na = read_internal(substr(cm.sometext, cm.ncf[nwordstest + 1], cm.ncl[nwordstest + 1])).int(na)
        if cm.nwords > nwordstest + 1:
            nb = read_internal(substr(cm.sometext, cm.ncf[nwordstest + 2], cm.ncl[nwordstest + 2])).int(nb)
            if cm.nwords > nwordstest + 2:
                nc = read_internal(substr(cm.sometext, cm.ncf[nwordstest + 3], cm.ncl[nwordstest + 3])).int(nc)
            else:
                units.write(6, fmt(_F606), advance=False)
                nc = units.list_read(5).int(nc)
        else:
            units.write(6, fmt(_F607), advance=False)
            r = units.list_read(5)
            nb = r.int(nb)
            nc = r.int(nc)
    else:
        units.write(6, fmt(_F608), advance=False)
        r = units.list_read(5)
        na = r.int(na)
        nb = r.int(nb)
        nc = r.int(nc)
    #
    # Read this file
    #
    symbols.fnameato = fstr(substr(dn.directory_name, 1, dn.nchardir) + substr(readfile, 1, ilenreadfile), 256)
    symbols.ilenfnameato = dn.nchardir + ilenreadfile
    try:
        units.open(10, symbols.fnameato, status='old')
        ierr = 0
    except FortranIOError as e:
        ierr = e.iostat
    if ierr == 0:
        # The lattice a, b and c axes
        r = units.list_read(10)
        lt.alat = r.real(lt.alat)
        lt.blat = r.real(lt.blat)
        lt.clat = r.real(lt.clat)
        # The lattice internal angles
        r = units.list_read(10)
        alpha = r.real(alpha)
        beta = r.real(beta)
        gam = r.real(gam)
        # Convert these to the internal EPSR18 angles phib,thetac,phic
        lt.phib, lt.thetac, lt.phic = auxiliary_routines.convertalphabetagamma(
            alpha, beta, gam, lt.phib, lt.thetac, lt.phic)
    else:
        units.write(6, fmt(_F602, substr(symbols.fnameato, 1, symbols.ilenfnameato)))
        return
    units.write(6, fmt(_F603, lt.alat, lt.blat, lt.clat, alpha, beta, gam, lt.phib, lt.thetac, lt.phic))
    # If the word 'absolute' is encountered here then absolute Cartesian coordinates are expected
    # If the word 'fractional' is encountered here then fractional unit cell coordinates are expected
    # Default is fractional.
    lt.absolute = False
    cm.sometext = fstr(units.read_line(10), cm.nstsize)
    nunitcell = 0
    if findex(cm.sometext, 'absolute') > 0:
        lt.absolute = True
        # Get the number of atoms in the unit cell
        nunitcell = units.list_read(10).int(nunitcell)
    elif findex(cm.sometext, 'fractional') > 0:
        lt.absolute = False
        # Get the number of atoms in the unit cell
        nunitcell = units.list_read(10).int(nunitcell)
    else:
        nunitcell = read_internal(cm.sometext).int(nunitcell)
    atsymbolunit = char_array(nunitcell)
    xunit = real_array(nunitcell)
    yunit = real_array(nunitcell)
    zunit = real_array(nunitcell)
    # Get the unit cell coordinates
    molecom.ntype = 0
    for i in range(1, nunitcell + 1):
        cm.sometext = fstr(units.read_line(10), cm.nstsize)
        cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
        if cm.nwords > 3:
            atsymbolunit[i] = fstr(substr(cm.sometext, cm.ncf[1], cm.ncl[1]), 3)
            # Store this temporarily to set up the list of atom types so that the list of potential
            # parameters can be checked.
            molecom.alabel[i] = atsymbolunit[i]
            ato_manipulation_rev4.ftype(i)
            r = read_internal(substr(cm.sometext, cm.ncf[2], cm.ncl[4]))
            xunit[i] = r.real(xunit[i])
            yunit[i] = r.real(yunit[i])
            zunit[i] = r.real(zunit[i])
        else:
            units.write(6, fmt(_F604))
            return
    units.write(6, fmt(_F605, molecom.ntype))
    for i in range(1, molecom.ntype + 1):
        units.write(6, fmt('(1x,i3,1x,a3)', i, Ch(molecom.type[i], 3)))
    #
    # read Lennard-Jones, mass and charge parameters for each atom type
    #
    for it in range(1, molecom.ntype + 1):
        cm.sometext, ierr = _read_text(10, cm.sometext, cm.nstsize)
        #
        # parse this text into words
        #
        cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
        a1 = fstr(substr(cm.sometext, cm.ncf[1], cm.ncl[1]), 3)
        if cm.nwords > 1:
            ats = fstr(substr(cm.sometext, cm.ncf[2], cm.ncl[2]), 3)
        else:
            ats = 'XXX'
        if cm.nwords > 2:
            iexch = read_internal(substr(cm.sometext, cm.ncf[3], cm.ncl[3])).int(iexch)
        else:
            iexch = -1
        # read(10,*,iostat=ierr) epsit,sigt,amass,qt,rqt  (items read before an error keep
        # their new values)
        ierr = 0
        try:
            r = units.list_read(10)
            epsit = r.real(epsit)
            sigt = r.real(sigt)
            amass = r.real(amass)
            qt = r.real(qt)
            rqt = r.real(rqt)
        except FortranEOF:
            ierr = -1
        except FortranIOError as e:
            ierr = e.iostat
        # NB: natomsph is not declared in makelattice.f, so it is an implicit INTEGER:
        # natomsph=1.0 stores the integer 1, and ctype's REAL dummy receives its bit
        # pattern, the denormal 1.4e-45 (Rust makelattice.rs does the same).
        natomsph = np.array([1], dtype=np.int32).view(np.float32)[0]
        radminsph = f32(0.0)
        radmaxsph = f32(0.5) * sigt
        pseudofacsph = f32(1.0)
        rhosph = f32(0.0)
        # assign type number
        ato_manipulation_rev4.ctype(a1, ats, iexch, epsit, sigt, amass, qt, rqt, natomsph, rhosph,
                                    radminsph, radmaxsph, pseudofacsph)
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
    if ierr != 0:
        molecom.ecore = f32(0.0)
        molecom.dcore = f32(0.0)
    #
    # initialise random number generator
    #
    auxiliary_routines.initran1(10)
    units.close(10)
    #
    # get the filename to put the mixture in
    #
    cm.nwords = 0
    while cm.nwords <= 0:
        units.write(6, fmt(_F35), advance=False)
        cm.sometext, ierr = _read_text(5, cm.sometext, cm.nstsize)
        cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
        if cm.nwords <= 0:
            units.write(6, fmt(_F499))
    readfile = fstr(substr(cm.sometext, cm.ncf[1], cm.ncl[1]), 256)
    ilenreadfile = int(cm.ncl[1]) - int(cm.ncf[1]) + 1
    #
    # strip out any extension that was supplied - extension must be .ato, otherwise
    # nothing is done.
    #
    ilenuptodot = findex(readfile, '.') - 1
    if ilenuptodot > 0:
        readfile = fstr(substr(readfile, 1, ilenuptodot), 256)
        ilenreadfile = ilenuptodot
    symbols.fnameato = fstr(substr(dn.directory_name, 1, dn.nchardir) + substr(readfile, 1, ilenreadfile)
                            + '.ato', 256)
    symbols.ilenfnameato = dn.nchardir + ilenreadfile + 4
    units.write(6, fmt(_F609, substr(symbols.fnameato, 1, symbols.ilenfnameato)))

    # OK Ready to build lattice with specified number of unit cells
    #
    # define the conversion from lattice vectors to Cartesian coordinates
    # (cosd/sind are gfortran's intrinsics here: makelattice.o calls _gfortran_cosd_r4)
    #
    lt.alatx = lt.alat
    lt.blatx = lt.blat * cosd_intrinsic(lt.phib)
    lt.blaty = lt.blat * sind_intrinsic(lt.phib)
    lt.clatx = lt.clat * sind_intrinsic(lt.thetac) * cosd_intrinsic(lt.phic)
    lt.claty = lt.clat * sind_intrinsic(lt.thetac) * sind_intrinsic(lt.phic)
    lt.clatz = lt.clat * cosd_intrinsic(lt.thetac)
    cc.nmol = 0
    molecom.simtk = f32(300.0)
    molecom.tolout = f32(0.1)
    molecom.stepmi = f32(0.01)
    molecom.stepri = f32(0.0)
    molecom.stepr = f32(0.001)
    molecom.stepm = f32(0.001)
    molecom.tempb = f32(65.0)
    molecom.tempangreq = f32(1.0)
    molecom.tempdihreq = f32(10.0)
    molecom.tempangwt = f32(1.0)
    molecom.tempdihwt = f32(10.0)
    ia = 0
    ig = 0  # noqa: F841
    im = 0
    xav = f32(0.0)
    yav = f32(0.0)
    zav = f32(0.0)
    molecom.natom = 0
    # ngrrot is a local here (makelattice does not use molecom_mod)
    ngrrot = 0  # noqa: F841
    xadd, yadd, zadd = cc.xadd, cc.yadd, cc.zadd
    #
    # set up position of centre of box (roughly)
    #
    acent = f32(0.5) * f32(na)
    bcent = f32(0.5) * f32(nb)
    ccent = f32(0.5) * f32(nc)
    for kc in range(1, nc + 1):
        fc = f32(kc - 1) - ccent
        zlat = fc * lt.clatz
        for kb in range(1, nb + 1):
            fb = f32(kb - 1) - bcent
            ylat = fb * lt.blaty + fc * lt.claty
            for ka in range(1, na + 1):
                fa = f32(ka - 1) - acent
                #
                # Cartesian coordinates of this lattice point
                #
                xlat = fa * lt.alatx + fb * lt.blatx + fc * lt.clatx
                #
                # now within this unit cell set up the atom positions
                #
                for jm in range(1, nunitcell + 1):
                    im = im + 1
                    cc.nato[im] = 1
                    if lt.absolute:
                        # Use absolute Cartesian coordinates
                        xadd[im] = xlat + xunit[jm]
                        yadd[im] = ylat + yunit[jm]
                        zadd[im] = zlat + zunit[jm]
                    else:
                        # Convert fractional coordinates to Cartesian values
                        xadd[im] = xlat + xunit[jm] * lt.alatx + yunit[jm] * lt.blatx + zunit[jm] * lt.clatx
                        yadd[im] = ylat + yunit[jm] * lt.blaty + zunit[jm] * lt.claty
                        zadd[im] = zlat + zunit[jm] * lt.clatz
                    #
                    # form mean coordinate of all the molecules
                    #
                    xav = xav + xadd[im]
                    yav = yav + yadd[im]
                    zav = zav + zadd[im]
                    cc.phix[im] = 0.0
                    cc.phiy[im] = 0.0
                    cc.phiz[im] = 0.0
                    cc.tethered[im] = True
                    # for each atom set its label
                    cc.naf[im] = ia + 1
                    for ib in range(1, int(cc.nato[im]) + 1):
                        ia = ia + 1
                        molecom.alabel[ia] = atsymbolunit[jm]
                        molecom.itype[ia] = molecom.itype[jm]
                        molecom.nbond[ia] = 0
                        molecom.xa[ia] = 0.0
                        molecom.ya[ia] = 0.0
                        molecom.za[ia] = 0.0
                    # number of internal rotation groups
                    cc.ngrot[im] = 0
                    cc.nal[im] = ia
    molecom.natom = ia
    cc.nmol = im
    #
    # coordinates of centre of box
    #
    xbox = xav / f32(cc.nmol)
    ybox = yav / f32(cc.nmol)
    zbox = zav / f32(cc.nmol)
    #
    # shift molecular coordinates so that centre of all molecules is at centre of
    # box
    #
    for im in range(1, cc.nmol + 1):
        xadd[im] = xadd[im] - xbox
        yadd[im] = yadd[im] - ybox
        zadd[im] = zadd[im] - zbox
        cc.phix[im] = xadd[im]
        cc.phiy[im] = yadd[im]
        cc.phiz[im] = zadd[im]
    units.write(6, fmt(_F610, cc.nmol, molecom.natom))
    units.write(6, fmt(_F611, xbox, ybox, zbox))
    #
    # make a list of the molecule types read in from the current file
    #
    cc.nmoltype = 0
    im = 0
    while im < cc.nmol:
        im = im + 1
        imoltype = 0
        ifound = 0
        while imoltype < cc.nmoltype and ifound == 0:
            imoltype = imoltype + 1
            if molecom.alabel[cc.naf[im]] == cc.moltype[imoltype]:
                ifound = imoltype
                cc.mtype[im] = ifound
        if ifound == 0:
            cc.nmoltype = cc.nmoltype + 1
            cc.moltype[cc.nmoltype] = molecom.alabel[cc.naf[im]]
            cc.moltethered[cc.nmoltype] = cc.tethered[im]
            #
            # initialise the template filename
            #
            cc.moltemplate[cc.nmoltype] = fstr(fmt('("moltype",i2.2)', cc.nmoltype), 256)
            cc.mtype[im] = cc.nmoltype

            # Initialise the rotational and translational step sizes

            cc.wholerot[cc.nmoltype] = molecom.stepr
            cc.wholetrans[cc.nmoltype] = molecom.stepm

    # Write these out as a check

    units.write(6, fmt(_F201))
    for imoltype in range(1, cc.nmoltype + 1):
        units.write(6, fmt(_F202, imoltype, Ch(cc.moltype[imoltype], 3), substr(cc.moltemplate[imoltype], 1, 10),
                           cc.wholerot[imoltype], cc.wholetrans[imoltype]))
    lt.alat = lt.alat * f32(na)
    lt.blat = lt.blat * f32(nb)
    lt.clat = lt.clat * f32(nc)
    # Write out ato file
    ato_read_and_save.moloutc(symbols.fnameato)

    return
