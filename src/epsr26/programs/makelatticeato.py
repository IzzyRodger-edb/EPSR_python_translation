"""Translation of fortran_src/makelatticeato.f: ``program standalone`` and ``subroutine
makelatticeato``.

Makes an EPSR .ato file from an existing .ato file by treating its box as a unit cell
and creating na x nb x nc copies of it along the a, b, c axes; optionally ("convert")
every atom then becomes a separate tethered one-atom molecule.
See docs/translation_notes/N2.md.
"""

from ..common import clustercom, commands, lattice, molecom, molecom_mod, symbols
from ..fortran import (Ch, FortranEOF, FortranIOError, PI4ATAN, f32, findex, fmt, fstr, lw, r4, read_internal,
                       real_array, substr, units)
from ..common.epsr_dimension import mmol
from ..library import ato_read_and_save, file_and_directory_commands, parse, reallocation_routines
from ..library import molecom_mod as molecom_mod_lib
from ..trace import routine

_F700 = "(/'makelatticeato> No .ato files found in the current directory - try another directory.'/)"
_F499 = "(/'makelaticeato> No file specified! - try again!'/)"
_F606 = "(/'makelattice> Give number of unit cells along c-axis: ')"
_F607 = "(/'makelattice> Give number of unit cells along b and c-axes: ')"
_F608 = "(/'makelattice> Give number of unit cells along a, b and c-axes: ')"
# 35 format(/'makelatticeato> Type name of file to put lattice in: ',$)
_F35 = "(/'makelatticeato> Type name of file to put lattice in: ')"
_F101 = "(a)"
_F613 = "('makelatticeato> Convert molecules to atoms? (y or n) ')"
_F599 = "(/'Problems with .ato file ',a/'Routine will exit now'/)"
# 603 format('makelatticeato> Lattice constants, angles and derived angles:'
#     13(/3(1x,f10.5)))       (the 1 is the continuation mark)
_F603 = "('makelatticeato> Lattice constants, angles and derived angles:',3(/3(1x,f10.5)))"
_F605 = "('makelatticeato> Number of atom types found: ',i3)"
_F609 = "('makelatticeato> Output being written to file: ',a)"
_F610 = "('makelatticeato> ',i8,' molecules and ',i8,' atoms created')"
_F611 = "('makelatticeato> Centre of box at: ',3(1x,e12.5))"
_F612 = "('makelatticeato> Value of convert = ',L)"


def main():
    """makelatticeato.f:standalone -- 2009/3/30 Attempt to make this a standalone program."""
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    makelatticeato()


@routine
def makelatticeato():
    """makelatticeato.f:makelatticeato -- Make an EPSR18 ato file from an exising ato file
    by creating the specified number of unit cells along the a,b,c axes."""
    # use reallocation_routines
    # use molecom_mod
    # character(len=256)  :: readfile
    # integer             :: ilenreadfile,ierr,i,nwordstest
    # real                :: pi,radconv,degconv
    # integer             :: nunitcell,na,nb,nc,ka,kb,kc,jb
    # character(len=4)    :: ext
    # character(len=1)    :: iaccept
    # real                :: acent,bcent,ccent,fa,fb,fc,xav,yav,zav,xbox,ybox,zbox
    # real  xadds(mmol),yadds(mmol),zadds(mmol),phixs(mmol),phiys(mmol),phizs(mmol)
    # character(len=256) :: text
    # logical :: convert
    cm = commands
    cc = clustercom
    lt = lattice
    mm = molecom_mod
    xadds = real_array(mmol)
    yadds = real_array(mmol)
    zadds = real_array(mmol)
    phixs = real_array(mmol)
    phiys = real_array(mmol)
    phizs = real_array(mmol)
    # NB: text is never assigned in the Fortran (parse is called on it, see below).  The
    # stack memory is zero, i.e. NUL characters, which are not blanks; parse is told the
    # length nstsize=3020 and so reads past the 256 characters of text into more zero
    # memory.  The reference output shows exactly that: one word, ncf(1)=1, ncl(1)=3020.
    # 3020 NULs reproduce it (the Rust port assumed blanks, which gives the same convert
    # answer in the test cases but not the same stdout).
    text = '\0' * commands.nstsize
    iaccept = ''
    na = nb = nc = 0
    # NB: nwordstest is set only when no file name is given; it is used only in the dead
    # store iread=nwordstest+3
    nwordstest = 0
    iread = 0  # noqa: F841
    pi = PI4ATAN
    radconv = pi / r4(180.0)  # noqa: F841
    degconv = r4(180.0) / pi  # noqa: F841
    #
    # generate the name of the file containing the unit cell coordinates
    #
    if cm.nwords < 2:
        cm.nwords = 0
        nwordstest = 1

        # generate the name of each .ato file in turn

        ext = '.ato'
        ilenext = 4
        while cm.nwords <= 0:
            iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)

            # if no .ato files are available in this directory, abandon this attempt at dockato

            if iaccept != 'y':
                units.write(6, fmt(_F700))
                return
            cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
            if cm.nwords <= 0:
                units.write(6, fmt(_F499))
        readfile = fstr(substr(cm.sometext, cm.ncf[1], cm.ncl[cm.nwords]), 256)
        ilenreadfile = int(cm.ncl[cm.nwords]) - int(cm.ncf[1]) + 1
    else:
        readfile = fstr(substr(cm.sometext, cm.ncf[2], cm.ncl[2]), 256)
        ilenreadfile = int(cm.ncl[2]) - int(cm.ncf[2]) + 1
    if cm.nwords > 2:
        na = read_internal(substr(cm.sometext, cm.ncf[3], cm.ncl[3])).int(na)
        if cm.nwords > 3:
            nb = read_internal(substr(cm.sometext, cm.ncf[4], cm.ncl[4])).int(nb)
            if cm.nwords > 4:
                nc = read_internal(substr(cm.sometext, cm.ncf[5], cm.ncl[5])).int(nc)
                iread = nwordstest + 3  # noqa: F841
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
    # create the full name of the input file by adding on the directory name if not present
    #
    symbols.fnameato, symbols.ilenfnameato = file_and_directory_commands.get_full_filename(
        readfile, ilenreadfile, symbols.fnameato, symbols.ilenfnameato)
    if findex(substr(symbols.fnameato, 1, symbols.ilenfnameato), '.ato') == 0:
        symbols.fnameato = fstr(substr(symbols.fnameato, 1, symbols.ilenfnameato) + '.ato', 256)
        symbols.ilenfnameato = symbols.ilenfnameato + 4
    #
    # get the filename to put the mixture in
    #
    if cm.nwords > 5:
        readfile = fstr(substr(cm.sometext, cm.ncf[6], cm.ncl[6]), 256)
        ilenreadfile = int(cm.ncl[6]) - int(cm.ncf[6]) + 1
    else:
        cm.nwords = 0
        iread = 0  # noqa: F841
        while cm.nwords <= 0:
            units.write(6, fmt(_F35), advance=False)
            try:
                cm.sometext = fstr(units.read_line(5), cm.nstsize)
            except (FortranEOF, FortranIOError):
                pass
            cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
            if cm.nwords <= 0:
                units.write(6, fmt(_F499))
        readfile = fstr(substr(cm.sometext, cm.ncf[1], cm.ncl[cm.nwords]), 256)
        ilenreadfile = int(cm.ncl[cm.nwords]) - int(cm.ncf[1]) + 1
    # Decide if to convert atoms to molecules
    if cm.nwords > 6:
        cm.sometext = fstr(substr(cm.sometext, cm.ncf[7], cm.ncl[7]), cm.nstsize)
    else:
        units.write(6, fmt(_F613), advance=False)
        cm.sometext = fstr(units.read_line(5), cm.nstsize)
    #      write(sometext,'(a)') 'y'
    # NB: the Fortran parses the uninitialised variable text, not sometext: nwords becomes 1
    # with ncf(1)=1, ncl(1)=nstsize (see above), so the whole of sometext is compared with
    # 'y' / 'Y' and written out (3020 characters)
    cm.nwords = parse.parse(text, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
    units.write(6, lw(cm.nwords, ' ', substr(cm.sometext, cm.ncf[1], cm.ncl[1])))
    convert = (substr(cm.sometext, cm.ncf[1], cm.ncl[1]).rstrip() == 'y'
               or substr(cm.sometext, cm.ncf[1], cm.ncl[1]).rstrip() == 'Y')
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
        units.write(6, fmt(_F599, substr(symbols.fnameato, 1, symbols.ilenfnameato)))
        return
    units.write(6, lw(cc.dmax, molecom.dmaxmol))
    units.write(6, fmt(_F603, lt.alat, lt.blat, lt.clat, lt.phib, lt.thetac, lt.phic))

    # Save the coordinates of these molecules

    for im in range(1, cc.nmol + 1):
        if cc.tethered[im]:  # For tethered we set coordinates to tethering positions to avoid problems with periodic boundary conditions
            xadds[im] = cc.phix[im]
            yadds[im] = cc.phiy[im]
            zadds[im] = cc.phiz[im]
        else:
            xadds[im] = cc.xadd[im]
            yadds[im] = cc.yadd[im]
            zadds[im] = cc.zadd[im]
        phixs[im] = cc.phix[im]
        phiys[im] = cc.phiy[im]
        phizs[im] = cc.phiz[im]
    # Get the number of molecules in the unit cell
    nunitcell = cc.nmol
    #      allocate(atsymbolunit(nunitcell),xunit(nunitcell),yunit(nunitcell),zunit(nunitcell))

    units.write(6, fmt(_F605, molecom.ntype))
    for i in range(1, molecom.ntype + 1):
        units.write(6, fmt('(1x,i3,1x,a3)', i, Ch(molecom.type[i], 3)))

    # Get full name of output file

    symbols.fnameato, symbols.ilenfnameato = file_and_directory_commands.get_full_filename(
        readfile, ilenreadfile, symbols.fnameato, symbols.ilenfnameato)
    if findex(substr(symbols.fnameato, 1, symbols.ilenfnameato), '.ato') == 0:
        symbols.fnameato = fstr(substr(symbols.fnameato, 1, symbols.ilenfnameato) + '.ato', 256)
        symbols.ilenfnameato = symbols.ilenfnameato + 4
    units.write(6, fmt(_F609, substr(symbols.fnameato, 1, symbols.ilenfnameato)))

    # OK Ready to build lattice with specified number of unit cells
    ia = 0
    ig = 0
    mm.ngrrot = 0
    #      ngrrotmax=0
    mm.nrottot = 0
    #      nrottotmax=0
    im = 0
    xav = f32(0.0)
    yav = f32(0.0)
    zav = f32(0.0)
    nadd = 0
    xa, ya, za = molecom.xa, molecom.ya, molecom.za
    alabel = molecom.alabel
    imageatom = molecom.imageatom
    nbond = molecom.nbond
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
                ja = 0
                irgref = 0
                for jm in range(1, nunitcell + 1):
                    im = im + 1
                    cc.nato[im] = cc.nato[jm]
                    cc.xadd[im] = xlat + xadds[jm]
                    cc.yadd[im] = ylat + yadds[jm]
                    cc.zadd[im] = zlat + zadds[jm]
                    cc.tethered[im] = cc.tethered[jm]
                    cc.mtype[im] = cc.mtype[jm]
                    if cc.tethered[im]:
                        cc.phix[im] = xlat + phixs[jm]
                        cc.phiy[im] = ylat + phiys[jm]
                        cc.phiz[im] = zlat + phizs[jm]
                    else:
                        cc.phix[im] = phixs[jm]
                        cc.phiy[im] = phiys[jm]
                        cc.phiz[im] = phizs[jm]
                    #
                    # form mean coordinate of all the molecules
                    #
                    xav = xav + cc.xadd[im]
                    yav = yav + cc.yadd[im]
                    zav = zav + cc.zadd[im]
                    # for each atom set its label
                    cc.naf[im] = ia + 1
                    for jb in range(1, int(cc.nato[jm]) + 1):
                        ja = ja + 1
                        ia = ia + 1
                        alabel[ia] = alabel[ja]
                        units.write(6, lw(im, jm, ia, ja, imageatom[ja]))
                        if imageatom[ja] > 0:
                            imageatom[ia] = imageatom[ja] + nadd
                        else:
                            imageatom[ia] = 0
                        nbond[ia] = nbond[ja]
                        xa[ia] = xa[ja]
                        ya[ia] = ya[ja]
                        za[ia] = za[ja]
                        molecom.itype[ia] = molecom.itype[ja]
                        molecom.mref[ia] = im
                        if nbond[ia] > 0:
                            for ib in range(1, int(nbond[ia]) + 1):
                                molecom.iabond[ia, ib] = molecom.iabond[ja, ib] + nadd
                                molecom.tbond[ia, ib] = molecom.tbond[ja, ib]
                    # number of internal rotation groups
                    cc.ngrot[im] = cc.ngrot[jm]
                    if cc.ngrot[im] > 0:
                        mm.ngrrot = mm.ngrrot + int(cc.ngrot[im])
                        molecom_mod_lib.testsize_rotgrps()
                        for ih in range(1, int(cc.ngrot[im]) + 1):
                            ig = ig + 1
                            irgref = irgref + 1
                            # save molecule in which this rotational group occurs
                            mm.mrot[ig] = im
                            mm.rotgrp[ig] = mm.rotgrp[irgref]
                            # for rotational groups read atom numbers of two atoms to form axis of rotation
                            for j in range(1, 2 + 1):
                                mm.labaxi[ig, j] = mm.labaxi[irgref, j] + nadd
                            # read number of atoms in group and their numbers
                            mm.nrot[ig] = mm.nrot[irgref]
                            mm.rotf[ig] = mm.nrottot + 1
                            mm.nrottot = mm.nrottot + int(mm.nrot[ig])
                            mm.rotl[ig] = mm.nrottot
                            molecom_mod_lib.testsize_labrot()
                            labrefold = int(mm.rotf[irgref]) - 1
                            labrefnew = int(mm.rotf[ig]) - 1
                            for j in range(1, int(mm.nrot[ig]) + 1):
                                labrefold = labrefold + 1
                                labrefnew = labrefnew + 1
                                mm.labrot[labrefnew] = mm.labrot[labrefold] + nadd
                            # Do any dihedral groups.
                            if mm.rotgrp[ig] == 'DIH':
                                molecom.ndihedralang[ig] = molecom.ndihedralang[irgref]
                                for k in range(1, int(molecom.ndihedralang[ig]) + 1):
                                    molecom.dihedralang[k, ig] = molecom.dihedralang[k, irgref]
                                    for j in range(1, 2 + 1):
                                        molecom.labrotdih[j, k, ig] = molecom.labrotdih[j, k, irgref] + nadd
                    cc.nal[im] = ia
                nadd = ia
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
        cc.xadd[im] = cc.xadd[im] - xbox
        cc.yadd[im] = cc.yadd[im] - ybox
        cc.zadd[im] = cc.zadd[im] - zbox
        if cc.tethered[im]:
            cc.phix[im] = cc.phix[im] - xbox
            cc.phiy[im] = cc.phiy[im] - ybox
            cc.phiz[im] = cc.phiz[im] - zbox
    molecom.natom = ia
    units.write(6, fmt(_F610, cc.nmol, molecom.natom))
    units.write(6, fmt(_F611, xbox, ybox, zbox))
    lt.alat = lt.alat * f32(na)
    lt.blat = lt.blat * f32(nb)
    lt.clat = lt.clat * f32(nc)
    # Retain the geometry if possible
    lt.cubic = lt.cubic and (na == nb and na == nc and nb == nc)
    if lt.cubic:
        cc.dmax = cc.dmax * f32(na)
    units.write(6, fmt(_F612, convert))
    if convert:
        # First make all the atom coordinates absolute
        for im in range(1, cc.nmol + 1):
            for ia in range(int(cc.naf[im]), int(cc.nal[im]) + 1):
                xa[ia] = xa[ia] + cc.xadd[im]
                ya[ia] = ya[ia] + cc.yadd[im]
                za[ia] = za[ia] + cc.zadd[im]
        # Set the number of molecules equal to number of atoms
        cc.nmol = molecom.natom
        # Define the coordinates of each molecule and tether it
        for im in range(1, cc.nmol + 1):
            ia = im
            cc.xadd[im] = xa[ia]
            cc.phix[im] = xa[ia]
            cc.yadd[im] = ya[ia]
            cc.phiy[im] = ya[ia]
            cc.zadd[im] = za[ia]
            cc.phiz[im] = za[ia]
            cc.tethered[im] = True
            cc.tetheredatom[im] = 0
            cc.nstepreject[im] = 0
            cc.nato[im] = 1
            cc.naf[im] = ia
            cc.nal[im] = ia
            nbond[im] = 0
            cc.ngrot[im] = 0
            xa[ia] = 0.0
            ya[ia] = 0.0
            za[ia] = 0.0
            imageatom[ia] = 0
        # Set the tolerance low so that atoms cannot move far from their starting positions
        molecom.tolout = f32(0.1)
        cc.nmoltype = 0
        #
        # make a list of the molecule types read in from the current file
        #
        nmoltypefirst = cc.nmoltype  # noqa: F841
        im = 0
        while im < cc.nmol:
            im = im + 1
            imoltype = 0
            ifound = 0
            while imoltype < cc.nmoltype and ifound == 0:
                imoltype = imoltype + 1
                if alabel[cc.naf[im]] == cc.moltype[imoltype]:
                    ifound = imoltype
                    cc.mtype[im] = ifound
            if ifound == 0:
                cc.nmoltype = cc.nmoltype + 1
                cc.moltype[cc.nmoltype] = alabel[cc.naf[im]]
                cc.moltethered[cc.nmoltype] = cc.tethered[im]
                cc.moltetheredatom[cc.nmoltype] = cc.tetheredatom[im]
                #
                # initialise the template filename
                #
                cc.moltemplate[cc.nmoltype] = fstr(fmt('("moltype",i2.2)', cc.nmoltype), 256)
                cc.mtype[im] = cc.nmoltype

                # Initialise the rotational and translational step sizes

                cc.wholerot[cc.nmoltype] = molecom.stepr
                cc.wholetrans[cc.nmoltype] = molecom.stepm

        # Initialise the relative broadening of each molecule type if coarse graining
        mm.rfractmol = reallocation_routines.reallocate1d_r(mm.rfractmol, cc.nmoltype)
        mm.rfractmol[:] = 0.0

    # Write out ato file
    ato_read_and_save.moloutc(symbols.fnameato)

    return
