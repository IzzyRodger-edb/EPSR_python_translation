"""Translation of fortran_src/writexyz.f: ``program standalone``, ``subroutine
run_writexyz`` and ``subroutine writexyz``.

Writes (part of) a configuration as an .xyz block: optionally centred on a molecule,
restricted to a slab, rotated by three Euler angles and with atom types excluded; the
block is appended to ``<file>append.xyz``, or written to ``<file>NNNNN.xyz`` with the
block counter kept in ``<file>.nblock``.  See docs/translation_notes/N2.md.
"""

from ..common import clustercom, commands, directory_name, lattice, molecom, plot_pars, plotato, symbols
from ..common.epsr_dimension import mcomp
from ..fortran import (Ch, FortranEOF, FortranIOError, cosd_intrinsic, f32, findex, fmt, fstr, int_array, lw,
                       read_internal, sind_intrinsic, substr, units)
from ..library import ato_manipulation_rev4, ato_read_and_save, file_and_directory_commands, minimagec, parse
from ..trace import routine

_F499 = "(/'writexyz> No .ato file specified! - routine will exit'/)"
_F599 = "(/'writexyz> Problems with .ato file ',a/'Routine will exit now'/)"
# 115 format(/'writexyz> Specify the chemical symbol for atom type ',a3,': ',$)
_F115 = "(/'writexyz> Specify the chemical symbol for atom type ',a3,': ')"
_F121 = "(/'writexyz> Do you want to append to an existing file? (y or n) ')"
_F122 = "(a)"
_F132 = "(/'plotato> Specify number of molecule at centre (0 to plot entire box): ')"
_F133 = ("(/'plotato> Give maximum distances along x and y axes,'"
         "/,'minimum and maximum distances along the z-axis'"
         "/' - use zero distances to show only the central molecule: ')")
_F134 = "(/'plotato> Give phi, theta and chi for rotation of coordinates (0 0 0) ')"
_F111 = "(/'writexyz> ',i3,' components available for plotting:-',/(1x,i3,1x,a3))"
_F131 = "(/'writexyz> Give component numbers to EXCLUDE from the plot - ',' type 0 to include all: ')"
_F101 = "(a)"
_F151 = "(/'writexyz> ',i3,' components to be printed:-',/(1x,i3,1x,a3))"
_F100 = "(a)"
_F113 = "(a3,4x,e10.3,2x,e10.3,2x,e10.3,2x,a3)"


def main():
    """writexyz.f:standalone -- 2009/3/30 Attempt to make this a standalone program."""
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_writexyz()


@routine
def run_writexyz():
    """writexyz.f:run_writexyz -- To generate a 3D plot of a .ato file.

    2006-11-15 Output options are (1) GNUplot, (2) PGPlot, (3) ATOMS.
    Also atoms to be plotted will be defined by their types, not atomic symbol.
    """
    # character*256 readfile
    # character*4 ext
    # character*1 iaccept
    # character*3020 sometextsave
    # character*4 yesorno
    # logical appendset
    # integer iread,nwordssave,ncfsave(nworddim),nclsave(nworddim)
    cm = commands
    cc = clustercom
    pa = plotato
    ncfsave = int_array(cm.nworddim)
    nclsave = int_array(cm.nworddim)
    iaccept = ''
    xaxismax = yaxismax = zaxismin = zaxismax = f32(0.0)
    zphi = ztheta = zchi = f32(0.0)
    itest = 0

    # First identify the number of arguments read in

    cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)

    # If number of words is greater than 2 then we need to store the remaining arguments in a text stream

    pa.appendxyz = True
    appendset = False

    # Expect to see the append type as word 3. If not found assume xyz file is to be appended

    if cm.nwords > 2:

        pa.appendxyz = not (substr(cm.sometext, cm.ncf[3], cm.ncl[3]).rstrip() == 'noappend')
        appendset = True

    if cm.nwords > 3:

        sometextsave = fstr(substr(cm.sometext, cm.ncf[4], cm.ncl[cm.nwords]), 3020)
        nwordssave = parse.parse(sometextsave, 0, ncfsave, nclsave, 5, cm.nstsize, cm.nworddim)
        units.write(6, lw(nwordssave, ' ', substr(sometextsave, ncfsave[1], nclsave[nwordssave])))
        iread = 1

    else:

        sometextsave = ''
        nwordssave = 0
        units.write(6, lw(nwordssave))
        iread = 1

    # generate the name of the .ato file

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
    #     call normalise_coordinates
    #
    # step through the atom types and make sure the chemical symbol has been specified
    #
    for it in range(1, molecom.ntype + 1):
        if molecom.atsymbol[it] == 'XXX' or molecom.atsymbol[it] == '':
            units.write(6, fmt(_F115, Ch(molecom.type[it], 3)), advance=False)
            molecom.atsymbol[it] = fstr(units.list_read(5).word(molecom.atsymbol[it]), 3)

    if not appendset:

        units.write(6, fmt(_F121), advance=False)
        yesorno = fstr(units.read_line(5), 4)
        pa.appendxyz = not (substr(yesorno, 1, 1) == 'n' or substr(yesorno, 1, 1) == 'N')
        appendset = True

    #
    # sort out which atoms will actually be plotted, based on whether all the atoms
    # are to be plotted, a group of them centred on one molecule, or just one molecule
    #
    if nwordssave < iread:
        units.write(6, fmt(_F132), advance=False)
        pa.nmolcentre = units.list_read(5).int(pa.nmolcentre)
    else:
        pa.nmolcentre = read_internal(substr(sometextsave, ncfsave[iread], nclsave[iread])).int(pa.nmolcentre)
    iread = iread + 1

    if pa.nmolcentre > 0:
        plot_pars.plottype = 2
    else:
        plot_pars.plottype = 1

    #
    # sort out which atoms will actually be plotted, based on whether all the atoms
    # are to be plotted, a group of them centred on one molecule, or just one molecule
    #
    #      if(nmolcentre.gt.0) then
    if nwordssave < iread + 3:
        units.write(6, fmt(_F133), advance=False)
        r = units.list_read(5)
    else:
        r = read_internal(substr(sometextsave, ncfsave[iread], nclsave[iread + 3]))
    xaxismax = r.real(xaxismax)
    yaxismax = r.real(yaxismax)
    zaxismin = r.real(zaxismin)
    zaxismax = r.real(zaxismax)
    iread = iread + 4

    if nwordssave < iread + 2:
        units.write(6, fmt(_F134), advance=False)
        r = units.list_read(5)
    else:
        r = read_internal(substr(sometextsave, ncfsave[iread], nclsave[iread + 2]))
    zphi = r.real(zphi)
    ztheta = r.real(ztheta)
    zchi = r.real(zchi)
    iread = iread + 3

    #      else
    #         xaxismax=dmax
    #         yaxismax=dmax
    #         zaxismin=-dmax
    #         zaxismax=dmax
    #         zphi=0.0
    #         ztheta=0.0
    #         zchi=0.0
    #      endif

    # Get sines and cosines
    # (gfortran's intrinsics: writexyz.o calls _gfortran_cosd_r4/_gfortran_sind_r4)

    czphi = cosd_intrinsic(zphi)
    szphi = sind_intrinsic(-zphi)
    cztheta = cosd_intrinsic(ztheta)
    sztheta = sind_intrinsic(-ztheta)
    czchi = cosd_intrinsic(zchi)
    szchi = sind_intrinsic(-zchi)

    # Ensure the coordinates of each atom in each molecule are relative to the centre of mass

    for im in range(1, cc.nmol + 1):

        ato_manipulation_rev4.centreofmass(im)

    #
    # step through molecules and determine which ones are to be plotted
    #
    iaplot = 0
    implot = 0

    # Ensure coordinates plotted are relative to molecule at centre

    if pa.nmolcentre > 0 and pa.nmolcentre <= cc.nmol:
        xaddcentre = cc.xadd[pa.nmolcentre]
        yaddcentre = cc.yadd[pa.nmolcentre]
        zaddcentre = cc.zadd[pa.nmolcentre]
    else:
        xaddcentre = f32(0.0)
        yaddcentre = f32(0.0)
        zaddcentre = f32(0.0)
    xa, ya, za = molecom.xa, molecom.ya, molecom.za
    xat, yat, zat = pa.xat, pa.yat, pa.zat
    iplot = pa.iplot
    for im in range(1, cc.nmol + 1):

        # Apply minimum image to centre of mass coordinates.

        xdif = cc.xadd[im] - xaddcentre
        ydif = cc.yadd[im] - yaddcentre
        zdif = cc.zadd[im] - zaddcentre
        if not lattice.cubic:

            xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)

        else:

            xdif = xdif / cc.dmax
            ydif = ydif / cc.dmax
            zdif = zdif / cc.dmax
            xdif, ydif, zdif = minimagec.minimage(xdif, ydif, zdif)
            xdif = xdif * cc.dmax
            ydif = ydif * cc.dmax
            zdif = zdif * cc.dmax

        #
        # Rotate these coordinates to the new axis
        #
        xdif, ydif, zdif = ato_manipulation_rev4.rotatecoords(czphi, szphi, cztheta, sztheta, czchi, szchi,
                                                              xdif, ydif, zdif)
        #
        # if one molecule is to be at centre, check to see if the current molecule
        # is within the correct distance
        #
        ifound = 0
        #         if(plottype.gt.1) then
        #            dtest=sqrt(xdif*xdif+ydif*ydif+zdif*zdif)
        #            if(dtest.le.radmax.or.im.eq.nmolcentre) then
        if ((abs(xdif) <= xaxismax and abs(ydif) <= yaxismax
             and zdif >= zaxismin and zdif <= zaxismax) or im == pa.nmolcentre):
            ifound = im
        #         else
        #            ifound=im
        #         endif
        if ifound > 0:
            implot = implot + 1
            for ia in range(int(cc.naf[im]), int(cc.nal[im]) + 1):
                iaplot = iaplot + 1
                iplot[ia] = 1
        else:
            for ia in range(int(cc.naf[im]), int(cc.nal[im]) + 1):
                iplot[ia] = 0

        # Rotate each molecule onto the new coordinate system

        for ia in range(int(cc.naf[im]), int(cc.nal[im]) + 1):

            xat[ia] = xa[ia]
            yat[ia] = ya[ia]
            zat[ia] = za[ia]
            xat[ia], yat[ia], zat[ia] = ato_manipulation_rev4.rotatecoords(czphi, szphi, cztheta, sztheta, czchi,
                                                                           szchi, xat[ia], yat[ia], zat[ia])

        #
        # generate the coordinates of each atom relative to the centre of the plot
        #
        for ia in range(int(cc.naf[im]), int(cc.nal[im]) + 1):

            xat[ia] = xat[ia] + xdif
            yat[ia] = yat[ia] + ydif
            zat[ia] = zat[ia] + zdif
    units.write(6, lw(implot, iaplot))

    # Make a list of atom types for those atoms that could be plotted.

    pa.ntypplot = 0
    for i in range(1, molecom.natom + 1):
        it = molecom.itype[i]
        if iplot[i] > 0:
            iflagt = 0
            j = 0
            while j < pa.ntypplot and iflagt == 0:
                j = j + 1
                if pa.itypplot[j] == it:
                    iflagt = j
            if iflagt == 0 and pa.ntypplot < mcomp:
                pa.ntypplot = pa.ntypplot + 1
                pa.itypplot[pa.ntypplot] = it
            else:
                pass

    if nwordssave < iread:
        units.write(6, fmt(_F111, pa.ntypplot, *[v for j in range(1, pa.ntypplot + 1)
                                                 for v in (pa.itypplot[j], Ch(molecom.type[pa.itypplot[j]], 3))]))
        units.write(6, fmt(_F131), advance=False)
        cm.sometext = fstr(units.read_line(5), cm.nstsize)
    else:
        # read(sometextsave(ncfsave(iread):nclsave(nwordssave)),101) sometext
        cm.sometext = fstr(substr(sometextsave, ncfsave[iread], nclsave[nwordssave]), cm.nstsize)
    cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
    #      write(6,*) iread,nwordssave,sometext(ncf(1):ncl(nwords))

    # Assume at the outset all atom types are to be plotted

    pa.nclplt = pa.ntypplot
    nplt = pa.nplt
    for i in range(1, pa.ntypplot + 1):

        nplt[i] = pa.itypplot[i]

    iw = 0
    ierr = 0
    while iw < cm.nwords and ierr == 0:

        iw = iw + 1
        try:
            itest = read_internal(substr(cm.sometext, cm.ncf[iw], cm.ncl[iw])).int(itest)
            ierr = 0
        except FortranEOF:
            ierr = -1
        except FortranIOError as e:
            ierr = e.iostat
        if ierr == 0:
            i = 0
            while i < pa.nclplt:
                i = i + 1
                if nplt[i] == itest:
                    j = i - 1
                    while j < pa.nclplt - 1:
                        j = j + 1
                        nplt[j] = nplt[j + 1]
                    pa.nclplt = j
    units.write(6, fmt(_F151, pa.nclplt, *[v for j in range(1, pa.nclplt + 1)
                                           for v in (nplt[j], Ch(molecom.type[nplt[j]], 3))]))

    # Based on this list, revise the values of iplot so that only the specified components
    # are plotted.

    for ia in range(1, molecom.natom + 1):
        it = molecom.itype[ia]
        if iplot[ia] > 0:
            iflagt = 0
            j = 0
            while j < pa.nclplt and iflagt == 0:
                j = j + 1
                if nplt[j] == it:
                    iflagt = j

            # If the type of this atom does not correspond to one of the specified types, we signal
            # we do not want to plot this atom

            if iflagt == 0 and molecom.mref[ia] != pa.nmolcentre:
                iplot[ia] = 0
        else:
            iplot[ia] = 0

    writexyz(readfile, ilenreadfile)

    return


@routine
def writexyz(readfile, ilenreadfile):
    """writexyz.f:writexyz -- write the selected atoms as one .xyz block (appended to
    <readfile>append.xyz, or in <readfile>NNNNN.xyz with the counter in <readfile>.nblock)."""
    # character*256 readfile
    # character*5 blocknumber
    # integer ilenreadfile,nblocks,ir
    # logical centralmolecule
    pa = plotato
    cc = clustercom
    dn = directory_name
    # NB: natplot is not initialised in the Fortran; it is only compared (ir.ne.natplot+2)
    # when no count was read, i.e. for an empty file, where any value but -2 rewinds; 0 is used.
    natplot = 0
    # NB: ntag is not initialised in the Fortran; it is only incremented, never used
    ntag = 0
    ir = 0
    if pa.appendxyz:

        #
        # generate name of appended xyz file.
        #
        pa.fnamesplot[0] = fstr(substr(dn.directory_name, 1, dn.nchardir) + substr(readfile, 1, ilenreadfile)
                                + 'append.xyz', 256)
        pa.ilenfnamesplot[0] = dn.nchardir + ilenreadfile + 4

        # Write xyz output file

        try:
            units.open(10, pa.fnamesplot[0], status='unknown')
            ierr = 0
        except FortranIOError as e:
            ierr = e.iostat
        if ierr == 0:

            # Read to the end of file

            nblocks = 0
            while ierr == 0:
                try:
                    natplot = units.list_read(10).int(natplot)
                    ierr = 0
                except FortranEOF:
                    ierr = -1
                except FortranIOError as e:
                    ierr = e.iostat
                #      write(6,*) natplot,nblocks,ierr
                if ierr == 0:
                    ir = 0
                while ierr == 0 and ir < natplot + 2:
                    ir = ir + 1
                    try:
                        units.skip_record(10)
                        ierr = 0
                    except FortranEOF:
                        ierr = -1
                    except FortranIOError as e:
                        ierr = e.iostat
                    #      write(6,*) ir,ierr
                if ierr == 0:
                    nblocks = nblocks + 1

            # Backspace the file to back before the end of file

            units.backspace(10)

            # If an error occurred before the correct number of records were read in, rewind the file to the beginning and start again

            #            write(6,*) ir,natplot+2,nblocks
            if ir != natplot + 2:
                units.rewind(10)
                nblocks = 0
            nblocks = nblocks + 1

            # First determine how many atoms to plot

            natplot = 0
            for ia in range(1, molecom.natom + 1):
                if pa.iplot[ia] > 0:
                    natplot = natplot + 1
            units.write(10, lw(natplot))
            #            write(6,*) natplot,nblocks,ierr
            blocknumber = fstr(fmt('(i5.5)', nblocks), 5)

            # Output some information text

            units.write(10, fmt(_F100, '#Filename: ' + substr(dn.directory_name, 1, dn.nchardir)
                                + substr(readfile, 1, ilenreadfile) + '.ato Block ' + blocknumber))
            for ia in range(1, molecom.natom + 1):
                j = 0
                iflagt = 0
                # Check if a central molecule has been defined.
                centralmolecule = pa.nmolcentre > 0 and pa.nmolcentre <= cc.nmol
                # Check if the current atom is in the central molecule
                if centralmolecule:
                    centralmolecule = ia >= cc.naf[pa.nmolcentre] and ia <= cc.nal[pa.nmolcentre]
                while j < pa.nclplt and iflagt == 0:
                    j = j + 1
                    if pa.iplot[ia] > 0 or centralmolecule:
                        iflagt = j
                        ntag = ntag + 1
                        units.write(10, fmt(_F113, Ch(molecom.atsymbol[molecom.itype[ia]], 3), pa.xat[ia],
                                            pa.yat[ia], pa.zat[ia], Ch(molecom.alabel[ia], 3)))

            # Write a blank line at end of file

            units.write(10, lw())

        units.close(10)

    else:

        # If separate files are being stored, then we need to get the last block number

        pa.fnamesplot[0] = fstr(substr(dn.directory_name, 1, dn.nchardir) + substr(readfile, 1, ilenreadfile)
                                + '.nblock', 256)

        # Write xyz output file

        try:
            units.open(10, pa.fnamesplot[0], status='unknown')
            ierr = 0
        except FortranIOError as e:
            ierr = e.iostat

        nblocks = 0
        if ierr == 0:
            try:
                ir = units.list_read(10).int(ir)
                ierr = 0
            except FortranEOF:
                ierr = -1
            except FortranIOError as e:
                ierr = e.iostat
        if ierr == 0:
            nblocks = ir
        #      write(6,*) ir,nblocks,ierr
        units.rewind(10)

        nblocks = nblocks + 1
        units.write(10, lw(nblocks))

        units.close(10)

        # Generate the name of the file to write the data to

        blocknumber = fstr(fmt('(i5.5)', nblocks), 5)
        pa.fnamesplot[0] = fstr(substr(dn.directory_name, 1, dn.nchardir) + substr(readfile, 1, ilenreadfile)
                                + blocknumber + '.xyz', 256)
        pa.ilenfnamesplot[0] = dn.nchardir + ilenreadfile + 5 + 4

        # Write xyz output file

        try:
            units.open(10, pa.fnamesplot[0], status='unknown')
            ierr = 0
        except FortranIOError as e:
            ierr = e.iostat

        if ierr == 0:

            # First determine how many atoms to plot

            natplot = 0
            for ia in range(1, molecom.natom + 1):
                if pa.iplot[ia] > 0:
                    natplot = natplot + 1
            units.write(10, lw(natplot))

            # Output some information text

            units.write(10, fmt(_F100, '#Filename: ' + substr(dn.directory_name, 1, dn.nchardir)
                                + substr(readfile, 1, ilenreadfile) + '.ato  Block ' + blocknumber))
            for ia in range(1, molecom.natom + 1):
                # Check if a central molecule has been defined.
                centralmolecule = pa.nmolcentre > 0 and pa.nmolcentre <= cc.nmol
                # Check if the current atom is in the central molecule
                if centralmolecule:
                    centralmolecule = ia >= cc.naf[pa.nmolcentre] and ia <= cc.nal[pa.nmolcentre]
                j = 0
                iflagt = 0
                while j < pa.nclplt and iflagt == 0:
                    j = j + 1
                    if pa.iplot[ia] > 0 or centralmolecule:
                        iflagt = j
                        ntag = ntag + 1
                        units.write(10, fmt(_F113, Ch(molecom.atsymbol[molecom.itype[ia]], 3), pa.xat[ia],
                                            pa.yat[ia], pa.zat[ia], Ch(molecom.alabel[ia], 3)))

            # Write a blank line at end of file

            units.write(10, lw())

        units.close(10)

    return
