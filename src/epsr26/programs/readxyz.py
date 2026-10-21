"""Translation of fortran_src/readxyz.f: ``program standalone``, ``subroutine run_readxyz``
and ``subroutine readxyz``.

Reads an .xyz file (atom count, a line with the cell a b c alpha beta gamma, then
``label x y z`` lines), makes every atom a one-atom molecule with default potential
parameters chosen from its chemical symbol, and writes the .ato file.

The reference build compiles readxyz.f with an explicit INTERFACE for ``ctype``
(tools/build_patched_readxyz_readpdb.sh): the original calls ctype with 11 of its 13
arguments and no interface and crashes.  This translation behaves like the patched
program: the two missing optional arguments are absent.  See docs/translation_notes/N2.md.
"""

from ..common import (chemical_symbols, clustercom, commands, lattice, molecom, molecom_mod, plot_pars,
                      symbols)
from ..common.epsr_dimension import matom, mmol
from ..fortran import (Ch, FortranEOF, FortranIOError, f32, fmt, fstr, int_array, lw, r4, read_internal,
                       sqrtf, substr, units)
from ..library import (ato_manipulation_rev4, ato_read_and_save, auxiliary_routines, file_and_directory_commands,
                       parse)
from ..trace import routine

_F499 = "(/'readxyz> No .xyz file specified! - routine will exit'/)"
_F599 = "(/'readxyz> Problems with .xyz file ',a/'Routine will exit now'/)"
# 115 format(/'readxyz> Specify the chemical symbol for atom type ',a3,': ',$)
_F115 = "(/'readxyz> Specify the chemical symbol for atom type ',a3,': ')"
_F201 = "(/'readxyz> Following molecule types defined initially:-'/)"
_F202 = "(1x,i3,1x,a3,1x,a10,2(1x,e12.5))"


def main():
    """readxyz.f:standalone -- 2009/3/30 Attempt to make this a standalone program."""
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_readxyz()


@routine
def run_readxyz():
    """readxyz.f:run_readxyz -- read an .xyz file and write it as a .ato file."""
    # character*256 readfile
    # character*4 ext
    # character*1 iaccept
    # character*3020 sometextsave
    # character*4 yesorno
    # integer iread,nwordssave,ncfsave(nworddim),nclsave(nworddim)
    cm = commands
    iaccept = ''

    # First identify the number of arguments read in

    cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)

    # generate the name of the .xyz file

    ext = '.xyz'
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
    symbols.fnameato = fstr(substr(symbols.fnameato, 1, symbols.ilenfnameato) + '.xyz', 256)
    symbols.ilenfnameato = symbols.ilenfnameato + 4
    #
    # initialise number of atom and molecule types to zero
    #
    molecom.ntype = 0
    clustercom.nmoltype = 0
    #
    # read .xyz file
    #
    molecom.sizefactor = f32(1.0)
    # rho (real) and ierr are locals of run_readxyz
    rho, ierr = readxyz(symbols.fnameato, f32(0.0), 0)
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

    # Ensure the coordinates of each atom in each molecule are relative to the centre of mass

    for im in range(1, clustercom.nmol + 1):

        ato_manipulation_rev4.centreofmass(im)

    # Replace the .xyz extension with the .ato extension

    ilendot = auxiliary_routines.indexlastdot(symbols.fnameato) - 1
    symbols.fnameato = fstr(substr(symbols.fnameato, 1, ilendot) + '.ato', 256)
    ato_read_and_save.moloutc(symbols.fnameato)

    return


@routine
def readxyz(readfile, rho, ierr):
    """readxyz.f:readxyz -- read the .xyz file readfile into single-atom molecules.

    Returns (rho, ierr).
    """
    # use molecom_mod
    # character*256 readfile
    # integer ilenreadfile,ierr,iline,natomexpect
    # real vanderWaalsFactor,rho,xdif,ydif,zdif,pi,piconv,alphaaxes,betaaxes,gammaaxes
    # real, dimension(melements) :: vanderWaalsRadius
    # integer, dimension(:), allocatable :: iatomicno
    # character*3 a1,ats
    # integer iexch,j,k,ir,ia,ja
    # external cosd,sind   -> EPSR's own cosd/sind (readxyz.o calls cosd_, sind_)
    # logical found
    cm = commands
    cc = clustercom
    lt = lattice
    mc = molecom
    cs = chemical_symbols
    # implicitly typed real locals (alpha, beta, gam are read before use; zmax, rmax, cell,
    # perpab... are assigned)
    alpha = beta = gam = f32(0.0)
    natomexpect = 0

    ilenreadfile = len(readfile)  # noqa: F841
    try:
        units.open(10, readfile, status='old')
        ierr = 0
    except FortranIOError as e:
        ierr = e.iostat
    #      write(6,*) 'readxyz> Open file ',ierr
    if ierr != 0:
        return rho, ierr  # Any errors opening file, return
    iline = 0
    natomexpect = 0
    mc.sizefactor = f32(1.0)
    lt.cubic = False
    cc.nmol = 0
    mc.simtk = f32(300.0)
    lt.phib = f32(90.0)
    lt.thetac = f32(0.0)
    lt.phic = f32(0.0)
    mc.tol = f32(0.0)
    mc.tolout = f32(0.0)
    # Ensure no atom moves are attempted
    mc.stepmi = f32(0.0)
    mc.stepri = f32(0.0)
    mc.stepr = f32(0.0)
    mc.stepm = f32(0.0)
    mc.tempangreq = f32(1.0)
    mc.tempdihreq = f32(10.0)
    ia = 0
    im = 0
    ig = 0  # noqa: F841
    mc.natom = 0
    mc.natomr = 0
    molecom_mod.ngrrot = 0
    molecom_mod.nrottot = 0
    mc.ntype = 0
    # Reset array size counters
    molecom_mod.ngrrotmax = 0
    molecom_mod.nrottotmax = 0
    while ierr == 0 and ia < matom and im < mmol:
        try:
            cm.sometext = fstr(units.read_line(10), cm.nstsize)
            ierr = 0
        except FortranEOF:
            ierr = -1
        except FortranIOError as e:
            ierr = e.iostat
        #         write(6,*) 'readxyz> Read line: ',iline,ierr
        if ierr == 0:
            cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 11, cm.nstsize, cm.nworddim)
            iline = iline + 1
            if iline == 1:
                natomexpect = read_internal(cm.sometext).int(natomexpect)
            elif iline == 2:
                if cm.nwords < 7:
                    ierr = 1
                else:  # Get the crystallographic axes and angles
                    try:
                        r = read_internal(substr(cm.sometext, cm.ncf[2], cm.ncl[7]))
                        lt.alat = r.real(lt.alat)
                        lt.blat = r.real(lt.blat)
                        lt.clat = r.real(lt.clat)
                        alpha = r.real(alpha)
                        beta = r.real(beta)
                        gam = r.real(gam)
                        ierr = 0
                    except FortranEOF:
                        ierr = -1
                    except FortranIOError as e:
                        ierr = e.iostat
                while ierr != 0:
                    units.write(6, lw('readxyz> Need a,b,c axes and 3 angles (alpha,beta,gamma):'))
                    try:
                        r = units.list_read(5)
                        lt.alat = r.real(lt.alat)
                        lt.blat = r.real(lt.blat)
                        lt.clat = r.real(lt.clat)
                        alpha = r.real(alpha)
                        beta = r.real(beta)
                        gam = r.real(gam)
                        ierr = 0
                    except FortranEOF:
                        ierr = -1
                    except FortranIOError as e:
                        ierr = e.iostat
                # Convert angles to internal spherical polar coordinates
                lt.phib, lt.thetac, lt.phic = auxiliary_routines.convertalphabetagamma(
                    alpha, beta, gam, lt.phib, lt.thetac, lt.phic)
                # set up lattice parameters for applying the minimum image convention
                lt.alatx = lt.alat
                lt.blatx = lt.blat * auxiliary_routines.cosd(lt.phib)
                lt.blaty = lt.blat * auxiliary_routines.sind(lt.phib)
                lt.clatx = lt.clat * auxiliary_routines.sind(lt.thetac) * auxiliary_routines.cosd(lt.phic)
                lt.claty = lt.clat * auxiliary_routines.sind(lt.thetac) * auxiliary_routines.sind(lt.phic)
                lt.clatz = lt.clat * auxiliary_routines.cosd(lt.thetac)
                # define the Cartesian coordinates of the lattice vectors
                lt.alatc[1] = lt.alatx
                lt.alatc[2] = 0.0
                lt.alatc[3] = 0.0
                lt.blatc[1] = lt.blatx
                lt.blatc[2] = lt.blaty
                lt.blatc[3] = 0.0
                lt.clatc[1] = lt.clatx
                lt.clatc[2] = lt.claty
                lt.clatc[3] = lt.clatz
                lt.bxax = lt.blatx / lt.alatx
                lt.cxax = lt.clatx / lt.alatx
                lt.cyby = lt.claty / lt.blaty
                # calculate the perpendicular distance between the faces of the lattice
                perpab = auxiliary_routines.perpd(lt.alatc, lt.blatc, lt.clatc)
                perpbc = auxiliary_routines.perpd(lt.blatc, lt.clatc, lt.alatc)
                perpca = auxiliary_routines.perpd(lt.clatc, lt.alatc, lt.blatc)
                units.write(6, lw(perpab, perpbc, perpca))
                # set up maximum radius to be used for Fourier transforms and to truncate
                # the potential. This is defined as the minimum perpendicular distance between
                # planes.
                #
                # G(r) and potentials will be calculated to rmax, or maximum allowed distance
                # based on number of points and step size, whichever is smaller.
                cc.dmax = min(perpab, perpbc, perpca)
                units.write(6, lw(cc.nmol, cc.dmax))
                cell = cc.dmax
                #
                # setup molecular size factor
                #
                mc.dmaxmol = cc.dmax * mc.sizefactor
                #
                # determine maximum distance across box
                #
                # G(r) and potentials will be calculated to rmax, or maximum allowed distance
                # based on number of points and step size, whichever is smaller.
                #
                plot_pars.xmax = abs(lt.alatx) + abs(lt.blatx) + abs(lt.clatx)
                plot_pars.ymax = abs(lt.blaty) + abs(lt.claty)
                zmax = abs(lt.clatz)
                rmax = f32(0.5) * sqrtf(plot_pars.xmax * plot_pars.xmax + plot_pars.ymax * plot_pars.ymax
                                        + zmax * zmax)
                units.write(6, lw(rmax, cell))
                #
                # convert temperature in K to kJ/mole
                #
                mc.simt = r4(0.008314) * mc.simtk
                #
                # calculate volume of box
                #
                cc.volume = lt.alatx * lt.blaty * lt.clatz
            else:  # Otherwise we read the list of atomic coordinates, creating single atom molecules
                if cm.nwords < 4:
                    ierr = 1
                else:
                    ia = ia + 1
                    im = im + 1
                    # NB: the Fortran reads yadd(ia), not yadd(im); ia = im here
                    try:
                        r = read_internal(substr(cm.sometext, cm.ncf[1], cm.ncl[4]))
                        mc.alabel[ia] = fstr(r.word(mc.alabel[ia]), 3)
                        cc.xadd[im] = r.real(cc.xadd[im])
                        cc.yadd[ia] = r.real(cc.yadd[ia])
                        cc.zadd[im] = r.real(cc.zadd[im])
                        ierr = 0
                    except FortranEOF:
                        ierr = -1
                    except FortranIOError as e:
                        ierr = e.iostat
                    mc.imageatom[ia] = 0
                    mc.xa[ia] = 0.0
                    mc.ya[ia] = 0.0
                    mc.za[ia] = 0.0
                    # set up 'type' of atom - this is used in assigning neutron weights
                    ato_manipulation_rev4.ftype(ia)
                    #
                    # do not count atoms which begin with 'Q' or 'q' - these are assumed
                    # invisible to the scattering radiation and experiment (although they
                    # can be used to calculate other structural quantities if necessary
                    #
                    if mc.qtype[ia] != 1:
                        mc.natomr = mc.natomr + 1
                    mc.natom = mc.natom + 1
                    mc.nbond[ia] = 0
                    cc.naf[im] = ia
                    cc.nal[im] = ia
                    cc.nato[im] = 1
                    cc.ngrot[im] = 0
                    molecom_mod.ngrrot = molecom_mod.ngrrot + int(cc.ngrot[im])
                    cc.tethered[im] = False
                    cc.phix[im] = 0.0
                    cc.phiy[im] = 0.0
                    cc.phiz[im] = 0.0
    units.close(10)
    cc.nmol = im
    #      do ia=1,natom
    #         write(6,*) natomexpect,ia,alabel(ia),xa(ia),ya(ia),za(ia)
    #      end do
    rho = f32(mc.natomr) / cc.volume
    units.write(6, lw('readxyz> rho = ', rho, mc.natomr, mc.natom, natomexpect))
    if mc.natom >= natomexpect:
        ierr = 0
    # Set up the atomic numbers of all the atoms
    iatomicno = None
    iatomicno = int_array(mc.natom)
    for ia in range(1, mc.natom + 1):
        ir = 0
        found = False
        while not found and ir < cs.melements:
            ir = ir + 1
            # It is assumed the atom labels are also the chemical symbols
            found = mc.alabel[ia] == cs.chsymbol[ir]
            if found:
                iatomicno[ia] = ir
        if not found:  # If a valid symbol is not found, then the file cannot be processed
            ierr = 1
            units.write(6, lw('readxyz> Cannot read this file: atom label not recognised.'))
            return rho, ierr
    # Use changeato to sort these atoms into molecules
    #
    # make a list of the atom types read in from the current file
    #
    mc.ntypefile = 0
    ia = 0
    while ia < mc.natom:
        ia = ia + 1
        itypefile = 0
        ifound = 0
        while itypefile < mc.ntypefile and ifound == 0:
            itypefile = itypefile + 1
            if mc.alabel[ia] == mc.typefile[itypefile]:
                ifound = itypefile
        if ifound == 0:
            mc.ntypefile = mc.ntypefile + 1
            mc.typefile[mc.ntypefile] = mc.alabel[ia]
            mc.iatomicnotype[mc.ntypefile] = iatomicno[ia]
    # Use these values to set the potential parameters
    for it in range(1, mc.ntypefile + 1):
        # epsilon, sigma, amass are implicit REAL locals
        if cs.chsymbol[mc.iatomicnotype[it]] == 'H':
            epsilon = f32(0)
            sigma = f32(0)
            amass = f32(2)
        elif cs.chsymbol[mc.iatomicnotype[it]] == 'O':
            epsilon = r4(0.65)
            sigma = r4(3.1)
            amass = cs.chatwt[mc.iatomicnotype[it]]
        elif cs.chsymbol[mc.iatomicnotype[it]] == 'N':
            epsilon = r4(0.7)
            sigma = r4(3.2)
            amass = cs.chatwt[mc.iatomicnotype[it]]
        elif cs.chsymbol[mc.iatomicnotype[it]] == 'C':
            epsilon = r4(0.8)
            sigma = r4(3.7)
            amass = cs.chatwt[mc.iatomicnotype[it]]
        else:
            epsilon = r4(0.8)
            sigma = r4(3.2)
            amass = cs.chatwt[mc.iatomicnotype[it]]
        # call ctype(typefile(it),chsymbol(iatomicnotype(it)),0,epsilon,sigma,amass,1.0,0.0,0.5*sigma,1.0,0.0)
        # 11 of 13 arguments: qt=1.0, rqt=0.0, natomsph=0.5*sigma, rhosph=1.0, radminsph=0.0;
        # radmaxsph and pseudofacsph are absent (patched reference build, see the module docstring)
        ato_manipulation_rev4.ctype(mc.typefile[it], fstr(cs.chsymbol[mc.iatomicnotype[it]], 2), 0, epsilon, sigma,
                                    amass, r4(1.0), r4(0.0), r4(0.5) * sigma, r4(1.0), r4(0.0))
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
            if mc.alabel[cc.naf[im]] == cc.moltype[imoltype]:
                ifound = imoltype
                cc.mtype[im] = ifound
        if ifound == 0:
            cc.nmoltype = cc.nmoltype + 1
            cc.moltype[cc.nmoltype] = mc.alabel[cc.naf[im]]
            cc.moltethered[cc.nmoltype] = cc.tethered[im]
            #
            # initialise the template filename
            #
            cc.moltemplate[cc.nmoltype] = fstr(fmt('("moltype",i2.2)', cc.nmoltype), 256)
            cc.mtype[im] = cc.nmoltype

            # Initialise the rotational and translational step sizes

            cc.interntrans[cc.nmoltype] = mc.stepmi
            cc.internrot[cc.nmoltype] = mc.stepri
            cc.wholerot[cc.nmoltype] = mc.stepr
            cc.wholetrans[cc.nmoltype] = mc.stepm

    # Write these out as a check

    units.write(6, fmt(_F201))

    for imoltype in range(1, cc.nmoltype + 1):

        units.write(6, fmt(_F202, imoltype, Ch(cc.moltype[imoltype], 3), substr(cc.moltemplate[imoltype], 1, 10),
                           cc.wholerot[imoltype], cc.wholetrans[imoltype]))

    return rho, ierr
