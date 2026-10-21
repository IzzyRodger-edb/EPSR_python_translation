"""Translation of fortran_src/plot3djmol.f: ``program standalone`` (the ``plot3djmol``
program), ``subroutine run_plot3djmol`` and ``subroutine map3djmol``.

Makes a Jmol cube file (and Jmol script) of an isosurface of a spatial density
function from the spherical harmonic coefficients written by SHARM/SDF, and starts
Jmol on it (through ``runtime.system`` with the Fortran's command string).

map3djmol declares ``external cosd,sind,acosd,asind,atan2d``: these are EPSR's own
(auxiliary_routines); plain sin/cos/acos/atan2/sqrt of reals are the libm float
functions.  The plotsharm_routines, setup_symbols_generic and plot_routines_g95 modules
are imported inside the functions.  See docs/translation_notes/N2.md.
"""

from ..common import SaveBlock, commands, directory_name, symbols, system_commands
from ..common import molecule_coordinates_mod as mcm
from ..fortran import (
    Ch,
    FortranIOError,
    FortranStop,
    acosf,
    atan2f,
    char_array,
    cosf,
    f32,
    findex,
    fmt,
    fstr,
    int_array,
    isign,
    lw,
    nint,
    r4,
    read_internal,
    real_array,
    runtime,
    sinf,
    sqrtf,
    substr,
    units,
)
from ..library import auxiliary_routines, file_and_directory_commands, parse
from ..trace import routine

_F51 = "(/'plot3djmol> Error code ',i5,' on attempting to open file ',a)"
_F102 = "(/a,1x,'Will use plot3djmol input file: \"',a,'\"')"
_F92 = "(1x,10(1x,i5))"
_F999 = "(1x,' No. of occupied voxels = ',i9,' out of ',i9,' total.')"
_F795 = (
    "('load \"',a,'\"; zoom 20; background {'"
    ",f4.1,1x,f4.1,1x,f4.1,'}; spacefill -15; rotate z -90; rotate x -90; set frank off; set defaulttranslucent '"
    ",f6.3,';center {0 0 0};set perspectiveDepth off ')"
)
_F807 = "('rotate y ',f10.5,'; rotate x ',f10.5,'; spin ',a3)"
_F797 = (
    "('isosurface ',a,1x,f10.5,' \"',a,'\"; color isosurface {'"
    ",f4.1,1x,f4.1,1x,f4.1,'}; isosurface ',a)"
)


def main():
    """plot3djmol.f:standalone -- 2009/3/30 Attempt to make this a standalone program."""
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_plot3djmol()


@routine
def run_plot3djmol():
    """plot3djmol.f:run_plot3djmol -- invokes the system command to perform 3d plotting
    surfaces using Jmol."""
    from ..library import plot_routines_g95, read_and_save_files
    from ..library import setup_symbols_generic

    # character*256 plot_output,fname_out,nextfile
    # integer ierrout,plotnumber
    sc = system_commands
    sy = symbols
    # NB: ierrout is set by setup_input_file before use; fname_out is set by the first
    # map3djmol call (newstart) unless it returns early (then '' here, garbage in the Fortran)
    ierrout = 0
    fname_out = ''
    #
    # setup the plot3d symbols
    #
    setup_symbols_generic.setup_symbols_plot3djmol()

    # Signal to start a new plot

    nextfile = '<newstart>'
    plotnumber = 0
    nin = 11
    nout = 6
    #
    # now setup the input file
    #
    while findex(nextfile, '<none>') == 0:
        nout, ierrout = file_and_directory_commands.setup_input_file(nin, nout, ierrout)
        #
        # read the input file
        #
        nout, ierrout = read_and_save_files.read_input_file(nin, nout, ierrout)
        if ierrout != 0:  # Error occurred - try reading using the old format
            ierrout = 0
            # Re-open the file that would have been closed by read_input_file
            try:
                units.open(nin, sy.fnameinp, status='old')
            except FortranIOError as e:
                ierrout = e.iostat
            nout, ierrout = plot_routines_g95.read_plot3djmol(nin, nout, ierrout)
        #
        # if an error was encountered on setup, exit now
        #
        if ierrout != 0:
            units.write(6, fmt(_F51, ierrout, substr(nextfile, 1, len(nextfile))))
            return
        units.write(6, fmt(_F102, substr(sy.prompt, 1, sy.ilenprompt), substr(sy.fname_inp, 1, sy.ilenfname_inp)))

        plotnumber = plotnumber + 1
        fname_out, nextfile = map3djmol(sy.fname_inp, fname_out, nextfile, plotnumber)
        commands.sometext = fstr(substr(sy.symbol_setup, 1, sy.ilensetup) + ' ' + nextfile, commands.nstsize)
        units.write(6, lw(substr(commands.sometext, 1, len(commands.sometext))))
    #
    # now generate the command that will actually run this file
    #
    plot_output = fstr(substr(sc.system_cd, 1, sc.ilensystem_cd)
                       + ' '
                       + substr(sc.system_quotes, 1, sc.ilensystem_quotes)
                       + substr(directory_name.directory_name, 1, directory_name.nchardir)
                       + substr(sc.system_quotes, 1, sc.ilensystem_quotes)
                       + ' ' + substr(sc.system_join, 1, sc.ilensystem_join)
                       + ' ' + substr(sc.system_jmol, 1, sc.ilensystem_jmol) + ' -s '
                       + substr(fname_out, 1, len(fname_out)), 256)
    units.write(6, lw(Ch(plot_output, 256)))
    ierr = runtime.system(plot_output)
    # call system(plot_output,ierr): the GNU SYSTEM subroutine returns the raw wait
    # status of system(3) (exit code * 256, or the signal number); runtime.system returns
    # the exit code (or -signal), so convert it.
    ierr = ierr * 256 if ierr >= 0 else -ierr
    units.write(6, lw(ierr))
    return


# character(len=256), save :: firstmoleculename
# integer, save :: lenfirstmoleculename,lennextfile,lenmolfile
# character(len=256), save :: molfile
_save = SaveBlock("map3djmol", firstmoleculename='', lenfirstmoleculename=0, lennextfile=0,
                  lenmolfile=0, molfile='')


def _pm1(n):
    """``(-1)**n`` for an integer n."""
    return 1 if n % 2 == 0 else -1


@routine
def map3djmol(infile, plotfile, nextfile, plotnumber):
    """plot3djmol.f:map3djmol -- one isosurface: g(r,omega) from the spherical harmonic
    coefficients, interpolated onto a cubic grid, written as a Gaussian cube file, plus
    the molecule cube file and the Jmol script.

    Returns (plotfile, nextfile).
    """
    # use molecule_coordinates_mod
    # use lngamma_mod
    from ..library import lngamma_mod, molecule_coordinates_mod, plotsharm_routines

    cm = commands
    sy = symbols
    ar = auxiliary_routines
    cosd = ar.cosd
    sind = ar.sind
    acosd = ar.acosd
    atan2d = ar.atan2d
    # external cosd,sind,acosd,asind,atan2d

    pi = r4(3.141592653)
    bohr = r4(0.52917721092)
    mexpan = 3000
    datac = real_array(500, mexpan)
    datas = real_array(500, mexpan)
    rd = real_array(500)
    zcoeffc = real_array(mexpan)
    zcoeffs = real_array(mexpan)
    grw = real_array(50, 50, 500)
    temp = real_array(500)
    rgbobj = real_array(3)
    rgbbak = real_array(3)
    l1ex = int_array(mexpan)
    l2ex = int_array(mexpan)
    lex = int_array(mexpan)
    n1ex = int_array(mexpan)
    n2ex = int_array(mexpan)
    titled = char_array(mexpan)  # character(len=29)
    titlex = ''  # character(len=14)
    # NB: locals that the Fortran leaves uninitialised on some paths start at 0 here:
    # phim, thetam, chim, thetal, phil (only some are read, depending on nvary),
    # phi, chi (if nvary is not 1-3), xcenterold/ycentreold/zcentreold (misspelt in the
    # Fortran: never assigned, but written to the -mol.cube file when a molecule is added).
    phim = thetam = chim = thetal = phil = f32(0.0)
    phi = chi = f32(0.0)
    xcenterold = ycentreold = zcentreold = f32(0.0)
    ierr = 0
    nr = 0
    nexpan = 0
    ident = 0
    kk0 = 0
    nsmoo = 0
    xmax = f32(0.0)
    mpixel = 0
    rmin = rmax = surfra = f32(0.0)
    l1only = l2only = lonly = n1only = n2only = m2fix = 0
    nvary = 0
    nspher = 0
    rot = elev = f32(0.0)
    scalefactor = translucency = f32(0.0)
    itrans = imesh = 0
    lennextfile = _save.lennextfile

    # Decide whether this is a new plot, or adding to an existing plot
    newstart = substr(nextfile, 1, len(nextfile)) == '<newstart>'
    # Generate an isosurface label
    surfacelabel = fstr(fmt("(a1,i5.5)", 's', plotnumber), 10)
    # Initialise random number generator
    units.open(10, 'ran1.dat', status='unknown')
    auxiliary_routines.initran1(10)
    units.close(10)
    nxdim = 500
    nfdim = mexpan
    nin = 10  # noqa: F841
    #
    # initialise log(factorial)
    #
    lngamma_mod.fctrln(250)
    fname = fstr(infile, 256)
    # Name of output file - remove any text after the last dot
    lenfirstdot = findex(fname, '.CUBE.txt')  # Make sure this is a .CUBE file
    if lenfirstdot > 0:
        outfile = fstr(substr(fname, 1, lenfirstdot) + 'cube', 256)
    else:
        # NB: the Fortran assigns indexlastdot to lendfirstdot (sic) and then tests
        # lenfirstdot, which is 0 here: '.cube' is always appended.
        lendfirstdot = ar.indexlastdot(fname)  # noqa: F841
        if lenfirstdot > 0:
            outfile = fstr(substr(fname, 1, lenfirstdot) + 'cube', 256)
        else:
            outfile = fstr(substr(fname, 1, len(fname)) + '.cube', 256)
    lenoutfile = len(outfile)
    if newstart:
        plotfile = fstr('plot' + substr(outfile, 1, ar.indexlastdot(outfile)) + 'txt', 256)
        _save.molfile = fstr(substr(outfile, 1, findex(outfile, '.cube') - 1) + '-mol.cube', 256)
        _save.lenmolfile = len(_save.molfile)
    molfile = _save.molfile
    iref = 0
    iref = iref + 1
    cm.sometext = fstr(sy.symbol_value[iref], cm.nstsize)
    cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
    _, cm.nwords, datafile, lendatafile = parse.concatenate(1, cm.nwords, '', 0, cm.nwords, cm.ncf, cm.ncl,
                                                            cm.sometext)
    iref = iref + 1
    nexpan = read_internal(sy.symbol_value[iref]).int(nexpan)
    # A value of ident .lt. 0 will be used to signify that we want molecule 2 to be
    # at the origin
    iref = iref + 1
    ident = read_internal(sy.symbol_value[iref]).int(ident)

    # save the sign of ident and keep separate from the absolute value

    if abs(ident) > 0:
        identsign = isign(1, ident)
        ident = abs(ident)
    else:
        identsign = 1
        ident = 0

    iref = iref + 1
    kk0 = read_internal(sy.symbol_value[iref]).int(kk0)
    if kk0 < 0:
        kk0 = 0
    if kk0 > 2:
        kk0 = 2
    # Get the titles from the input file
    datafile, nexpan, titlex, ierr = ar.r_int_titles(datafile, nexpan, titlex, titled, nfdim, ierr)
    for k in range(1, nexpan + 1):
        temptitle = fstr(titled[k], 29)
        lentemptitle = 29
        #         write(6,*) k,temptitle
        ifirst = 1
        ifirst = findex(substr(temptitle, ifirst, lentemptitle), '(') + 1
        ilast = ifirst + findex(substr(temptitle, ifirst, lentemptitle), ',') - 2
        l1ex[k] = read_internal(substr(temptitle, ifirst, ilast)).int(l1ex[k])
        ifirst = ilast + 2
        ilast = ifirst + findex(substr(temptitle, ifirst, lentemptitle), ',') - 2
        l2ex[k] = read_internal(substr(temptitle, ifirst, ilast)).int(l2ex[k])
        ifirst = ilast + 2
        ilast = ifirst + findex(substr(temptitle, ifirst, lentemptitle), ';') - 2
        lex[k] = read_internal(substr(temptitle, ifirst, ilast)).int(lex[k])
        ifirst = ilast + 2
        ilast = ifirst + findex(substr(temptitle, ifirst, lentemptitle), ',') - 2
        n1ex[k] = read_internal(substr(temptitle, ifirst, ilast)).int(n1ex[k])
        ifirst = ilast + 2
        ilast = ifirst + findex(substr(temptitle, ifirst, lentemptitle), ')') - 2
        n2ex[k] = read_internal(substr(temptitle, ifirst, ilast)).int(n2ex[k])
        units.write(6, fmt(_F92, k, l1ex[k], l2ex[k], lex[k], n1ex[k], n2ex[k]))
    #
    # read in spherical harmonics data
    #
    nr, datafile, ierr = ar.r_int(nr, datafile, nexpan, rd, datac, datas, nxdim, nfdim, ierr)
    rstep = rd[3] - rd[2]
    iref = iref + 1
    nsmoo = read_internal(sy.symbol_value[iref]).int(nsmoo)
    iref = iref + 1
    xmax = read_internal(sy.symbol_value[iref]).real(xmax)
    iref = iref + 1
    mpixel = read_internal(sy.symbol_value[iref]).int(mpixel)
    iref = iref + 1
    r = read_internal(sy.symbol_value[iref])
    rmin = r.real(rmin)
    rmax = r.real(rmax)
    iref = iref + 1
    surfra = read_internal(sy.symbol_value[iref]).real(surfra)
    iref = iref + 1
    cm.sometext = fstr(sy.symbol_value[iref], cm.nstsize)
    cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
    r = read_internal(substr(cm.sometext, cm.ncf[1], cm.ncl[cm.nwords]))
    if cm.nwords > 2:
        l1only = r.int(l1only)
        l2only = r.int(l2only)
        lonly = r.int(lonly)
    else:
        l1only = r.int(l1only)
        l2only = r.int(l2only)
        lonly = 1
    iref = iref + 1
    r = read_internal(sy.symbol_value[iref])
    n1only = r.int(n1only)
    n2only = r.int(n2only)
    iref = iref + 1
    m2fix = read_internal(sy.symbol_value[iref]).int(m2fix)
    iref = iref + 1
    forcemirrorsym = sy.symbol_value[iref] == 'y'
    iref = iref + 1
    nvary = read_internal(sy.symbol_value[iref]).int(nvary)
    iref = iref + 1
    r = read_internal(sy.symbol_value[iref])
    if nvary == 1:
        phim = r.real(phim)
        thetam = r.real(thetam)
        chim = r.real(chim)
    elif nvary == 2:
        thetal = r.real(thetal)
        phil = r.real(phil)
        chim = r.real(chim)
    elif nvary == 3:
        thetal = r.real(thetal)
        phil = r.real(phil)
        phim = r.real(phim)
    iref = iref + 1
    nspher = read_internal(sy.symbol_value[iref]).int(nspher)
    nsize = nspher
    rdcoor = real_array(nsize)
    thcoor = real_array(nsize)
    phcoor = real_array(nsize)
    radbal = real_array(nsize)  # noqa: F841
    rgbbal = real_array(3, nsize)  # noqa: F841
    chsymbol = char_array(nsize)  # character(len=3)
    piconv = r4(3.141592653) / r4(180.0)
    lenmoleculename = 0
    moleculename = ''
    for is_ in range(1, nspher + 1):
        iref = iref + 1
        chsymbol[is_] = fstr(substr(sy.symbol_value[iref], 1, 2), 3)  # '(a2)'
        iref = iref + 1
        r = read_internal(sy.symbol_value[iref])
        rdcoor[is_] = r.real(rdcoor[is_])
        thc = r.real(f32(0.0))
        phc = r.real(f32(0.0))
        thcoor[is_] = thc * piconv
        phcoor[is_] = phc * piconv
        if is_ == 1:
            moleculename = fstr(chsymbol[is_], 256)
            lenmoleculename = len(chsymbol[is_])
        elif lenmoleculename < 256 - 3:
            moleculename = fstr(substr(moleculename, 1, lenmoleculename) + chsymbol[is_], 256)
            lenmoleculename = lenmoleculename + len(chsymbol[is_])
    if kk0 < 2:  # Use atoms read in from xyz file
        ic = 0
        # Get the name of the central molecule
        moleculefile = fstr(substr(datafile, 1, ar.indexlastdot(datafile) - 1) + '-first.xyz', 256)
        lenmoleculefile = len(moleculefile)
        molecule_coordinates_mod.get_coords_from_xyz(moleculefile, lenmoleculefile)
        moleculename = mcm.atomsequence
        lenmoleculename = mcm.lenatomsequence
        if nvary > 1:  # In this case read the second molecule
            moleculefile = fstr(substr(datafile, 1, ar.indexlastdot(datafile) - 1) + '-second.xyz', 256)
            lenmoleculefile = len(moleculefile)
            molecule_coordinates_mod.get_coords_from_xyz(moleculefile, lenmoleculefile)
        nsize = mcm.nofatoms
        rdcoor = real_array(nsize)
        thcoor = real_array(nsize)
        phcoor = real_array(nsize)
        radbal = real_array(nsize)  # noqa: F841
        rgbbal = real_array(3, nsize)  # noqa: F841
        chsymbol = char_array(nsize)
        is_ = 0
        while is_ < mcm.nofatoms and ic < nsize:
            is_ = is_ + 1
            if mcm.centralmolecule[is_] >= kk0:
                ic = ic + 1
                chsymbol[ic] = mcm.chsymbolt[is_]
                rdcoor[ic] = mcm.rt[is_]
                thcoor[ic] = mcm.thetat[is_]
                phcoor[ic] = mcm.phit[is_]
        nspher = ic
    # Save molecule name if this the first plot of the series
    if newstart:
        _save.firstmoleculename = moleculename
        _save.lenfirstmoleculename = lenmoleculename
    # Need to check that the same atom sequence is used for central molecule as for previous contours
    if (substr(moleculename, 1, lenmoleculename).rstrip(' ')
            != substr(_save.firstmoleculename, 1, _save.lenfirstmoleculename).rstrip(' ')):
        units.write(6, fmt('(a,a,a,a)', 'map3djmol> Molecule symbols: ', substr(moleculename, 1, lenmoleculename),
                           ':', substr(_save.firstmoleculename, 1, _save.lenfirstmoleculename)))
        units.write(6, fmt('(a/a/a)', ' map3djmol> You need to ensure the central molecule for file ',
                           substr(fname, 1, len(fname)), ' is the same as for the first plot.'))
        return plotfile, nextfile
    iref = iref + 1
    r = read_internal(sy.symbol_value[iref])
    rgbbak[1] = r.real(rgbbak[1])
    rgbbak[2] = r.real(rgbbak[2])
    rgbbak[3] = r.real(rgbbak[3])
    iref = iref + 1
    r = read_internal(sy.symbol_value[iref])
    rgbobj[1] = r.real(rgbobj[1])
    rgbobj[2] = r.real(rgbobj[2])
    rgbobj[3] = r.real(rgbobj[3])
    iref = iref + 1
    itrans = read_internal(sy.symbol_value[iref]).int(itrans)
    if itrans == 0:
        transparency = 'fill'
    elif itrans == 1:
        transparency = 'translucent'
    else:
        transparency = 'nofill'
    lentransparency = len(transparency)
    iref = iref + 1
    imesh = read_internal(sy.symbol_value[iref]).int(imesh)
    #      if(itrans.gt.1.and.imesh.eq.0) imesh=1 !Ensure something is plotted!
    if imesh == 0:
        mesh = 'nomesh nodots'
    elif imesh == 1:
        mesh = 'mesh'
    else:
        mesh = 'dots'
    lenmesh = len(mesh)
    surfacetype = fstr(substr(mesh, 1, lenmesh) + ' ' + substr(transparency, 1, lentransparency), 80)
    iref = iref + 1
    r = read_internal(sy.symbol_value[iref])
    rot = r.real(rot)
    elev = r.real(elev)
    iref = iref + 1
    scalefactor = read_internal(sy.symbol_value[iref]).real(scalefactor)
    iref = iref + 1
    translucency = read_internal(sy.symbol_value[iref]).real(translucency)
    iref = iref + 1
    rotatecentralmolecule = sy.symbol_value[iref] == 'y'
    if nvary == 1:
        rotatecentralmolecule = False  # Cannot rotate the central molecule for a spatial density function.
    iref = iref + 1
    cm.sometext = fstr(sy.symbol_value[iref], cm.nstsize)
    cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
    # NB: concatenate assigns its nlast dummy, here the literal 1, when nwords is 0 (blank
    # spin): gfortran then dies with SIGSEGV (constant in read-only memory); reproduced.
    if 1 > cm.nwords:
        raise FortranStop(139, "Program received signal SIGSEGV: Segmentation fault - invalid memory reference.")
    _, _, spin, ilen = parse.concatenate(1, 1, '', 0, cm.nwords, cm.ncf, cm.ncl, cm.sometext, len_return_string=3)
    if spin != 'x' and spin != 'y' and spin != 'z':
        # Ensure a sensible value of spin has been set. Actual behaviour will be determined
        # by last plot of a sequence.
        spin = 'off'
    iref = iref + 1
    cm.sometext = fstr(sy.symbol_value[iref], cm.nstsize)
    cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
    _, cm.nwords, nextfile, lennextfile = parse.concatenate(1, cm.nwords, nextfile, lennextfile, cm.nwords,
                                                            cm.ncf, cm.ncl, cm.sometext)
    _save.lennextfile = lennextfile
    # Set up pixels and voxels. Currently a cubic plot area is assumed.
    nx = mpixel
    xstep = xmax / f32(nx)
    ny = nx
    nz = nx
    ystep = xstep
    zstep = xstep
    nx0 = 2 * nx + 1
    ny0 = 2 * ny + 1
    nz0 = 2 * nz + 1
    dens = real_array(nx0, ny0, nz0)
    nvoxel = nx0 * ny0 * nz0
    # calculate g(r,w) at 10 deg intervals
    nt0 = 19
    tstep = f32(10)
    np0 = 36
    pstep = f32(10)
    #
    # initial maximum and minimum values
    #
    dhigh = f32(1.0)
    dlow = f32(1.0)
    phimmax = f32(0.0)
    thetammax = f32(0.0)
    chimmax = f32(0.0)
    irmax = 0  # Reference to distance at which maximum intensity is found
    for it in range(1, nt0 + 1):
        t = f32(it - 1) * tstep
        for ip in range(1, np0 + 1):
            p = f32(ip - 1) * pstep
            if nvary == 1:
                thetal = t
                phil = p
                phi = phim
                chi = chim
            elif nvary == 2:
                phi = p
                thetam = t
                chi = chim
            elif nvary == 3:
                phi = phim
                thetam = t
                chi = p
            ctl = cosd(thetal)
            cz = cosd(thetam)
            for k in range(1, nexpan + 1):
                sumc = f32(0.0)
                sums = f32(0.0)
                ll1 = int(l1ex[k])
                ll2 = int(l2ex[k])
                ll = int(lex[k])
                n1 = int(n1ex[k])
                n2 = int(n2ex[k])
                l1tmp = ll1
                l2tmp = ll2
                n1tmp = n1
                n2tmp = n2
                test = (l1only == 0 and l1tmp == 0) or l1only != 0
                test = test and ((l2only == 0 and l2tmp == 0) or l2only != 0)
                test = test and ((lonly == 0 and ll == 0) or lonly != 0)
                test = test and ((n1only == 0 and n1tmp == 0) or n1only != 0)
                test = test and ((n2only == 0 and n2tmp == 0) or n2only != 0)
                if test:
                    phs = f32(1.0)
                    phsn = f32(_pm1(l1tmp + l2tmp + ll + n1tmp + n2tmp))

                    # Normal algorithm for molecule 1 at the origin

                    if identsign >= 0:
                        if m2fix > 0:
                            m2f = -l2tmp
                            m2l = l2tmp
                        else:
                            m2f = 0
                            m2l = 0
                        sumc, sums = plotsharm_routines.sumoverm2(l1tmp, l2tmp, ll, n1tmp, m2f, m2l, n2tmp, phil, ctl,
                                                                  phi, cz, chi, phs, sumc, sums)
                        # now reverse signs of n1,n2 if either is non-zero
                        if n1 > 0 or (n1 == 0 and n2 > 0):
                            n1tmp = -n1
                            n2tmp = -n2
                            phs = phs * phsn
                            sumc, sums = plotsharm_routines.sumoverm2(l1tmp, l2tmp, ll, n1tmp, m2f, m2l, n2tmp, phil,
                                                                      ctl, phi, cz, chi, phs, sumc, sums)
                    else:
                        if m2fix > 0:
                            m2f = -l1tmp
                            m2l = l1tmp
                        else:
                            m2f = 0
                            m2l = 0

                        # Phase factor to reverse parity of centre to centre vector

                        phs = f32(_pm1(ll))
                        sumc, sums = plotsharm_routines.sumoverm2(l2tmp, l1tmp, ll, n2tmp, m2f, m2l, n1tmp, phil, ctl,
                                                                  phi, cz, chi, phs, sumc, sums)
                        # now reverse signs of n1,n2 if either is non-zero
                        if n1 > 0 or (n1 == 0 and n2 > 0):
                            n1tmp = -n1
                            n2tmp = -n2
                            phs = phsn * phs
                            sumc, sums = plotsharm_routines.sumoverm2(l2tmp, l1tmp, ll, n2tmp, m2f, m2l, n1tmp, phil,
                                                                      ctl, phi, cz, chi, phs, sumc, sums)
                # now repeat for interchange of l1,n1,with l2,n2 for identical molecules
                # bearing in mind the phase convention on the h's - the case when l1=l2 AND
                # n1=n2 doesn't count here.
                #
                if ident == 0 and (ll1 != ll2):  # l1 will always be .ge. l2 (from the way they are calculated).
                    phs = f32(_pm1(ll1 + ll2))
                    l1tmp = ll2
                    l2tmp = ll1
                    n1tmp = n2
                    n2tmp = n1
                    test = (l1only == 0 and l1tmp == 0) or l1only != 0
                    test = test and ((l2only == 0 and l2tmp == 0) or l2only != 0)
                    test = test and ((lonly == 0 and ll == 0) or lonly != 0)
                    test = test and ((n1only == 0 and n1tmp == 0) or n1only != 0)
                    test = test and ((n2only == 0 and n2tmp == 0) or n2only != 0)
                    if test:
                        phsn = f32(_pm1(l1tmp + l2tmp + ll + n1tmp + n2tmp))
                        if m2fix > 0:
                            m2f = -l2tmp
                            m2l = l2tmp
                        else:
                            m2f = 0
                            m2l = 0
                        sumc, sums = plotsharm_routines.sumoverm2(l1tmp, l2tmp, ll, n1tmp, m2f, m2l, n2tmp, phil, ctl,
                                                                  phi, cz, chi, phs, sumc, sums)
                        # now reverse signs of n1,n2 if either is non-zero
                        if n1 > 0 or (n1 == 0 and n2 > 0):
                            n1tmp = -n2
                            n2tmp = -n1
                            phs = phs * phsn
                            sumc, sums = plotsharm_routines.sumoverm2(l1tmp, l2tmp, ll, n1tmp, m2f, m2l, n2tmp, phil,
                                                                      ctl, phi, cz, chi, phs, sumc, sums)
                zcoeffc[k] = sumc
                if forcemirrorsym:
                    zcoeffs[k] = 0.0
                else:
                    zcoeffs[k] = sums
                #               write(6,*) ll1,ll2,ll,n1,n2,phil,ctl,phi,cz,sumc,sums
            for i in range(1, nr + 1):
                sumt = f32(0)
                for k in range(1, nexpan + 1):
                    sumt = sumt + datac[i, k] * zcoeffc[k] + datas[i, k] * zcoeffs[k]
                #
                # test for maximum and minimum in selected r range
                #
                r = rd[i]
                temp[i] = 0.0
                if r >= rmin:
                    if r <= rmax:
                        temp[i] = sumt + f32(1.0)
                    else:
                        temp[i] = 1.0
            #
            # perform smoothing on coefficients if needed - this is primarily to remove
            # hard edges and so avoid fringes in the maps
            #
            if nsmoo > 0:
                plotsharm_routines.smoosharm(1, nr, nsmoo, temp)
            for i in range(1, nr + 1):
                grw[it, ip, i] = temp[i]
                if grw[it, ip, i] < dlow:
                    dlow = grw[it, ip, i]
                if grw[it, ip, i] > dhigh:
                    dhigh = grw[it, ip, i]
                    if rd[i] >= rmin:
                        irmax = i
                        phimmax = phi
                        thetammax = thetam
                        chimmax = chi
    units.write(6, lw('Max. and min. values in this function: ', dhigh, dlow))
    if surfra >= 0.0:
        dsurf = dhigh + r4(0.9) * (dlow - dhigh)
    else:
        dsurf = -surfra
    #
    # form the sums of r**2 * g(r,w) between rmin and rmax, above DSURF, and in
    # total
    #
    fracold = f32(1.0)
    frac = f32(0.9)
    ntimes = 0
    while True:
        # 200 continue
        ntimes = ntimes + 1
        sumt = f32(0.0)
        sump = f32(0.0)
        delr = r4(4.0) * pi * (rd[3] - rd[2])
        for i in range(1, nr + 1):
            r = rd[i]
            if r >= rmin and r <= rmax:
                r2 = r * r * delr
                for ip in range(1, np0 + 1):
                    for it in range(1, nt0 + 1):
                        tempgr = grw[it, ip, i]
                        sumt = sumt + r2 * tempgr
                        if tempgr >= dsurf:
                            sump = sump + r2 * tempgr
        sumt = sumt / f32(nt0) / f32(np0)
        sump = sump / f32(nt0) / f32(np0)
        units.write(6, lw('Total integral between ', rmin, ' and ', rmax, sumt))
        units.write(6, lw('Integral for this map above ', dsurf, 'between ', rmin, ' and ', rmax, sump))
        if sumt > 0.0:
            srat = sump / sumt
            units.write(6, lw('Ratio = ', f32(100) * srat, '%'))
            #
            # if surfra is +ve it represents the fraction of molecules
            # at this distance
            # range which are to be plotted in the contour plot. Therefore adjust the
            # value of dsurf until the fraction of molecules plotted is at the specified
            # level
            #
            if surfra > 0.0:
                fact = srat / surfra
                if fact > r4(1.01) or fact < r4(0.99):
                    #
                    # estimate next value
                    #
                    if srat < surfra:
                        #
                        # increase frac
                        #
                        frac = r4(0.5) * (frac + fracold)
                    else:
                        fracold = frac
                        frac = r4(0.9) * fracold
                    dsurf = dhigh + frac * (dlow - dhigh)
                    #
                    # dsurf cannot be above dhigh or below dlow
                    #
                    if dsurf > dhigh:
                        dsurf = dhigh
                    if dsurf < dlow:
                        dsurf = dlow
                    if ntimes < 500:
                        continue  # go to 200
        break
    units.write(6, lw('Phi, theta, chi at maximum density, and maximum density: ', phimmax, thetammax, chimmax,
                      dhigh))
    phirot = phimmax
    thetarot = thetammax
    chirot = chimmax

    # For NVARY = 3 we need to rotate the central molecule about the z-axis by chirot degrees
    # because that is the way the orientational spatial density function is plotted, namely the
    # distribution out of the z-x plane represents the variation with chim.

    if nvary == 3:
        phirot = phirot + chirot
    #      write(6,*) phirot,thetarot,chirot,dhigh
    phirot = phirot * piconv
    thetarot = thetarot * piconv
    chirot = chirot * piconv
    # In case it is needed define the translation for the plot
    xtrans = f32(0.0)
    ytrans = f32(0.0)
    ztrans = f32(0.0)
    #
    # apply the rotation. Do this in reverse order about the lab. axes
    #
    if nvary > 1 and rotatecentralmolecule:
        # NB: irmax stays 0 if no maximum is found above rmin; rd(0) is then read out of
        # bounds in the Fortran (the unused element rd[0] = 0 here).
        ztrans = rd[irmax] * cosd(thetal)
        rtrans = rd[irmax] * sind(thetal)
        xtrans = rtrans * cosd(phil)
        ytrans = rtrans * sind(phil)
        for is_ in range(1, nspher + 1):
            # rotate about lab z axis
            phcoor[is_] = phcoor[is_] + chirot
            # rotate about lab y axis, then about z axis again
            z = cosf(thcoor[is_])
            x = sinf(thcoor[is_]) * cosf(phcoor[is_])
            y = sinf(thcoor[is_]) * sinf(phcoor[is_])
            d = sqrtf(x * x + z * z)
            alpha = atan2f(x, z) + thetarot
            znew = d * cosf(alpha)
            xnew = d * sinf(alpha)
            thcoor[is_] = acosf(znew)
            phcoor[is_] = atan2f(y, xnew) + phirot
    #
    # now interpolate onto a cubic grid and determine whether pixel is occupied or not
    #
    quadinterp = plotsharm_routines.quadinterp
    suminterp = plotsharm_routines.suminterp
    rdnr = rd[nr]
    iref = 0
    for k in range(1, nz0 + 1):
        z = f32(k - nz - 1) * zstep
        z2 = z * z
        for j in range(1, ny0 + 1):
            y = f32(j - ny - 1) * ystep
            #            yy(j)=y
            y2 = y * y
            for i in range(1, nx0 + 1):
                x = f32(i - nx - 1) * xstep
                #               xx(i)=x
                x2 = x * x
                rxy = sqrtf(x2 + y2)
                r = sqrtf(x2 + y2 + z2)
                # determine pair of rd values
                sum_ = f32(1.0)
                if r <= rdnr:
                    ir = nint(r / rstep) + 1
                    ir = max(2, ir)  # This is to prevent irm referencing undefined grw at low r
                    ir = min(nr - 1, ir)  # This is to prevent irp referencing undefined grw at high r
                    irp = ir + 1
                    irm = ir - 1
                    rrat = r / rstep - f32(ir - 1)
                    # determine nearest theta value
                    theta = f32(0.0)
                    if r > 0.0:
                        costh = z / r
                        costh = min(f32(1.0), costh)
                        costh = max(f32(-1.0), costh)
                        theta = acosd(costh)
                    nth = nint(theta / tstep) + 1
                    nth = max(2, nth)
                    nth = min(nt0 - 1, nth)
                    nthp = nth + 1
                    nthm = nth - 1
                    thrat = theta / tstep - f32(nth - 1)
                    # determine nearest phi value
                    phi = f32(0)
                    if rxy > 0.:
                        cosphi = x / rxy
                        sinphi = y / rxy
                        phi = atan2d(sinphi, cosphi)
                        if phi < 0.0:
                            phi = phi + f32(360.0)
                    nph = nint(phi / pstep) + 1
                    phrat = phi / pstep - f32(nph - 1)
                    while nph > np0:
                        nph = nph - np0
                    while nph < 1:
                        nph = nph + np0
                    # Try averaging over relevant phi values. We will try to keep the solid angle of this averaging constant
                    # for different theta.
                    thetatrial = f32(nth) * tstep
                    sintheta = sind(thetatrial)
                    if sintheta > 0.0:
                        nphtrial = nint(r4(0.5) / sintheta)
                    else:
                        nphtrial = np0 // 2
                    # if(nptrial.lt.1) nptrial=1   (nptrial, sic: an unused variable)
                    nphp = nph + nphtrial
                    nphm = nph - nphtrial
                    if nphtrial == 1:
                        while nphp > np0:
                            nphp = nphp - np0
                        while nphm < 1:
                            nphm = nphm + np0
                        # Perform quadratic interpolation
                        sum_ = quadinterp(grw, rrat, thrat, phrat, ir, irp, irm, nth, nthp, nthm, nph, nphp, nphm,
                                          50, 50, 500)
                    else:
                        # Form mean over phi values
                        sum_, nphp, nphm = suminterp(grw, rrat, thrat, ir, irp, irm, nth, nthp, nthm, np0, nphp, nphm,
                                                     50, 50, 500)
                dens[i, j, k] = sum_
                if sum_ >= dsurf:
                    iref = iref + 1
    nref = iref
    units.write(6, fmt(_F999, nref, nvoxel))
    nout = 10 + 1  # nout=nin+1
    # Translate the centre of the plot if requested
    xcentre = (xtrans - f32(nx) * xstep * scalefactor) / bohr
    ycentre = (ytrans - f32(ny) * ystep * scalefactor) / bohr
    zcentre = (ztrans - f32(nz) * zstep * scalefactor) / bohr
    # The molecular coordinates will be stored in a separate file and added to if more than one surface is to be plotted
    try:
        units.open(nout, molfile, status='old')
        ierr = 0
    except FortranIOError as e:
        ierr = e.iostat
    if ierr == 0 and not newstart and rotatecentralmolecule:
        # This file exists, this is not a newstart, we are not plotting another SDF, and ARE
        # rotating/translating the central molecule, so we add the new molecule to the
        # existing central molecule.
        molfileline1 = fstr(units.read_line(nout), 256)
        molfileline2 = fstr(units.read_line(nout), 256)
        r = units.list_read(nout)
        nspherold = r.int(0)
        xcentremol = r.real(f32(0.0))  # noqa: F841
        ycentremol = r.real(f32(0.0))  # noqa: F841
        zcentremol = r.real(f32(0.0))  # noqa: F841
        molfileline = char_array(nspherold + 3)
        for i in range(1, nspherold + 3 + 1):
            molfileline[i] = fstr(units.read_line(nout), 256)
        units.rewind(nout)
        units.write(nout, fmt('(a)', substr(molfileline1, 1, len(molfileline1))))
        units.write(nout, fmt('(a)', substr(molfileline2, 1, len(molfileline2))))
        units.write(nout, lw(nspherold + nspher, xcenterold, ycentreold, zcentreold))
        for i in range(1, nspherold + 3 + 1):
            units.write(nout, fmt('(a)', substr(molfileline[i], 1, len(molfileline[i]))))
        # Add in the new lines
        for is_ in range(1, nspher + 1):
            xbal = rdcoor[is_] * sinf(thcoor[is_]) * cosf(phcoor[is_]) + xtrans
            ybal = rdcoor[is_] * sinf(thcoor[is_]) * sinf(phcoor[is_]) + ytrans
            zbal = rdcoor[is_] * cosf(thcoor[is_]) + ztrans
            iatomicno = ar.get_atomicnumber(chsymbol[is_])
            units.write(nout, lw(iatomicno, f32(iatomicno), xbal / bohr, ybal / bohr, zbal / bohr))
        units.close(nout)
        molfileline = None
    elif newstart:
        try:
            units.open(nout, molfile, status='unknown')
            ierr = 0
        except FortranIOError as e:
            ierr = e.iostat
        _write_cube_header(nout, nspher, xcentre, ycentre, zcentre, 0, xstep, ystep, zstep, scalefactor, bohr)
        _write_atoms(nout, nspher, rdcoor, thcoor, phcoor, chsymbol, xtrans, ytrans, ztrans, bohr)
        units.close(nout)
    # Open the Cube output file. Units need to be in Bohr atomic radii.
    units.open(nout, outfile, status='unknown')
    if scalefactor > 0.0:
        _write_cube_header(nout, nspher, xcentre, ycentre, zcentre, (nx0, ny0, nz0), xstep, ystep, zstep,
                           scalefactor, bohr)
    else:
        _write_cube_header(nout, nspher, xcentre, ycentre, zcentre, 0, xstep, ystep, zstep, scalefactor, bohr)
    _write_atoms(nout, nspher, rdcoor, thcoor, phcoor, chsymbol, xtrans, ytrans, ztrans, bohr)
    if scalefactor > 0.0:
        for ix in range(1, nx0 + 1):
            for iy in range(1, ny0 + 1):
                for iz in range(1, nz0 + 1, 10):
                    itmax = min(9, nz0 - iz)
                    units.write(nout, fmt('(10(1x,e13.6))', *[dens[ix, iy, iz + it] for it in range(0, itmax + 1)]))
    units.close(nout)
    try:
        if newstart:
            units.open(nout, plotfile, status='unknown')
        else:
            units.open(nout, plotfile, status='old', position='append')
        ierr = 0
    except FortranIOError as e:
        ierr = e.iostat
    if ierr == 0:
        if newstart:
            units.write(nout, fmt(_F795, substr(molfile, 1, _save.lenmolfile), rgbbak[1], rgbbak[2], rgbbak[3],
                                  translucency))
            units.write(nout, fmt(_F807, rot, elev, Ch(spin, 3)))
        if scalefactor > 0.0:
            units.write(nout, fmt(_F797, substr(surfacelabel, 1, len(surfacelabel)), dsurf,
                                  substr(outfile, 1, lenoutfile), rgbobj[1], rgbobj[2], rgbobj[3],
                                  substr(surfacetype, 1, len(surfacetype))))
    units.close(nout)

    # Save random number generator

    units.open(10, 'ran1.dat', status='unknown')
    auxiliary_routines.saveran1(10)
    units.close(10)
    return plotfile, nextfile


def _write_cube_header(nout, nspher, xcentre, ycentre, zcentre, ngrid, xstep, ystep, zstep, scalefactor, bohr):
    """The six header records of a cube file, as map3djmol writes them (three times in
    the Fortran): ngrid is 0 or (nx0, ny0, nz0)."""
    n1, n2, n3 = (0, 0, 0) if ngrid == 0 else ngrid
    units.write(nout, lw('# Cube file output by map3djmol'))
    units.write(nout, lw('Outer loop x, middle loop y, inner loop z'))
    units.write(nout, lw(nspher, xcentre, ycentre, zcentre))
    units.write(nout, lw(n1, xstep * scalefactor / bohr, f32(0.0), f32(0.0)))
    units.write(nout, lw(n2, f32(0.0), ystep * scalefactor / bohr, f32(0.0)))
    units.write(nout, lw(n3, f32(0.0), f32(0.0), zstep * scalefactor / bohr))


def _write_atoms(nout, nspher, rdcoor, thcoor, phcoor, chsymbol, xtrans, ytrans, ztrans, bohr):
    """The atom records of a cube file (the ``do is=1,nspher`` loop of map3djmol)."""
    for is_ in range(1, nspher + 1):
        xbal = rdcoor[is_] * sinf(thcoor[is_]) * cosf(phcoor[is_]) + xtrans
        ybal = rdcoor[is_] * sinf(thcoor[is_]) * sinf(phcoor[is_]) + ytrans
        zbal = rdcoor[is_] * cosf(thcoor[is_]) + ztrans
        iatomicno = auxiliary_routines.get_atomicnumber(chsymbol[is_])
        units.write(nout, lw(iatomicno, f32(iatomicno), xbal / bohr, ybal / bohr, zbal / bohr))
