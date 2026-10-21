"""Translation of fortran_src/sharm_rev2.f: ``program standalone`` (sharm),
``subroutine run_sharm`` and ``subroutine sharm`` -- spherical harmonic coefficients of
the orientational pair correlation function of two molecule types, from a .ato file.

The molecular centre is the mass-weighted centre of the listed atom types.  The shared
routines are in library/sharm_routines_rev1.py, the data of ``module sharm_mod`` in
common/sharm_mod.py.  ``setup_symbols_sharm`` is in setup_symbols_generic.f
(library/setup_symbols_generic.py, imported inside the routine: another work package).

The run-time FORMAT of the column titles is the one of the reference build, which
patches ``'(1h('//il1fmt...`` into ``'(1h(,'//il1fmt...`` (tools/build_fortran_reference.sh;
current gfortran rejects the original at run time).  See docs/translation_notes/Q.md.
"""

from ..common import clustercom, commands, lattice, molecom, sharm_mod, sharm_sdf, symbols
from ..common import sharm as sharm_inc  # sharm.inc (the name sharm is the routine below)
from ..fortran import (
    PI4ATAN,
    Ch,
    FortranEOF,
    FortranIOError,
    acosf,
    char_array,
    cosf,
    f32,
    fmt,
    fstr,
    int_array,
    iostat_of,
    lw,
    r4,
    read_internal,
    real_array,
    substr,
    units,
)
from ..library import (
    auxiliary_routines,
    file_and_directory_commands,
    parse,
    read_and_save_files,
    reallocation_routines,
    sharm_routines_rev1,
)
from ..library import ato_read_and_save
from ..trace import routine

_F201 = "(/'sharm> Atom type ',a,' does not exist in .ato file ',a)"
_F333 = "(1x,'sharm> Number of type 1 molecules = ',i10/1x, 'sharm> Number of type 2 molecules = ',i10)"
_FXYZ = "(a3,1x,f10.5,1x,f10.5,1x,f10.5,1x,a3)"


def main():
    """sharm_rev2.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_sharm()


def _parse(iref, nspc):
    """call parse(symbol_value(iref),nwords,ncf,ncl,nspc,256,nworddim)"""
    commands.nwords = parse.parse(symbols.symbol_value[iref], commands.nwords, commands.ncf, commands.ncl, nspc,
                                  256, commands.nworddim)


def _word(iref, ic):
    """symbol_value(iref)(ncf(ic):ncl(ic))"""
    return substr(symbols.symbol_value[iref], commands.ncf[ic], commands.ncl[ic])


@routine
def run_sharm():
    """sharm_rev2.f:run_sharm -- performs calculation of spherical harmonic coefficients
    from .ato file.

    uses arithmetic centre of molecule rather than a particular atom to define
    molecular centre.
    """
    from ..library import setup_symbols_generic

    S = sharm_mod
    sv = symbols.symbol_value
    # character*50 formatout
    # character*2 ifmt1,ifmt2,ifmt3,ifmt4,il1fmt,il2fmt,ilfmt,in1fmt,in2fmt
    il1fmt = il2fmt = ilfmt = in1fmt = in2fmt = ''
    #
    # setup the SHARM symbols
    #
    setup_symbols_generic.setup_symbols_sharm()
    #
    # now setup the input file
    #
    nin = 11
    nout = 6
    # NB: ierrout is not initialised in the Fortran; setup_input_file sets it
    ierrout = 0
    nout, ierrout = file_and_directory_commands.setup_input_file(nin, nout, ierrout)
    #
    # if an error occurred on setup, exit now
    #
    if ierrout != 0:
        return

    sharm_sdf.pi = PI4ATAN  # 4.0*atan(1.0)
    sharm_sdf.piconv = sharm_sdf.pi / r4(180.0)
    #
    # read the input file
    #
    nout, ierrout = read_and_save_files.read_input_file(nin, nout, ierrout)
    #
    # if an error occurred on read, exit now
    #
    if ierrout != 0:
        return
    # NB: ierr is only set by the iostat= reads below; 0 if none of them is executed
    ierr = 0
    iref = 1
    #
    # set up number of radius values
    #
    iref = iref + 1
    sharm_inc.nz = read_internal(sv[iref]).int(sharm_inc.nz)
    #
    # set up rmax of calculation
    #
    iref = iref + 1
    sharm_inc.zmax = read_internal(sv[iref]).real(sharm_inc.zmax)
    #
    # set up number of configurations already stored
    #
    iref = iref + 1
    sharm_inc.nsumt = read_internal(sv[iref]).int(sharm_inc.nsumt)
    #
    # set up number of coefficients
    #
    iref = iref + 1
    S.ndist = read_internal(sv[iref]).int(S.ndist)
    #
    # set up l1 values
    #
    iref = iref + 1
    _parse(iref, 1)
    S.nll1v = commands.nwords
    if S.ll1v is not None:
        S.ll1v = None
    S.ll1v = int_array(S.nll1v)
    for i in range(1, S.nll1v + 1):
        try:
            S.ll1v[i] = read_internal(_word(iref, i)).int(S.ll1v[i])
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = iostat_of(e)
    #
    # set up l2 values
    #
    iref = iref + 1
    _parse(iref, 1)
    S.nll2v = commands.nwords
    if S.ll2v is not None:
        S.ll2v = None
    S.ll2v = int_array(S.nll2v)
    for i in range(1, S.nll2v + 1):
        try:
            S.ll2v[i] = read_internal(_word(iref, i)).int(S.ll2v[i])
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = iostat_of(e)
    #
    # set up l values
    #
    iref = iref + 1
    _parse(iref, 5)
    S.nllv = commands.nwords
    if S.llv is not None:
        S.llv = None
    S.llv = int_array(S.nllv)
    for i in range(1, S.nllv + 1):
        S.llv[i] = read_internal(_word(iref, i)).int(S.llv[i])
    #
    # set up n1step
    #
    iref = iref + 1
    S.n1step = read_internal(sv[iref]).int(S.n1step)
    #
    # set up n2step
    #
    iref = iref + 1
    S.n2step = read_internal(sv[iref]).int(S.n2step)
    #
    # for molecule 1 read number and type of atoms used to define the molecular
    # centre
    #
    iref = iref + 1
    _parse(iref, 5)
    sharm_inc.ncent1 = 0
    #
    # skip through the atom types to find out the type number of the specified atom
    #
    for ic in range(1, commands.nwords + 1):
        it = 0
        ifound = 0
        while it < molecom.ntype and ifound == 0:
            it = it + 1
            if molecom.type[it] == _word(iref, ic):
                ifound = it
        if ifound <= 0:
            units.write(6, fmt(_F201, _word(iref, ic), substr(symbols.fname_ato, 1, symbols.ilenfname_ato)))
            return
        else:
            if sharm_inc.ncent1 < sharm_inc.mcentre:
                sharm_inc.ncent1 = sharm_inc.ncent1 + 1
                sharm_inc.lctyp1[sharm_inc.ncent1] = ifound
    #
    # for molecule 1 read primary axis type (1=X, 2=Y, 3=Z), and
    # the number and atom numbers of atoms used to define the molecule's nzt1
    # axis
    #
    iref = iref + 1
    _parse(iref, 5)
    sharm_inc.nzt1 = fstr(_word(iref, 1), 2)
    sharm_inc.ncz1 = commands.nwords - 1
    if sharm_inc.ncz1 > sharm_inc.mcentre:
        sharm_inc.ncz1 = sharm_inc.mcentre
    for i in range(1, sharm_inc.ncz1 + 1):
        ic = i + 1
        sharm_inc.lztyp1[i] = read_internal(_word(iref, ic)).int(sharm_inc.lztyp1[i])
    #
    # for molecule 1 read secondary axis type (1=X, 2=Y, 3=Z), and
    # the number and atom numbers of atoms used to define the
    # molecule's nyt1 axis
    #
    iref = iref + 1
    _parse(iref, 5)
    sharm_inc.nyt1 = fstr(_word(iref, 1), 2)
    sharm_inc.ncy1 = commands.nwords - 1
    if sharm_inc.ncy1 > sharm_inc.mcentre:
        sharm_inc.ncy1 = sharm_inc.mcentre
    for i in range(1, sharm_inc.ncy1 + 1):
        ic = i + 1
        sharm_inc.lytyp1[i] = read_internal(_word(iref, ic)).int(sharm_inc.lytyp1[i])

    # for molecule 2 read number and type of atoms used to define the molecular
    # centre. Also check whether this molecule is identical to the first

    iref = iref + 1
    _parse(iref, 5)
    sharm_inc.ncent2 = 0
    #
    # skip through the atom types to find out the type number of the specified atom
    #
    for ic in range(1, commands.nwords + 1):
        it = 0
        ifound = 0
        while it < molecom.ntype and ifound == 0:
            it = it + 1
            if molecom.type[it] == _word(iref, ic):
                ifound = it
        if ifound == 0:
            units.write(6, fmt(_F201, _word(iref, ic), substr(symbols.fname_ato, 1, symbols.ilenfname_ato)))
            return
        else:
            if sharm_inc.ncent2 < sharm_inc.mcentre:
                sharm_inc.ncent2 = sharm_inc.ncent2 + 1
                sharm_inc.lctyp2[sharm_inc.ncent2] = ifound

    # test for identity of central and second molecules

    S.ldistinct = 0
    if sharm_inc.ncent2 == sharm_inc.ncent1:
        S.ldistinct = 0
        for it in range(1, sharm_inc.ncent1 + 1):
            if sharm_inc.lctyp2[it] != sharm_inc.lctyp1[it]:
                S.ldistinct = 1
    else:
        S.ldistinct = 1

    units.write(6, lw(sharm_inc.ncent1, sharm_inc.ncent2, S.ldistinct))

    if S.ldistinct != 0:

        # if molecules are different, then repeat this for the second molecule

        # for molecule 2 read primary axis type (1=X, 2=Y, 3=Z), and
        # the number and atom numbers of atoms used to define the molecule's nzt1
        # axis

        iref = iref + 1
        _parse(iref, 5)
        sharm_inc.nzt2 = fstr(_word(iref, 1), 2)
        sharm_inc.ncz2 = commands.nwords - 1
        if sharm_inc.ncz2 > sharm_inc.mcentre:
            sharm_inc.ncz2 = sharm_inc.mcentre
        for i in range(1, sharm_inc.ncz2 + 1):
            ic = i + 1
            sharm_inc.lztyp2[i] = read_internal(_word(iref, ic)).int(sharm_inc.lztyp2[i])

        # for molecule 2 read secondary axis type (1=X, 2=Y, 3=Z), and
        # the number and atom numbers of atoms used to define the
        # molecule's nyt1 axis

        iref = iref + 1
        _parse(iref, 5)
        sharm_inc.nyt2 = fstr(_word(iref, 1), 2)
        sharm_inc.ncy2 = commands.nwords - 1
        if sharm_inc.ncy2 > sharm_inc.mcentre:
            sharm_inc.ncy2 = sharm_inc.mcentre
        for i in range(1, sharm_inc.ncy2 + 1):
            ic = i + 1
            sharm_inc.lytyp2[i] = read_internal(_word(iref, ic)).int(sharm_inc.lytyp2[i])
    else:

        # otherwise set all the values to the same as the central molecule

        sharm_inc.ncent2 = sharm_inc.ncent1
        sharm_inc.nzt2 = sharm_inc.nzt1
        sharm_inc.nyt2 = sharm_inc.nyt1
        if sharm_inc.ncent2 > 0:
            for ic in range(1, sharm_inc.ncent1 + 1):
                sharm_inc.lctyp2[ic] = sharm_inc.lctyp1[ic]
        sharm_inc.ncz2 = sharm_inc.ncz1
        for i in range(1, sharm_inc.ncz1 + 1):
            sharm_inc.lztyp2[i] = sharm_inc.lztyp1[i]
        sharm_inc.ncy2 = sharm_inc.ncy1
        for i in range(1, abs(sharm_inc.ncy1) + 1):
            sharm_inc.lytyp2[i] = sharm_inc.lytyp1[i]
        iref = iref + 2
    # c     write(6,*) ncent1,ncent2,ncz1,ncy1,ncz2,ncy2
    # Decide whether to use simulated phim, thetam, chim
    iref = iref + 1
    _parse(iref, 5)
    S.performtest = _word(iref, 1) != 'notest' and commands.nwords > 7
    if S.performtest:
        ierr = _read_test_angles(substr(sv[iref], commands.ncf[1], commands.ncl[commands.nwords]))
        p = sharm_sdf.piconv
        sharm_sdf.phim1 = sharm_sdf.phim1 * p
        sharm_sdf.thetam1 = sharm_sdf.thetam1 * p
        sharm_sdf.chim1 = sharm_sdf.chim1 * p
        sharm_sdf.phim2 = sharm_sdf.phim2 * p
        sharm_sdf.thetam2 = sharm_sdf.thetam2 * p
        sharm_sdf.chim2 = sharm_sdf.chim2 * p
        sharm_sdf.thetallab = sharm_sdf.thetallab * p
        sharm_sdf.phillab = sharm_sdf.phillab * p
    else:
        sharm_sdf.phim1 = f32(0.0)
        sharm_sdf.thetam1 = f32(0.0)
        sharm_sdf.chim1 = f32(0.0)
        sharm_sdf.phim2 = f32(0.0)
        sharm_sdf.thetam2 = f32(0.0)
        sharm_sdf.chim2 = f32(0.0)
        sharm_sdf.thetallab = f32(0.0)
        sharm_sdf.phillab = f32(0.0)
    if ierr != 0:
        S.performtest = False
    iref = iref + 1
    S.outputaveragemol = sv[iref] == 'y'

    # Only perform the calculation when nsum is zero or a multiple of nmult

    nmult = 5
    nsum = sharm_inc.nsumt // nmult
    nt = nsum * nmult
    if nt == sharm_inc.nsumt:
        #
        # normalise atomic coordinates to box size
        #
        ato_read_and_save.normalise_coordinates()
        sharm_inc.zmaxd = sharm_inc.zmax / clustercom.dmax
        if sharm_inc.zmaxd > 0.5:
            sharm_inc.zmaxd = f32(0.5)
            sharm_inc.zmax = sharm_inc.zmaxd * clustercom.dmax
        # Set up l and n values and allocate some needed arrays
        sharm_routines_rev1.get_l_and_n_values()
        S.mz = sharm_inc.nz
        S.mdist = S.ndist
        if S.z is not None:
            S.z = S.gofz = S.gofzs = S.gzt = S.gzt2 = S.titledat = None
        S.z = real_array(S.mz)
        S.gofz = real_array(S.mz, S.mdist)
        S.gofzs = real_array(S.mz, S.mdist)
        S.gzt = real_array(S.mz, S.mdist)
        S.gzt2 = real_array(S.mz, S.mdist)
        S.titledat = char_array(S.mdist)
        # set up step size and z values
        sharm_inc.zstep = sharm_inc.zmaxd / f32(sharm_inc.nz)
        sharm_inc.zhalf = r4(0.5) * sharm_inc.zstep
        for j in range(1, sharm_inc.nz + 1):
            S.z[j] = sharm_inc.zstep * f32(j) - sharm_inc.zhalf
        # zero gofz
        read_and_save_files.initgofzs(sharm_inc.nz, S.ndist, S.gofz, S.gofzs, S.mz, S.mdist)
        sharm(nin, symbols.fname_inp)
        #
        # set the z values
        #
        for i in range(1, sharm_inc.nz + 1):
            S.z[i] = S.z[i] * clustercom.dmax
        #
        # axis title
        #
        sharm_inc.titler = '   r'
        #
        # possible formats for integers in title
        #
        ifmt1 = 'i1'
        ifmt2 = 'i2'
        ifmt3 = 'i3'
        ifmt4 = 'i4'
        #
        # set the column titles
        #
        for id_ in range(1, S.ndist + 1):
            il1 = int(S.ll1[id_])
            il2 = int(S.ll2[id_])
            il = int(S.ll[id_])
            in1 = int(S.n1[id_])
            in2 = int(S.n2[id_])
            #
            # set up the format for the title
            #
            if il1 < 10:
                il1fmt = ifmt1
            elif il1 < 100:
                il1fmt = ifmt2
            elif il1 < 1000:
                il1fmt = ifmt3
            if il2 < 10:
                il2fmt = ifmt1
            elif il2 < 100:
                il2fmt = ifmt2
            elif il2 < 1000:
                il2fmt = ifmt3
            if il < 10:
                ilfmt = ifmt1
            elif il < 100:
                ilfmt = ifmt2
            elif il < 1000:
                ilfmt = ifmt3
            if in1 < 10:
                in1fmt = ifmt1
            elif in1 < 100:
                in1fmt = ifmt2
            elif in1 < 1000:
                in1fmt = ifmt3
            if abs(in2) < 10:
                if in2 < 0:
                    in2fmt = ifmt2
                else:
                    in2fmt = ifmt1
            elif abs(in2) < 100:
                if in2 < 0:
                    in2fmt = ifmt3
                else:
                    in2fmt = ifmt2
            elif abs(in2) < 1000:
                if in2 < 0:
                    in2fmt = ifmt4
                else:
                    in2fmt = ifmt3
            # Patched as in the reference build: '(1h(' -> '(1h(,' (tools/build_fortran_reference.sh)
            formatout = fstr('(1h(,' + il1fmt + ',1h,,' + il2fmt + ',1h,,' + ilfmt + ',1h;,' + in1fmt + ',1h,,'
                             + in2fmt + ',1h))', 50)
            #
            # write the title for this distribution
            #
            S.titledat[id_] = fstr(fmt(formatout, il1, il2, il, in1, in2), 29)
        sharm_inc.ext = 'h'
        nsum = read_and_save_files.outgofzs_sharm(symbols.fnameinp, sharm_inc.ext, sharm_inc.nz, S.ndist, nsum, S.z, S.gofz,
                                                  S.gofzs, S.gzt, S.gzt2, sharm_inc.titler, S.titledat, S.mz, S.mdist)
    if not S.performtest:
        sharm_inc.nsumt = sharm_inc.nsumt + 1
    #
    # write input file
    #
    #
    # make sure all current values are stored
    #
    sv[4] = fstr(lw(sharm_inc.nsumt), 256)
    sv[5] = fstr(lw(S.ndist), 256)
    units.write(6, lw(S.ndist, ' ', Ch(sv[5], 256)))
    read_and_save_files.save_files_others(14)
    units.close(14)
    if S.performtest:
        units.write(6, lw('sharm> Test performed.'))
    return


def _read_test_angles(text):
    """read(text,*,iostat=ierr) phim1,thetam1,chim1,phim2,thetam2,chim2,thetallab,phillab
    (the sharm_sdf.inc variables); the items read before an error keep their new values.
    Returns ierr."""
    names = ('phim1', 'thetam1', 'chim1', 'phim2', 'thetam2', 'chim2', 'thetallab', 'phillab')
    r = read_internal(text)
    try:
        for name in names:
            setattr(sharm_sdf, name, r.real(getattr(sharm_sdf, name)))
    except (FortranEOF, FortranIOError) as e:
        return iostat_of(e)
    return 0


def _write_xyz(nout, fname_xyz, natocentral, atsym, xs, ys, zs, attype, centraldef):
    """The coordinate-file block of sharm (repeated three times in the Fortran):
    open(nout,...,iostat=ierr); if(ierr.eq.0) then ... close(nout)."""
    try:
        units.open(nout, fname_xyz, status='unknown')
        ierr = 0
    except FortranIOError as e:
        ierr = iostat_of(e)
    if ierr == 0:
        # Write number of atoms in molecule
        units.write(nout, fmt('(i6)', natocentral))
        # Write information line - not used for plotting
        units.write(nout, fmt('(a)', 'Central molecule coordinate file: ' + fname_xyz))
        for ia in range(1, natocentral + 1):
            textline = fstr(fmt(_FXYZ, Ch(atsym[ia], 3), xs[ia], ys[ia], zs[ia], Ch(attype[ia], 3)), 256)
            # Append a comment if this atom is part of the central molecule definition
            if centraldef[ia]:
                textline = fstr(textline + '  #CENTRALMOLECULE', 256)
            units.write(nout, fmt('(a)', textline))
        units.close(nout)
    return ierr


@routine
def sharm(nout, fname_inp):
    """sharm_rev2.f:sharm -- set up the molecular centres and Euler angles, write the
    averaged first and second molecules, and sum the coefficients
    (form_sum_prod_dlmn_cgcoeff).

    Returns nout.
    """
    S = sharm_mod
    M = molecom
    C = clustercom
    # real amass,amasssum
    # logical :: centralmolecule,found,foundatom
    # character(len=256) :: fname_xyz,textline
    xacentral = yacentral = zacentral = None  # real, dimension(:), allocatable
    atsymbolcentral = attypecentral = None  # character(len=3), dimension(:), allocatable
    centraldefinition = None  # logical, dimension(:), allocatable
    xat = yat = zat = None  # real, dimension(:), allocatable
    natocentral = 0

    sharm_sdf.pi = PI4ATAN  # 4.0*atan(1.0)
    sharm_sdf.piconv = sharm_sdf.pi / r4(180.0)
    degconv = r4(180.0) / sharm_sdf.pi  # noqa: F841
    sharm_inc.zstep = sharm_inc.zmaxd / f32(sharm_inc.nz)
    gfac = r4(4.0) * sharm_sdf.pi / r4(3.0)  # noqa: F841
    sharm_inc.vol = lattice.alatx * (lattice.blaty * lattice.clatz)
    # set up counter for second atoms
    nsec = 0  # noqa: F841
    # Arrays to store Euler angles
    if S.nm is not None:
        S.nm = S.phim = S.cthetm = S.chim = S.xcent = S.ycent = S.zcent = None
        S.dlmnmc = S.chiphi = S.dlmnc = S.philab = None
    sharm_routines_rev1.get_max_size_dlmn()
    nmol = C.nmol
    S.nm = int_array(nmol)
    S.phim = real_array(nmol)
    S.cthetm = real_array(nmol)
    S.chim = real_array(nmol)
    S.xcent = real_array(nmol)
    S.ycent = real_array(nmol)
    S.zcent = real_array(nmol)
    S.dlmnmc = real_array(2, nmol, S.mdlmn)
    S.chiphi = real_array(2, nmol, S.mdlmn)
    S.dlmnc = real_array(2, S.mdlmn)
    S.philab = real_array(S.mdlmn)
    if S.performtest:
        # Allocate suitable arrays if performing a test
        if S.phimt1 is not None:
            S.phimt1 = S.cthetmt1 = S.chimt1 = S.phimt2 = S.cthetmt2 = S.chimt2 = S.dlmnmct = S.chiphit = None
        S.phimt1 = real_array(nmol)
        S.cthetmt1 = real_array(nmol)
        S.chimt1 = real_array(nmol)
        S.phimt2 = real_array(nmol)
        S.cthetmt2 = real_array(nmol)
        S.chimt2 = real_array(nmol)
        S.dlmnmct = real_array(2, nmol, S.mdlmn)
        S.chiphit = real_array(2, nmol, S.mdlmn)
    #
    # set up coordinates of Laboratory X, Y and Z axes
    #
    alab = sharm_sdf.alab
    blab = sharm_sdf.blab
    clab = sharm_sdf.clab
    alab[1] = 1.0
    alab[2] = 0.0
    alab[3] = 0.0
    blab[1] = 0.0
    blab[2] = 1.0
    blab[3] = 0.0
    clab[1] = 0.0
    clab[2] = 0.0
    clab[3] = 1.0
    #
    # set up coordinates of molecular centres and Euler angles for each molecule
    # in laboratory coordinate frame
    #
    nsmmo1 = 0
    nsmmo2 = 0
    for imol1 in range(1, nmol + 1):
        # define origin of this molecule
        xo1 = f32(0.0)
        yo1 = f32(0.0)
        zo1 = f32(0.0)
        amasssum = f32(0.0)
        ncref = 0  # noqa: F841
        S.nm[imol1] = 0
        #
        # test to see if this is molecule 1 or molecule 2. Need to check that all the specified
        # types occur - not just one
        #
        found = True
        ic = 0
        while ic < sharm_inc.ncent1 and found:
            ic = ic + 1
            foundatom = False
            ia1 = C.naf[imol1] - 1
            while ia1 < C.nal[imol1] and not foundatom:
                ia1 = ia1 + 1
                foundatom = M.itype[ia1] == sharm_inc.lctyp1[ic]
            if ic == 1:
                found = foundatom
            else:
                found = found and foundatom
        if found:
            S.nm[imol1] = 1
        else:  # molecule 1 not found, try molecule 2
            found = True
            ic = 0
            while ic < sharm_inc.ncent2 and found:
                ic = ic + 1
                foundatom = False
                ia1 = C.naf[imol1] - 1
                while ia1 < C.nal[imol1] and not foundatom:
                    ia1 = ia1 + 1
                    foundatom = M.itype[ia1] == sharm_inc.lctyp2[ic]
                if ic == 1:
                    found = foundatom
                else:
                    found = found and foundatom
            if found:
                S.nm[imol1] = 2
        if S.nm[imol1] == 1:
            nsmmo1 = nsmmo1 + 1
        if S.nm[imol1] == 2:
            nsmmo2 = nsmmo2 + 1
        if S.nm[imol1] == 1 and sharm_inc.ncent1 > 0:
            for ia1 in range(C.naf[imol1], C.nal[imol1] + 1):
                for ic in range(1, sharm_inc.ncent1 + 1):
                    if M.itype[ia1] == sharm_inc.lctyp1[ic]:
                        amass = M.amassa[M.itype[ia1]]
                        xo1 = xo1 + M.xa[ia1] * amass
                        yo1 = yo1 + M.ya[ia1] * amass
                        zo1 = zo1 + M.za[ia1] * amass
                        amasssum = amasssum + amass
        elif S.nm[imol1] == 2 and sharm_inc.ncent2 > 0:
            for ia1 in range(C.naf[imol1], C.nal[imol1] + 1):
                for ic in range(1, sharm_inc.ncent2 + 1):
                    if M.itype[ia1] == sharm_inc.lctyp2[ic]:
                        amass = M.amassa[M.itype[ia1]]
                        xo1 = xo1 + M.xa[ia1] * amass
                        yo1 = yo1 + M.ya[ia1] * amass
                        zo1 = zo1 + M.za[ia1] * amass
                        amasssum = amasssum + amass
        if amasssum > 0.0 and S.nm[imol1] != 0:
            xo1 = xo1 / amasssum
            yo1 = yo1 / amasssum
            zo1 = zo1 / amasssum
            # c     write(6,991) xo1,yo1,zo1
            S.xcent[imol1] = C.xadd[imol1] + xo1
            S.ycent[imol1] = C.yadd[imol1] + yo1
            S.zcent[imol1] = C.zadd[imol1] + zo1
            # Form atomic coordinates relative to this centre
            xat = reallocation_routines.reallocate1d_r(xat, C.nato[imol1])
            yat = reallocation_routines.reallocate1d_r(yat, C.nato[imol1])
            zat = reallocation_routines.reallocate1d_r(zat, C.nato[imol1])
            for ia in range(1, C.nato[imol1] + 1):
                iaref = ia + C.naf[imol1] - 1
                xat[ia] = M.xa[iaref] - xo1
                yat[ia] = M.ya[iaref] - yo1
                zat[ia] = M.za[iaref] - zo1
            if S.nm[imol1] == 1:
                S.phim[imol1], S.cthetm[imol1], S.chim[imol1] = sharm_routines_rev1.form_eulerangles(
                    imol1, xat, yat, zat, sharm_inc.nzt1, sharm_inc.nyt1, sharm_inc.ncz1, sharm_inc.ncy1, sharm_inc.lztyp1, sharm_inc.lytyp1,
                    alab, blab, clab, S.phim[imol1], S.cthetm[imol1], S.chim[imol1])
            elif S.nm[imol1] == 2:
                S.phim[imol1], S.cthetm[imol1], S.chim[imol1] = sharm_routines_rev1.form_eulerangles(
                    imol1, xat, yat, zat, sharm_inc.nzt2, sharm_inc.nyt2, sharm_inc.ncz2, sharm_inc.ncy2, sharm_inc.lztyp2, sharm_inc.lytyp2,
                    alab, blab, clab, S.phim[imol1], S.cthetm[imol1], S.chim[imol1])
            # If testing orientations, supply a fictitious but known distribution of angles
            if S.performtest:
                #  write(6,499) imol1,phim1,thetam1,chim1,phim2,thetam2,chim2
                #  call populate_disk(pi,phim1,thetam1,chim1,phinew,thetanew,chinew,thetaspread,phispread)
                # Rotate this molecule onto desired orientation
                S.phimt1[imol1] = sharm_sdf.phim1
                S.cthetmt1[imol1] = cosf(sharm_sdf.thetam1)
                S.chimt1[imol1] = sharm_sdf.chim1
                #  call populate_disk(pi,phim2,thetam2,chim2,phinew,thetanew,chinew,thetaspread,phispread)
                S.phimt2[imol1] = sharm_sdf.phim2
                S.cthetmt2[imol1] = cosf(sharm_sdf.thetam2)
                S.chimt2[imol1] = sharm_sdf.chim2
                # 499   format(1x,i5,6(1x,e12.5))
            ierr = 0
            ierr = sharm_routines_rev1.form_dlmn(imol1, ierr)
            if ierr != 0:
                units.write(6, lw('sharm> Error when calculating dlm coefficients for molecule ', imol1))
                return nout
            #
            # step through all the atoms and ignore those which do not match the required
            # ltype's
        else:
            S.nm[imol1] = 0
    units.write(6, fmt(_F333, nsmmo1, nsmmo2))
    #      stop
    # Write out coordinates of centralmolecule. We will use the average coordinates of all
    # these molecules.
    im = 0
    ncentral = 0
    centralmolecule = False
    while im < nmol:
        im = im + 1
        centralmolecule = S.nm[im] == 1
        if centralmolecule:  # Write out coordinates of this molecule
            if ncentral == 0 or S.outputaveragemol:
                ncentral = ncentral + 1
                if ncentral == 1:  # First time allocate the arrays to store summed coordinates
                    natocentral = int(C.nato[im])
                    xacentral = real_array(natocentral)
                    yacentral = real_array(natocentral)
                    zacentral = real_array(natocentral)
                    atsymbolcentral = char_array(natocentral)
                    attypecentral = char_array(natocentral)
                    centraldefinition = [False] * (natocentral + 1)
                # Shift coordinates to centre used by SHARM/SDF
                xshift = C.xadd[im] - S.xcent[im]
                yshift = C.yadd[im] - S.ycent[im]
                zshift = C.zadd[im] - S.zcent[im]
                for ia in range(C.naf[im], C.nal[im] + 1):
                    # Temporary coordinates
                    xt = M.xa[ia] + xshift
                    yt = M.ya[ia] + yshift
                    zt = M.za[ia] + zshift
                    # Need to rotate these to align with the laboratory axes - performed in reverse
                    # order about the laboratory axes
                    if C.nato[im] > 1:
                        xr, yr, zr = auxiliary_routines.rotateeuler(xt, yt, zt, -S.chim[im], -acosf(S.cthetm[im]),
                                                                    -S.phim[im], None, None, None)
                    else:
                        xr = xt
                        yr = yt
                        zr = zt
                    iaref = ia - C.naf[im] + 1
                    if ncentral == 1:
                        xacentral[iaref] = xr
                        yacentral[iaref] = yr
                        zacentral[iaref] = zr
                        atsymbolcentral[iaref] = M.atsymbol[M.itype[ia]]
                        attypecentral[iaref] = M.type[M.itype[ia]]
                        # Check whether this atom is part of the central molecule definition
                        centraldefinition[iaref] = _central_definition(ia, iaref)
                    elif S.outputaveragemol:
                        xacentral[iaref] = xacentral[iaref] + xr
                        yacentral[iaref] = yacentral[iaref] + yr
                        zacentral[iaref] = zacentral[iaref] + zr
    # Form average coordinates for central molecule
    if ncentral > 0:
        for ia in range(1, natocentral + 1):
            xacentral[ia] = xacentral[ia] * C.dmax / f32(ncentral)
            yacentral[ia] = yacentral[ia] * C.dmax / f32(ncentral)
            zacentral[ia] = zacentral[ia] * C.dmax / f32(ncentral)
        fname_xyz = fstr(substr(fname_inp, 1, auxiliary_routines.indexlastdot(fname_inp) - 1) + '-first.xyz', 256)
        _write_xyz(nout, fname_xyz, natocentral, atsymbolcentral, xacentral, yacentral, zacentral, attypecentral,
                   centraldefinition)
        if S.ldistinct == 0:  # Also write this out as the second molecule
            fname_xyz = fstr(substr(fname_inp, 1, auxiliary_routines.indexlastdot(fname_inp) - 1) + '-second.xyz',
                             256)
            _write_xyz(nout, fname_xyz, natocentral, atsymbolcentral, xacentral, yacentral, zacentral, attypecentral,
                       centraldefinition)
    # Now repeat this for second molecules if they are different from first molecules
    # Write out coordinates of centralmolecule. We will use the average coordinates of all
    # these molecules.
    im = 0
    ncentral = 0
    centralmolecule = False
    while im < nmol:
        im = im + 1
        centralmolecule = S.nm[im] == 2
        if centralmolecule:  # Write out coordinates of this molecule
            if ncentral == 0 or S.outputaveragemol:
                ncentral = ncentral + 1
                if ncentral == 1:  # First time allocate the arrays to store summed coordinates
                    natocentral = int(C.nato[im])
                    xacentral = real_array(natocentral)
                    yacentral = real_array(natocentral)
                    zacentral = real_array(natocentral)
                    atsymbolcentral = char_array(natocentral)
                    attypecentral = char_array(natocentral)
                    centraldefinition = [False] * (natocentral + 1)
                # Shift coordinates to centre used by SHARM/SDF
                xshift = C.xadd[im] - S.xcent[im]
                yshift = C.yadd[im] - S.ycent[im]
                zshift = C.zadd[im] - S.zcent[im]
                for ia in range(C.naf[im], C.nal[im] + 1):
                    # Temporary coordinates
                    xt = M.xa[ia] + xshift
                    yt = M.ya[ia] + yshift
                    zt = M.za[ia] + zshift
                    # Need to rotate these to align with the laboratory axes - performed in reverse
                    # order about the laboratory axes
                    if C.nato[im] > 1:
                        xr, yr, zr = auxiliary_routines.rotateeuler(xt, yt, zt, -S.chim[im], -acosf(S.cthetm[im]),
                                                                    -S.phim[im], None, None, None)
                    else:
                        xr = xt
                        yr = yt
                        zr = zt
                    iaref = ia - C.naf[im] + 1
                    if ncentral == 1:
                        xacentral[iaref] = xr
                        yacentral[iaref] = yr
                        zacentral[iaref] = zr
                        atsymbolcentral[iaref] = M.atsymbol[M.itype[ia]]
                        attypecentral[iaref] = M.type[M.itype[ia]]
                        # Check whether this atom is part of the central molecule definition
                        # (the Fortran tests the first-molecule definition here too)
                        centraldefinition[iaref] = _central_definition(ia, iaref)
                    else:
                        xacentral[iaref] = xacentral[iaref] + xr
                        yacentral[iaref] = yacentral[iaref] + yr
                        zacentral[iaref] = zacentral[iaref] + zr
    # Form average coordinates for central molecule
    if ncentral > 0:
        for ia in range(1, natocentral + 1):
            xacentral[ia] = xacentral[ia] * C.dmax / f32(ncentral)
            yacentral[ia] = yacentral[ia] * C.dmax / f32(ncentral)
            zacentral[ia] = zacentral[ia] * C.dmax / f32(ncentral)
        fname_xyz = fstr(substr(fname_inp, 1, auxiliary_routines.indexlastdot(fname_inp) - 1) + '-second.xyz', 256)
        _write_xyz(nout, fname_xyz, natocentral, atsymbolcentral, xacentral, yacentral, zacentral, attypecentral,
                   centraldefinition)
    nmolsharm = nmol
    ierr = 0
    ierr = sharm_routines_rev1.form_sum_prod_dlmn_cgcoeff(nmolsharm, ierr)
    return nout


def _central_definition(ia, iaref):
    """The test 'Check whether this atom is part of the central molecule definition' of
    sharm (always with the first molecule's ncent1, lctyp1, ncz1, lztyp1, ncy1, lytyp1).
    Returns found."""
    found = False
    ic = 0
    while not found and ic < sharm_inc.ncent1:
        ic = ic + 1
        found = molecom.itype[ia] == sharm_inc.lctyp1[ic]
    iz = 0
    while not found and iz < sharm_inc.ncz1:
        iz = iz + 1
        found = iaref == sharm_inc.lztyp1[iz]
    iy = 0
    while not found and iy < sharm_inc.ncy1:
        iy = iy + 1
        found = iaref == sharm_inc.lytyp1[iy]
    return bool(found)
