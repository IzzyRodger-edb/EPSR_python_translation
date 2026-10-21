"""Translation of fortran_src/sdf_rev2.f: ``program standalone`` (sdf),
``subroutine run_sdf`` and ``subroutine sdf`` -- spherical harmonic coefficients of the
spatial density function of a second "molecule" around a central "molecule", each
defined by a set of atom types and the distances between them.

The shared routines are in library/sharm_routines_rev1.py, the data of
``module sharm_mod`` in common/sharm_mod.py.  ``setup_symbols_sdf`` is in
setup_symbols_generic.f (library/setup_symbols_generic.py, imported inside the routine:
another work package).

Quirks reproduced (see docs/translation_notes/Q.md): run_sdf does not include
sharm_sdf.inc, so its pi, piconv and the test angles phim1 ... phillab are local
variables (sdf and form_dlmn see the sharm_sdf.inc ones, which stay zero); the check of
the remaining intra-group distances uses ``zdiff=ytemp(i2)-ztemp(i1)``.  The run-time
FORMAT of the column titles is the patched one of the reference build
(``'(1h(,'//il1fmt...``, tools/build_fortran_reference.sh).
"""

from ..common import clustercom, commands, lattice, molecom, sharm_mod, sharm_sdf, symbols
from ..common import sharm as sharm_inc  # sharm.inc
from ..fortran import (
    PI4ATAN,
    Ch,
    FortranEOF,
    FortranIOError,
    FortranStop,
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
    sqrtf,
    substr,
    units,
)
from ..library import (
    auxiliary_routines,
    file_and_directory_commands,
    minimagec,
    parse,
    read_and_save_files,
    reallocation_routines,
    sharm_routines_rev1,
)
from ..library import ato_read_and_save
from ..trace import routine

_F201 = "(/'sdf> Atom type ',a,' does not exist in .ato file ',a)"
_F202 = "(/'sdf> None of the specified atom types were recognised. '/' - returning to main program')"
_F203 = "(/'sdf> No. of atoms used to define centre molecule = ',i2)"
_F204 = "(/'sdf> Not enough pair distances defined ','for central molecule. Need: ',i5)"
_F205 = "(/'sdf> Not enough pair distances defined ','for second molecule. Need: ',i5)"
_F206 = "(/'sdf> No. of atoms used to define second molecule = ',i2)"
_F333 = "(1x,'sdf> Number of type 1 molecules = ',i10/1x, 'sdf> Number of type 2 molecules = ',i10)"
_FXYZ = "(a3,1x,f10.5,1x,f10.5,1x,f10.5,1x,a3)"


def main():
    """sdf_rev2.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_sdf()


def _parse(iref, nspc):
    """call parse(symbol_value(iref),nwords,ncf,ncl,nspc,256,nworddim)"""
    commands.nwords = parse.parse(symbols.symbol_value[iref], commands.nwords, commands.ncf, commands.ncl, nspc,
                                  256, commands.nworddim)


def _word(iref, ic):
    """symbol_value(iref)(ncf(ic):ncl(ic))"""
    return substr(symbols.symbol_value[iref], commands.ncf[ic], commands.ncl[ic])


@routine
def run_sdf():
    """sdf_rev2.f:run_sdf -- performs calculation of spherical harmonic coefficients
    from .ato file.

    uses a specified set of atoms to define the central "molecule", then calculates the
    spatial density function (SDF) of a specified second "molecule" around this.
    """
    from ..library import setup_symbols_generic

    S = sharm_mod
    H = sharm_inc
    sv = symbols.symbol_value
    # character*50 formatout
    # character*2 ifmt1,ifmt2,ifmt3,ifmt4,il1fmt,il2fmt,ilfmt,in1fmt,in2fmt
    il1fmt = il2fmt = ilfmt = in1fmt = in2fmt = ''
    # real pi,piconv (locals: run_sdf does not include sharm_sdf.inc, so the test angles
    # phim1 ... phillab below are local implicit reals too)
    phim1 = thetam1 = chim1 = phim2 = thetam2 = chim2 = thetallab = phillab = f32(0.0)
    #
    # setup the SDF symbols
    #
    setup_symbols_generic.setup_symbols_sdf()
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

    pi = PI4ATAN  # 4.0*atan(1.0)
    piconv = pi / r4(180.0)
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
    H.nz = read_internal(sv[iref]).int(H.nz)
    #
    # set up rmax of calculation
    #
    iref = iref + 1
    H.zmax = read_internal(sv[iref]).real(H.zmax)
    #
    # set up number of configurations already stored
    #
    iref = iref + 1
    H.nsumt = read_internal(sv[iref]).int(H.nsumt)
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
    # For the central "molecule" read the list of atom types to be used to define the
    # unit. The first in the list will be assigned to the centre of the molecule.
    #
    iref = iref + 1
    _parse(iref, 5)
    H.natofirst = 0

    # skip through the atom types to find out the type number of the specified atom

    ic = 0
    while ic < commands.nwords and H.natofirst < H.mcentre:
        ic = ic + 1
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
            H.natofirst = H.natofirst + 1
            H.lctyp1[H.natofirst] = ifound
            H.typefirst[H.natofirst] = ifound
    if H.natofirst <= 0:
        units.write(6, fmt(_F202))
        return
    else:
        units.write(6, fmt(_F203, H.natofirst))
    H.ncent1 = 1

    # For central molecule read the array of maximum distances between atoms in the order
    # (1,2) (1,3) (1,4) ... (2,3) (2,4) ... (3,4) etc.
    # There should be natofirst*(natofirst-1)/2 of them

    iref = iref + 1
    _parse(iref, 5)
    npairfirstexpec = H.natofirst * (H.natofirst - 1) // 2
    if commands.nwords < npairfirstexpec:
        units.write(6, fmt(_F204, npairfirstexpec))
        return
    i = 0
    while i < npairfirstexpec:
        i = i + 1
        H.distfirst[i] = read_internal(_word(iref, i)).real(H.distfirst[i])
    H.npairfirst = i

    # Variation of distance (half width) about the specified distance for an atom to be accepted.
    # This is expressed as a fraction of the specified distance. For example fraction of 0.2
    # for a bond distance of 2A means atoms within 2 +- 0.4A will be accepted

    iref = iref + 1
    H.fractfirst = read_internal(sv[iref]).real(H.fractfirst)

    #
    # for molecule 1 read primary axis type (1=X, 2=Y, 3=Z), and
    # the number and atom numbers of atoms used to define the molecule's nzt1
    # axis
    #
    iref = iref + 1
    _parse(iref, 5)
    H.nzt1 = fstr(_word(iref, 1), 2)
    H.ncz1 = commands.nwords - 1
    if H.ncz1 > H.mcentre:
        H.ncz1 = H.mcentre
    for i in range(1, H.ncz1 + 1):
        ic = i + 1
        H.lztyp1[i] = read_internal(_word(iref, ic)).int(H.lztyp1[i])
    #
    # for molecule 1 read secondary axis type (1=X, 2=Y, 3=Z), and
    # the number and atom numbers of atoms used to define the
    # molecule's nyt1 axis
    #
    iref = iref + 1
    _parse(iref, 5)
    H.nyt1 = fstr(_word(iref, 1), 2)
    H.ncy1 = commands.nwords - 1
    if H.ncy1 > H.mcentre:
        H.ncy1 = H.mcentre
    for i in range(1, H.ncy1 + 1):
        ic = i + 1
        H.lytyp1[i] = read_internal(_word(iref, ic)).int(H.lytyp1[i])

    # For second "molecule" read type of atom used to define the molecular
    # centre. Only the first in the list is used.

    iref = iref + 1
    _parse(iref, 5)

    # Skip through the specified atom types to find out the type number of the specified atom
    # We also test to see if this atom is the same as the corresponding atom on the central
    # molecule. If not the molecules are treated as distinct.

    H.natosecond = 0
    S.ldistinct = 0
    for ic in range(1, commands.nwords + 1):
        it = 0
        ifound = 0
        while it < molecom.ntype and ifound == 0 and H.natosecond < H.mcentre:
            it = it + 1
            if molecom.type[it] == _word(iref, ic):
                ifound = it
        if ifound == 0:
            units.write(6, fmt(_F201, _word(iref, ic), substr(symbols.fname_ato, 1, symbols.ilenfname_ato)))
            return
        else:
            H.natosecond = H.natosecond + 1
            H.lctyp2[H.natosecond] = ifound

    if H.natosecond <= 0:
        units.write(6, fmt(_F202))
        return
    else:
        units.write(6, fmt(_F206, H.natosecond))
    H.ncent2 = 1

    # test for identity of central and second molecules

    S.ldistinct = 0
    if H.natosecond == H.natofirst:
        S.ldistinct = 0
        for it in range(1, H.natofirst + 1):
            if H.lctyp2[it] != H.lctyp1[it]:
                S.ldistinct = 1
    else:
        S.ldistinct = 1

    # If the second and central molecules are identical, the rest of the input for the
    # second molecule is identical to the first.

    if S.ldistinct != 0:

        # For second molecule read the array of maximum distances between atoms in the order
        # (1,2) (1,3) (1,4) ... (2,3) (2,4) ... (3,4) etc.
        # There should be natosecond*(natosecond-1)/2 of them

        iref = iref + 1
        _parse(iref, 5)
        npairsecondexpec = H.natosecond * (H.natosecond - 1) // 2
        if commands.nwords < npairsecondexpec:
            units.write(6, fmt(_F205, npairsecondexpec))
            return
        i = 0
        while i < npairsecondexpec:
            i = i + 1
            H.distsecond[i] = read_internal(_word(iref, i)).real(H.distsecond[i])
        H.npairsecond = i

        # Variation of distance (half width) about the specified distance for an atom to be
        # accepted.  This is expressed as a fraction of the specified distance. For example
        # fraction of 0.2 for a bond distance of 2A means atoms within 2 +- 0.4A will be accepted

        iref = iref + 1
        H.fractsecond = read_internal(sv[iref]).real(H.fractsecond)

        # for molecule 2 read primary axis type (1=X, 2=Y, 3=Z), and
        # the number and atom numbers of atoms used to define the molecule's nzt1
        # axis

        iref = iref + 1
        _parse(iref, 5)
        H.nzt2 = fstr(_word(iref, 1), 2)
        H.ncz2 = commands.nwords - 1
        for i in range(1, H.ncz2 + 1):
            ic = i + 1
            H.lztyp2[i] = read_internal(_word(iref, ic)).int(H.lztyp2[i])

        # for molecule 2 read secondary axis type (1=X, 2=Y, 3=Z), and
        # the number and atom numbers of atoms used to define the
        # molecule's nyt1 axis

        iref = iref + 1
        _parse(iref, 5)
        H.nyt2 = fstr(_word(iref, 1), 2)
        H.ncy2 = commands.nwords - 1
        for i in range(1, H.ncy2 + 1):
            ic = i + 1
            H.lytyp2[i] = read_internal(_word(iref, ic)).int(H.lytyp2[i])

    else:

        # For identical molecules set all the parameters for the second molecule
        # equal to the first molecule.

        H.ncent2 = H.ncent1
        H.natosecond = H.natofirst
        H.npairsecond = H.npairfirst
        for i in range(1, H.npairsecond + 1):
            H.distsecond[i] = H.distfirst[i]
        H.fractsecond = H.fractfirst
        H.nzt2 = H.nzt1
        H.ncz2 = H.ncz1
        for i in range(1, H.ncz2 + 1):
            H.lztyp2[i] = H.lztyp1[i]
        H.nyt2 = H.nyt1
        H.ncy2 = H.ncy1
        for i in range(1, H.ncy2 + 1):
            H.lytyp2[i] = H.lytyp1[i]
        iref = iref + 4
    # Decide whether to use simulated phim, thetam, chim
    iref = iref + 1
    _parse(iref, 5)
    S.performtest = _word(iref, 1) != 'notest' and commands.nwords > 7
    if S.performtest:
        # read(...,*,iostat=ierr) phim1,thetam1,chim1,phim2,thetam2,chim2,thetallab,phillab
        vals = [phim1, thetam1, chim1, phim2, thetam2, chim2, thetallab, phillab]
        r = read_internal(substr(sv[iref], commands.ncf[1], commands.ncl[commands.nwords]))
        ierr = 0
        try:
            for k in range(8):
                vals[k] = r.real(vals[k])
        except (FortranEOF, FortranIOError) as e:
            ierr = iostat_of(e)
        phim1, thetam1, chim1, phim2, thetam2, chim2, thetallab, phillab = vals
        phim1 = phim1 * piconv
        thetam1 = thetam1 * piconv
        chim1 = chim1 * piconv
        phim2 = phim2 * piconv
        thetam2 = thetam2 * piconv
        chim2 = chim2 * piconv
        thetallab = thetallab * piconv
        phillab = phillab * piconv
    else:
        phim1 = f32(0.0)
        thetam1 = f32(0.0)
        chim1 = f32(0.0)
        phim2 = f32(0.0)
        thetam2 = f32(0.0)
        chim2 = f32(0.0)
        thetallab = f32(0.0)
        phillab = f32(0.0)
    if ierr != 0:
        S.performtest = False
    iref = iref + 1
    S.outputaveragemol = sv[iref] == 'y'

    # Only perform the calculation when nsum is zero or a multiple of nmult

    nmult = 5
    nsum = H.nsumt // nmult
    nt = nsum * nmult
    if nt == H.nsumt:
        #
        # normalise atomic coordinates to box size
        #
        ato_read_and_save.normalise_coordinates()
        H.zmaxd = H.zmax / clustercom.dmax
        if H.zmaxd > 0.5:
            H.zmaxd = f32(0.5)
            H.zmax = H.zmaxd * clustercom.dmax
        # Set up l and n values and allocate some needed arrays
        sharm_routines_rev1.get_l_and_n_values()
        S.mz = H.nz
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
        H.zstep = H.zmaxd / f32(H.nz)
        H.zhalf = r4(0.5) * H.zstep
        for j in range(1, H.nz + 1):
            S.z[j] = H.zstep * f32(j) - H.zhalf
        units.write(6, lw('sdf> zmaxd: ', H.zmax, clustercom.dmax, H.zmaxd, H.zstep))
        # zero gofz
        read_and_save_files.initgofzs(H.nz, S.ndist, S.gofz, S.gofzs, S.mz, S.mdist)
        sdf(nin, symbols.fname_inp)
        #
        # set the z values
        #
        for i in range(1, H.nz + 1):
            S.z[i] = S.z[i] * clustercom.dmax
        #
        # axis title
        #
        H.titler = '   r'
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
        H.ext = 'h'
        nsum = read_and_save_files.outgofzs_sharm(symbols.fnameinp, H.ext, H.nz, S.ndist, nsum, S.z, S.gofz,
                                                  S.gofzs, S.gzt, S.gzt2, H.titler, S.titledat, S.mz, S.mdist)
    if not S.performtest:
        H.nsumt = H.nsumt + 1
    #
    # write input file
    #
    #
    # make sure all current values are stored
    #
    sv[4] = fstr(lw(H.nsumt), 256)
    sv[5] = fstr(lw(S.ndist), 256)
    units.write(6, lw(S.ndist, ' ', Ch(sv[5], 256)))
    read_and_save_files.save_files_others(14)
    units.close(14)
    if S.performtest:
        units.write(6, lw('sdf> Test performed.'))
    return


@routine
def sdf(nout, fname_inp):
    """sdf_rev2.f:sdf -- find the central and second "molecules" (groups of atoms of
    the specified types at the specified distances), their Euler angles, write the averaged
    molecules, and sum the coefficients (form_sum_prod_dlmn_cgcoeff).

    Returns nout.
    """
    S = sharm_mod
    H = sharm_inc
    M = molecom
    C = clustercom
    xat = yat = zat = None  # real, dimension(:), allocatable
    iasave = None  # integer, dimension(:,:), allocatable: atom numbers of atoms used to define
    # molecules, so as not to be used in calculating coefficients
    # character(len=256) :: fname_xyz,textline
    mcentre = H.mcentre
    xfirstc = real_array(mcentre)
    yfirstc = real_array(mcentre)
    zfirstc = real_array(mcentre)
    xsecondc = real_array(mcentre)
    ysecondc = real_array(mcentre)
    zsecondc = real_array(mcentre)
    atsymbolfirstc = char_array(mcentre)
    atsymbolsecondc = char_array(mcentre)
    attypefirstc = char_array(mcentre)
    attypesecondc = char_array(mcentre)
    atomno = sharm_sdf.atomno
    xtemp = H.xtemp
    ytemp = H.ytemp
    ztemp = H.ztemp
    distcentre = H.distcentre
    lctyp = H.lctyp

    sharm_sdf.pi = PI4ATAN  # 4.0*atan(1.0)
    sharm_sdf.piconv = sharm_sdf.pi / r4(180.0)
    degconv = r4(180.0) / sharm_sdf.pi  # noqa: F841
    H.zstep = H.zmaxd / f32(H.nz)
    # volume of sphere of radius 0.5
    H.vol = lattice.alatx * (lattice.blaty * lattice.clatz)
    # set up counter for second atoms
    nsec = 0  # noqa: F841
    # Arrays to store Euler angles
    if S.nm is not None:
        S.nm = S.phim = S.cthetm = S.chim = S.xcent = S.ycent = S.zcent = None
        S.dlmnmc = S.chiphi = S.dlmnc = S.philab = None
    sharm_routines_rev1.get_max_size_dlmn()
    nmol = C.nmol
    units.write(6, lw('sdf> ', nmol, S.mdlmn))
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
    # Allocate arrays to store atom numbers used in molecule definitions
    iasave = int_array(max(H.natofirst, H.natosecond), nmol)

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
    # Save the label of the first centre molecule
    imolcentre = 0  # noqa: F841

    # Look for centre molecules. When found calculate their Euler angles.

    imol1 = 0
    ncentral = 0
    nsecond = 0
    dmax = C.dmax

    # step through the atoms and find those which correspond to ltyp1(1)
    for ia1 in range(1, M.natom + 1):

        # First set up the test conditions

        # Signal that so far a true molecule has not been found

        iflag = 0

        itmol = 0

        while itmol <= S.ldistinct:
            itmol = itmol + 1
            if itmol == 1:
                H.natomol = H.natofirst
                H.npairmol = H.npairfirst
                for i in range(1, H.natomol + 1):
                    lctyp[i] = H.lctyp1[i]
                for i in range(1, H.npairmol + 1):
                    distcentre[i] = H.distfirst[i] / dmax
                H.fractcentre = H.fractfirst
            else:
                H.natomol = H.natosecond
                H.npairmol = H.npairsecond
                for i in range(1, H.natomol + 1):
                    lctyp[i] = H.lctyp2[i]
                for i in range(1, H.npairmol + 1):
                    distcentre[i] = H.distsecond[i] / dmax
                H.fractcentre = H.fractsecond
            natomol = H.natomol
            npairmol = H.npairmol
            fractcentre = H.fractcentre

            # Test to see if this atom corresponds to the required centre atom.

            if M.itype[ia1] == lctyp[1]:

                # Increment the molecule counter in case this is a valid molecule. It will later be
                # decremented if this turns out to not to be valid.

                if imol1 < nmol:
                    imol1 = imol1 + 1

                # Set the label whether this is a first or second molecule

                S.nm[imol1] = itmol

                # Now search for other atoms which satisfy the specified label and distance criteria

                nfound = 1
                atomno[1] = ia1
                im1 = M.mref[ia1]
                xo1 = C.xadd[im1] + M.xa[ia1]
                yo1 = C.yadd[im1] + M.ya[ia1]
                zo1 = C.zadd[im1] + M.za[ia1]
                # Temporary coordinates of this atom relative to the centre of the molecule
                xtemp[1] = 0.0
                ytemp[1] = 0.0
                ztemp[1] = 0.0
                S.xcent[imol1] = xo1
                S.ycent[imol1] = yo1
                S.zcent[imol1] = zo1

                # At this point there if there is only 1 atom in the molecule and it must be of the
                # correct type, otherwise we wouldn't have got here, we signal immediately that a
                # molecule has been found, and there is no need to check intra-molecular distances.

                if natomol > nfound:

                    # Set the other atom numbers for this molecule to zero. This will be used to
                    # signal whether or not an atom has been found at each index. If the correct set
                    # of atom types with the specified distances from neighbours cannot be found,
                    # this molecule is abandoned as a possible centre molecule

                    for i in range(2, natomol + 1):
                        atomno[i] = 0

                    # Step through second atoms and see which if any satisfy the type and distance
                    # constraints

                    ia2 = 0
                    while ia2 < M.natom and nfound < natomol:

                        ia2 = ia2 + 1

                        # First check whether this atom already exists for this molecule

                        itest = 1
                        while itest <= natomol and atomno[itest] != ia2:
                            itest = itest + 1

                        # itest .gt. natomol means it was not found so far

                        if itest > natomol:

                            # Now see whether it fits the type number for one of the atoms of this
                            # molecule.

                            i = 1
                            ifound = 0
                            while i < natomol and ifound == 0:
                                i = i + 1

                                # Check to see if ia2 matches the current label and that this slot
                                # has not already been assigned

                                if lctyp[i] == M.itype[ia2] and atomno[i] == 0:

                                    # now measure its distance from the centre atom

                                    im2 = M.mref[ia2]
                                    xo2 = C.xadd[im2] + M.xa[ia2]
                                    yo2 = C.yadd[im2] + M.ya[ia2]
                                    zo2 = C.zadd[im2] + M.za[ia2]
                                    xdiff = xo2 - xo1
                                    ydiff = yo2 - yo1
                                    zdiff = zo2 - zo1
                                    xdiff, ydiff, zdiff = minimagec.minimagec(xdiff, ydiff, zdiff)
                                    dmeas = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff
                                    dmeas = sqrtf(dmeas)
                                    dtestlow = distcentre[i - 1] * (r4(1.0) - fractcentre)
                                    dtesthigh = distcentre[i - 1] * (r4(1.0) + fractcentre)
                                    if dmeas >= dtestlow and dmeas <= dtesthigh:

                                        # atom found at correct distance from centre - save the atom
                                        # number in atomno.

                                        ifound = 1
                                        atomno[i] = ia2
                                        xtemp[i] = xdiff
                                        ytemp[i] = ydiff
                                        ztemp[i] = zdiff
                                        nfound = nfound + 1

                    # If atoms of the correct number and type and at the correct distance
                    # from the central atom have been found, then this group is a possible candidate
                    # for a molecule. Now check that the other distances are correct.

                    if nfound == natomol:
                        iflag = 1
                        iref = natomol - 1
                        i1 = 2
                        i2 = 3
                        while iref < npairmol and iflag == 1:
                            xdiff = xtemp[i2] - xtemp[i1]
                            ydiff = ytemp[i2] - ytemp[i1]
                            zdiff = ytemp[i2] - ztemp[i1]  # sic (the Fortran uses ytemp here)
                            xdiff, ydiff, zdiff = minimagec.minimagec(xdiff, ydiff, zdiff)
                            dmeas = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff
                            dmeas = sqrtf(dmeas)
                            iref = iref + 1
                            dtestlow = distcentre[iref] * (r4(1.0) - fractcentre)
                            dtesthigh = distcentre[iref] * (r4(1.0) + fractcentre)
                            if dmeas < dtestlow or dmeas > dtesthigh:
                                iflag = 0
                            i2 = i2 + 1
                            if i2 > natomol:
                                i1 = i1 + 1
                                i2 = i1 + 1
                    else:
                        iflag = 0
                else:
                    iflag = 1
                # If iflag is still 1 after this it means a valid "molecule" has been found
                # so we determine its Euler angles
                if iflag == 1:
                    # Save the atom numbers that form this molecule
                    for i in range(1, natomol + 1):
                        iasave[i, imol1] = atomno[i]
                    if S.nm[imol1] == 1:
                        nsmmo1 = nsmmo1 + 1
                    if S.nm[imol1] == 2:
                        nsmmo2 = nsmmo2 + 1
                    # Form atomic coordinates relative to this centre
                    xat = reallocation_routines.reallocate1d_r(xat, natomol)
                    yat = reallocation_routines.reallocate1d_r(yat, natomol)
                    zat = reallocation_routines.reallocate1d_r(zat, natomol)
                    for ia in range(1, natomol + 1):
                        xat[ia] = xtemp[ia]
                        yat[ia] = ytemp[ia]
                        zat[ia] = ztemp[ia]
                    if S.nm[imol1] == 1:
                        S.phim[imol1], S.cthetm[imol1], S.chim[imol1] = sharm_routines_rev1.form_eulerangles(
                            imol1, xat, yat, zat, H.nzt1, H.nyt1, H.ncz1, H.ncy1, H.lztyp1, H.lytyp1, alab, blab,
                            clab, S.phim[imol1], S.cthetm[imol1], S.chim[imol1])
                    elif S.nm[imol1] == 2:
                        S.phim[imol1], S.cthetm[imol1], S.chim[imol1] = sharm_routines_rev1.form_eulerangles(
                            imol1, xat, yat, zat, H.nzt2, H.nyt2, H.ncz2, H.ncy2, H.lztyp2, H.lytyp2, alab, blab,
                            clab, S.phim[imol1], S.cthetm[imol1], S.chim[imol1])
                    # If testing orientations, supply a fictitious but known distribution of angles
                    if S.performtest:
                        #  write(6,499) imol1,phim1,thetam1,chim1,phim2,thetam2,chim2
                        #  call populate_disk(pi,phim1,thetam1,chim1,phinew,thetanew,chinew,thetaspread,phispread)
                        # Rotate this molecule onto desired orientation
                        # (phim1 ... chim2 are the sharm_sdf.inc variables here, not run_sdf's)
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
                        units.write(6, lw('sdf> Error when calculating dlm coefficients for molecule ', imol1))
                        return nout
                    # If this is a centre molecule, we form its average rotated coordinates for later
                    # output
                    if S.nm[imol1] == 1:
                        if ncentral == 0 or S.outputaveragemol:
                            ncentral = ncentral + 1
                            for ia in range(1, natomol + 1):
                                # Rotate coordinates onto lab axes
                                # Temporary coordinates
                                xt = xtemp[ia]
                                yt = ytemp[ia]
                                zt = ztemp[ia]
                                # Need to rotate these to align with the laboratory axes - performed
                                # in reverse order about the laboratory axes
                                if H.natofirst > 1:
                                    xr, yr, zr = auxiliary_routines.rotateeuler(
                                        xt, yt, zt, -S.chim[imol1], -acosf(S.cthetm[imol1]), -S.phim[imol1],
                                        None, None, None)
                                else:
                                    xr = xt
                                    yr = yt
                                    zr = zt
                                units.write(6, lw(atomno[ia], xr * dmax, yr * dmax, zr * dmax))
                                if ncentral == 1:
                                    xfirstc[ia] = xr
                                    yfirstc[ia] = yr
                                    zfirstc[ia] = zr
                                    atsymbolfirstc[ia] = M.atsymbol[M.itype[atomno[ia]]]
                                    attypefirstc[ia] = M.type[M.itype[atomno[ia]]]
                                else:
                                    xfirstc[ia] = xfirstc[ia] + xr
                                    yfirstc[ia] = yfirstc[ia] + yr
                                    zfirstc[ia] = zfirstc[ia] + zr
                    # If this is a centre molecule, we form its average rotated coordinates for later
                    # output
                    if S.nm[imol1] == 2:
                        if nsecond == 0 or S.outputaveragemol:
                            nsecond = nsecond + 1
                            for ia in range(1, natomol + 1):
                                # Rotate coordinates onto lab axes
                                # Temporary coordinates
                                xt = xtemp[ia]
                                yt = ytemp[ia]
                                zt = ztemp[ia]
                                # Need to rotate these to align with the laboratory axes - performed
                                # in reverse order about the laboratory axes
                                if H.natosecond > 1:
                                    xr, yr, zr = auxiliary_routines.rotateeuler(
                                        xt, yt, zt, -S.chim[imol1], -acosf(S.cthetm[imol1]), -S.phim[imol1],
                                        None, None, None)
                                else:
                                    xr = xt
                                    yr = yt
                                    zr = zt
                                #  write(6,*) atomno(ia),xr*dmax,yr*dmax,zr*dmax
                                if nsecond == 1:
                                    xsecondc[ia] = xr
                                    ysecondc[ia] = yr
                                    zsecondc[ia] = zr
                                    atsymbolsecondc[ia] = M.atsymbol[M.itype[atomno[ia]]]
                                    attypesecondc[ia] = M.type[M.itype[atomno[ia]]]
                                else:
                                    xsecondc[ia] = xsecondc[ia] + xr
                                    ysecondc[ia] = ysecondc[ia] + yr
                                    zsecondc[ia] = zsecondc[ia] + zr
                else:
                    # If this is not a valid molecule then we ignore this candidate
                    if imol1 > 0:
                        imol1 = imol1 - 1
    nmolsharm = imol1
    # If ldistinct = 1 (molecules are distinct) we need to ensure that none of the atoms that
    # occur in second molecules do not also occur in first molecules
    if S.ldistinct > 0:
        imol1 = 0
        while imol1 < nmolsharm:
            imol1 = imol1 + 1
            if S.nm[imol1] == 1:
                # Test second molecules to see if there any atoms common with first molecule
                imol2 = 0
                while imol2 < nmolsharm:
                    imol2 = imol2 + 1
                    iflag = 0
                    if imol1 != imol2 and S.nm[imol2] == 2:
                        ia1 = 0
                        while ia1 < H.natofirst and iflag == 0:
                            ia1 = ia1 + 1
                            ia2 = 0
                            while ia2 < H.natosecond and iflag == 0:
                                ia2 = ia2 + 1
                                if iasave[ia1, imol1] == iasave[ia2, imol2]:
                                    iflag = 1
                    if iflag == 1:  # This second molecule must be removed from list
                        units.write(6, fmt('(a,i5)', 'sdf> Removing molecule: ', imol2))
                        imol = imol2
                        while imol < nmolsharm:
                            imolt = imol + 1
                            S.nm[imol] = S.nm[imolt]
                            S.phim[imol] = S.phim[imolt]
                            S.cthetm[imol] = S.cthetm[imolt]
                            S.chim[imol] = S.chim[imolt]
                            if S.nm[imol] == 1:
                                for ia in range(1, H.natofirst + 1):
                                    iasave[ia, imol] = iasave[ia, imolt]
                            else:
                                for ia in range(1, H.natosecond + 1):
                                    iasave[ia, imol] = iasave[ia, imolt]
                            imol = imolt
                        nmolsharm = nmolsharm - 1
                        nsmmo2 = nsmmo2 - 1
                        if imol1 >= imol2:
                            # Need to decrement imol1 if it is .gt. imol2, since imol2 has been
                            # removed from the list
                            imol1 = imol1 - 1
                        imol2 = imol2 - 1
        # Write out lists of atoms to check that algorithm works
        for imol in range(1, nmolsharm + 1):
            if S.nm[imol] == 1:
                units.write(6, lw(imol, S.nm[imol], *[iasave[ia, imol] for ia in range(1, H.natofirst + 1)]))
        for imol in range(1, nmolsharm + 1):
            if S.nm[imol] == 2:
                units.write(6, lw(imol, S.nm[imol], *[iasave[ia, imol] for ia in range(1, H.natosecond + 1)]))
    units.write(6, fmt(_F333, nsmmo1, nsmmo2))

    # Write out coordinates of central molecule

    fname_xyz = fstr(substr(fname_inp, 1, auxiliary_routines.indexlastdot(fname_inp) - 1) + '-first.xyz', 256)
    ierr = _open(nout, fname_xyz)
    if ierr == 0 and ncentral > 0:
        # Write number of atoms in molecule
        units.write(nout, fmt('(i6)', H.natofirst))
        # Write information line - not used for plotting
        units.write(nout, fmt('(a)', 'Central molecule coordinate file: ' + fname_xyz))
        # Shift coordinates to centre used by SHARM/SDF
        for ia in range(1, H.natofirst + 1):
            xr = xfirstc[ia] * dmax / f32(ncentral)
            yr = yfirstc[ia] * dmax / f32(ncentral)
            zr = zfirstc[ia] * dmax / f32(ncentral)
            textline = fstr(fmt(_FXYZ, Ch(atsymbolfirstc[ia], 3), xr, yr, zr, Ch(attypefirstc[ia], 3)), 256)
            textline = fstr(textline + '  #CENTRALMOLECULE', 256)
            units.write(nout, fmt('(a)', textline))
        units.close(nout)
    else:
        units.write(6, lw('sdf> Error writing first molecule file'))
        raise FortranStop(0)
    # Write out coordinates of second molecule
    fname_xyz = fstr(substr(fname_inp, 1, auxiliary_routines.indexlastdot(fname_inp) - 1) + '-second.xyz', 256)
    ierr = _open(nout, fname_xyz)
    if ierr == 0 and S.ldistinct == 0 and ncentral > 0:
        # Write number of atoms in molecule
        units.write(nout, fmt('(i6)', H.natofirst))
        # Write information line - not used for plotting
        units.write(nout, fmt('(a)', 'Central molecule coordinate file: ' + fname_xyz))
        # Shift coordinates to centre used by SHARM/SDF
        for ia in range(1, H.natofirst + 1):
            xr = xfirstc[ia] * dmax / f32(ncentral)
            yr = yfirstc[ia] * dmax / f32(ncentral)
            zr = zfirstc[ia] * dmax / f32(ncentral)
            textline = fstr(fmt(_FXYZ, Ch(atsymbolfirstc[ia], 3), xr, yr, zr, Ch(attypefirstc[ia], 3)), 256)
            textline = fstr(textline + '  #CENTRALMOLECULE', 256)
            units.write(nout, fmt('(a)', textline))
        units.close(nout)
    elif ierr == 0 and S.ldistinct == 1 and nsecond > 0:  # Second molecule is different from central molecule
        # Write number of atoms in molecule
        units.write(nout, fmt('(i6)', H.natosecond))
        # Write information line - not used for plotting
        units.write(nout, fmt('(a)', 'Central molecule coordinate file: ' + fname_xyz))
        # Shift coordinates to centre used by SHARM/SDF
        for ia in range(1, H.natosecond + 1):
            xr = xsecondc[ia] * dmax / f32(nsecond)
            yr = ysecondc[ia] * dmax / f32(nsecond)
            zr = zsecondc[ia] * dmax / f32(nsecond)
            textline = fstr(fmt(_FXYZ, Ch(atsymbolsecondc[ia], 3), xr, yr, zr, Ch(attypesecondc[ia], 3)), 256)
            textline = fstr(textline + '  #CENTRALMOLECULE', 256)
            units.write(nout, fmt('(a)', textline))
        units.close(nout)
    else:
        units.write(6, lw('sdf> Error writing second molecule file'))
        raise FortranStop(0)
    ierr = sharm_routines_rev1.form_sum_prod_dlmn_cgcoeff(nmolsharm, ierr)
    return nout


def _open(nout, fname_xyz):
    """open(nout,file=fname_xyz,status='unknown',iostat=ierr).  Returns ierr."""
    try:
        units.open(nout, fname_xyz, status='unknown')
    except FortranIOError as e:
        return iostat_of(e)
    return 0
