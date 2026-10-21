"""Translation of fortran_src/lsi.f: ``program standalone``, ``subroutine run_lsi`` and
``subroutine find_lsi`` -- distribution of local structure indices (LSI) for a set of
bonded/non-bonded molecules.

``setup_symbols_lsi`` is in setup_symbols_generic.f (library/setup_symbols_generic.py,
imported inside the routine: it is translated by another work package).
``match`` is in chains_routines.f.  The COMMON of lsi.inc is ``common.lsi``, imported as
``lsicom`` because find_lsi has a local array called ``lsi``.
"""

from ..common import clustercom, commands, molecom, symbols
from ..common import lsi as lsicom
from ..common.epsr_dimension import matom, mmol
from ..common.lsi import mbond, mdist, mz
from ..fortran import PI4ATAN, Ch, acosf, cosf, f32, findex_padded, fmt, fstr, logical_array, lw, nint, r4
from ..fortran import read_internal, real_array, sinf, sqrtf, substr, units
from ..library import ato_read_and_save
from ..library import chains_routines
from ..library import file_and_directory_commands
from ..library import minimagec
from ..library import parse
from ..library import read_and_save_files
from ..trace import routine

_F200 = "('lsi> Not found enough items in distribution ',i2,', symbol: ',a11)"
_F100 = "(4x,a,'-',a,'-',a)"
_F101 = "(2x,a,'-',a,'-',a)"
_F102 = "(4x,a,'-',a,'-',a2)"
_F994 = "(1x,' Average of ',e12.5,' bonds/mol and ',e12.5,' non bonds/mol (in first shell)')"


def main():
    """lsi.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_lsi()


def _parse5():
    """call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)"""
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)


def _word(iw):
    """sometext(ncf(iw):ncl(iw))"""
    return substr(commands.sometext, commands.ncf[iw], commands.ncl[iw])


@routine
def run_lsi():
    """lsi.f:run_lsi -- 2019-01-04 Derived from chains_sp_criterion.

    This program calculates the distribution of local structure indices (LSI) for a set
    of bonded/non-bonded molecules.  The method is given by Shi, Russo, and Tanaka,
    J. Chem. Phys. 149, 224502 (2018); and Russo and Tanaka, 2014 NATURE COMMUNICATIONS
    DOI: 10.1038/ncomms4556.  In this method a list of bonded and non-bonded near
    neighbours for each molecule is set up.  The LSI is defined as the difference in
    distance from the central molecule between the nearest non-bonded and the furthest
    bonded molecule.
    """
    from ..library import setup_symbols_generic

    c = lsicom
    ireff = file_and_directory_commands.ireff
    # integer*4 ierrout
    # character*20 text
    # real*8 fc,fc1

    # real, allocatable, dimension(:) :: avbondsave,navbondedsave,navnonbondedsave
    # real, allocatable, dimension(:) :: avbondsave2,navbondedsave2,navnonbondedsave2
    pi = PI4ATAN  # 4.0*atan(1.0)
    piconv = pi / r4(180.0)
    #
    # setup the clusters symbols
    #
    setup_symbols_generic.setup_symbols_lsi()
    #
    # now setup the input file
    #
    nin = 11
    nout = 6

    # NB: ierrout is not initialised in the Fortran; setup_input_file sets it
    ierrout = 0
    nout, ierrout = file_and_directory_commands.setup_input_file(nin, nout, ierrout)
    #
    # if an error was encountered on setup, exit now
    #
    if ierrout != 0:
        return

    #
    # read the input file
    #

    # c$$$  RH: Actually read the .dat parameters into symbols arrays
    units.write(6, lw("Ran setup_input_file ok."))

    nout, ierrout = read_and_save_files.read_input_file(nin, nout, ierrout)
    units.write(6, lw("Ran read_input_file ok."))
    #
    # if an error was encountered on read, exit now
    #
    if ierrout != 0:
        return

    # c$$$     RH:
    # c$$$     Using the variables from the input file (*CLUSTERS.dat), which
    # c$$$     a restored in the char array symbol_value.
    # c$$$
    # c$$$     Values read:
    # c$$$       nz -- the max cluster size
    # c$$$       nsumt -- number of accumulations
    # c$$$  and the values for each cluster distribuition

    #
    # set up the number of lsi bins
    #
    iref = 2
    c.nz = read_internal(symbols.symbol_value[iref]).int(c.nz)
    if c.nz > mz:
        c.nz = mz
    #
    # set up min and max values of lsi
    #
    iref = iref + 1
    r = read_internal(symbols.symbol_value[iref])
    c.lsimin = r.real(c.lsimin)
    c.lsimax = r.real(c.lsimax)
    #
    # set up current accumulator
    #
    iref = iref + 1
    c.nsumt = read_internal(symbols.symbol_value[iref]).int(c.nsumt)
    #
    # check the number of distributions does not exceed the array size
    #

    # c$$$  RH: nfiles is declared in symbols.inc, and set in
    # c$$$  setup_symbols_clusters() in setup_symbols.f, and it is the number
    # c$$$  of distributions
    c.ndist = symbols.nfiles
    if c.ndist > mdist // 6:
        c.ndist = mdist // 6
    navbondedsave = real_array(c.ndist)
    navnonbondedsave = real_array(c.ndist)
    avbondsave = real_array(c.ndist)
    navbondedsave2 = real_array(c.ndist)
    navnonbondedsave2 = real_array(c.ndist)
    avbondsave2 = real_array(c.ndist)
    #
    # Two sets of distributions are calculated:
    # 1) a histogram of bonded neighbours
    # 2) a histogram of lsis
    #
    c.ndist3 = 10 * c.ndist
    #
    # setup the type numbers for the three atoms used to define a bond
    # also the distance ranges 1 - 3 and 2 - 3 needed to define the bond
    #
    for id_ in range(1, c.ndist + 1):
        iref = ireff(symbols.nsymboluptodata + 1, id_)
        commands.sometext = fstr(symbols.symbol_value[iref], 3020)
        _parse5()
        # Save the number of words read to test the other values
        c.nwordsave[id_] = commands.nwords
        for iw in range(1, commands.nwords + 1):
            c.atom1[id_, iw] = fstr(_word(iw), 4)
        iref = iref + 1
        commands.sometext = fstr(symbols.symbol_value[iref], 3020)
        _parse5()
        if commands.nwords >= c.nwordsave[id_]:
            for iw in range(1, commands.nwords + 1):
                c.atom2[id_, iw] = fstr(_word(iw), 4)
        else:
            units.write(6, fmt(_F200, id_, Ch(symbols.symbol[iref], 11)))
            return
        iref = iref + 1
        commands.sometext = fstr(symbols.symbol_value[iref], 3020)
        _parse5()
        if commands.nwords >= c.nwordsave[id_]:
            for iw in range(1, commands.nwords + 1):
                c.atom3[id_, iw] = fstr(_word(iw), 4)
        else:
            units.write(6, fmt(_F200, id_, Ch(symbols.symbol[iref], 11)))
            return
        # dff1, dff2, dfs1, dfs2, dft1, dft2 (the same block six times in the Fortran)
        for arr in (c.dff1, c.dff2, c.dfs1, c.dfs2, c.dft1, c.dft2):
            iref = iref + 1
            commands.sometext = fstr(symbols.symbol_value[iref], 3020)
            _parse5()
            if commands.nwords >= c.nwordsave[id_]:
                for iw in range(1, commands.nwords + 1):
                    arr[id_, iw] = read_internal(_word(iw)).real(arr[id_, iw])
            else:
                units.write(6, fmt(_F200, id_, Ch(symbols.symbol[iref], 11)))
                return
        # Read the minimum and maximum angles for defining a bond. Please note that since we
        # use the cosine rule, the minimum angle is actually the maximum in cosine, and the maximum angle
        # is the minimum in cosine. The ranges are therefore read in reverse order to that specified in the symbols table.
        for arr in (c.ang213max, c.ang213min):
            iref = iref + 1
            commands.sometext = fstr(symbols.symbol_value[iref], 3020)
            _parse5()
            if commands.nwords >= c.nwordsave[id_]:
                for iw in range(1, commands.nwords + 1):
                    arr[id_, iw] = read_internal(_word(iw)).real(arr[id_, iw])
                    arr[id_, iw] = cosf(arr[id_, iw] * piconv)
            else:
                units.write(6, fmt(_F200, id_, Ch(symbols.symbol[iref], 11)))
                return
        # Read in the maximum inter-atomic distance to be used when calculating triangle or q values.
        # 0.0 signals to use all distances.
        iref = iref + 1
        commands.sometext = fstr(symbols.symbol_value[iref], 3020)
        _parse5()
        c.rtrimax = read_internal(_word(1)).real(c.rtrimax)
        # Read current average number of bonded molecules per atom1, average number of non-bonded
        # molecules per atom1, and average bond length for this distribution (same block three times)
        for arr, arr2 in ((navbondedsave, navbondedsave2), (navnonbondedsave, navnonbondedsave2),
                          (avbondsave, avbondsave2)):
            iref = iref + 1
            commands.sometext = fstr(symbols.symbol_value[iref], 3020)
            _parse5()
            if commands.nwords >= 2:
                r = read_internal(symbols.symbol_value[iref])
                arr[id_] = r.real(arr[id_])
                arr2[id_] = r.real(arr2[id_])
                arr2[id_] = arr2[id_] * arr2[id_]
            else:
                arr[id_] = read_internal(symbols.symbol_value[iref]).real(arr[id_])
                arr2[id_] = 0.0
    units.write(6, lw(c.nz, c.ndist))
    c.zstep = (c.lsimax - c.lsimin) / f32(c.nz - 1)
    c.zhalf = r4(0.5) * c.zstep
    for j in range(1, c.nz + 1):
        c.z[j] = c.lsimin + f32(j - 1) * c.zstep

    read_and_save_files.initgofz(c.nz, c.ndist3, c.gofz, mz, mdist)

    # normalise atomic coordinates to box size

    ato_read_and_save.normalise_coordinates()
    #
    # calculate the average value and std dev for each quantity
    #
    if c.nsumt > 0:
        fc = 1.0 / float(c.nsumt)
    else:
        fc = 1.0
    fc1 = 1.0 - fc

    # set up the table of lsi indices

    find_lsi(fc, fc1, navbondedsave, navnonbondedsave, avbondsave, navbondedsave2, navnonbondedsave2, avbondsave2)

    # axis title

    c.titler = "lsi"
    idref10 = 0
    for id_ in range(1, c.ndist + 1):
        for iw in range(1, c.nwordsave[id_] + 1):
            units.write(6, lw(id_, Ch(c.atom1[id_, iw], 4), " ", Ch(c.atom2[id_, iw], 4), " ", Ch(c.atom3[id_, iw], 4)))
            units.write(6, lw(c.dff1[id_, iw], c.dff2[id_, iw], c.dfs1[id_, iw], c.dfs2[id_, iw]))
        idref1 = idref10 + 1
        idref2 = idref1 + 1
        idref3 = idref2 + 1
        idref4 = idref3 + 1
        idref5 = idref4 + 1
        idref6 = idref5 + 1
        idref7 = idref6 + 1
        idref8 = idref7 + 1
        idref9 = idref8 + 1
        idref10 = idref9 + 1

        # set the column title for this pair of atoms

        ilent1 = findex_padded(c.atom1[id_, 1], " ", 4) - 1
        ilent2 = findex_padded(c.atom2[id_, 1], " ", 4) - 1
        ilent3 = findex_padded(c.atom3[id_, 1], " ", 4) - 1
        a1 = substr(c.atom1[id_, 1], 1, ilent1)
        a2 = substr(c.atom2[id_, 1], 1, ilent2)
        a3 = substr(c.atom3[id_, 1], 1, ilent3)
        a3l1 = substr(c.atom3[id_, 1], 1, ilent1)  # atom3(id,1)(1:ilent1), as in the Fortran
        c.titledat[idref1] = fstr(fmt(_F100, a1, a2, a3), 29)
        c.titledat[idref2] = c.titledat[idref1]
        c.titledat[idref3] = c.titledat[idref1]
        c.titledat[idref4] = fstr(fmt(_F101, a3, a1, a3 + "_1"), 29)
        c.titledat[idref5] = fstr(fmt(_F101, a3, a1, a3 + "_2"), 29)
        c.titledat[idref6] = fstr(fmt(_F101, a3, a1, a3 + "_3"), 29)
        c.titledat[idref7] = fstr(fmt(_F102, a1, a3l1, "r3"), 29)
        c.titledat[idref8] = fstr(fmt(_F102, a1, a3l1, "r4"), 29)
        c.titledat[idref9] = fstr(fmt(_F102, a1, a3l1, "r5"), 29)
        c.titledat[idref10] = fstr(fmt(_F102, a1, a3l1, "r6"), 29)

    units.write(6, lw(*[x for id_ in range(1, idref10 + 1) for x in (Ch(c.titledat[id_], 29), " ")]))

    # output the data

    c.ext = "n"

    c.nsumt = read_and_save_files.outgofz(symbols.fnameinp, c.ext, c.nz, c.ndist3, c.nsumt, c.z, c.gofz, c.gzt,
                                          c.gzt2, c.titler, c.titledat, mz, mdist)
    for id_ in range(1, c.ndist + 1):
        # Calculate average bond length, average number of bonded molecules, an average number of non-bonded molecules
        pass
    # Calculate uncertainties on these numbers
    for id_ in range(1, c.ndist + 1):
        navbondedsave2[id_] = sqrtf(navbondedsave2[id_])
        navnonbondedsave2[id_] = sqrtf(navnonbondedsave2[id_])
        avbondsave2[id_] = sqrtf(avbondsave2[id_])
    #
    # write input file
    #
    #
    # update the current value of nsumt
    #
    icount = 0
    text = "nsize"
    icount = file_and_directory_commands.search_for_symbol(text, icount)
    iref = ireff(icount, 1)
    symbols.symbol_value[iref] = fstr(lw(c.nz), 256)
    iref = iref + 2
    symbols.symbol_value[iref] = fstr(lw(c.nsumt), 256)
    iref = iref + 1
    symbols.symbol_value[iref] = fstr(lw(c.ndist), 256)
    for id_ in range(1, c.ndist + 1):
        iref = ireff(symbols.nsymboluptodata + 13, id_)
        symbols.symbol_value[iref] = fstr(fmt("(2(1x,e12.5))", navbondedsave[id_], navbondedsave2[id_]), 256)
        iref = iref + 1
        symbols.symbol_value[iref] = fstr(fmt("(2(1x,e12.5))", navnonbondedsave[id_], navnonbondedsave2[id_]), 256)
        iref = iref + 1
        symbols.symbol_value[iref] = fstr(fmt("(2(1x,e12.5))", avbondsave[id_], avbondsave2[id_]), 256)

    read_and_save_files.save_files_others(14)

    return


def _dfs_range(id_, iw):
    """dfs1(id,iw), dfs2(id,iw), including iw = 0.

    NB: iwmin(im1,icount) is 0 for a neighbour for which no atom pair matched; the Fortran
    then reads dfs1(id,0) and dfs2(id,0), which in COMMON /lsi_r/ (dff1,dff2,dfs1,dfs2,...)
    are dff2(id,mdist) and dfs1(id,mdist) (Rust lsi.rs dfs_range).
    """
    c = lsicom
    if iw >= 1:
        return c.dfs1[id_, iw], c.dfs2[id_, iw]
    return c.dff2[id_, mdist], c.dfs1[id_, mdist]


def _atom_xyz(ia):
    """xa(ia), ya(ia), za(ia), including ia = 0.

    NB: ia2min(im1,icount) is 0 for a neighbour for which no atom pair matched; the Fortran
    then reads xa(0), ya(0), za(0), which in COMMON /molecom_r/ are the words before xa(1),
    ya(1), za(1): stepmi, xa(matom), ya(matom) (Rust lsi.rs atom_xyz).
    """
    mol = molecom
    if ia >= 1:
        return mol.xa[ia], mol.ya[ia], mol.za[ia]
    return mol.stepmi, mol.xa[matom], mol.ya[matom]


def _amax(a, b):
    """gfortran max(a,b) for reals: b if a is NaN (Rust lsi.rs fmax)."""
    if a != a or b > a:
        return b
    return a


def _amin(a, b):
    """gfortran min(a,b) for reals: b if a is NaN (Rust lsi.rs fmin)."""
    if a != a or b < a:
        return b
    return a


@routine
def find_lsi(fc, fc1, navbondedsave, navnonbondedsave, avbondsave, navbondedsave2, navnonbondedsave2, avbondsave2):
    """lsi.f:find_lsi -- calculates the distribution of local structure indices (LSI) for
    a set of bonded/non-bonded molecules (Shi, Russo, and Tanaka, J. Chem. Phys. 149,
    224502 (2018); Russo and Tanaka, 2014 NATURE COMMUNICATIONS DOI: 10.1038/ncomms4556).

    The save arrays are updated in place; fc, fc1 are real*8.
    """
    c = lsicom
    mol = molecom
    cl = clustercom
    match = chains_routines.match

    # logical bneighb(mmol,mbond) (lsi.inc, not in COMMON)
    bneighb = logical_array(mmol, mbond)

    # integer icount,icount1,icountmax,ia1save,ia2save,iwsave,nbonded,nbondedsum
    # integer lenbondsum,lennonbondsum,idref1,idref2,idref3,idref4,idref5,idref6

    # real mindist,anglefound,dnonbondedmin,dbondedav,dbondedsum
    # real lsi(mmol),r3(mmol),r4(mmol),r5(mmol),r6(mmol)
    lsi = real_array(mmol)
    r3 = real_array(mmol)
    r4_ = real_array(mmol)  # r4(mmol): renamed, r4 is the single-precision literal helper
    r5 = real_array(mmol)
    r6 = real_array(mmol)
    # real avbondmol(mmol),nneighbbond(mmol),nneighbnonbond(mmol)
    avbondmol = real_array(mmol)
    nneighbbond = real_array(mmol)
    nneighbnonbond = real_array(mmol)
    # real navbonded,navnonbonded,avbond
    # real navbondedsave(*),navnonbondedsave(*),avbondsave(*)
    # real navbondedsave2(*),navnonbondedsave2(*),avbondsave2(*)
    # real*8 fc,fc1

    c.zstep = (c.lsimax - c.lsimin) / f32(c.nz - 1)
    pi = PI4ATAN  # 4.0*atan(1.0)
    piconv = pi / r4(180.0)  # noqa: F841
    dmaxsqr = cl.dmax * cl.dmax
    idref10 = 0
    for id_ in range(1, c.ndist + 1):
        idref1 = idref10 + 1
        idref2 = idref1 + 1
        idref3 = idref2 + 1
        idref4 = idref3 + 1
        idref5 = idref4 + 1
        idref6 = idref5 + 1
        idref7 = idref6 + 1
        idref8 = idref7 + 1
        idref9 = idref8 + 1
        idref10 = idref9 + 1
        c.gofz[:] = 0.0  # gofz=0.0 (the whole array, for every distribution)
        #
        # First set up for each molecule a list of bonded neighbours.
        #
        for i in range(1, cl.nmol + 1):
            c.nneighb[i] = 0
            c.dneighb[i, 1] = 0
            nneighbbond[i] = 0
        # step through all origin molecules

        nfsum = 0  # noqa: F841
        for im1 in range(1, cl.nmol + 1):
            xc1 = cl.xadd[im1]
            yc1 = cl.yadd[im1]
            zc1 = cl.zadd[im1]
            for im2 in range(1, cl.nmol + 1):
                # Need to find the minimum separation between atoms atom1 on molecule 1 and atom 3 on molecule 2
                # The two atoms cannot be on the same molecule
                if im1 != im2:
                    xc2 = cl.xadd[im2]
                    yc2 = cl.yadd[im2]
                    zc2 = cl.zadd[im2]
                    mindist = cl.dmax
                    ia1save = 0
                    ia2save = 0
                    iwsave = 0
                    iw = 0
                    while iw < c.nwordsave[id_]:
                        iw = iw + 1
                        ia1 = int(cl.naf[im1]) - 1
                        while ia1 < cl.nal[im1]:
                            ia1 = ia1 + 1
                            foundia1 = match(mol.type[mol.itype[ia1]], c.atom1[id_, iw])
                            if foundia1:
                                xa1 = mol.xa[ia1] + xc1
                                ya1 = mol.ya[ia1] + yc1
                                za1 = mol.za[ia1] + zc1
                                ia2 = int(cl.naf[im2]) - 1
                                while ia2 < cl.nal[im2]:
                                    ia2 = ia2 + 1
                                    foundia2 = match(mol.type[mol.itype[ia2]], c.atom3[id_, iw])
                                    if foundia2:
                                        # Get the coordinates of the third atom
                                        xa2 = mol.xa[ia2] + xc2
                                        ya2 = mol.ya[ia2] + yc2
                                        za2 = mol.za[ia2] + zc2
                                        # Calculate the distance between ia1 and ia2
                                        xdif = xa2 - xa1
                                        ydif = ya2 - ya1
                                        zdif = za2 - za1
                                        # apply periodic boundary conditions
                                        xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                                        dsq13 = (xdif * xdif + ydif * ydif + zdif * zdif) * dmaxsqr
                                        da13 = sqrtf(dsq13)
                                        if da13 < mindist:
                                            mindist = da13
                                            ia1save = ia1
                                            ia2save = ia2
                                            iwsave = iw

                    # Step through existing distances for im1 and if this distance is less than any existing distances,
                    # or the total number of distances is less than mbond, insert the relevant value in order
                    icount = 0
                    foundneighb = False
                    while icount < c.nneighb[im1] and not foundneighb:
                        icount = icount + 1
                        if mindist < c.dneighb[im1, icount]:
                            # Move all the existing values up by one, up to the limit mbond, but only if im2 is
                            # different from im1 - otherwise simply replace the current value
                            icountmax = min(mbond, int(c.nneighb[im1]) + 1)
                            icount2 = icountmax
                            while icount2 > icount:
                                icount1 = icount2 - 1
                                c.neighb[im1, icount2] = c.neighb[im1, icount1]
                                c.dneighb[im1, icount2] = c.dneighb[im1, icount1]
                                c.ia1min[im1, icount2] = c.ia1min[im1, icount1]
                                c.ia2min[im1, icount2] = c.ia2min[im1, icount1]
                                c.iwmin[im1, icount2] = c.iwmin[im1, icount1]
                                icount2 = icount1
                            # Finally insert the new value
                            c.nneighb[im1] = icountmax
                            c.neighb[im1, icount] = im2
                            c.dneighb[im1, icount] = mindist
                            c.ia1min[im1, icount] = ia1save
                            c.ia2min[im1, icount] = ia2save
                            c.iwmin[im1, icount] = iwsave
                            foundneighb = True
                            #                     icount=icountmax !Make sure it doesn't add
                    # If no values were inserted, and number of bonds is less that maximum, then insert this value at
                    # the end.
                    if not foundneighb and c.nneighb[im1] < mbond:
                        c.nneighb[im1] = c.nneighb[im1] + 1
                        c.neighb[im1, c.nneighb[im1]] = im2
                        c.dneighb[im1, c.nneighb[im1]] = mindist
                        c.ia1min[im1, c.nneighb[im1]] = ia1save
                        c.ia2min[im1, c.nneighb[im1]] = ia2save
                        c.iwmin[im1, c.nneighb[im1]] = iwsave
            # Now find out which of the nearest neighbours of im1 are bonded to it
            for icount in range(1, c.nneighb[im1] + 1):
                iw = int(c.iwmin[im1, icount])
                da13 = c.dneighb[im1, icount]
                # First check that the distance of this pair is within the specified distance
                dfs1, dfs2 = _dfs_range(id_, iw)
                found = da13 >= dfs1 and da13 <= dfs2
                foundia2 = False
                if found:
                    dsq13 = da13 * da13
                    ia1 = int(c.ia1min[im1, icount])
                    im2 = int(c.neighb[im1, icount])
                    ia3 = int(c.ia2min[im1, icount])
                    xc3 = cl.xadd[im2]
                    yc3 = cl.yadd[im2]
                    zc3 = cl.zadd[im2]
                    xa3 = xc3 + mol.xa[ia3]
                    ya3 = yc3 + mol.ya[ia3]
                    za3 = zc3 + mol.za[ia3]
                    ia2 = int(cl.naf[im1]) - 1
                    while ia2 < cl.nal[im1] and not foundia2:
                        ia2 = ia2 + 1
                        # Make sure atom2 is the right type
                        foundia2 = match(mol.type[mol.itype[ia2]], c.atom2[id_, iw])
                        if foundia2:
                            # Check the separation of ia1 and ia2
                            xdif = mol.xa[ia2] - mol.xa[ia1]
                            ydif = mol.ya[ia2] - mol.ya[ia1]
                            zdif = mol.za[ia2] - mol.za[ia1]
                            dsq12 = (xdif * xdif + ydif * ydif + zdif * zdif) * dmaxsqr
                            da12 = sqrtf(dsq12)
                            foundia2 = da12 >= c.dff1[id_, iw] and da12 <= c.dff2[id_, iw]
                            if foundia2:
                                # Check the separation of ia2 and ia3
                                xa2 = xc1 + mol.xa[ia2]
                                ya2 = yc1 + mol.ya[ia2]
                                za2 = zc1 + mol.za[ia2]
                                xdif = xa3 - xa2
                                ydif = ya3 - ya2
                                zdif = za3 - za2
                                xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                                dsq23 = (xdif * xdif + ydif * ydif + zdif * zdif) * dmaxsqr
                                da23 = sqrtf(dsq23)
                                foundia2 = da23 >= c.dft1[id_, iw] and da23 <= c.dft2[id_, iw]
                                if foundia2:
                                    # Finally need to check the angle between da12 and da13
                                    anglefound = (dsq12 + dsq13 - dsq23) / (r4(2.0) * da12 * da13)
                                    foundia2 = anglefound >= c.ang213min[id_, iw] and anglefound <= c.ang213max[id_, iw]
                bneighb[im1, icount] = found and foundia2
        # For each non-bonded neighbour of a molecule, we need to check that that molecule is not bonded to the first
        # molecule
        for im1 in range(1, cl.nmol + 1):
            for icount1 in range(1, c.nneighb[im1] + 1):
                if not bneighb[im1, icount1]:
                    im2 = int(c.neighb[im1, icount1])
                    # Is im2 bonded to im1?
                    icount2 = 0
                    found = False
                    while icount2 < c.nneighb[im2] and not found:
                        icount2 = icount2 + 1
                        found = c.neighb[im2, icount2] == im1 and bneighb[im2, icount2]
                    bneighb[im1, icount1] = found
            # Calculate the lsi for this molecule
            dnonbondedmin = cl.dmax
            dbondedav = f32(0.0)
            nbonded = 0
            for icount in range(1, c.nneighb[im1] + 1):
                if bneighb[im1, icount]:
                    nbonded = nbonded + 1
                    dbondedav = (f32(nbonded - 1) * dbondedav + c.dneighb[im1, icount]) / f32(nbonded)
                else:
                    if c.dneighb[im1, icount] < dnonbondedmin:
                        dnonbondedmin = c.dneighb[im1, icount]
            lsi[im1] = dnonbondedmin - dbondedav
            r3[im1] = c.dneighb[im1, 3] - dbondedav
            r4_[im1] = c.dneighb[im1, 4] - dbondedav
            r5[im1] = c.dneighb[im1, 5] - dbondedav
            r6[im1] = c.dneighb[im1, 6] - dbondedav
            units.write(6, fmt("(1x,i5,1x,i5,1x,f10.4)", im1, c.nneighb[im1], lsi[im1]))
            units.write(6, fmt("(5(1x,i5,1x,f10.4,1x,l))",
                               *[x for j in range(1, c.nneighb[im1] + 1)
                                 for x in (c.neighb[im1, j], c.dneighb[im1, j], bool(bneighb[im1, j]))]))

        # Calculate the number of bonds, non-bonds, and average bond length for each molecule

        for im1 in range(1, cl.nmol + 1):
            nneighbbond[im1] = 0
            nneighbnonbond[im1] = 0
            avbondmol[im1] = 0.0
            if c.nneighb[im1] >= 0:
                for icount in range(1, c.nneighb[im1] + 1):
                    if bneighb[im1, icount]:
                        nneighbbond[im1] = nneighbbond[im1] + r4(1.0)
                        avbondmol[im1] = avbondmol[im1] + c.dneighb[im1, icount]
                    else:
                        iw = int(c.iwmin[im1, icount])
                        da13 = c.dneighb[im1, icount]
                        dfs1, dfs2 = _dfs_range(id_, iw)
                        if da13 >= dfs1 and da13 <= dfs2:
                            nneighbnonbond[im1] = nneighbnonbond[im1] + r4(1.0)
            if nneighbbond[im1] > 0:
                avbondmol[im1] = avbondmol[im1] / nneighbbond[im1]
        # Form averages
        navbonded = f32(0)
        navnonbonded = f32(0)
        avbond = f32(0)
        lensum = 0
        for im1 in range(1, cl.nmol + 1):
            navbonded = navbonded + nneighbbond[im1]
            navnonbonded = navnonbonded + nneighbnonbond[im1]
            if avbondmol[im1] > 0.0:
                avbond = avbond + avbondmol[im1]
                lensum = lensum + 1
        navbonded = navbonded / f32(cl.nmol)
        navnonbonded = navnonbonded / f32(cl.nmol)
        if lensum > 0:
            avbond = avbond / f32(lensum)
        units.write(6, fmt(_F994, navbonded, navnonbonded))
        units.write(6, fmt('("Average distance of bonded molecules: ",e12.5)', avbond))
        navbondedsave[id_] = f32(fc * float(navbonded) + fc1 * float(navbondedsave[id_]))
        navnonbondedsave[id_] = f32(fc * float(navnonbonded) + fc1 * float(navnonbondedsave[id_]))
        avbondsave[id_] = f32(fc * float(avbond) + fc1 * float(avbondsave[id_]))
        # Form average deviations
        navbonded = f32(0)
        navnonbonded = f32(0)
        avbond = f32(0)
        lensum = 0
        for im1 in range(1, cl.nmol + 1):
            diff = nneighbbond[im1] - navbondedsave[id_]
            navbonded = navbonded + diff * diff
            diff = nneighbnonbond[im1] - navnonbondedsave[id_]
            navnonbonded = navnonbonded + diff * diff
            if avbondmol[im1] > 0.0:
                diff = avbondmol[im1] - avbondsave[id_]
                avbond = avbond + diff * diff
                lensum = lensum + 1
        navbonded = navbonded / f32(cl.nmol)
        navnonbonded = navnonbonded / f32(cl.nmol)
        if lensum > 0:
            avbond = avbond / f32(lensum)
        navbondedsave2[id_] = f32(fc * float(navbonded) + fc1 * float(navbondedsave2[id_]))
        navnonbondedsave2[id_] = f32(fc * float(navnonbonded) + fc1 * float(navnonbondedsave2[id_]))
        avbondsave2[id_] = f32(fc * float(avbond) + fc1 * float(avbondsave2[id_]))

        # Histogram number of bonds, lsi and r5 for each molecule
        for im1 in range(1, cl.nmol + 1):
            izref = 1 + nint(nneighbbond[im1])
            c.gofz[izref, idref1] = c.gofz[izref, idref1] + r4(1.0)
            izref = 1 + nint(nneighbnonbond[im1])
            c.gofz[izref, idref2] = c.gofz[izref, idref2] + r4(1.0)
            for val, idr in ((lsi[im1], idref3), (r3[im1], idref7), (r4_[im1], idref8), (r5[im1], idref9),
                             (r6[im1], idref10)):
                izref = nint((val - c.lsimin + c.zhalf) / c.zstep) + 1
                if izref > 0 and izref < c.nz + 1:
                    c.gofz[izref, idr] = c.gofz[izref, idr] + r4(1.0)

        # Calculate the triangle distributions

        astep = pi / f32(c.nz)
        ahalf = r4(0.5) * astep
        qstep = r4(4.0) / f32(c.nz)
        qhalf = r4(0.5) * qstep  # noqa: F841
        for im1 in range(1, cl.nmol + 1):
            # q2found, q3found, q4found, n2found ... n5found are only used in commented-out code
            for icount2 in range(1, min(3, int(c.nneighb[im1]) - 1) + 1):
                ia1 = int(c.ia1min[im1, icount2])  # noqa: F841
                da12 = c.dneighb[im1, icount2]
                dsq12 = da12 * da12
                im2 = int(c.neighb[im1, icount2])
                ia2 = int(c.ia2min[im1, icount2])
                xc2 = cl.xadd[im2]
                yc2 = cl.yadd[im2]
                zc2 = cl.zadd[im2]
                xq, yq, zq = _atom_xyz(ia2)
                xa2 = xc2 + xq
                ya2 = yc2 + yq
                za2 = zc2 + zq
                for icount3 in range(icount2 + 1, min(4, int(c.nneighb[im1])) + 1):
                    da13 = c.dneighb[im1, icount3]
                    if (c.rtrimax > 0.0 and da12 <= c.rtrimax and da13 <= c.rtrimax) or c.rtrimax == 0.0:
                        dsq13 = da13 * da13
                        im3 = int(c.neighb[im1, icount3])
                        ia3 = int(c.ia2min[im1, icount3])
                        xc3 = cl.xadd[im3]
                        yc3 = cl.yadd[im3]
                        zc3 = cl.zadd[im3]
                        xq, yq, zq = _atom_xyz(ia3)
                        xa3 = xc3 + xq
                        ya3 = yc3 + yq
                        za3 = zc3 + zq
                        # Distance between ia2 and ia3
                        xdif = xa3 - xa2
                        ydif = ya3 - ya2
                        zdif = za3 - za2
                        xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                        dsq23 = (xdif * xdif + ydif * ydif + zdif * zdif) * dmaxsqr
                        da23 = sqrtf(dsq23)  # noqa: F841
                        anglefound = (dsq12 + dsq13 - dsq23) / (r4(2.0) * da12 * da13)
                        anglefound = _amax(anglefound, r4(-1.0))
                        anglefound = _amin(anglefound, r4(+1.0))
                        # (the q2found/q3found/q4found accumulation is commented out in the Fortran)
                        izref = nint((acosf(anglefound) + ahalf) / astep)
                        izref = min(c.nz, izref)
                        izref = max(1, izref)
                        if icount3 == 2:
                            c.gofz[izref, idref4] = c.gofz[izref, idref4] + r4(1.0)
                        if icount3 == 3:
                            c.gofz[izref, idref5] = c.gofz[izref, idref5] + r4(1.0)
                        if icount3 == 4:
                            c.gofz[izref, idref6] = c.gofz[izref, idref6] + r4(1.0)
            # (the q-value histograms are commented out in the Fortran)
        # Normalise the distributions
        sum1 = f32(0.0)
        sum2 = f32(0.0)
        sum3 = f32(0.0)
        sum4 = f32(0.0)
        sum5 = f32(0.0)
        sum6 = f32(0.0)
        sum7 = f32(0.0)
        sum8 = f32(0.0)
        sum9 = f32(0.0)
        sum10 = f32(0.0)
        for iz in range(1, c.nz + 1):
            betar = sinf(pi * (f32(iz) - r4(0.5)) / f32(c.nz))
            c.gofz[iz, idref4] = c.gofz[iz, idref4] / betar
            c.gofz[iz, idref5] = c.gofz[iz, idref5] / betar
            c.gofz[iz, idref6] = c.gofz[iz, idref6] / betar
            sum1 = sum1 + c.gofz[iz, idref1]
            sum2 = sum2 + c.gofz[iz, idref2]
            sum3 = sum3 + c.gofz[iz, idref3]
            sum4 = sum4 + c.gofz[iz, idref4]
            sum5 = sum5 + c.gofz[iz, idref5]
            sum6 = sum6 + c.gofz[iz, idref6]
            sum7 = sum7 + c.gofz[iz, idref7]
            sum8 = sum8 + c.gofz[iz, idref8]
            sum9 = sum9 + c.gofz[iz, idref9]
            sum10 = sum10 + c.gofz[iz, idref10]
        sum3 = sum3 * c.zstep
        #      sum4=sum4*4.0/real(nz)
        sum4 = sum4 * r4(180.0) / f32(c.nz)
        sum5 = sum5 * r4(180.0) / f32(c.nz)
        sum6 = sum6 * r4(180.0) / f32(c.nz)
        sum7 = sum7 * c.zstep
        sum8 = sum8 * c.zstep
        sum9 = sum9 * c.zstep
        sum10 = sum10 * c.zstep
        sums = (None, sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10)
        idrefs = (None, idref1, idref2, idref3, idref4, idref5, idref6, idref7, idref8, idref9, idref10)
        for iz in range(1, c.nz + 1):
            for k in range(1, 11):
                if sums[k] > 0.0:
                    c.gofz[iz, idrefs[k]] = c.gofz[iz, idrefs[k]] / sums[k]

    # end of id loop

    return
