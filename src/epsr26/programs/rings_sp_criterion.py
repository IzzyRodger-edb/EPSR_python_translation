"""Translation of fortran_src/rings_sp_criterion.f: ``program standalone``,
``subroutine run_rings`` and ``subroutine find_rings`` -- estimate the number of rings
in a distribution of molecules according to the 'shortest path' method.

``setup_symbols_rings`` is in setup_symbols_generic.f (library/setup_symbols_generic.py,
imported inside the routine: it is translated by another work package).
``match`` and ``shortest_path`` are in chains_routines.f.
"""

from ..common import chains, clustercom, commands, molecom, symbols
from ..common.chains import mbond, mdist, mring, mz
from ..common.epsr_dimension import mmol
from ..fortran import PI4ATAN, Ch, cosf, f32, findex_padded, fmt, fstr, idiv, int_array, lw, r4
from ..fortran import read_internal, sqrtf, substr, units
from ..library import ato_read_and_save
from ..library import chains_routines
from ..library import file_and_directory_commands
from ..library import minimagec
from ..library import parse
from ..library import read_and_save_files
from ..trace import routine

_F200 = "('rings> Not found enough items in distribution ',i2,', symbol: ',a11)"
_F100 = "(4x,a,'-',a,'-',a)"
_F102 = "(1x,i3,6(1x,e12.5))"
_F994 = "(1x,i5,' molecules, with average of ',e12.5,' bonds/mol.')"
_F996 = "(1x,i5,' rings of average size = ',e12.5,' molecules')"


def main():
    """rings_sp_criterion.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_rings()


def _parse5():
    """call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)"""
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)


def _word(iw):
    """sometext(ncf(iw):ncl(iw))"""
    return substr(commands.sometext, commands.ncf[iw], commands.ncl[iw])


@routine
def run_rings():
    """rings_sp_criterion.f:run_rings -- program to estimate number of rings in a
    distribution of molecules according to the 'shortest path' method."""
    from ..library import setup_symbols_generic

    c = chains
    ireff = file_and_directory_commands.ireff
    # integer*4 ierrout
    # character*20 text

    pi = PI4ATAN  # 4.0*atan(1.0)
    piconv = pi / r4(180.0)
    #
    # setup the rings symbols
    #
    setup_symbols_generic.setup_symbols_rings()
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
    nout, ierrout = read_and_save_files.read_input_file(nin, nout, ierrout)
    #
    # if an error was encountered on read, exit now
    #
    if ierrout != 0:
        return
    #
    # set up the number of ring length number bins
    #
    iref = 2
    c.nz = read_internal(symbols.symbol_value[iref]).int(c.nz)
    if c.nz > mring:
        c.nz = mring
    #
    # set up current accumulator
    #
    iref = iref + 1
    c.nsumt = read_internal(symbols.symbol_value[iref]).int(c.nsumt)
    #
    # check that the number of distributions is not too big for array sizes
    #
    nfilesmax = mdist // 2
    if symbols.nfiles > nfilesmax:
        symbols.nfiles = nfilesmax
    c.ndist = symbols.nfiles
    #
    # three sets of distributions are calculated:
    # 1) a histogram of bonded neighbours
    # 2) a histogram of ring lengths
    #
    c.ndist3 = 2 * c.ndist
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
    units.write(6, lw(c.nz, c.ndist))
    c.zstep = r4(0.5) / f32(c.nz)
    c.zhalf = r4(0.5) * c.zstep
    for j in range(1, c.nz + 1):
        c.z[j] = f32(j - 1)

    read_and_save_files.initgofz(c.nz, c.ndist3, c.gofz, mz, mdist)

    # normalise atomic coordinates to box size

    ato_read_and_save.normalise_coordinates()

    # set up the table of molecular bonds

    find_rings()

    # axis title

    c.titler = "Ring_size"
    idref2 = 0
    for id_ in range(1, c.ndist + 1):
        for iw in range(1, c.nwordsave[id_] + 1):
            units.write(6, lw(id_, Ch(c.atom1[id_, iw], 4), " ", Ch(c.atom2[id_, iw], 4), " ", Ch(c.atom3[id_, iw], 4)))
            units.write(6, lw(c.dff1[id_, iw], c.dff2[id_, iw], c.dfs1[id_, iw], c.dfs2[id_, iw]))
        idref1 = idref2 + 1
        idref2 = idref1 + 1

        # set the column title for this pair of atoms

        ilent1 = findex_padded(c.atom1[id_, 1], " ", 4) - 1
        ilent2 = findex_padded(c.atom2[id_, 1], " ", 4) - 1
        ilent3 = findex_padded(c.atom3[id_, 1], " ", 4) - 1
        c.titledat[idref1] = fstr(fmt(_F100, substr(c.atom1[id_, 1], 1, ilent1), substr(c.atom2[id_, 1], 1, ilent2),
                                      substr(c.atom3[id_, 1], 1, ilent3)), 29)
        c.titledat[idref2] = c.titledat[idref1]

    units.write(6, lw(*[x for id_ in range(1, idref2 + 1) for x in (Ch(c.titledat[id_], 29), " ")]))

    # output the data

    c.ext = "n"

    c.nsumt = read_and_save_files.outgofz(symbols.fnameinp, c.ext, c.nz, c.ndist3, c.nsumt, c.z, c.gofz, c.gzt,
                                          c.gzt2, c.titler, c.titledat, mz, mdist)

    #
    # calculate the average value and std dev for each quantity
    #
    for id_ in range(1, c.ndist3 + 1):
        # Calculate normalising constant

        anorm = f32(0.0)
        for i in range(1, c.nz + 1):
            anorm = anorm + c.gzt[i, id_]

        if anorm > 0:
            sum1 = f32(0.0)
            sum2 = f32(0.0)
            for i in range(1, c.nz + 1):
                term = c.z[i] * c.gzt[i, id_] / anorm
                sum1 = sum1 + term
                sum2 = sum2 + c.z[i] * term
            c.chainlen[id_, 1] = sum1
            std = sum2 - sum1 * sum1
            if std > 0.0:
                c.chainlen[id_, 2] = sqrtf(std)
            else:
                c.chainlen[id_, 2] = 0.0

        else:

            c.chainlen[id_, 1] = 0.0
            c.chainlen[id_, 2] = 0.0
    #
    # write input file
    #
    #
    # update the current value of nsumt
    #
    c.icount = 0
    text = "nsize"
    c.icount = file_and_directory_commands.search_for_symbol(text, c.icount)
    iref = ireff(c.icount, 1)
    symbols.symbol_value[iref] = fstr(lw(c.nz), 256)
    iref = iref + 1
    symbols.symbol_value[iref] = fstr(lw(c.nsumt), 256)
    iref = iref + 1
    symbols.symbol_value[iref] = fstr(lw(c.ndist), 256)

    read_and_save_files.save_files_others(14)

    # Output mean values (unit 14 was closed by save_files_others: these go to fort.14)

    for id_ in range(1, c.ndist + 1):
        idfirst = 2 * id_ - 1
        idlast = idfirst + 1
        units.write(14, fmt(_F102, id_, *[x for idr in range(idfirst, idlast + 1)
                                          for x in (c.chainlen[idr, 1], c.chainlen[idr, 2])]))

    units.close(14)

    return


@routine
def find_rings():
    """rings_sp_criterion.f:find_rings -- method used is based on the "Shortest Path"
    criterion of D S Franzblau, Phys Rev B, vol 44, p4925-4930 (1991). This seems a very
    reasonable criterion, applicable in a variety of situations."""
    c = chains
    mol = molecom
    cl = clustercom
    match = chains_routines.match

    # nring is the total number of rings found
    # lenring(*) gives the number of molecules in each ring

    # integer*4 nring,lenring(mmol),id,iw
    lenring = int_array(mmol)
    # integer, dimension(:,:), allocatable :: nspvaluet
    nspvaluet = None

    c.zstep = r4(0.5) / f32(c.nz)
    pi = r4(3.141592653)  # noqa: F841
    dmaxsqr = cl.dmax * cl.dmax
    idref2 = 0
    for id_ in range(1, c.ndist + 1):
        idref1 = idref2 + 1
        idref2 = idref1 + 1
        nring = 0
        #
        # First set up for each molecule a list of bonded neighbours.
        #
        for i in range(1, cl.nmol + 1):
            c.nneighb[i] = 0
        # step through all origin molecules

        nfsum = 0
        nssum = 0  # noqa: F841
        for im1 in range(1, cl.nmol + 1):
            ia1 = int(cl.naf[im1])
            ib1 = int(cl.naf[im1])
            for im2 in range(1, cl.nmol + 1):
                ia2 = int(cl.naf[im2])
                # Need to check that these two molecules are not the same and that they
                # are not already included in each others lists

                found = im1 == im2
                imt = 0
                while imt < c.nneighb[im1] and not found:
                    imt = imt + 1
                    found = im2 == c.neighb[im1, imt]

                # Also need to check there is room for more neighbours in the arrays

                found = found or c.nneighb[im1] == mbond or c.nneighb[im2] == mbond
                # Start the loop through sets of atoms
                foundiw = False or found
                iw = 0
                while not foundiw:
                    iw = iw + 1
                    foundia1 = False or found
                    ia1 = int(cl.naf[im1]) - 1
                    while not foundia1:
                        ia1 = ia1 + 1
                        foundia1 = match(mol.type[mol.itype[ia1]], c.atom1[id_, iw])
                        foundib1 = False or found
                        ib1 = int(cl.naf[im1]) - 1
                        while foundia1 and not foundib1:
                            ib1 = ib1 + 1
                            foundib1 = match(mol.type[mol.itype[ib1]], c.atom2[id_, iw])
                            if foundib1:
                                nfsum = nfsum + 1
                            # Get the coordinates of the first two atoms
                            xc1 = cl.xadd[im1]
                            yc1 = cl.yadd[im1]
                            zc1 = cl.zadd[im1]
                            xa1 = mol.xa[ia1] + xc1
                            ya1 = mol.ya[ia1] + yc1
                            za1 = mol.za[ia1] + zc1
                            xb1 = mol.xa[ib1] + xc1
                            yb1 = mol.ya[ib1] + yc1
                            zb1 = mol.za[ib1] + zc1
                            foundia2 = False or found
                            ia2 = int(cl.naf[im2]) - 1
                            while foundib1 and not foundia2:
                                ia2 = ia2 + 1
                                foundia2 = match(mol.type[mol.itype[ia2]], c.atom3[id_, iw])
                                if foundia2:
                                    # Get the coordinates of the third atom
                                    xc2 = cl.xadd[im2]
                                    yc2 = cl.yadd[im2]
                                    zc2 = cl.zadd[im2]
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
                                    # check that distance is within bounds prescribed for this distribution
                                    foundd13 = da13 >= c.dfs1[id_, iw] and da13 <= c.dfs2[id_, iw]
                                    if foundd13:
                                        # Distance between atoms 1 and 2
                                        xdif = xa1 - xb1
                                        ydif = ya1 - yb1
                                        zdif = za1 - zb1
                                        xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                                        dsq12 = (xdif * xdif + ydif * ydif + zdif * zdif) * dmaxsqr
                                        da12 = sqrtf(dsq12)
                                        foundd12 = da12 >= c.dff1[id_, iw] and da12 <= c.dff2[id_, iw]
                                        if foundd12:
                                            # Calculate the distance between ib1 and ia2
                                            xdif = xa2 - xb1
                                            ydif = ya2 - yb1
                                            zdif = za2 - zb1
                                            # apply periodic boundary conditions
                                            xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                                            dsq23 = (xdif * xdif + ydif * ydif + zdif * zdif) * dmaxsqr
                                            da23 = sqrtf(dsq23)
                                            foundd23 = da23 >= c.dft1[id_, iw] and da23 <= c.dft2[id_, iw]
                                            # Don't do the angle check if the angle range is silly or approximately =  pi
                                            if (foundd23 and c.ang213min[id_, iw] <= c.ang213max[id_, iw]
                                                    and (c.ang213max[id_, iw] - c.ang213min[id_, iw]) < r4(1.9)):
                                                anglefound = (dsq12 + dsq13 - dsq23) / (r4(2.0) * da12 * da13)
                                                found = (anglefound >= c.ang213min[id_, iw]
                                                         and anglefound <= c.ang213max[id_, iw])
                                            else:
                                                found = foundd23
                                        if found:
                                            # O.K. found a bond - save molecule numbers
                                            c.nneighb[im1] = c.nneighb[im1] + 1
                                            c.nneighb[im2] = c.nneighb[im2] + 1
                                            c.neighb[im1, c.nneighb[im1]] = im2
                                            c.neighb[im2, c.nneighb[im2]] = im1
                                foundia2 = ia2 >= cl.nal[im2] or found
                            foundib1 = ib1 >= cl.nal[im1] or found

                        # Ensure the loop will continue to the end in case no bonds have been found so far.

                        foundia1 = ia1 >= cl.nal[im1] or found
                    foundiw = iw >= c.nwordsave[id_] or found
                found = found and iw > 0
                #          if(found) write(6,*) id,im1,im2

        for im1 in range(1, cl.nmol + 1):
            units.write(6, lw("rings> Neighbours: ", im1, c.nneighb[im1],
                              *[c.neighb[im1, j] for j in range(1, c.nneighb[im1] + 1)]))

        # Calculate the average number of bonds per molecule

        lensum = 0
        nsum = 0
        for im1 in range(1, cl.nmol + 1):
            if c.nneighb[im1] >= 0:
                lensum = lensum + int(c.nneighb[im1])
                nsum = nsum + 1
        if nsum > 0:
            avbond = f32(lensum) / f32(nsum)
        else:
            avbond = f32(0.0)
        units.write(6, fmt(_F994, nsum, avbond))
        units.write(6, lw(nfsum))

        # maximum distance to calculate shortest path

        nspmax = idiv(c.nz, 2) + 1
        nspvaluet = int_array(cl.nmol, nspmax)

        # Starting from each molecule in turn look for rings of length not exceeding nz.
        # Cannot form rings containing less than three molecules.

        nring = 0
        im1 = 0

        while im1 < cl.nmol and nring < mmol:
            im1 = im1 + 1
            # !      do while(nring.lt.mmol)
            # !         im1=1

            # No chance of finding a ring if this molecule has less than 2 neighbours

            if c.nneighb[im1] > 1:

                # Initially allow all the molecules .ge. im1 to be included in the search

                for im in range(1, cl.nmol + 1):
                    if im < im1:
                        c.carryon[im] = False
                    else:
                        c.carryon[im] = True

                # Set up the shortest paths to this molecule

                chains_routines.shortest_path(im1, nspmax, cl.nmol, c.nneighb, c.neighb, c.rneighb1, c.rneighb2,
                                              c.spvaluet, nspvaluet, mmol, mbond)

                for im in range(1, cl.nmol + 1):
                    c.spvalue[im, 1] = c.spvaluet[im]

                # Now begin the search for rings proper.

                # Place molecule im1 in the stack

                c.nrlab[0] = im1
                c.nrlab[1] = im1
                c.spseq[0] = -1
                c.spseq[1] = 0
                c.rneighb[1] = 0
                c.spmin = int(c.spvalue[c.neighb[im1, 1], 1])
                # Determine the minimum shortest path for the neighbours of im1 - this should be unity for all neighbours!
                for itest in range(2, c.nneighb[im1] + 1):
                    if c.spvalue[c.neighb[im1, itest], 1] < c.spmin:
                        c.spmin = int(c.spvalue[c.neighb[im1, itest], 1])

                # Now search the neighbours of molecule nrlab(icount) if there are any left

                c.icount = 1

                # ibase will count the base to be used in the 'antipodal' shortest path check

                ibase = 1

                while c.icount > 0 and nring < mmol:

                    # If we have searched all the neighbours of the current stack pointer need to backtrack
                    # and try another molecule

                    # For the origin molecule, we do not need to use the last neighbour as all possible rings
                    # with this molecule will already have been found.

                    isub = 0
                    if c.icount == 1:
                        isub = 1
                    if c.rneighb[c.icount] >= (c.nneighb[c.nrlab[c.icount]] - isub):

                        # If this molecule is a neighbour of im1, ensure it will not be used again, so as to
                        # avoid double counting rings.

                        if c.icount == 2:
                            c.carryon[c.nrlab[c.icount]] = False
                        c.icount = c.icount - 1
                        if ibase > 1:
                            ibase = ibase - 1
                        # Get the shortest shortest path to im1 from the neighbours of nrlab(icount), but ignoring the
                        # previous member of the ring.
                        if c.icount > 1:
                            c.spmin = cl.nmol + 1
                            for itest in range(1, c.nneighb[c.nrlab[c.icount]] + 1):
                                if (c.spvalue[c.neighb[c.nrlab[c.icount], itest], 1] < c.spmin
                                        and c.neighb[c.nrlab[c.icount], itest] != c.nrlab[c.icount - 1]):
                                    c.spmin = int(c.spvalue[c.neighb[c.nrlab[c.icount], itest], 1])
                        else:
                            c.spmin = 1
                    else:

                        # Increment the neighbour counter for this member of the ring

                        c.rneighb[c.icount] = c.rneighb[c.icount] + 1

                        # Get the molecule number of the new neighbour

                        im2 = int(c.neighb[c.nrlab[c.icount], c.rneighb[c.icount]])

                        # im2 must be allowed, otherwise the rings associated with it have already been counted.
                        # To get the shortest path: we only consider neighbours that have the shortest path back to im1
                        # Molecules with undefined shortest paths are automatically eliminated
                        # Also the number of bonds associated with im2 must be greater than 1
                        # Also im2 cannot be equal to the one immediately preceeding icount in the ring
                        # Also don't waste time looking at rings longer than nz
                        # Also don't waste time if the ring counter has reached its maximum value

                        if (c.carryon[im2]
                                and c.spvalue[im2, 1] == c.spmin
                                and c.spvalue[im2, 1] < cl.nmol
                                and c.nneighb[im2] > 1
                                and im2 != c.nrlab[c.icount - 1]
                                and c.icount < c.nz
                                and nring < mmol):

                            # Check this molecule does not already exist in the ring (apart from the first one
                            # of course, in which case the ring is closed!!)

                            ifound = 0
                            ic = 1
                            while ifound == 0 and ic < c.icount:

                                ic = ic + 1
                                if im2 == c.nrlab[ic]:
                                    ifound = im2

                            if ifound == 0:

                                # Appear to have found a valid neighbour - use it.

                                # If im2 .eq. im1, the ring is closed.

                                if im2 == im1:

                                    # Before this ring can be accepted, need to prove that each member is on the shortest
                                    # path around the ring to every other member (bearing in mind there are two paths
                                    # around the ring to each member), i.e. there are no possible short cuts.

                                    found = False
                                    # Don't need to check interactions with im1 since this has already been tested. Also
                                    # only need to check 2nd neighbours and beyond.
                                    ic1 = 1
                                    while ic1 < c.icount - 2 and not found:
                                        ic1 = ic1 + 1
                                        ic2 = ic1 + 1
                                        # Set up the shortest paths to this ring molecule
                                        chains_routines.shortest_path(int(c.nrlab[ic1]), nspmax, cl.nmol, c.nneighb,
                                                                      c.neighb, c.rneighb1, c.rneighb2, c.spvaluet,
                                                                      nspvaluet, mmol, mbond)
                                        while ic2 < c.icount and not found:
                                            ic2 = ic2 + 1
                                            # Get the shortest path around the ring between these two molecules
                                            c.spmin = min((ic2 - ic1), (c.icount + ic1 - ic2))
                                            c.spmax = max((ic2 - ic1), (c.icount + ic1 - ic2))
                                            found = c.spmin > c.spvaluet[c.nrlab[ic2]]
                                            # Need to check that a short circuit does not occur between these two atoms.
                                            # This is manifest if spmax is .gt. spmin and there is more than one route of
                                            # length spmin between ic1 and ic2.
                                            # If spmin and spmax are the same it means the two molecules are opposite each
                                            # other in the ring. In that case there is nothing to test since a short
                                            # circuit will already have been tested for above.
                                            if not found and c.spmax > c.spmin:
                                                ic = c.spmin
                                                found = nspvaluet[c.nrlab[ic2], ic] > 1
                                                # Also need to check if there are other short circuits less than spmax.
                                                # In this case any number greater than 0 is a short circuit.
                                                while not found and ic < c.spmax - 1 and ic < nspmax:
                                                    ic = ic + 1
                                                    found = nspvaluet[c.nrlab[ic2], ic] > 0

                                    if not found and nring < mmol:

                                        nring = nring + 1
                                        lenring[nring] = c.icount
                                        units.write(6, lw())
                                        units.write(6, lw(im1, nring, lenring[nring],
                                                          *[c.nrlab[ic] for ic in range(1, c.icount + 1)]))
                                        units.write(6, lw(im1, nring, lenring[nring],
                                                          *[c.spseq[ic] for ic in range(1, c.icount + 1)]))
                                        units.write(6, lw())

                                else:
                                    c.icount = c.icount + 1
                                    c.nrlab[c.icount] = im2
                                    c.rneighb[c.icount] = 0
                                    c.spseq[c.icount] = c.spvalue[im2, 1]
                                    # Get the shortest shortest path to im1 from the neighbours of im2, but ignoring the
                                    # previous member of the ring which has already been included.
                                    if c.icount > 1:
                                        c.spmin = cl.nmol + 1
                                        for itest in range(1, c.nneighb[im2] + 1):
                                            if (c.spvalue[c.neighb[im2, itest], 1] < c.spmin
                                                    and c.neighb[im2, itest] != c.nrlab[c.icount - 1]):
                                                c.spmin = int(c.spvalue[c.neighb[im2, itest], 1])
                                    else:
                                        c.spmin = 1

            # As all possible rings associated with this molecule have been counted, it cannot
            # be used again

            c.carryon[im1] = False

        if nring > 0:
            lensum = 0
            for in_ in range(1, nring + 1):
                lensum = lensum + int(lenring[in_])
            avlen = f32(lensum) / f32(nring)
        else:
            avlen = f32(0.0)
        units.write(6, fmt(_F996, nring, avlen))

        # form histogram of number of bonded neighbours

        for im1 in range(1, cl.nmol + 1):
            iref = int(c.nneighb[im1]) + 1
            if iref <= c.nz:
                c.gofz[iref, idref1] = c.gofz[iref, idref1] + r4(1.0)

        # form histogram of ring lengths

        for ir in range(1, nring + 1):
            iref = int(lenring[ir]) + 1
            if iref <= c.nz:
                c.gofz[iref, idref2] = c.gofz[iref, idref2] + r4(1.0)

        # normalise the distributions to probabilities
        # (commented out in the Fortran)
    return
