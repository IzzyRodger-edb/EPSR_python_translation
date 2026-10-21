"""Translation of fortran_src/voids.f: ``program standalone`` (the ``voids`` program),
``run_voids``, ``find_voids``, ``ranmine``, ``occupypix``, ``markpixel``,
``search_voids``, ``get_void_list`` and ``spheregr``.

Void-size distribution: the (cubic) box is pixellated, pixels within a radius of an
atom are occupied, adjacent unoccupied pixels are grouped into voids; then the void
size histogram and the void-pixel pair correlation function (sampled from 100 random
origins, with random broadening) are accumulated.

The ``interface`` block for find_voids in run_voids (SIGNATURES lists it as a second
find_voids at line 24) only declares the routine; there is one find_voids.

``voidpix`` (integer*1) and ``voidnum`` (integer*2) are held in Python lists (index 0
unused): their values (0/1, and void numbers <= maxvoid = 32767) cannot overflow those
kinds, and list indexing is much faster than NumPy element access in the pixel loops.
The library modules are imported inside the routines.  See docs/translation_notes/O.md.
"""

from ..common import clustercom
from ..common import clusters
from ..common import commands
from ..common import lattice
from ..common import molecom
from ..common import ran1com
from ..common import symbols
from ..fortran import Ch
from ..fortran import FortranStop
from ..fortran import f32
from ..fortran import fmt
from ..fortran import fstr
from ..fortran import idiv
from ..fortran import ifix
from ..fortran import int_array
from ..fortran import lw
from ..fortran import nint
from ..fortran import powf
from ..fortran import r4
from ..fortran import read_internal
from ..fortran import sqrtf
from ..fortran import units
from ..trace import routine

mpix3 = clusters.mpix * clusters.mpix * clusters.mpix  # parameter (mpix3=mpix*mpix*mpix)
maxvoid = 2**15 - 1  # parameter (maxvoid=2**15-1)


def main():
    """voids.f:standalone -- program to perform the EPSR voids command.

    2009/3/30 Attempt to make this a standalone program
    """
    from ..library import file_and_directory_commands

    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_voids()

    raise FortranStop(0)


_F102_run = "('   Radius_',f5.3)"
_F101_run = "(/'Distribution ',i3,': fraction of pixels used = ',f10.5)"


@routine
def run_voids():
    """voids.f:run_voids -- programme to estimate the distribution of void sizes
    within a structure.

    Voids are estimated by pixellating the box into a grid of small cubes, then
    estimating the distance from each pixel to the nearest atom. If this distance is
    greater than a specified distance - a nominal radius of the atoms - then the pixel
    is marked as being in a void.  Then look at each pixel and check the surrounding 26
    neighbouring pixels. If any of them are occupied, then both pixels are given the
    same (normally the lowest non-zero) void number of the two of them.  Pixels that are
    not in a void are put in a void of zero volume.  Finally a histogram is built up of
    the number pixels in each void. These numbers are divided by the total number of
    pixels in the box, so that the effect of different pixel sizes can be checked.
    """
    from ..library import (
        ato_read_and_save,
        auxiliary_routines,
        file_and_directory_commands,
        parse,
        read_and_save_files,
        reallocation_routines,
        setup_symbols_generic,
    )

    # use reallocation_routines
    p = clusters
    ireff = file_and_directory_commands.ireff
    # integer*4 ierrout,nexcludemax
    # character*20 text
    nexclude = int_array(p.mdist)  # integer(kind=4), dimension(mdist) :: nexclude
    exclude = None  # character(len=3), dimension(:,:), allocatable :: exclude
    # character*256 texclude
    texclude = ""
    ilentexclude = 0
    # logical found

    nexcludemax = 0
    #
    # setup the clusters symbols
    #
    setup_symbols_generic.setup_symbols_voids()

    #
    # now setup the input file
    #
    nin = 11
    nout = 6
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
    if not lattice.cubic:
        units.write(6, lw("voids> Cannot currently process non-cubic boxes!! Apologies!!"))
        return
    #
    # set up the number of voids size bins
    #
    iref = 2
    p.nz = read_internal(symbols.symbol_value[iref]).int(p.nz)
    p.nz = min(p.nz, p.mz)
    #
    # set up number of pixels
    #
    iref = iref + 1
    p.npix = read_internal(symbols.symbol_value[iref]).int(p.npix)

    # Ensure the total number of pixels across the box is an odd number, so there is
    # always a pixel at the centre of the box.

    p.npix = idiv(p.npix, 2)
    npixtotal = 2 * p.npix + 1
    if npixtotal > p.mpix:
        p.npix = (p.mpix - 1) // 2
    units.write(6, lw(p.npix, npixtotal))

    # Setup maximum void volume in units of Angstroms^3

    iref = iref + 1
    p.voidmax = read_internal(symbols.symbol_value[iref]).real(p.voidmax)
    #
    # set up current accumulator
    #
    iref = iref + 1
    p.nsumt = read_internal(symbols.symbol_value[iref]).int(p.nsumt)
    #
    # check the number of distributions does not exceed the array size
    #
    p.ndist = symbols.nfiles
    #
    # Two sets of distributions are calculated:
    # 1) a histogram of void sizes
    # 2) a void-pixel pair correlation function.
    #
    ndist3 = 3 * p.ndist
    if ndist3 > p.mdist:
        ndist3 = p.mdist
    p.ndist = ndist3 // 3

    for id_ in range(1, p.ndist + 1):
        iref = ireff(symbols.nsymboluptodata + 1, id_)
        #
        # setup the distance used to define whether a pixel is within an atom
        #
        p.damax[id_] = read_internal(symbols.symbol_value[iref]).real(p.damax[id_])

        iref = iref + 1
        commands.sometext = fstr(symbols.symbol_value[iref], 3020)
        commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                      commands.nstsize, commands.nworddim)
        if commands.nwords > nexcludemax:
            nexcludemax = commands.nwords
            exclude = reallocation_routines.reallocate2d_c(exclude, 3, nexcludemax, p.ndist)
        nexclude[id_] = 0
        ie = 0
        while ie < commands.nwords:
            ie = ie + 1
            # call concatenate(ie,ie,texclude,ilentexclude,nwords,ncf,ncl,sometext)
            # (ie is passed for both nfirst and nlast; it is within 1..nwords, so not clipped)
            _, ie, texclude, ilentexclude = parse.concatenate(ie, ie, texclude, ilentexclude, commands.nwords,
                                                              commands.ncf, commands.ncl, commands.sometext)
            # Check whether this atom type exists
            t3 = fstr(texclude, 3)  # texclude(1:3)
            found = False
            it = 0
            while it < molecom.ntype and not found:
                it = it + 1
                found = molecom.type[it] == t3
            if found:
                nexclude[id_] = nexclude[id_] + 1
                exclude[nexclude[id_], id_] = t3
        if nexclude[id_] > 0:
            units.write(6, lw(id_, p.damax[id_], *[Ch(exclude[ie, id_], 3) for ie in range(1, nexclude[id_] + 1)]))
        else:
            units.write(6, lw(id_, p.damax[id_]))

    units.write(6, lw(p.nz, p.ndist))

    # setup step size and z values
    # voidmax represents the maximum void volume sampled in units of the mean volume
    # per atom, i.e. 1/rho

    p.zstep = clustercom.dmax / f32(p.nz)
    p.zhalf = r4(0.5) * p.zstep
    for j in range(1, p.nz + 1):
        p.z[j] = f32(j) * p.zstep - p.zhalf

    read_and_save_files.initgofz(p.nz, ndist3, p.gofz, p.mz, p.mdist)

    #
    # normalise atomic coordinates to box size
    #
    ato_read_and_save.normalise_coordinates()

    # Initialise random number generator - this must be done here, after the .ato
    # file has been read, otherwise the random generator is always set to the values
    # read in from the .ato file. At the end of voids, the .ato file is not saved.

    units.open(14, "randomno.dat", status="unknown")
    auxiliary_routines.initran1(14)
    units.write(6, lw(ran1com.mcall, ran1com.idum, ran1com.iy))
    units.close(14)
    #
    # set up the table of molecular bonds
    #
    find_voids(nexclude, exclude)
    #
    # axis title
    #
    p.titler = "Void_radius"
    idref3 = 0
    for id_ in range(1, p.ndist + 1):
        idref1 = idref3 + 1
        idref2 = idref1 + 1
        idref3 = idref2 + 1
        #
        # set the column title for this pair of atoms
        #
        p.titledat[idref1] = fstr(fmt(_F102_run, p.damax[id_]), 29)
        p.titledat[idref2] = fstr(fmt(_F102_run, p.damax[id_]), 29)
        p.titledat[idref3] = fstr(fmt(_F102_run, p.damax[id_]), 29)

    # output the data

    p.ext = "n"

    p.nsumt = read_and_save_files.outgofz(symbols.fnameinp, p.ext, p.nz, ndist3, p.nsumt, p.z, p.gofz, p.gzt,
                                          p.gzt2, p.titler, p.titledat, p.mz, p.mdist)

    # For each distribution calculate the fraction of box volume occupied by voids

    dmax = clustercom.dmax
    boxv = dmax * dmax * dmax
    idref3 = 0
    for id_ in range(1, p.ndist + 1):

        idref1 = idref3 + 1
        idref2 = idref1 + 1
        idref3 = idref2 + 1
        sum1 = f32(0.0)
        for i in range(1, p.nz + 1):

            # Void linear dimension

            voidl = p.z[i]

            # Void volume

            voidv = voidl * voidl * voidl
            sum1 = sum1 + voidv * p.gzt[i, idref1]

        p.fract2[id_] = sum1 / boxv

    #
    # write input file
    #

    # update the current value of nsumt and ndist

    icount = 0
    text = "nsize"
    icount = file_and_directory_commands.search_for_symbol(text, icount)
    iref = ireff(icount, 1)
    symbols.symbol_value[iref] = fstr(lw(p.nz), 256)
    iref = iref + 1
    symbols.symbol_value[iref] = fstr(lw(p.nsumt), 256)
    iref = iref + 1
    symbols.symbol_value[iref] = fstr(lw(p.ndist), 256)

    read_and_save_files.save_files_others(14)

    for id_ in range(1, p.ndist + 1):
        #
        # for each void distribution indicate the fraction of unoccupied pixels
        #
        units.write(14, fmt(_F101_run, id_, p.fract2[id_]))

    units.close(14)

    # Save random number parameters.

    units.open(14, "randomno.dat", status="unknown")
    auxiliary_routines.saveran1(14)
    units.close(14)

    return


_F103 = ("(//'find_voids> Distribution ',i2"
         "//'find_voids> Adjusting number of pixels to give optimum'"
         "/ '            fit to specified void distance.',"
         "//'find_voids> New number of pixels = ', i8,"
         "//'find_voids> Specified radius = ',f10.5,"
         "//'find_voids> Actual radius used = ',f10.5)")
_F102 = "(/'find_voids> No. of pixels in sphere = ',i8)"
_F101 = ("(/'find_voids> ',i10,' pixels used out of ',i10,' pixels in total.'"
         ",//'find_voids> Fraction used = ',f10.5)")


@routine
def find_voids(nexclude, exclude):
    """voids.f:find_voids -- pixellate the box, mark the pixels within damax of the
    (selected) atoms, group the unoccupied pixels into voids, and accumulate the
    void-size histogram and the void-pixel pair correlation functions."""
    from ..library import minimagec

    p = clusters
    m = molecom
    c = clustercom
    # parameter (mpix3=mpix*mpix*mpix,maxvoid=2**15-1)
    voidpix = [0] * (mpix3 + 1)  # integer*1 voidpix(mpix3)
    voidnum = [0] * (mpix3 + 1)  # integer*2 voidnum(mpix3)
    # integer*4 ih,ik,il,npixtot,npixall,nvoidpix
    labvoid = int_array(maxvoid)
    lenvoid = int_array(maxvoid)
    # integer(kind=4), dimension(:) :: nexclude
    # character(len=3), dimension(:,:) :: exclude
    mi = minimagec.minimagec
    gofz = p.gofz
    dmax = c.dmax
    nz = p.nz
    # NB: npixfound is not initialised in the Fortran; it is set by occupypix whenever an
    # atom occupies pixels (always the case in the tests); 0 otherwise.
    npixfound = 0
    # NB: iscount is incremented before it is first set to 0 (the increment is then
    # overwritten); it is initialised to 0 here.
    iscount = 0

    nfdim = p.mdist  # noqa: F841  (not used)

    # Mean volume per atom

    volatom = r4(1.0) / f32(m.natom)

    p.zstep = r4(1.0) / f32(nz)
    zstep = p.zstep

    # Keep track of the two distributions

    idref3 = 0

    #
    # look at each requested atom radius
    #
    for id_ in range(1, p.ndist + 1):

        # Keep track of the three distributions

        idref1 = idref3 + 1
        idref2 = idref1 + 1
        idref3 = idref2 + 1

        # set up the test radius

        radtest = p.damax[id_] / dmax

        # Parameters to define the pixel positions

        pixstep = r4(1.0) / f32(2 * p.npix + 1)

        # No. of pixels needed to define sphere.

        npixsphere = nint(radtest / pixstep)

        # Need to ensure pixel size is adjusted to give closest possible match to
        # specified radius

        # Give the ideal pixel size

        pixstep = radtest / f32(npixsphere)

        # Convert this to the ideal number of pixels. The final number must be the
        # the nearest ODD integer.

        ratio = r4(1.0) / pixstep
        npixt = ifix(ratio / r4(2.0))

        # Check the remainder to be sure this the nearest ODD integer

        remain = ratio - f32(2 * npixt + 1)
        units.write(6, lw(ratio, npixt, remain))
        if remain > 1.0:
            npixt = npixt + 1

        if (2 * npixt + 1) > p.mpix:
            npixt = (p.mpix - 1) // 2

        # Get the correct pixel size

        pixstep = r4(1.0) / f32(2 * npixt + 1)

        # No. of pixels needed to define sphere.

        npixsphere = nint(radtest / pixstep)
        radactual = f32(npixsphere) * pixstep

        npixtot = 2 * npixt + 1

        # Total number of pixels

        npixall = npixtot * npixtot * npixtot

        units.write(6, fmt(_F103, id_, npixtot, radtest * dmax, radactual * dmax))

        # Mark all pixels as unoccupied initially

        for ipix in range(1, npixall + 1):

            voidpix[ipix] = 0
            voidnum[ipix] = 0

        # Step through all the atoms

        ia = 0
        while ia < m.natom:
            ia = ia + 1
            ifound = 0
            if nexclude[id_] > 0:  # If excluded atoms have been defined, the occupied pixels are defined ONLY for those atoms
                ie = 0
                while ie < nexclude[id_] and ifound == 0:
                    ie = ie + 1
                    if m.type[m.itype[ia]] == exclude[ie, id_]:
                        ifound = int(m.itype[ia])
                #      write(6,*) ia,ie,type(itype(ia)),exclude(ie,id),ifound
            else:  # Otherwise all atoms create occupied pixels.
                ifound = int(m.itype[ia])
            #            write(6,*) ia,ifound
            if ifound != 0:
                im = m.mref[ia]
                xatom = m.xa[ia] + c.xadd[im]
                yatom = m.ya[ia] + c.yadd[im]
                zatom = m.za[ia] + c.zadd[im]

                # Just as a precaution apply minimum image to make sure this atom is definitely
                # in the range +- 0.5 in each direction

                xatom, yatom, zatom = mi(xatom, yatom, zatom)

                # Find the closest pixel

                ih = nint(xatom / pixstep)
                ik = nint(yatom / pixstep)
                il = nint(zatom / pixstep)

                # Mark this pixel as occupied, plus all the surrounding pixels within the
                # sphere.

                #     write(6,*) xatom,yatom,zatom,ih,ik,il,nvoidpix
                npixfound = occupypix(npixt, npixsphere, npixfound, ih, ik, il, voidpix, mpix3)

        # If nexclude(id) > 0, we are interested in the way the 'excluded' atoms are distributed,
        # so instead of calculating the distribution of unoccuped pixels, we calculate the distribution of occupied pixels.
        # To use the same code we now need to define the currently occupied pixels as unoccupied, and the currently unoccupied ones as occupied.

        if nexclude[id_] > 0:
            for ipix in range(1, npixall + 1):
                if voidpix[ipix] != 0:
                    voidpix[ipix] = 0
                else:
                    voidpix[ipix] = 1

        units.write(6, fmt(_F102, npixfound))

        # Set the initial number of unoccupied pixels

        nvoidpix = 0

        # Count the number of unoccupied pixels

        for ipix in range(1, npixall + 1):

            if voidpix[ipix] == 0:
                nvoidpix = nvoidpix + 1

        # Write the total number of unoccupied pixels

        p.fract2[id_] = f32(nvoidpix) / f32(npixall)

        units.write(6, fmt(_F101, nvoidpix, npixall, p.fract2[id_]))

        # Now correlate neighbouring unoccupied pixels into voids.
        # This is done in three stages. First we make a list of the unoccupied pixels
        # that are immediately adjacent to each unoccupied pixel. Then we assign all
        # of these unoccupied, adjacent pixels the same (lowest) void number. Finally we
        # count up the total number of voids. The latter two processes are repeated
        # until there is no further reduction in the number of voids.

        # Initially we assume each pixel is its own void

        nvoidsold = nvoidpix
        nvoidsnew = 0
        it = 0

        while nvoidsold > nvoidsnew:

            if it > 0:
                nvoidsold = nvoidsnew

            nvoidsnew = search_voids(nvoidsnew, npixt, npixall, voidpix, voidnum, mpix3)

            it = it + 1
            units.write(6, lw(it, nvoidsold, nvoidsnew))

        # Finally the size of each void is tallied and formed into a histogram.

        # First zero the tally and label

        for ivoid in range(1, nvoidsnew + 1):

            labvoid[ivoid] = 0
            lenvoid[ivoid] = 0

        nvoid = 0

        for ipix in range(1, npixall + 1):

            # Get an unoccupied void with finite void number

            if voidpix[ipix] == 0 and voidnum[ipix] > 0:

                # Check that this void number does not already exist

                iflag = 0
                ivoid = 0
                while ivoid < nvoid and iflag == 0:

                    ivoid = ivoid + 1
                    if labvoid[ivoid] == voidnum[ipix]:
                        iflag = ivoid

                if iflag == 0:

                    nvoid = nvoid + 1
                    labvoid[nvoid] = voidnum[ipix]
                    lenvoid[nvoid] = 1

                else:

                    lenvoid[iflag] = lenvoid[iflag] + 1

        # Form a histogram of number of voids versus volume.

        ivoid = 0
        volpix = pixstep * pixstep * pixstep
        volpixstar = volpix / volatom
        third = r4(1.0) / r4(3.0)
        while ivoid < nvoid:

            ivoid = ivoid + 1

            if ivoid < 51:
                units.write(6, lw(ivoid, f32(lenvoid[ivoid]) * volpixstar))

            # The size of a void is expressed as the cube root of its volume

            iref = ifix(powf(f32(lenvoid[ivoid]) * volpix, third) / zstep) + 1
            if iref > nz:
                iref = nz
            gofz[iref, idref1] = gofz[iref, idref1] + r4(1.0)

        # Now calculate the radial distribution function of void pixels. We ensure we
        # only count each pair once.

        npixtot = 2 * npixt + 1
        npixtot2 = npixtot * npixtot
        ivcount = 0
        nvcount = 100
        ivoidsum = 0

        # In the event that the step in z is much larger than the pixel size,
        # we will broaden the bin that the data are entered into

        nbroad = nint(f32(nz) / f32(npixtot))
        # Choose an amount to broaden each point
        pixhwidth = r4(0.5) / f32(npixtot)
        if nbroad < 0:
            nbroad = 0

        if nvoidpix == 0:
            # The Fortran loops forever below (no unoccupied pixel to choose as origin).
            raise FortranStop(1, "voids: no unoccupied pixels (the Fortran loops forever here)")

        fnpixtot = f32(npixtot)
        fnz = f32(nz)
        fnpixall = f32(npixall)
        one = r4(1.0)
        two = r4(2.0)
        while ivcount < nvcount:

            ivcount = ivcount + 1

            # Choose the origin pixel at random

            itest = 0

            while itest == 0:

                iv1 = ifix(fnpixall * ranmine() + one)
                if voidpix[iv1] == 0:
                    itest = 1

            if voidpix[iv1] == 0:
                # Include the correlation of this pixel with itself.
                iscount = iscount + 1
                # Add the broadening component
                xdif = pixhwidth * (two * ranmine() - one)
                ydif = pixhwidth * (two * ranmine() - one)
                zdif = pixhwidth * (two * ranmine() - one)
                dist = sqrtf(xdif * xdif + ydif * ydif + zdif * zdif)
                iref = ifix(dist * fnz) + 1
                if iref < 1:
                    iref = 1
                if iref > nz:
                    iref = nz
                gofz[iref, idref2] = gofz[iref, idref2] + one
                ivoidsum = ivoidsum + 1

                # Get the coordinates of this pixel - convert these into fractional coordinates

                ilv1 = (iv1 - 1) // npixtot2
                ikv1 = (iv1 - 1 - ilv1 * npixtot2) // npixtot
                ihv1 = iv1 - 1 - npixtot * ikv1 - npixtot2 * ilv1
                xv1 = f32(ihv1) / fnpixtot
                yv1 = f32(ikv1) / fnpixtot
                zv1 = f32(ilv1) / fnpixtot

                iscount, ivoidsum = _pair_loop(iv1, npixall, npixtot, npixtot2, xv1, yv1, zv1, ihv1, ikv1, ilv1,
                                               pixhwidth, voidpix, voidnum, idref2, idref3, ivoidsum)
                units.write(6, lw(ivcount, iscount))

        # Normalise this g(r) to 4*pi*rho*r**2delta-r

        rhovoid = f32(nvcount * iscount)
        rvstep = r4(1.0) / f32(nz)
        pi = r4(3.141592653)
        fact = r4(3.0) / (r4(4.0) * pi * rhovoid)
        units.write(6, lw(nvcount, rhovoid, rvstep, pi, fact))
        fracfac1 = r4(1.0) - p.fract2[id_]
        fracfac2 = p.fract2[id_]  # noqa: F841  (used only in commented-out code)
        fracfac3 = fracfac1 * fracfac1  # noqa: F841
        rv1 = f32(0.0)
        rv13 = rv1 * (rv1 * rv1)
        for iz in range(1, nz + 1):

            rv2 = f32(iz) * rvstep
            rv23 = rv2 * (rv2 * rv2)
            rvfact = fact / (rv23 - rv13)
            gofz[iz, idref2] = gofz[iz, idref2] * rvfact
            gofz[iz, idref3] = gofz[iz, idref3] * rvfact

            # c Subtract the 'uninteresting' spherical g(r) term
            #
            #           rstar=real(iz)/(real(nz)*radactual)
            #     write(6,*) iz,rstar,spheregr(rstar),fracfac1,fracfac2,fracfac3
            #           if(fract2(id).gt.0.0005) then
            #              gofz(iz,idref3)=spheregr(rstar)/frac2(id)
            #           else
            #              gofz(iz,idref3)=0.0
            #           endif
            rv1 = rv2
            rv13 = rv23

        # Repeat for other distributions

    return


def _pair_loop(iv1, npixall, npixtot, npixtot2, xv1, yv1, zv1, ihv1, ikv1, ilv1, pixhwidth, voidpix, voidnum,
               idref2, idref3, ivoidsum):
    """The ``do while(iv2.lt.npixall)`` loop of find_voids (the inner loop over the
    second void pixel), split out only to keep the local variables of the hot loop
    local.  Returns (iscount, ivoidsum)."""
    from ..library import minimagec

    mi = minimagec.minimagec
    gofz = clusters.gofz
    nz = clusters.nz
    fnpixtot = f32(npixtot)
    fnz = f32(nz)
    one = r4(1.0)
    two = r4(2.0)
    vn1 = voidnum[iv1]
    iscount = 0
    iv2 = 0
    while iv2 < npixall:

        iv2 = iv2 + 1
        if voidpix[iv2] == 0 and iv1 != iv2:

            iscount = iscount + 1

            # Get the coordinates of this pixel - convert these into fractional coordinates

            #                    ivt=int(ranmine()*npixall)+1
            ivt = iv2
            ilv2 = (ivt - 1) // npixtot2
            ikv2 = (ivt - 1 - ilv2 * npixtot2) // npixtot
            ihv2 = ivt - 1 - npixtot * ikv2 - npixtot2 * ilv2
            #     write(6,*) ivt,ilv2,ikv2,ihv2
            xv2 = f32(ihv2) / fnpixtot
            yv2 = f32(ikv2) / fnpixtot
            zv2 = f32(ilv2) / fnpixtot

            # Calculate separation of pixels, applying the minimum image convention

            xdif = xv2 - xv1
            ydif = yv2 - yv1
            zdif = zv2 - zv1

            # Add the broadening component

            xdif = xdif + pixhwidth * (two * ranmine() - one)
            ydif = ydif + pixhwidth * (two * ranmine() - one)
            zdif = zdif + pixhwidth * (two * ranmine() - one)
            #     xdif=ranmine()-0.5
            #     ydif=ranmine()-0.5
            #     zdif=ranmine()-0.5
            xdif, ydif, zdif = mi(xdif, ydif, zdif)
            dist = sqrtf(xdif * xdif + ydif * ydif + zdif * zdif)
            #                       write(6,*) ihv1,ikv1,ilv1,ihv2,ikv2,ilv2,npixtot
            if xdif > 0.5 or ydif > 0.5 or zdif > 0.5:
                units.write(6, lw(ihv1, ikv1, ilv1, ihv2, ikv2, ilv2, npixtot))
            else:
                iref = ifix(dist * fnz) + 1
                if iref < 1:
                    iref = 1
                if iref > nz:
                    iref = nz
                gofz[iref, idref2] = gofz[iref, idref2] + one
                if vn1 != voidnum[iv2]:
                    gofz[iref, idref3] = gofz[iref, idref3] + one
                ivoidsum = ivoidsum + 1

    return iscount, ivoidsum


@routine
def ranmine():
    """voids.f:ranmine -- Get a random number (ran1).  Returns ranmine."""
    from ..library import auxiliary_routines

    # include 'ran1com.inc'

    ranmine = auxiliary_routines.ran1()

    return ranmine


@routine
def occupypix(npix, npixsphere, npixfound, ih, ik, il, voidpix, mpix3):
    """voids.f:occupypix -- Sets the occupation of pixels within radius npixsphere from
    the point (ih,ik,il).

    ih,ik and il are specified relative to the centre of the box, whereas pixels are
    counted from a corner. If the calculated pixel position is outside the box then the
    image position is used.  Returns npixfound.
    """
    # nvoidpix keeps track of the number of unoccupied pixels

    # integer*1 voidpix(mpix3)
    # integer*4 npix,npixsphere,npixsphere2,npixfound
    # integer*4 ih,ik,il,ils,iks,ihs

    # Square no. of pixels in sphere to save square roots later

    npixsphere2 = npixsphere * npixsphere
    ils = 0
    iks = 0
    ihs = 0

    # Count the number of pixels found in the sphere

    npixfound = 0

    while ils <= npixsphere:
        ils2 = ils * ils
        iks2pils2 = ils2
        iks = 0
        while iks <= npixsphere:
            iks2 = iks * iks
            iks2pils2 = ils2 + iks2
            ihs = 0
            while ihs <= npixsphere:

                ihs2 = ihs * ihs
                ihs2piks2pils2 = iks2pils2 + ihs2

                # Test for inside sphere

                if ihs2piks2pils2 <= npixsphere2:

                    # Found an occupied pixel

                    ihabs = ih + ihs
                    ikabs = ik + iks
                    ilabs = il + ils
                    npixfound = markpixel(npix, ihabs, ikabs, ilabs, npixfound, voidpix, mpix3)

                    # Reverse the sign of ihs, if appropriate

                    if ihs > 0:
                        ihabs = ih - ihs
                        npixfound = markpixel(npix, ihabs, ikabs, ilabs, npixfound, voidpix, mpix3)

                    # Repeat for the reversed sign of iks, if appropriate

                    if iks > 0:

                        ihabs = ih + ihs
                        ikabs = ik - iks
                        ilabs = il + ils
                        npixfound = markpixel(npix, ihabs, ikabs, ilabs, npixfound, voidpix, mpix3)

                        if ihs > 0:
                            ihabs = ih - ihs
                            npixfound = markpixel(npix, ihabs, ikabs, ilabs, npixfound, voidpix, mpix3)

                    # Repeat for reversed sign of ils, if appropriate.

                    if ils > 0:

                        ihabs = ih + ihs
                        ikabs = ik + iks
                        ilabs = il - ils
                        npixfound = markpixel(npix, ihabs, ikabs, ilabs, npixfound, voidpix, mpix3)

                        if ihs > 0:
                            ihabs = ih - ihs
                            npixfound = markpixel(npix, ihabs, ikabs, ilabs, npixfound, voidpix, mpix3)

                        if iks > 0:

                            ihabs = ih + ihs
                            ikabs = ik - iks
                            npixfound = markpixel(npix, ihabs, ikabs, ilabs, npixfound, voidpix, mpix3)

                            if ihs > 0:
                                ihabs = ih - ihs
                                npixfound = markpixel(npix, ihabs, ikabs, ilabs, npixfound, voidpix, mpix3)
                ihs = ihs + 1
            iks = iks + 1
        ils = ils + 1
    return npixfound


_F100_markpixel = "('markpixel> Error generating pixel reference numbers:',/6(1x,i4),1x,i10)"


@routine
def markpixel(npix, ihabs, ikabs, ilabs, npixfound, voidpix, mpix3):
    """voids.f:markpixel -- Routine to mark a pixel as occupied.  Returns npixfound."""
    # integer*1 voidpix(mpix3)
    # integer*4 npix,offset,npixtot,ihabs,ikabs,ilabs,npixref
    #     1,ihref,ikref,ilref,npixfound

    offset = npix
    npixtot = 2 * npix + 1
    fnpixtot = f32(npixtot)

    # Apply periodic boundary conditions, and shift coordinates of pixel so that
    # origin is at one corner of the box instead of at the centre.

    ihref = ihabs - npixtot * nint(f32(ihabs) / fnpixtot) + offset
    ikref = ikabs - npixtot * nint(f32(ikabs) / fnpixtot) + offset
    ilref = ilabs - npixtot * nint(f32(ilabs) / fnpixtot) + offset

    # Generate sequence number of pixel

    npixref = 1 + ihref + npixtot * (ikref + npixtot * ilref)
    if npixref <= mpix3:
        #     write(6,*) ihabs,ikabs,ilabs,ihref,ikref,ilref,npixref
        voidpix[npixref] = 1

        # Increment pixel counter

        npixfound = npixfound + 1

    else:

        units.write(6, fmt(_F100_markpixel, ihabs, ikabs, ilabs, ihref, ikref, ilref, npixref))
        raise FortranStop(0)

    return npixfound


_F100_search_voids = ("(//'search_voids> Unable to allocate any more void numbers!'"
                      "//'search_voids> Try using a larger void distance?')")


@routine
def search_voids(nvoids, npix, npixall, voidpix, voidnum, mpix3):
    """voids.f:search_voids -- give every unoccupied pixel the lowest non-zero void
    number of its unoccupied neighbours (or a new one), then count the voids.

    Set maximum number of voids to maximum word size for voidnum.  Returns nvoids.
    """
    # parameter (maxvoid=2**15-1)
    # integer*1 voidpix(mpix3)
    # integer*2 voidnum(mpix3)
    voidlist = [0] * 28  # integer*4 voidlist(27)
    # integer*4 nvoid,npix,npixall
    # integer*4 npixtot,nvoidlist,voidtally(maxvoid)
    voidtally = [0] * (maxvoid + 1)

    npixtot = 2 * npix + 1  # noqa: F841  (not used)
    nvoid = nvoids
    nvoidlist = 0

    # Step through the unoccupied pixels

    for ipix in range(1, npixall + 1):

        if voidpix[ipix] == 0:

            # Get the list of unoccupied pixel numbers which are adjacent to ipix

            nvoidlist = get_void_list(ipix, npix, nvoidlist, voidlist, voidpix, mpix3)

            # Search for the lowest non-zero void number within this list

            ilowest = maxvoid
            for it in range(1, nvoidlist + 1):

                itest = voidnum[voidlist[it]]
                if itest < ilowest and itest > 0:
                    ilowest = itest

            # If ilowest is maxvoid, it means none of the unoccupied pixels have been
            # given a void number so create one.

            if ilowest == maxvoid and nvoid < maxvoid:

                nvoid = nvoid + 1
                ilowest = nvoid

            elif nvoid >= maxvoid:

                units.write(6, fmt(_F100_search_voids))
                raise FortranStop(0)

            # Set the void number of all the unoccupied pixels to this (lowest)
            # void number.

            for it in range(1, nvoidlist + 1):

                voidnum[voidlist[it]] = ilowest

    # Now count the total number of voids

    nvoids = 0

    # voidtally will store the list of void numbers

    # Step through the unoccupied pixels

    for ipix in range(1, npixall + 1):

        if voidpix[ipix] == 0:

            ifind = 0
            itest = 0
            vn = voidnum[ipix]
            while itest < nvoids and ifind == 0:

                itest = itest + 1
                if vn == voidtally[itest]:
                    ifind = vn

            if ifind == 0:

                nvoids = nvoids + 1
                voidtally[nvoids] = vn

    return nvoids


@routine
def get_void_list(ipix, npix, nvoidlist, voidlist, voidpix, mpix3):
    """voids.f:get_void_list -- Makes a list of the unoccupied pixels adjacent to ipix
    (including ipix).  Returns nvoidlist."""
    # integer*1 voidpix(mpix3)
    # list of local void pixel numbers
    # integer*4 voidlist(27)
    # integer*4 ipix,npix,nvoidlist
    # integer*4 offset,npixtot,npixtot2
    # integer*4 ihc,ikc,ilc,ih,ik,il,pixref

    offset = npix
    npixtot = 2 * npix + 1
    npixtot2 = npixtot * npixtot
    nvoidlist = 0
    wrap = _wrap_table(npixtot)

    # Split ipix into its component indices relative to corner of box

    ilc = (ipix - 1) // npixtot2
    ikc = (ipix - 1 - ilc * npixtot2) // npixtot
    ihc = ipix - 1 - npixtot * ikc - npixtot2 * ilc

    # Shift these to centre of box

    ilc = ilc - offset
    ikc = ikc - offset
    ihc = ihc - offset
    #     write(6,*) ipix,ihc,ikc,ilc

    # Step through the 27 pixels around ipix (including ipix) and determine the
    # smallest, non-zero void number

    ilt = ilc - 2
    while ilt < ilc + 1:

        ilt = ilt + 1

        # Apply periodic boundary condition to il

        il = wrap[ilt]  # ilt-npixtot*nint(real(ilt)/real(npixtot))
        iladd = (il + offset) * npixtot2

        ikt = ikc - 2
        while ikt < ikc + 1:

            ikt = ikt + 1
            ik = wrap[ikt]  # ikt-npixtot*nint(real(ikt)/real(npixtot))
            ikadd = iladd + (ik + offset) * npixtot

            iht = ihc - 2
            while iht < ihc + 1:

                iht = iht + 1
                ih = wrap[iht]  # iht-npixtot*nint(real(iht)/real(npixtot))

                pixref = 1 + ih + offset + ikadd
                #     write(6,*) ih,ik,il,pixref
                # Test this pixel for occupation and void number

                if voidpix[pixref] == 0:

                    nvoidlist = nvoidlist + 1
                    voidlist[nvoidlist] = pixref

    return nvoidlist


_wrap_cache = {}


def _wrap_table(npixtot):
    """{it: it - npixtot*nint(real(it)/real(npixtot))} for the indices get_void_list
    can meet (|it| <= npixtot): the same float32 expression, evaluated once per index
    instead of 27 times per pixel."""
    t = _wrap_cache.get(npixtot)
    if t is None:
        fn = f32(npixtot)
        t = {it: it - npixtot * nint(f32(it) / fn) for it in range(-npixtot - 1, npixtot + 2)}
        _wrap_cache.clear()
        _wrap_cache[npixtot] = t
    return t


@routine
def spheregr(rstar):
    """voids.f:spheregr -- Gets the form factor (g(r)) for a solid uniform sphere.

    (Only referenced in commented-out code.)  Returns spheregr.
    """
    cons1 = r4(0.75)  # parameter (cons1=0.75,cons2=0.0625)
    cons2 = r4(0.0625)
    # real spheregr

    spheregr = f32(0.0)
    if rstar < 2.0:

        spheregr = f32(1) - rstar * (cons1 - rstar * rstar * cons2)

    return spheregr
