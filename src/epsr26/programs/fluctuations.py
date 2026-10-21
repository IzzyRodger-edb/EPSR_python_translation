"""Translation of fortran_src/fluctuations.f: ``program standalone`` (the
``fluctuations`` program), ``run_fluctuations`` and ``find_fluctuations``.

Number-density fluctuations: the number of molecule centres within spheres of given
radii placed at random in the box, its histogram, and running averages of N and N^2
(S(0) from the number fluctuations).  The random-number state is read from and saved to
randomno.dat.  The library modules are imported inside the routines.  See
docs/translation_notes/O.md.
"""

from ..common import clustercom
from ..common import clusters
from ..common import lattice
from ..common import ran1com
from ..common import symbols
from ..fortran import PI4ATAN
from ..fortran import adjustl
from ..fortran import dble_array
from ..fortran import f32
from ..fortran import findex_padded
from ..fortran import fmt
from ..fortran import fstr
from ..fortran import int_array
from ..fortran import lw
from ..fortran import nint
from ..fortran import r4
from ..fortran import read_internal
from ..fortran import real_array
from ..fortran import sqrtf
from ..fortran import substr
from ..fortran import units
from ..trace import routine


def main():
    """fluctuations.f:standalone -- program to perform the EPSR fluctuations command.

    2009/3/30 Attempt to make this a standalone program
    """
    from ..library import file_and_directory_commands

    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_fluctuations()


_F100 = "(e12.5)"


@routine
def run_fluctuations():
    """fluctuations.f:run_fluctuations -- programme to estimate the density
    fluctuations within a specified radius.

    (The rest of the Fortran header comment is copied from voids.f.)
    """
    from ..library import (
        ato_read_and_save,
        auxiliary_routines,
        file_and_directory_commands,
        read_and_save_files,
        setup_symbols_generic,
    )

    p = clusters
    ireff = file_and_directory_commands.ireff
    # integer*4 ierrout
    # character*20 text
    sumnumber = dble_array(p.mdist)
    sumnumbersqr = dble_array(p.mdist)
    localnumber = dble_array(p.mdist)
    localnumbersqr = dble_array(p.mdist)
    szero = dble_array(p.mdist)
    # real*8 frac,frac1
    # real damaxmin
    #
    # setup the clusters symbols
    #
    setup_symbols_generic.setup_symbols_fluctuations()

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
    #
    # set up the number of density bins
    #
    pi = PI4ATAN  # noqa: F841  (not used)
    iref = 2
    p.nz = read_internal(symbols.symbol_value[iref]).int(p.nz)
    p.nz = min(p.nz, p.mz)
    #
    # Get the number of samples already probed
    #
    iref = iref + 1
    nsumts = 0
    nsumts = read_internal(symbols.symbol_value[iref]).int(nsumts)
    if nsumts < 0:
        p.nsumt = 0
    else:
        p.nsumt = nsumts

    # Set the number of samples to take from the current configuration

    iref = iref + 1
    p.npix = read_internal(symbols.symbol_value[iref]).int(p.npix)

    # Get the number of different volumes to sample

    p.ndist = symbols.nfiles
    if p.ndist > p.mdist:
        p.ndist = p.mdist

    # Find the largest value of damax

    damaxmin = f32(0.0)
    for id_ in range(1, p.ndist + 1):
        #
        # Set up the radius of the volume to test
        #
        iref = ireff(symbols.nsymboluptodata + 1, id_)
        p.damax[id_] = read_internal(symbols.symbol_value[iref]).real(p.damax[id_])
        if id_ == 1:
            damaxmin = p.damax[id_]
        else:
            damaxmin = min(damaxmin, p.damax[id_])

        # Average of N so far

        iref = iref + 1
        sumnumber[id_] = read_internal(symbols.symbol_value[iref]).dble(sumnumber[id_])

        # Average of N^2 so far

        iref = iref + 1
        sumnumbersqr[id_] = read_internal(symbols.symbol_value[iref]).dble(sumnumbersqr[id_])

        # S(0) from number fluctuations

        iref = iref + 1
        szero[id_] = read_internal(symbols.symbol_value[iref]).dble(szero[id_])

        if p.nsumt == 0:

            sumnumber[id_] = 0.0
            sumnumbersqr[id_] = 0.0

    #
    # setup step size and z values
    #
    rho = f32(clustercom.nmol) / (lattice.alatx * lattice.blaty * lattice.clatz)
    zstep1 = r4(2.0) * rho / f32(p.nz)  # noqa: F841  (not used)
    #      zstep2=1.0/(damaxmin*damaxmin*damaxmin)
    #      zstep=max(zstep1,zstep2)

    #      zstep=1.0/sqrt(4.0*pi*rho*damaxmin*damaxmin*damaxmin/3.0)
    p.zstep = f32(1)
    units.write(6, lw(p.nz, p.ndist, p.zstep))

    for j in range(1, p.nz + 1):
        p.z[j] = f32(j - 1 - p.nz // 2) * p.zstep

    read_and_save_files.initgofz(p.nz, p.ndist, p.gofz, p.mz, p.mdist)
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

    # set up the table of numbers of molecules in each sample volume

    find_fluctuations(rho, localnumber, localnumbersqr)

    #
    # axis title
    #
    p.titler = "Density"
    for id_ in range(1, p.ndist + 1):
        #
        # set the column title for this distribution
        #
        text = fstr(fmt(_F100, p.damax[id_]), 20)
        text = adjustl(text, 20)
        ilent = findex_padded(text, " ", 20) - 1
        p.titledat[id_] = fstr("r_max:" + substr(text, 1, ilent), 29)

    # output the data

    p.ext = "n"

    p.nsumt = read_and_save_files.outgofz(symbols.fnameinp, p.ext, p.nz, p.ndist, p.nsumt, p.z, p.gofz, p.gzt,
                                          p.gzt2, p.titler, p.titledat, p.mz, p.mdist)

    # Update the sums

    frac = float(r4(1.0) / f32(p.nsumt))  # 1.0/real(nsumt) is single precision
    frac1 = 1.0 - frac

    for id_ in range(1, p.ndist + 1):

        sumnumber[id_] = sumnumber[id_] * frac1 + localnumber[id_] * frac
        sumnumbersqr[id_] = sumnumbersqr[id_] * frac1 + localnumbersqr[id_] * frac
        if sumnumber[id_] > 0.0:
            szero[id_] = sumnumbersqr[id_] / sumnumber[id_] - sumnumber[id_]
        else:
            szero[id_] = 0.0

    # update the current value of nsumt and ndist

    icount = 0
    text = "nsize"
    icount = file_and_directory_commands.search_for_symbol(text, icount)
    iref = ireff(icount, 1)
    symbols.symbol_value[iref] = fstr(lw(p.nz), 256)
    iref = iref + 1
    if nsumts >= 0:
        nsumts = p.nsumt
    symbols.symbol_value[iref] = fstr(lw(nsumts), 256)
    iref = iref + 1
    symbols.symbol_value[iref] = fstr(lw(p.npix), 256)
    iref = iref + 1
    symbols.symbol_value[iref] = fstr(lw(p.ndist), 256)

    for id_ in range(1, p.ndist + 1):
        iref = iref + 1
        symbols.symbol_value[iref] = fstr(lw(p.damax[id_]), 256)
        iref = iref + 1
        symbols.symbol_value[iref] = fstr(lw(float(sumnumber[id_])), 256)
        iref = iref + 1
        symbols.symbol_value[iref] = fstr(lw(float(sumnumbersqr[id_])), 256)
        iref = iref + 1
        symbols.symbol_value[iref] = fstr(lw(float(szero[id_])), 256)

    # write input file

    read_and_save_files.save_files_others(14)

    # Save random number parameters.

    units.open(14, "randomno.dat", status="unknown")
    auxiliary_routines.saveran1(14)
    units.close(14)

    return


@routine
def find_fluctuations(rho, sumnumber, sumnumbersqr):
    """fluctuations.f:find_fluctuations -- count the molecule centres within each
    test radius around npix random points; histogram the deviation from the mean."""
    from ..library import auxiliary_routines, minimagec

    p = clusters
    c = clustercom
    L = lattice
    # real*4 rho,voidtest(mdist),xcent,ycent,zcent,pi
    # real*8 sumnumber(mdist),sumnumbersqr(mdist)
    voidtest = real_array(p.mdist)
    rhotest = int_array(p.mdist)
    # logical test
    ran1 = auxiliary_routines.ran1
    mi = minimagec.minimagec

    pi = PI4ATAN

    # Renormalise maximum radius for each distribution

    for id_ in range(1, p.ndist + 1):

        voidtest[id_] = p.damax[id_] / c.dmax

        # Setup average number of atoms in this volume

        rhotest[id_] = nint(r4(4.0) * pi * rho * p.damax[id_] * p.damax[id_] * p.damax[id_] / r4(3.0))

        # zero local sums

        sumnumber[id_] = 0.0
        sumnumbersqr[id_] = 0.0

    ndist = p.ndist
    nz = p.nz
    nmol = c.nmol
    xadd, yadd, zadd = c.xadd, c.yadd, c.zadd
    ltyp = p.ltyp
    gofz = p.gofz
    vt = [voidtest[id_] for id_ in range(0, ndist + 1)]  # local copy of voidtest (same float32 values)
    i = 0
    while i < p.npix:

        # Choose a position at random within the box

        fa = (ran1() - r4(0.5))
        fb = (ran1() - r4(0.5))
        fc = (ran1() - r4(0.5))

        # Convert to Cartesian coordinates

        #
        # calculate corrected coordinates
        #
        zcent = fc * L.clatz
        ycent = fb * L.blaty + fc * L.claty
        xcent = fa * L.alatx + fb * L.blatx + fc * L.clatx

        test = False
        if not test:

            i = i + 1

            # Zero accumulators

            for id_ in range(1, ndist + 1):
                ltyp[id_] = 0

            # Look at all the other molecules and pick those that are within the box of dimenstion damax(id)/dmax
            # look at each requested atom radius

            for im in range(1, nmol + 1):

                # Calculate the separation of these molecules

                xdif = xadd[im] - xcent
                ydif = yadd[im] - ycent
                zdif = zadd[im] - zcent

                # Apply minimum image convention

                xdif, ydif, zdif = mi(xdif, ydif, zdif)
                diftest = sqrtf(xdif * xdif + ydif * ydif + zdif * zdif)

                for id_ in range(1, ndist + 1):
                    if diftest <= vt[id_]:

                        ltyp[id_] = ltyp[id_] + 1

            # Accumulate number sums

            for id_ in range(1, ndist + 1):

                lt = int(ltyp[id_])
                sumnumber[id_] = sumnumber[id_] + float(lt)
                sumnumbersqr[id_] = sumnumbersqr[id_] + float(lt * lt)

                # density of this box

                size = p.damax[id_]  # noqa: F841  (not used)
                iref = lt - int(rhotest[id_]) + 1 + nz // 2
                #                  iref=int((real(ltyp(id))-rhotest(id))/(zstep*sqrt(rhotest(id))))+1+nz/2
                if iref < 1:
                    iref = 1
                if iref > nz:
                    iref = nz
                gofz[iref, id_] = gofz[iref, id_] + r4(1.0)

    # Normalise sums to number of terms

    for id_ in range(1, ndist + 1):

        sumnumber[id_] = sumnumber[id_] / float(p.npix)
        sumnumbersqr[id_] = sumnumbersqr[id_] / float(p.npix)

    # Normalise the density distributions to probabilities

    for id_ in range(1, ndist + 1):

        sum_ = f32(0.0)
        for iz in range(1, nz + 1):
            sum_ = sum_ + gofz[iz, id_]
        sum_ = sum_ * p.zstep
        if sum_ > 0.0:
            for iz in range(1, nz + 1):
                gofz[iz, id_] = gofz[iz, id_] / sum_
        else:
            for iz in range(1, nz + 1):
                gofz[iz, id_] = 0.0
    return
