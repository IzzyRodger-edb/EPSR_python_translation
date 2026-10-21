"""Translation of fortran_src/partials.f: ``program standalone`` (the ``partials``
program), ``run_partials`` and ``find_partials``.

Site-site (intermolecular) radial distribution functions between pairs of atom types
of the .ato configuration, optionally restricted to a centre atom type (atom-c), or
computed as density profiles along z (atom-z) or in the xy plane (atom-r).

The library modules are imported inside the routines (setup_symbols_generic is being
translated at the same time as this file).  See docs/translation_notes/O.md.
"""

from ..common import clustercom
from ..common import commands
from ..common import lattice
from ..common import molecom
from ..common import partials
from ..common import symbols
from ..common.epsr_dimension import mpairs
from ..fortran import PI4ATAN
from ..fortran import Ch
from ..fortran import f32
from ..fortran import findex_padded
from ..fortran import fmt
from ..fortran import fstr
from ..fortran import ifix
from ..fortran import int_array
from ..fortran import lw
from ..fortran import r4
from ..fortran import read_internal
from ..fortran import sqrtf
from ..fortran import substr
from ..fortran import units
from ..trace import routine


def main():
    """partials.f:standalone -- program to perform the EPSR partials command.

    2009/3/30 Attempt to make this a standalone program
    """
    from ..library import file_and_directory_commands

    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_partials()


_F104 = "(1x,i5,5x,a3,1x,a3)"


@routine
def run_partials():
    """partials.f:run_partials -- site-site radial distribution functions.

    (The Fortran header comment, copied from clusters.f, says: programme to estimate
    number of clusters in a distribution of molecules.)
    """
    from ..library import (
        ato_read_and_save,
        file_and_directory_commands,
        parse,
        read_and_save_files,
        setup_symbols_generic,
    )

    p = partials
    # integer*4 ierrout
    # character*20 text
    #
    # setup the paritals symbols
    #
    setup_symbols_generic.setup_symbols_partials()
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
    # set up the number of cluster size bins
    #
    iref = 2
    p.nz = read_internal(symbols.symbol_value[iref]).int(p.nz)
    p.nz = min(p.nz, p.mz)
    # Set up range of calculation
    iref = iref + 1
    p.rnormalise = f32(1.0)
    p.rnormalise = read_internal(symbols.symbol_value[iref]).real(p.rnormalise)
    if p.rnormalise > 1.0:
        p.rnormalise = f32(1.0)
    #
    # set up current accumulator
    #
    iref = iref + 1
    p.nsumt = read_internal(symbols.symbol_value[iref]).int(p.nsumt)
    #
    # Set up the centre atom, atom-c, atom-z or atom-r (only 1 can be specified).
    # This atom will be used to preselect the centre atom about which densities will be calculated
    # specified distance range of this atom.
    iref = iref + 1
    p.ltypc = 0
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    p.atomc = fstr(symbols.symbol_value[iref], 3)  # read(symbol_value(iref),'(a)',iostat=ierr) atomc
    # Check this symbol exists in the current file
    it = 0
    while it < molecom.ntype and p.ltypc == 0:
        it = it + 1
        if p.atomc == molecom.type[it]:
            p.ltypc = it
    #
    # set up whether rho-z is to be calculated
    #
    iref = iref + 1
    p.dorhoz = False
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    # For atomz to work we need to specify an atom type to act as the origin.
    p.atomz = fstr(symbols.symbol_value[iref], 3)  # read(symbol_value(iref),'(a)',iostat=ierr) atomz
    # Check this symbol exists in the current file
    it = 0
    while it < molecom.ntype and p.ltypc == 0:
        it = it + 1
        if p.atomz == molecom.type[it]:
            p.ltypc = it
            p.dorhoz = True
    #
    # set up whether rho-r is to be calculated. Note that if both atomz and atomr are set,
    # only the first will be used.
    #
    iref = iref + 1
    p.dorhor = False
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    # For atomr to work we need to specify an atom type to act as the origin.
    p.atomr = fstr(symbols.symbol_value[iref], 3)  # read(symbol_value(iref),'(a)',iostat=ierr) atomr
    # Check this symbol exists in the current file
    it = 0
    while it < molecom.ntype and p.ltypc == 0:
        it = it + 1
        if p.atomr == molecom.type[it]:
            p.ltypc = it
            p.dorhor = True
    units.write(6, fmt("(l2,l2,1x,i4)", p.dorhoz, p.dorhor, p.ltypc))
    #
    # check the number of distributions does not exceed the array size
    #
    p.ndist = symbols.nfiles
    if p.ndist > p.mdist:
        p.ndist = p.mdist
    units.write(6, lw(p.nz, p.ndist))
    #
    # set up number of partial g(r)'s to calculate
    #
    p.ndist = 0
    for it in range(1, molecom.ntype + 1):
        for jt in range(1, molecom.ntype + 1):
            if p.ndist < mpairs:
                # Restricted to specified atom at origin
                if p.ltypc > 0:
                    if it == p.ltypc:
                        p.ndist = p.ndist + 1
                        p.ltyp1[p.ndist] = it
                        p.ltyp2[p.ndist] = jt
                        p.atom1[p.ndist] = fstr(molecom.type[it], 4)
                        p.atom2[p.ndist] = fstr(molecom.type[jt], 4)
                elif jt >= it:
                    # Use all atom pairs
                    p.ndist = p.ndist + 1
                    p.ltyp1[p.ndist] = it
                    p.ltyp2[p.ndist] = jt
                    p.atom1[p.ndist] = fstr(molecom.type[it], 4)
                    p.atom2[p.ndist] = fstr(molecom.type[jt], 4)
    for id_ in range(1, p.ndist + 1):
        units.write(6, fmt("(1x,i4,2(1x,a3))", id_, Ch(p.atom1[id_], 4), Ch(p.atom2[id_], 4)))
    #
    # setup step size and z values
    #
    # For dorhoz and dorhor set up factor to capture true distance range in each case
    if p.dorhoz:
        p.rnormalise = p.rnormalise * lattice.clatz / clustercom.dmax
    if p.dorhor:
        p.rnormalise = p.rnormalise * min(lattice.alatx, lattice.blaty) / clustercom.dmax
    p.zstep = r4(0.5) * p.rnormalise / f32(p.nz)
    p.zhalf = r4(0.5) * p.zstep
    for j in range(1, p.nz + 1):
        p.z[j] = (p.zstep * f32(j) - p.zhalf) * clustercom.dmax

    read_and_save_files.initgofz(p.nz, p.ndist, p.gofz, p.mz, p.mdist)

    #
    # normalise atomic coordinates to box size
    #
    ato_read_and_save.normalise_coordinates()
    #
    # set up the table of molecular bonds
    #
    find_partials()
    #
    # axis title
    #
    p.titler = "r"
    for id_ in range(1, p.ndist + 1):
        #
        # set the column title for this pair of atoms
        #
        ilent1 = findex_padded(p.atom1[id_], " ", 4) - 1
        ilent2 = findex_padded(p.atom2[id_], " ", 4) - 1
        p.titledat[id_] = fstr("   " + substr(p.atom1[id_], 1, ilent1) + "-" + substr(p.atom2[id_], 1, ilent2), 29)

    # output the data

    p.ext = "g"

    p.nsumt = read_and_save_files.outgofz(symbols.fnameinp, p.ext, p.nz, p.ndist, p.nsumt, p.z, p.gofz, p.gzt,
                                          p.gzt2, p.titler, p.titledat, p.mz, p.mdist)

    #
    # write input file
    #

    # update the current value of nsumt and ndist

    icount = 0
    text = "nz"
    icount = file_and_directory_commands.search_for_symbol(text, icount)
    iref = file_and_directory_commands.ireff(icount, 1)
    symbols.symbol_value[iref] = fstr(lw(p.nz), 256)
    iref = iref + 2
    symbols.symbol_value[iref] = fstr(lw(p.nsumt), 256)
    iref = iref + 4
    symbols.symbol_value[iref] = fstr(lw(p.ndist), 256)

    read_and_save_files.save_files_others(14)

    for id_ in range(1, p.ndist + 1):
        #
        # for each partial correlation function write the two atom type numbers needed.
        #
        units.write(14, fmt(_F104, id_, Ch(p.atom1[id_], 4), Ch(p.atom2[id_], 4)))

    units.close(14)

    return


_F193 = "(1x,i5,5(1x,e13.6))"


@routine
def find_partials():
    """partials.f:find_partials -- accumulate the site-site histograms and normalise
    them to g(r) (or rho(z), rho(r))."""
    from ..library import minimagec

    p = partials
    m = molecom
    c = clustercom
    lsum1 = int_array(mpairs)
    lsum2 = int_array(mpairs)
    # integer, dimension(:), allocatable     :: xlist,ylist,zlist,xlistsave,ylistsave,zlistsave
    pi = PI4ATAN
    dmax = c.dmax
    # volume of sphere of radius 0.5
    if p.dorhoz:
        gfac = f32(1.0)
        vol = r4(2.0) * lattice.alatx * lattice.blaty * dmax * dmax  # The factor of 2 is needed because we symmetrise about z = 0
    elif p.dorhor:
        gfac = pi
        vol = lattice.clatz * dmax
    else:
        gfac = r4(4.0) * pi / r4(3.0)
        vol = lattice.alatx * lattice.blaty * lattice.clatz * dmax * dmax * dmax
    #
    # step through all the atoms and ignore those which do not match the required
    # ltyp1
    #
    # set up counters for number of origin and non-origin atoms for each partial
    #
    natom = m.natom
    itype = m.itype
    mref = m.mref
    xa, ya, za = m.xa, m.ya, m.za
    xadd, yadd, zadd = c.xadd, c.yadd, c.zadd
    gofz = p.gofz
    nz = p.nz
    zstep = p.zstep
    dorhoz, dorhor = p.dorhoz, p.dorhor
    mi = minimagec.minimagec
    for id_ in range(1, p.ndist + 1):
        lsum1[id_] = 0
        lsum2[id_] = 0
        lt1 = p.ltyp1[id_]
        lt2 = p.ltyp2[id_]
        for ia1 in range(1, natom + 1):
            if itype[ia1] == lt2:
                lsum2[id_] = lsum2[id_] + 1
            if itype[ia1] == lt1:
                im1 = mref[ia1]
                xc1 = xadd[im1]
                yc1 = yadd[im1]
                zc1 = zadd[im1]
                x1 = xa[ia1] + xc1
                y1 = ya[ia1] + yc1
                z1 = za[ia1] + zc1
                doloop = True
                #            if(lsum1(id).gt.0) then
                #               doloop=.not.((rhoz.eq.1.and.(abs(z1-z1s)).lt.0.05/dmax)
                #     1.or.(rhoxy.eq.1.and.(abs(x1-x1s)+abs(y1-y1s)).lt.0.05/dmax))
                #            end if
                if doloop:
                    lsum1[id_] = lsum1[id_] + 1
                    for ia2 in range(1, natom + 1):
                        im2 = mref[ia2]
                        if itype[ia2] == lt2 and im1 != im2:
                            #
                            # calculate atom1 to atom2 distance, applying periodic boundary conditions
                            #
                            xc2 = xadd[im2]
                            yc2 = yadd[im2]
                            zc2 = zadd[im2]
                            x2 = xa[ia2] + xc2
                            y2 = ya[ia2] + yc2
                            z2 = za[ia2] + zc2
                            xdif = x2 - x1
                            ydif = y2 - y1
                            zdif = z2 - z1
                            # apply periodic boundary conditions
                            xdif, ydif, zdif = mi(xdif, ydif, zdif)
                            if dorhoz:
                                d12 = abs(zdif)
                            elif dorhor:
                                dsq12 = xdif * xdif + ydif * ydif
                                d12 = sqrtf(dsq12)
                            else:
                                dsq12 = xdif * xdif + ydif * ydif + zdif * zdif
                                d12 = sqrtf(dsq12)
                            iref = ifix(d12 / zstep) + 1
                            # check that distance is within bounds prescribed for this distribution
                            if iref > 0 and iref <= nz:
                                gofz[iref, id_] = gofz[iref, id_] + r4(1.0)
                # Save the coordinates of this atom. For rhoz=1 or rhoxy=1, if the next atom of this type
                # has barely changed values of z (rhoz=1) or (x,y) (rhoxy=1) then we don't need to include it
                # since the second atoms will already have been included in these cases
                # (x1s, y1s, z1s are saved but not used)
    # renormalize g(r)'s to total number of atoms
    for id_ in range(1, p.ndist + 1):
        # density of second atoms in sphere
        if p.dorhoz or p.dorhor:
            rho = f32(lsum1[id_]) * vol
        elif p.ltypc > 0:
            #               rho=1.0
            rho = f32(lsum1[id_])
        else:
            rho = f32(lsum1[id_]) * f32(lsum2[id_]) / vol
        units.write(6, fmt(_F193, id_, f32(lsum1[id_]), f32(lsum2[id_]), rho))
        #           do iz=100,nz
        #                 write(6,193) iz,z(iz),gofz(iz,id)
        #           end do
        gfac1 = gfac * rho
        z0 = f32(0.0)
        z03 = f32(0.0)
        for iz in range(1, nz + 1):
            zm = p.z[iz]
            zdiff = zm - z0
            z1 = zm + zdiff
            if p.dorhoz:
                z13 = z1
            elif p.dorhor:
                z13 = z1 * z1
            else:
                z13 = z1 * z1 * z1
            # average number of particles in this volume element
            rat = gfac1 * (z13 - z03)
            gofz[iz, id_] = gofz[iz, id_] / rat
            z0 = z1
            z03 = z13
    return
