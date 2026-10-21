"""Translation of fortran_src/triangles.f: ``program standalone`` (the ``triangles``
program), ``run_triangles`` and ``find_triangles``.

Distribution of the included angles of triplets of atoms at specified distance
ranges, optionally restricted to middle atoms near an atom-c / atom-z / atom-r type.

``sind``/``cosd`` in run_triangles are not declared external, so gfortran resolves them
to its intrinsics (``sind_intrinsic``/``cosd_intrinsic``; Rust anal_common.rs).
The library modules are imported inside the routines.  See docs/translation_notes/O.md.
"""

from ..common import clustercom
from ..common import commands
from ..common import molecom
from ..common import symbols
from ..common import triangles
from ..fortran import FortranStop
from ..fortran import acosf
from ..fortran import cosd_intrinsic
from ..fortran import f32
from ..fortran import findex_padded
from ..fortran import fmt
from ..fortran import fstr
from ..fortran import ifix
from ..fortran import int_array
from ..fortran import lw
from ..fortran import r4
from ..fortran import read_internal
from ..fortran import sind_intrinsic
from ..fortran import sinf
from ..fortran import sqrtf
from ..fortran import substr
from ..fortran import units
from ..trace import routine


def main():
    """triangles.f:standalone -- program to perform the EPSR triangles command.

    2009/3/30 Attempt to make this a standalone program
    """
    from ..library import file_and_directory_commands

    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_triangles()

    raise FortranStop(0)


_F100 = "(/'triangles> atoms not matched to ','.ato file for distribution ',i3)"


@routine
def run_triangles():
    """triangles.f:run_triangles -- programme to estimate the distribution of included
    angles between specified triplets of atoms at specified distance ranges."""
    from ..library import (
        ato_read_and_save,
        file_and_directory_commands,
        parse,
        read_and_save_files,
        setup_symbols_generic,
    )

    p = triangles
    ireff = file_and_directory_commands.ireff
    dmax = clustercom.dmax
    # integer*4 ifound1,ifound2,ierrout
    # character*20 text
    #
    # setup the triangles symbols
    #
    setup_symbols_generic.setup_symbols_triangles()
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
    dmax = clustercom.dmax
    #
    # set up the number of cos(theta) bins
    #
    iref = 2
    p.nz = read_internal(symbols.symbol_value[iref]).int(p.nz)
    p.nz = min(p.nz, p.mz)
    #
    # set up current accumulator
    #
    iref = iref + 1
    p.nsumt = read_internal(symbols.symbol_value[iref]).int(p.nsumt)
    #
    # Set up the centre atom, atom-c, atom-z or atom-r (only 1 can be specified).
    # This atom will be used to preselect the middle atom of the triangle to be within a
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
    # set up the number of triangle distributions
    #
    p.ndist = symbols.nfiles
    if p.ndist > p.mdist:
        p.ndist = p.mdist
    #
    # find the symbol numbers corresponding to atom1, atom2 and atom3
    #
    icount1 = 0
    text = "atom1"
    icount1 = file_and_directory_commands.search_for_symbol(text, icount1)
    icount2 = icount1 + 1
    icount3 = icount2 + 1
    #
    # setup the type numbers of the three atoms used to define
    # whether two molecules are bonded, and the maximum and minimum distances
    #
    for id_ in range(1, symbols.nfiles + 1):
        iref1 = ireff(icount1, id_)
        iref2 = ireff(icount2, id_)
        iref3 = ireff(icount3, id_)
        ifound1 = 0
        ifound2 = 0
        ifound3 = 0
        it = 0
        while it < molecom.ntype:
            it = it + 1
            if molecom.type[it] == symbols.symbol_value[iref1]:
                ifound1 = it
            if molecom.type[it] == symbols.symbol_value[iref2]:
                ifound2 = it
            if molecom.type[it] == symbols.symbol_value[iref3]:
                ifound3 = it
        if ifound1 == 0 or ifound2 == 0 or ifound3 == 0:
            units.write(6, fmt(_F100, id_))
            return
        p.ltyp1[id_] = ifound1
        p.ltyp2[id_] = ifound2
        p.ltyp3[id_] = ifound3
        p.atom1[id_] = fstr(molecom.type[ifound1], 4)
        p.atom2[id_] = fstr(molecom.type[ifound2], 4)
        p.atom3[id_] = fstr(molecom.type[ifound3], 4)
        iref = iref3 + 1
        p.ltyp[id_] = read_internal(symbols.symbol_value[iref]).int(p.ltyp[id_])
        iref = iref + 1
        p.damin[id_, 1] = read_internal(symbols.symbol_value[iref]).real(p.damin[id_, 1])
        iref = iref + 1
        p.damax[id_, 1] = read_internal(symbols.symbol_value[iref]).real(p.damax[id_, 1])
        iref = iref + 1
        p.damin[id_, 2] = read_internal(symbols.symbol_value[iref]).real(p.damin[id_, 2])
        iref = iref + 1
        p.damax[id_, 2] = read_internal(symbols.symbol_value[iref]).real(p.damax[id_, 2])
        iref = iref + 1
        r = read_internal(symbols.symbol_value[iref])
        p.rminc[id_] = r.real(p.rminc[id_])
        p.rmaxc[id_] = r.real(p.rmaxc[id_])
        #         write(6,*) id,ltyp(id),ltyp1(id),ltyp2(id),ltyp3(id),ltypc
        #         write(6,*) (damin(id,i),damax(id,i),i=1,2),rminc(id),rmaxc(id)
        for i in range(1, 3):
            p.damin[id_, i] = p.damin[id_, i] / dmax
            p.damax[id_, i] = p.damax[id_, i] / dmax
        p.rminc[id_] = p.rminc[id_] / dmax
        p.rmaxc[id_] = p.rmaxc[id_] / dmax
    #
    # define the z values
    #
    p.zstep = r4(180.0) / f32(p.nz)
    p.zhalf = r4(0.5) * p.zstep
    for j in range(1, p.nz + 1):
        p.z[j] = p.zstep * f32(j) - p.zhalf

    # set gofz to zero

    read_and_save_files.initgofz(p.nz, p.ndist, p.gofz, p.mz, p.mdist)

    # normalise atomic coordinates to box size

    ato_read_and_save.normalise_coordinates()

    # set up the distribution of tri-angles

    find_triangles()

    # set the axis title

    p.titler = "theta[deg.]"

    # set the column titles

    for id_ in range(1, p.ndist + 1):
        ilent1 = findex_padded(p.atom1[id_], " ", 4) - 1
        ilent2 = findex_padded(p.atom2[id_], " ", 4) - 1
        ilent3 = findex_padded(p.atom3[id_], " ", 4) - 1
        p.titledat[id_] = fstr("   " + substr(p.atom1[id_], 1, ilent1) + "-" + substr(p.atom2[id_], 1, ilent2)
                               + "-" + substr(p.atom3[id_], 1, ilent3), 29)

    # set the file extension

    p.ext = "c"

    # write the accumulated results

    p.nsumt = read_and_save_files.outgofz(symbols.fnameinp, p.ext, p.nz, p.ndist, p.nsumt, p.z, p.gofz, p.gzt,
                                          p.gzt2, p.titler, p.titledat, p.mz, p.mdist)

    # update the current value of nsumt

    icount = 0
    text = "nsize"
    icount = file_and_directory_commands.search_for_symbol(text, icount)
    iref = ireff(icount, 1)
    symbols.symbol_value[iref] = fstr(lw(p.nz), 256)
    iref = iref + 1
    symbols.symbol_value[iref] = fstr(lw(p.nsumt), 256)
    iref = iref + 4
    symbols.symbol_value[iref] = fstr(lw(p.ndist), 256)

    # write input file

    read_and_save_files.save_files_others(14)

    # Calculate q value for each distribution

    for id_ in range(1, p.ndist + 1):

        sum1 = f32(0.0)
        sum2 = f32(0.0)
        for iz in range(1, p.nz + 1):
            term = p.gzt[iz, id_] * sind_intrinsic(p.z[iz])
            q = cosd_intrinsic(p.z[iz])
            if p.ltyp[id_] == 1:
                q = q + r4(1.0) / r4(3.0)
                q = q * q
            sum1 = sum1 + term
            sum2 = sum2 + q * term
        if sum1 > 0.0:
            if p.ltyp[id_] == 1:
                units.write(14, lw(id_, r4(0.5) * (p.rminc[id_] + p.rmaxc[id_]) * dmax,
                                   (r4(1.0) - (r4(9.0) * sum2) / (r4(4.0) * sum1))))
            else:
                units.write(14, lw(id_, r4(0.5) * (p.rminc[id_] + p.rmaxc[id_]) * dmax, (sum2 / sum1)))
        else:
            units.write(14, lw(id_, r4(0.5) * (p.rminc[id_] + p.rmaxc[id_]) * dmax, r4(0.0)))

    units.close(14)
    units.write(6, lw("Triangles completed successfully"))

    return


def _add_flat(gofz, iref, id_, v):
    """``gofz(iref,id) = gofz(iref,id) + v`` with the Fortran's column-major address:
    a bin index beyond the first dimension (an angle of exactly 180 degrees with
    nz = mz) addresses the next column, as the unchecked Fortran does (Rust add_flat)."""
    n1 = gofz.shape[0] - 1
    off = (iref - 1) + (id_ - 1) * n1
    j = off // n1 + 1
    i = off - (j - 1) * n1 + 1
    gofz[i, j] = gofz[i, j] + v


@routine
def find_triangles():
    """triangles.f:find_triangles -- histogram of the included angles of the triplets
    atom1-atom2-atom3, normalised to unit area (sin(theta) weighted)."""
    from ..library import minimagec

    p = triangles
    m = molecom
    c = clustercom
    # logical carryon
    atom1list = None  # integer, dimension(:), allocatable
    atom2list = None
    atom3list = None
    p.zstep = r4(180.0) / f32(p.nz)
    piconv = r4(180.0) / r4(3.141592653)
    natom = m.natom
    itype = m.itype
    mref = m.mref
    nbond = m.nbond
    iabond = m.iabond
    xa, ya, za = m.xa, m.ya, m.za
    xadd, yadd, zadd = c.xadd, c.yadd, c.zadd
    gofz = p.gofz
    zstep = p.zstep
    ltypc = p.ltypc
    dorhoz, dorhor = p.dorhoz, p.dorhor
    mi = minimagec.minimagec
    two = r4(2.0)
    one = r4(1.0)
    #
    # step through all the atoms and ignore those which do not match the required
    # ltyp1
    #
    for id_ in range(1, p.ndist + 1):
        lt1 = p.ltyp1[id_]
        lt2 = p.ltyp2[id_]
        lt3 = p.ltyp3[id_]
        # Make a list of all the atoms of type ltyp2(id) which are within the specified distance
        # of atoms of type ltypc
        natom1list = 0
        natom2list = 0
        natom3list = 0
        for ia2 in range(1, natom + 1):
            if itype[ia2] == lt1:
                natom1list = natom1list + 1
            if itype[ia2] == lt2:
                natom2list = natom2list + 1
            if itype[ia2] == lt3:
                natom3list = natom3list + 1
        #         write(6,*) id,natom1list,natom2list,natom3list
        if atom1list is not None:
            atom1list = None
        if atom2list is not None:
            atom2list = None
        if atom3list is not None:
            atom3list = None
        if natom1list > 0:
            atom1list = int_array(natom1list)
        if natom2list > 0:
            atom2list = int_array(natom2list)
        if natom3list > 0:
            atom3list = int_array(natom3list)
        # Repeat this saving the atom numbers
        iref1 = 0
        iref2 = 0
        iref3 = 0
        for ia2 in range(1, natom + 1):
            if itype[ia2] == lt1:
                iref1 = iref1 + 1
                atom1list[iref1] = ia2
            if itype[ia2] == lt3:
                iref3 = iref3 + 1
                atom3list[iref3] = ia2
            carryon = False
            if itype[ia2] == lt2:
                # Need to select those middle atoms which are within specified distance range of atomc, if defined
                carryon = True
                if ltypc > 0:
                    carryon = False
                    im2 = mref[ia2]
                    xc2 = xadd[im2]
                    yc2 = yadd[im2]
                    zc2 = zadd[im2]
                    x2 = xa[ia2] + xc2
                    y2 = ya[ia2] + yc2
                    z2 = za[ia2] + zc2
                    ia4 = 0
                    while ia4 < natom and not carryon:
                        ia4 = ia4 + 1
                        if itype[ia4] == ltypc:
                            im4 = mref[ia4]
                            xdif = x2 - (xa[ia4] + xadd[im4])
                            ydif = y2 - (ya[ia4] + yadd[im4])
                            zdif = z2 - (za[ia4] + zadd[im4])
                            # apply periodic boundary conditions
                            xdif, ydif, zdif = mi(xdif, ydif, zdif)
                            if dorhoz:
                                dsq12 = zdif * zdif
                            elif dorhor:
                                dsq12 = xdif * xdif + ydif * ydif
                            else:
                                dsq12 = xdif * xdif + ydif * ydif + zdif * zdif
                            d12 = sqrtf(dsq12)
                            carryon = d12 >= p.rminc[id_] and d12 < p.rmaxc[id_]
                            #                    write(6,*) ia2,ia4,itype(ia4),carryon
                if carryon:
                    iref2 = iref2 + 1
                    atom2list[iref2] = ia2
        # If the type of the first atom is the same as that of the centre atom (ltypc) then we use the minimum and maximum distance for this atom as specified for this pair
        if lt1 == ltypc:
            p.damin[id_, 1] = p.rminc[id_]
            p.damax[id_, 1] = p.rmaxc[id_]
        natom1list = iref1
        natom2list = iref2
        natom3list = iref3
        units.write(6, lw(id_, natom1list, natom2list, natom3list))
        damin1 = p.damin[id_, 1]
        damax1 = p.damax[id_, 1]
        damin2 = p.damin[id_, 2]
        damax2 = p.damax[id_, 2]
        centre1 = lt1 == ltypc
        ltyp_id = p.ltyp[id_]
        # Now do main calculation using these lists of ia1, ia2 and ia3 atoms
        for iref1 in range(1, natom1list + 1):
            ia1 = int(atom1list[iref1])
            xc1 = xadd[mref[ia1]]
            yc1 = yadd[mref[ia1]]
            zc1 = zadd[mref[ia1]]
            x1 = xa[ia1] + xc1
            y1 = ya[ia1] + yc1
            z1 = za[ia1] + zc1
            for iref2 in range(1, natom2list + 1):
                ia2 = int(atom2list[iref2])
                im2 = mref[ia2]
                xc2 = xadd[im2]
                yc2 = yadd[im2]
                zc2 = zadd[im2]
                x2 = xa[ia2] + xc2
                y2 = ya[ia2] + yc2
                z2 = za[ia2] + zc2
                #
                # calculate atom1 to atom2 distance, applying periodic boundary conditions
                #
                xdif = x2 - x1
                ydif = y2 - y1
                zdif = z2 - z1
                # apply periodic boundary conditions
                xdif, ydif, zdif = mi(xdif, ydif, zdif)
                if dorhoz and centre1:
                    dsq12 = zdif * zdif
                elif dorhor and centre1:
                    dsq12 = xdif * xdif + ydif * ydif
                else:
                    dsq12 = xdif * xdif + ydif * ydif + zdif * zdif
                d12 = sqrtf(dsq12)
                # check that distance is within bounds prescribed for this distribution
                if d12 >= damin1 and d12 <= damax1:
                    # search through all the other atoms for the third atom to define
                    # the bond angle
                    if ltyp_id == 1:
                        for iref3 in range(1, natom3list + 1):
                            ia3 = int(atom3list[iref3])
                            if ia3 != ia2 and ia3 != ia1:
                                im3 = mref[ia3]
                                xc3 = xadd[im3]
                                yc3 = yadd[im3]
                                zc3 = zadd[im3]
                                x3 = xa[ia3] + xc3
                                y3 = ya[ia3] + yc3
                                z3 = za[ia3] + zc3
                                xdif = x3 - x2
                                ydif = y3 - y2
                                zdif = z3 - z2
                                xdif, ydif, zdif = mi(xdif, ydif, zdif)
                                dsq23 = xdif * xdif + ydif * ydif + zdif * zdif
                                d23 = sqrtf(dsq23)
                                if d23 >= damin2 and d23 <= damax2:
                                    # If centre atom is defined, ignore respective coordinates for dorhoz or dorhor
                                    if dorhoz and centre1:
                                        xdif = x3 - x2
                                        ydif = y3 - y2
                                        zdif = z3 - z1
                                    elif dorhor and centre1:
                                        xdif = x3 - x1
                                        ydif = y3 - y1
                                        zdif = z3 - z2
                                    else:
                                        xdif = x3 - x1
                                        ydif = y3 - y1
                                        zdif = z3 - z1
                                    xdif, ydif, zdif = mi(xdif, ydif, zdif)
                                    dsq13 = xdif * xdif + ydif * ydif + zdif * zdif
                                    coss = (dsq12 + dsq23 - dsq13) / (two * d12 * d23)
                                    if coss > 1.0:
                                        coss = one
                                    if coss < -1.0:
                                        coss = -one
                                    deg = piconv * acosf(coss)
                                    iref = ifix(deg / zstep) + 1
                                    _add_flat(gofz, iref, id_, one)
                    else:
                        x3 = f32(0.0)
                        y3 = f32(0.0)
                        z3 = f32(0.0)
                        nf = 0
                        for ib in range(1, int(nbond[ia2]) + 1):
                            ia3 = int(iabond[ia2, ib])
                            if itype[ia3] == lt3 and ia3 != ia2:
                                x3 = x3 + xa[ia3]
                                y3 = y3 + ya[ia3]
                                z3 = z3 + za[ia3]
                                nf = nf + 1
                        if nf > 0:
                            x3 = x3 / f32(nf) + xc2
                            y3 = y3 / f32(nf) + yc2
                            z3 = z3 / f32(nf) + zc2
                            xdif = x3 - x2
                            ydif = y3 - y2
                            zdif = z3 - z2
                            xdif, ydif, zdif = mi(xdif, ydif, zdif)
                            dsq23 = xdif * xdif + ydif * ydif + zdif * zdif
                            d23 = sqrtf(dsq23)
                            # If centre atom is defined, ignore respective coordinates for dorhoz or dorhor
                            if dorhoz and centre1:
                                xdif = x3 - x2
                                ydif = y3 - y2
                                zdif = z3 - z1
                            elif dorhor and centre1:
                                xdif = x3 - x1
                                ydif = y3 - y1
                                zdif = z3 - z2
                            else:
                                xdif = x3 - x1
                                ydif = y3 - y1
                                zdif = z3 - z1
                            xdif, ydif, zdif = mi(xdif, ydif, zdif)
                            dsq13 = xdif * xdif + ydif * ydif + zdif * zdif
                            coss = (dsq12 + dsq23 - dsq13) / (two * d12 * d23)
                            if coss > 1.0:
                                coss = one
                            if coss < -1.0:
                                coss = -one
                            deg = piconv * acosf(coss)
                            iref = ifix(deg / zstep) + 1
                            _add_flat(gofz, iref, id_, one)
    if atom1list is not None:
        atom1list = None
    if atom2list is not None:
        atom2list = None
    if atom3list is not None:
        atom3list = None
    #
    # normalise each distribution function to have an area of unity
    #
    for id_ in range(1, p.ndist + 1):
        sumg = f32(0.0)
        beta = r4(-0.5) * zstep
        for iz in range(1, p.nz + 1):
            beta = beta + zstep
            betar = beta / piconv
            gofz[iz, id_] = gofz[iz, id_] / sinf(betar)
            sumg = sumg + gofz[iz, id_]
        if sumg > 0.0:
            for iz in range(1, p.nz + 1):
                gofz[iz, id_] = gofz[iz, id_] / (sumg * zstep)
    return
