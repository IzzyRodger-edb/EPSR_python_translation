"""Translation of fortran_src/coord.f: ``program standalone`` (the ``coord`` program),
``run_coord`` and ``find_coord``.

Coordination number histograms (and the mean bond length for each coordination number)
for specified pairs of atom types, optionally restricted to centre atoms within a
distance range of an atom-c / atom-z / atom-r type.

The library modules are imported inside the routines.  See docs/translation_notes/O.md.
"""

from ..common import clustercom
from ..common import commands
from ..common import coord
from ..common import molecom
from ..common import symbols
from ..common.epsr_dimension import matom
from ..fortran import Ch
from ..fortran import f32
from ..fortran import findex_padded
from ..fortran import fmt
from ..fortran import fstr
from ..fortran import int_array
from ..fortran import lw
from ..fortran import r4
from ..fortran import read_internal
from ..fortran import real_array
from ..fortran import sqrtf
from ..fortran import substr
from ..fortran import units
from ..trace import routine


def main():
    """coord.f:standalone -- program to perform the EPSR coord command.

    2009/3/30 Attempt to make this a standalone program
    """
    from ..library import file_and_directory_commands

    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_coord()


_F100 = "(1x,a,'-',a,'(',e9.3,',',e9.3,')')"
_F101 = "(1x,a3,1x,a3,4(1x,f10.5))"


@routine
def run_coord():
    """coord.f:run_coord -- programme to estimate coordination number histogram for
    specified atom type pairs."""
    from ..library import (
        ato_read_and_save,
        file_and_directory_commands,
        parse,
        read_and_save_files,
        setup_symbols_generic,
    )

    p = coord
    ireff = file_and_directory_commands.ireff
    # integer*4 ifound1,ifound2,ierrout
    # character*20 text
    #
    # setup the coord symbols
    #
    setup_symbols_generic.setup_symbols_coord()
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
    # set up the number of coordination number size bins
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
    # check the number of distributions does not exceed the array size
    #
    p.ndist = symbols.nfiles
    if p.ndist > p.mdist // 2:
        p.ndist = p.mdist // 2
    p.ndist2 = 2 * p.ndist
    #
    # setup the type number of the first and second atoms used to define
    # the distance range over which the coordination number is to be calculated
    #
    for id_ in range(1, p.ndist + 1):
        it = 0
        ifound1 = 0
        ifound2 = 0
        iref1 = ireff(symbols.nsymboluptodata + 1, id_)
        iref2 = iref1 + 1
        iref3 = iref2 + 1
        iref4 = iref3 + 1
        while it < molecom.ntype and (ifound1 == 0 or ifound2 == 0):
            it = it + 1
            if molecom.type[it] == symbols.symbol_value[iref1]:
                ifound1 = it
            if molecom.type[it] == symbols.symbol_value[iref2]:
                ifound2 = it
        p.ltyp1[id_] = ifound1
        p.ltyp2[id_] = ifound2
        p.atom1[id_] = fstr(molecom.type[ifound1], 4)
        p.atom2[id_] = fstr(molecom.type[ifound2], 4)
        p.damin[id_] = read_internal(symbols.symbol_value[iref3]).real(p.damin[id_])
        p.damax[id_] = read_internal(symbols.symbol_value[iref4]).real(p.damax[id_])
        units.write(6, lw(id_, p.ltyp1[id_], p.ltyp2[id_], p.damin[id_], p.damax[id_]))
        iref = iref4 + 1
        r = read_internal(symbols.symbol_value[iref])
        p.rminc[id_] = r.real(p.rminc[id_])
        p.rmaxc[id_] = r.real(p.rmaxc[id_])
    units.write(6, lw(p.nz, p.ndist))
    p.zstep = r4(0.5) / f32(p.nz)
    p.zhalf = r4(0.5) * p.zstep
    #
    # two sets of distributions are calculated: one of coordination number distribution
    # the other of mean bond length for each coordination number
    #
    for j in range(1, p.nz + 1):
        p.z[j] = f32(j - 1)

    read_and_save_files.initgofz(p.nz, p.ndist2, p.gofz, p.mz, p.mdist)

    #
    # normalise atomic coordinates to box size
    #
    ato_read_and_save.normalise_coordinates()
    #
    # set up the table of molecular bonds
    #
    find_coord()
    #
    # axis title
    #
    p.titler = "Coordination_number"
    idref2 = 0
    for id_ in range(1, symbols.nfiles + 1):
        idref1 = idref2 + 1
        idref2 = idref1 + 1
        #
        # set the column title for this pair of atoms
        #
        ilent1 = findex_padded(p.atom1[id_], " ", 4) - 1
        ilent2 = findex_padded(p.atom2[id_], " ", 4) - 1
        units.write(6, lw(ilent1, ilent2, p.ndist, p.ndist2))
        units.write(6, fmt(_F100, substr(p.atom1[id_], 1, ilent1), substr(p.atom2[id_], 1, ilent2),
                           p.damin[id_], p.damax[id_]))
        p.titledat[idref1] = fstr(fmt(_F100, substr(p.atom1[id_], 1, ilent1), substr(p.atom2[id_], 1, ilent2),
                                      p.damin[id_], p.damax[id_]), 29)
        #
        # title of corresponding column for bond length distribution
        #
        p.titledat[idref2] = p.titledat[idref1]
        units.write(6, lw(Ch(p.titledat[id_], 29)))

    # output the data

    p.ext = "n"

    p.nsumt = read_and_save_files.outgofz(symbols.fnameinp, p.ext, p.nz, p.ndist2, p.nsumt, p.z, p.gofz, p.gzt,
                                          p.gzt2, p.titler, p.titledat, p.mz, p.mdist)

    #
    # calculate average coordination number for each distribution
    #
    idref2 = 0
    for id_ in range(1, p.ndist + 1):
        idref1 = idref2 + 1
        idref2 = idref1 + 1
        sum1 = f32(0.0)
        sum2 = f32(0.0)
        for i in range(1, p.nz + 1):
            term = p.z[i] * p.gzt[i, idref1]
            sum1 = sum1 + term
            sum2 = sum2 + p.z[i] * term
        p.coordn[id_, 1] = sum1
        std = sum2 - sum1 * sum1
        if std > 0.0:
            p.coordn[id_, 2] = sqrtf(std)
        else:
            p.coordn[id_, 2] = 0.0

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
    symbols.symbol_value[iref] = fstr(lw(p.nz), 256)
    iref = iref + 1
    symbols.symbol_value[iref] = fstr(lw(p.nsumt), 256)
    iref = iref + 4
    symbols.symbol_value[iref] = fstr(lw(p.ndist), 256)

    read_and_save_files.save_files_others(14)

    for id_ in range(1, p.ndist + 1):
        #
        # for each partial correlation function write the two atom types,
        # the average coordination number, and it rms deviation
        #
        if p.ltypc == 0:
            units.write(14, fmt(_F101, Ch(p.atom1[id_], 4), Ch(p.atom2[id_], 4), p.damin[id_], p.damax[id_],
                                p.coordn[id_, 1], p.coordn[id_, 2]))
        else:
            units.write(14, fmt(_F101, Ch(p.atom1[id_], 4), Ch(p.atom2[id_], 4), p.damin[id_], p.damax[id_],
                                p.rminc[id_], p.rmaxc[id_], p.coordn[id_, 1], p.coordn[id_, 2]))

    units.close(14)

    return


_F994 = "(1x,'Average coord. number for type ',i2,' around type ',i2,' = ',f10.5)"


@routine
def find_coord():
    """coord.f:find_coord -- coordination number histogram and mean bond length per
    coordination number for each distribution."""
    from ..library import minimagec

    p = coord
    m = molecom
    c = clustercom
    bav = real_array(matom)
    mbond = int_array(matom)
    # logical carryon
    atom1list = None  # integer, dimension(:), allocatable
    atom2list = None

    p.zstep = f32(1)
    pi = r4(3.141592653)
    # volume of sphere of radius 0.5
    gfac = r4(4.0) * pi / r4(3.0)
    vol = gfac * r4(0.125) * c.dmax * c.dmax * c.dmax  # noqa: F841  (not used)
    natom = m.natom
    itype = m.itype
    mref = m.mref
    xa, ya, za = m.xa, m.ya, m.za
    xadd, yadd, zadd = c.xadd, c.yadd, c.zadd
    dmax = c.dmax
    gofz = p.gofz
    nz = p.nz
    mi = minimagec.minimagec
    for id_ in range(1, p.ndist + 1):
        #
        # define the reference integers for gofz
        #
        idref1 = 2 * id_ - 1
        idref2 = 2 * id_
        lt1 = p.ltyp1[id_]
        lt2 = p.ltyp2[id_]
        # Make a list of all the atoms of type ltyp2(id) which are within the specified distance
        # of atoms of type ltypc
        natom1list = 0
        natom2list = 0
        for ia2 in range(1, natom + 1):
            if itype[ia2] == lt1:
                natom1list = natom1list + 1
            if itype[ia2] == lt2:
                natom2list = natom2list + 1
        #         write(6,*) id,natom1list,natom2list,natom3list
        if atom1list is not None:
            atom1list = None
        if atom2list is not None:
            atom2list = None
        if natom1list > 0:
            atom1list = int_array(natom1list)
        if natom2list > 0:
            atom2list = int_array(natom2list)
        # Repeat this saving the atom numbers
        iref1 = 0
        iref2 = 0
        for ia2 in range(1, natom + 1):
            if itype[ia2] == lt1:
                # Need to select those middle atoms which are within specified distance range of atomc, if defined
                carryon = True
                if p.ltypc > 0:
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
                        if itype[ia4] == p.ltypc:
                            im4 = mref[ia4]
                            xdif = x2 - (xa[ia4] + xadd[im4])
                            ydif = y2 - (ya[ia4] + yadd[im4])
                            zdif = z2 - (za[ia4] + zadd[im4])
                            # apply periodic boundary conditions
                            xdif, ydif, zdif = mi(xdif, ydif, zdif)
                            if p.dorhoz:
                                dsq12 = zdif * zdif
                            elif p.dorhor:
                                dsq12 = xdif * xdif + ydif * ydif
                            else:
                                dsq12 = xdif * xdif + ydif * ydif + zdif * zdif
                            d12 = sqrtf(dsq12) * dmax
                            carryon = d12 >= p.rminc[id_] and d12 < p.rmaxc[id_]
                            #                    write(6,*) ia2,ia4,itype(ia4),carryon
                if carryon:
                    iref1 = iref1 + 1
                    atom1list[iref1] = ia2
            if itype[ia2] == lt2:
                iref2 = iref2 + 1
                atom2list[iref2] = ia2
        natom1list = iref1
        natom2list = iref2
        #
        # set up for each molecule the number of bonded neighbours.
        #
        for i in range(1, natom + 1):
            mbond[i] = 0
            bav[i] = 0.0
        #
        # number of origin atoms
        #
        momol = 0
        #
        # number of bonded atoms
        #
        mbmol = 0
        damin_id = p.damin[id_]
        damax_id = p.damax[id_]
        #
        # step through all origin molecules
        #
        for iref1 in range(1, natom1list + 1):
            ia1 = int(atom1list[iref1])
            im1 = mref[ia1]
            momol = momol + 1
            xc1 = xadd[im1]
            yc1 = yadd[im1]
            zc1 = zadd[im1]
            xa1 = xa[ia1] + xc1
            ya1 = ya[ia1] + yc1
            za1 = za[ia1] + zc1
            # now look at all the other molecules and determine which ones are
            # bonded
            for iref2 in range(1, natom2list + 1):
                ia2 = atom2list[iref2]
                im2 = mref[ia2]
                #
                # calculate atom1 to atom2 distance, applying periodic boundary conditions
                #
                xc2 = xadd[im2]
                yc2 = yadd[im2]
                zc2 = zadd[im2]
                x2 = xa[ia2] + xc2
                y2 = ya[ia2] + yc2
                z2 = za[ia2] + zc2
                xdif = x2 - xa1
                ydif = y2 - ya1
                zdif = z2 - za1
                # apply periodic boundary conditions
                xdif, ydif, zdif = mi(xdif, ydif, zdif)
                dsq12 = xdif * xdif + ydif * ydif + zdif * zdif
                da12 = sqrtf(dsq12) * dmax
                # check that distance is within bounds prescribed for this distribution
                if da12 >= damin_id and da12 <= damax_id:
                    mbond[ia1] = mbond[ia1] + 1
                    bav[ia1] = bav[ia1] + da12
                    mbmol = mbmol + 1
        #
        # accumulate the histogram for this molecule
        #
        for iref1 in range(1, natom1list + 1):
            ia1 = atom1list[iref1]
            if mbond[ia1] > 0:
                bav[ia1] = bav[ia1] / f32(mbond[ia1])
            else:
                bav[ia1] = 0.0
            iref = int(mbond[ia1]) + 1
            if iref > 0 and iref < nz:
                gofz[iref, idref1] = gofz[iref, idref1] + r4(1.0)
                gofz[iref, idref2] = gofz[iref, idref2] + bav[ia1]
        units.write(6, lw(mbmol, momol))
        #
        # calculate the average number of bonds per molecule and normalise histograms
        #
        if momol > 0:
            avbond = f32(mbmol) / f32(momol)
            for iref in range(1, nz + 1):
                if gofz[iref, idref1] > 0.0:
                    gofz[iref, idref2] = gofz[iref, idref2] / gofz[iref, idref1]
                else:
                    gofz[iref, idref2] = 0.0
                gofz[iref, idref1] = gofz[iref, idref1] / f32(momol)
        else:
            avbond = f32(0.0)
        units.write(6, fmt(_F994, p.ltyp2[id_], p.ltyp1[id_], avbond))
    return
