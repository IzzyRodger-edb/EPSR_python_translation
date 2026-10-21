"""Translation of fortran_src/torangles.f: ``program standalone`` (the ``torangles``
program), ``run_torangles`` and ``find_torangles``.

Distribution of torsion angles within molecules, defined by two base atoms and, at
each end, the midpoint of two atoms.  ``vecpro``/``scapro`` are the routines of
auxiliary_routines.f (the copies at the end of torangles.f are commented out).
The library modules are imported inside the routines.  See docs/translation_notes/O.md.
"""

from ..common import clustercom
from ..common import molecom
from ..common import symbols
from ..common import torangles
from ..fortran import Ch
from ..fortran import FortranEOF
from ..fortran import FortranIOError
from ..fortran import FortranStop
from ..fortran import atan2f
from ..fortran import f32
from ..fortran import findex_padded
from ..fortran import fmt
from ..fortran import fstr
from ..fortran import ifix
from ..fortran import lw
from ..fortran import r4
from ..fortran import read_internal
from ..fortran import real_array
from ..fortran import sqrtf
from ..fortran import substr
from ..fortran import units
from ..trace import routine


def main():
    """torangles.f:standalone -- program to perform the EPSR torangles command.

    2009/3/30 Attempt to make this a standalone program
    """
    from ..library import file_and_directory_commands

    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_torangles()

    raise FortranStop(0)


_F100 = "(/'torangles> atom not matched to ','.ato file for distribution ',i3)"
_F332_run = "('(3x,a,2h:[,a,',a2,',1h-,a,',a2,',1h])')"


@routine
def run_torangles():
    """torangles.f:run_torangles -- programme to estimate the distribution of torsion
    angles (the Fortran header comment is copied from triangles.f)."""
    from ..library import (
        ato_read_and_save,
        file_and_directory_commands,
        read_and_save_files,
        setup_symbols_generic,
    )

    p = torangles
    m = molecom
    c = clustercom
    ireff = file_and_directory_commands.ireff
    # integer*4 ifound,ierrout
    # character*20 text
    # character*50 formatout
    # character*2 fmti1,fmti2,fmti3,ltypfmt,ltyp1fmt
    # NB: ltypfmt/ltyp1fmt are not initialised in the Fortran; they are only left unset
    # when ltyp >= 1000 (not reached by the tests); blank is used.
    ltypfmt = ""
    ltyp1fmt = ""

    #
    # setup the triangles symbols
    #
    setup_symbols_generic.setup_symbols_torangles()
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
    # set up the number of phi bins (between -180 and 180)
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
    # set up the number of triangle distributions
    #
    p.ndist = symbols.nfiles
    if p.ndist > p.mdist:
        p.ndist = p.mdist
    #
    # find the symbol numbers corresponding to atom1
    #
    icount = 0
    text = "atom1"
    icount = file_and_directory_commands.search_for_symbol(text, icount)
    #
    # setup the type numbers of the atom used to define
    # which molecules to use
    #
    for id_ in range(1, symbols.nfiles + 1):
        iref = ireff(icount, id_)
        ifound = 0
        it = 0
        while it < m.ntype and ifound == 0:
            it = it + 1
            if m.type[it] == symbols.symbol_value[iref]:
                ifound = it
        if ifound == 0:
            units.write(6, fmt(_F100, id_))
            return
        p.nfirst[id_] = ifound
        p.atom1[id_] = fstr(m.type[ifound], 4)
        #
        # for each torsion angle read the two base atom numbers needed to define the
        # bond
        #
        iref = iref + 1
        p.ltyp[id_] = read_internal(symbols.symbol_value[iref]).int(p.ltyp[id_])
        iref = iref + 1
        p.ltyp1[id_] = read_internal(symbols.symbol_value[iref]).int(p.ltyp1[id_])
        #
        # search through the .ato file to find out the types of atoms ltyp and ltyp1
        #
        im = 0
        ifound = 0
        while im < c.nmol and ifound == 0:
            im = im + 1
            if p.atom1[id_] == m.type[m.itype[c.naf[im]]]:
                ifound = im
                p.atom2[id_] = fstr(m.type[m.itype[c.naf[im] + p.ltyp[id_] - 1]], 4)
                p.atom3[id_] = fstr(m.type[m.itype[c.naf[im] + p.ltyp1[id_] - 1]], 4)

        #
        # and for each base atom the two atom numbers needed to define the torsion angle.
        # (These can be the same atom if desired so that only 1 atom at each end is needed
        # to define the torsion angle.)
        #
        iref = iref + 1
        r = read_internal(symbols.symbol_value[iref])
        for j in range(1, 3):
            p.ltyp2[id_, j] = r.int(p.ltyp2[id_, j])
        iref = iref + 1
        r = read_internal(symbols.symbol_value[iref])
        for j in range(1, 3):
            p.ltyp3[id_, j] = r.int(p.ltyp3[id_, j])
    #
    # define the z values
    #
    p.zstep = f32(360) / f32(p.nz)
    p.zhalf = r4(0.5) * p.zstep
    for j in range(1, p.nz + 1):
        p.z[j] = p.zstep * f32(j) - p.zhalf - r4(180.0)

    # set gofz to zero

    read_and_save_files.initgofz(p.nz, p.ndist, p.gofz, p.mz, p.mdist)

    # normalise atomic coordinates to box size

    ato_read_and_save.normalise_coordinates()

    # set up the distribution of torsional angles

    find_torangles()

    # set the axis title

    p.titler = "Phi[deg]"

    # set the column titles

    fmti1 = "i1"
    fmti2 = "i2"
    fmti3 = "i3"
    for id_ in range(1, p.ndist + 1):
        if p.ltyp[id_] < 10:
            ltypfmt = fmti1
        elif p.ltyp[id_] < 100:
            ltypfmt = fmti2
        elif p.ltyp[id_] < 1000:
            ltypfmt = fmti3
        if p.ltyp1[id_] < 10:
            ltyp1fmt = fmti1
        elif p.ltyp1[id_] < 100:
            ltyp1fmt = fmti2
        elif p.ltyp1[id_] < 1000:
            ltyp1fmt = fmti3
        formatout = fstr(fmt(_F332_run, Ch(ltypfmt, 2), Ch(ltyp1fmt, 2)), 50)
        units.write(6, lw(Ch(formatout, 50)))
        ilent1 = findex_padded(p.atom1[id_], " ", 4) - 1
        ilent2 = findex_padded(p.atom2[id_], " ", 4) - 1
        ilent3 = findex_padded(p.atom3[id_], " ", 4) - 1
        p.titledat[id_] = fstr(fmt(formatout, substr(p.atom1[id_], 1, ilent1), substr(p.atom2[id_], 1, ilent2),
                                   p.ltyp[id_], substr(p.atom3[id_], 1, ilent3), p.ltyp1[id_]), 29)
        units.write(6, lw(Ch(p.titledat[id_], 29)))

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
    iref = iref + 1
    symbols.symbol_value[iref] = fstr(lw(p.ndist), 256)

    # write input file

    read_and_save_files.save_files_others(14)

    units.close(14)
    units.write(6, lw("Torangles completed successfully"))

    return


_F332 = "(a7,i4,a11)"
_F945 = "(1x,'No of molecules found = ',i5)"


def _add_flat(gofz, iref, id_, v):
    """``gofz(iref,id) = gofz(iref,id) + v`` with the Fortran's column-major address
    (phi = +180 exactly gives iref = nz+1; Rust anal_common add_flat)."""
    n1 = gofz.shape[0] - 1
    off = (iref - 1) + (id_ - 1) * n1
    j = off // n1 + 1
    i = off - (j - 1) * n1 + 1
    gofz[i, j] = gofz[i, j] + v


@routine
def find_torangles():
    """torangles.f:find_torangles -- histogram of the torsion angles of every
    molecule whose first atom matches nfirst(id); the angles of each molecule are also
    written to toranglesvalues.txt."""
    from ..library import auxiliary_routines

    p = torangles
    m = molecom
    c = clustercom
    xa, ya, za = m.xa, m.ya, m.za
    naf = c.naf
    cd = real_array(3)
    co1 = real_array(3)
    co2 = real_array(3)
    cv1 = real_array(3)
    cv2 = real_array(3)
    cpar = real_array(3)
    phisave = real_array(p.mdist)
    # character*80 formatout
    scapro = auxiliary_routines.scapro
    vecpro = auxiliary_routines.vecpro
    half = r4(0.5)
    pi = r4(3.141592653)
    p.zstep = r4(360.0) / f32(p.nz)
    #
    # step through all the molecules
    #
    ifindm = 0
    ierr = 0
    try:
        units.open(29, "toranglesvalues.txt", status="unknown")
    except (FortranEOF, FortranIOError) as e:
        ierr = e  # noqa: F841  (iostat=ierr, not tested)

    # Setup the output format

    formatout = fstr(fmt(_F332, "(1x,i6,", p.ndist, "(1x,e14.7))"), 80)
    for im1 in range(1, c.nmol + 1):
        for id_ in range(1, p.ndist + 1):
            phisave[id_] = 0.0
        for id_ in range(1, p.ndist + 1):
            if m.itype[naf[im1]] == p.nfirst[id_]:
                ifindm = ifindm + 1
                ia1 = int(naf[im1]) + int(p.ltyp[id_]) - 1
                ia2 = int(naf[im1]) + int(p.ltyp1[id_]) - 1
                x1 = xa[ia1]
                y1 = ya[ia1]
                z1 = za[ia1]
                x2 = xa[ia2]
                y2 = ya[ia2]
                z2 = za[ia2]
                #
                # define the separation vector cd between the two base atoms and normalize it to
                # unit value
                #
                cd[1] = (x2 - x1)
                cd[2] = (y2 - y1)
                cd[3] = (z2 - z1)
                dmag = sqrtf(scapro(cd, cd))
                cd[1] = cd[1] / dmag
                cd[2] = cd[2] / dmag
                cd[3] = cd[3] / dmag
                #
                # now determine the orientation vectors co for base atom 1 and base atom 2, using
                # the values of ltyp2 and ltyp3
                #
                io1a = int(naf[im1]) + int(p.ltyp2[id_, 1]) - 1
                io1b = int(naf[im1]) + int(p.ltyp2[id_, 2]) - 1
                io2a = int(naf[im1]) + int(p.ltyp3[id_, 1]) - 1
                io2b = int(naf[im1]) + int(p.ltyp3[id_, 2]) - 1
                co1[1] = half * (xa[io1a] + xa[io1b]) - x1
                co1[2] = half * (ya[io1a] + ya[io1b]) - y1
                co1[3] = half * (za[io1a] + za[io1b]) - z1
                co2[1] = half * (xa[io2a] + xa[io2b]) - x2
                co2[2] = half * (ya[io2a] + ya[io2b]) - y2
                co2[3] = half * (za[io2a] + za[io2b]) - z2
                #
                # define new vectors normal to both cd and co
                #
                vecpro(co1, cd, cv1)
                vecpro(co2, cd, cv2)
                #
                # renormalize these to unit vectors
                #
                dmag = sqrtf(scapro(cv1, cv1))
                cv1[1] = cv1[1] / dmag
                cv1[2] = cv1[2] / dmag
                cv1[3] = cv1[3] / dmag
                dmag = sqrtf(scapro(cv2, cv2))
                cv2[1] = cv2[1] / dmag
                cv2[2] = cv2[2] / dmag
                cv2[3] = cv2[3] / dmag
                # calculate the cosine of the angle between cv1 and cv2
                cosphi = scapro(cv1, cv2)
                if cosphi > 1.0:
                    cosphi = r4(1.0)
                if cosphi < -1.0:
                    cosphi = r4(-1.0)
                #
                # calculate the vector product of cv1 and cv2 - this should be either parallel
                # (phi positive) or antiparallel to cd (phi negative)
                #
                vecpro(cv1, cv2, cpar)
                sinphi = scapro(cpar, cd)
                phis = atan2f(sinphi, cosphi) * f32(180) / pi
                phisave[id_] = phis
                #           write(6,991) ia1,ia2,io1a,io1b,io2a,io2b,phi,phis
                iref = ifix((phis + r4(180.0)) / p.zstep) + 1
                _add_flat(p.gofz, iref, id_, r4(1.0))
                # 991   format(1x,6(i4,1x),2(e13.6,1x))
        units.write(29, fmt(formatout, im1, *[phisave[id_] for id_ in range(1, p.ndist + 1)]))
    #
    # ensure final distributions are normalised to the sum of all terms and calculate average
    # chain length
    #
    for id_ in range(1, p.ndist + 1):
        sum0 = f32(0.0)
        for i in range(1, p.nz + 1):
            sum0 = sum0 + p.gofz[i, id_]
        sum0 = sum0 / f32(p.nz)
        if sum0 != 0.0:
            for i in range(1, p.nz + 1):
                p.gofz[i, id_] = p.gofz[i, id_] / sum0
    units.write(6, fmt(_F945, ifindm))
    return

#     subroutine vecpro(a,b,c)       (commented out in the Fortran: auxiliary_routines.f's is used)
#     function scapro(a,b)
