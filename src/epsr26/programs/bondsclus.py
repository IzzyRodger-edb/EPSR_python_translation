"""Translation of fortran_src/bondsclus.f: ``program standalone`` (the EPSR ``bonds``
program) and ``subroutine bondsclus``.

Calculates the average distance (upper triangle) and its standard deviation (lower
triangle) between each atom pair within each run of molecules of the same type, and
writes them to ``<file>.bonds``.  See docs/translation_notes/N2.md.
"""

from ..common import clustercom, commands, directory_name, molecom, symbols
from ..fortran import Ch, f32, fmt, fstr, lw, findex, real_array, char_array, sqrtf, substr, units
from ..library import ato_read_and_save, file_and_directory_commands, parse
from ..trace import routine

_F499 = "(/'bonds> No .ato file specified! - routine will exit'/)"
_F599 = "(/'Problems with .ato file ',a/'Routine will exit now'/)"
_F104 = "(//7x,40(1x,i5,a3,2x))"
_F105 = "(i3,a3,1x,(40(1x,e10.3)))"


def main():
    """bondsclus.f:standalone -- 2009/3/30 Attempt to make this a standalone program."""
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    bondsclus()


@routine
def bondsclus():
    """bondsclus.f:bondsclus -- calculates average distance between each atom pair within
    each molecule."""
    # use reallocation_routines
    # character*256 readfile
    # character*4 ext
    # character*1 iaccept
    # character*3 chtype
    # character*256 fnamebonds
    # integer natommol
    # character(len=3), dimension(:), allocatable   :: alabelt
    # real, dimension(:,:), allocatable             :: bond
    cm = commands
    alabelt = None
    bond = None
    iaccept = ''
    # NB: natommol and nsum are set before use whenever nmol > 0 (Fortran: uninitialised)
    natommol = 0
    nsum = 0
    #
    # generate the name of the .ato file
    #
    ext = '.ato'
    ilenext = 4
    if cm.nwords < 2:
        cm.nwords = 0
        while cm.nwords <= 0:
            iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)
            cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
            if cm.nwords <= 0:
                units.write(6, fmt(_F499))
                return
        readfile = fstr(substr(cm.sometext, cm.ncf[1], cm.ncl[1]), 256)
        ilenreadfile = int(cm.ncl[1]) - int(cm.ncf[1]) + 1
    else:
        readfile = fstr(substr(cm.sometext, cm.ncf[2], cm.ncl[2]), 256)
        ilenreadfile = int(cm.ncl[2]) - int(cm.ncf[2]) + 1
    #
    # strip out any extension that was supplied - extension must be .ato, otherwise
    # nothing is done.
    #
    ilendot = findex(readfile, '.') - 1
    if ilendot > 0:
        readfile = fstr(substr(readfile, 1, ilendot), 256)
        ilenreadfile = ilendot
    #
    # create the full name of the input file by adding on the directory name if not present
    #
    symbols.fnameato, symbols.ilenfnameato = file_and_directory_commands.get_full_filename(
        readfile, ilenreadfile, symbols.fnameato, symbols.ilenfnameato)
    symbols.fnameato = fstr(substr(symbols.fnameato, 1, symbols.ilenfnameato) + '.ato', 256)
    symbols.ilenfnameato = symbols.ilenfnameato + 4
    #
    # generate name of output file  - .bonds
    #
    fnamebonds = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                      + substr(readfile, 1, ilenreadfile) + '.bonds', 256)
    symbols.ilenfnameato = directory_name.nchardir + ilenreadfile + 6
    #
    # initialise number of atom and molecule types to zero
    #
    molecom.ntype = 0
    clustercom.nmoltype = 0
    #
    # read .ato file
    #
    molecom.sizefactor = f32(1.0)
    # rmaxcell, cellmax, upol, rho, amubar, amubar2 (real) and ndist (integer) are locals
    rmaxcell, cellmax, upol, rho, amubar, amubar2, ndist, ierr = ato_read_and_save.minitc(
        symbols.fnameato, f32(0.0), f32(0.0), f32(0.0), f32(0.0), f32(0.0), f32(0.0), 0, 0)
    if rho == 0.0 or ierr != 0:
        units.write(6, fmt(_F599, Ch(symbols.fnameato, 256)))
        return
    #
    # open output file
    #
    units.open(11, fnamebonds, status='unknown')
    #
    # zero bond array
    #
    ibcnt = 0
    chtype = ''
    iaref = 0  # noqa: F841 (unused in the Fortran too)
    im = 0
    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za
    naf = clustercom.naf
    while im < clustercom.nmol:
        im = im + 1
        if molecom.alabel[naf[im]] != chtype:
            #
            # A new molecule has been found
            #
            # if ibcnt.ne.0 then write out previous bond lengths
            #
            if ibcnt > 0:
                units.write(6, lw(ibcnt, natommol, nsum))
                _write_block(natommol, nsum, alabelt, bond)
            #
            # save new values
            #
            ibcnt = ibcnt + 1
            chtype = fstr(molecom.alabel[naf[im]], 3)
            natommol = int(clustercom.nato[im])
            if alabelt is not None:
                alabelt = None
            if bond is not None:
                bond = None
            alabelt = char_array(natommol)
            # NB: allocate gives uninitialised storage; every element written out is
            # assigned below except bond(1,1) of a one-atom molecule, which is 0 here.
            bond = real_array(natommol, natommol)
            for i in range(1, natommol + 1):
                iref = i + int(naf[im]) - 1
                alabelt[i] = molecom.alabel[iref]
            #
            # start new bond values
            #
            nsum = 1
            for i in range(1, natommol - 1 + 1):
                bond[i, i] = 0.0
                iref = i + int(naf[im]) - 1
                x1 = xa[iref]
                y1 = ya[iref]
                z1 = za[iref]
                for j in range(i + 1, natommol + 1):
                    bond[j, j] = 0.0
                    jref = j + int(naf[im]) - 1
                    x2 = xa[jref]
                    y2 = ya[jref]
                    z2 = za[jref]
                    xd = x2 - x1
                    yd = y2 - y1
                    zd = z2 - z1
                    dist = xd * xd + yd * yd + zd * zd
                    bond[i, j] = sqrtf(dist)
                    bond[j, i] = dist
        else:
            nsum = nsum + 1
            for i in range(1, natommol - 1 + 1):
                iref = i + int(naf[im]) - 1
                x1 = xa[iref]
                y1 = ya[iref]
                z1 = za[iref]
                for j in range(i + 1, natommol + 1):
                    jref = j + int(naf[im]) - 1
                    x2 = xa[jref]
                    y2 = ya[jref]
                    z2 = za[jref]
                    xd = x2 - x1
                    yd = y2 - y1
                    zd = z2 - z1
                    dist = xd * xd + yd * yd + zd * zd
                    bond[i, j] = bond[i, j] + sqrtf(dist)
                    bond[j, i] = bond[j, i] + dist
    #
    # write out the results from the last molecule
    #
    # (the Fortran repeats the block-writing code here; it uses float(nsum) for the
    # first division and real(nsum) for the second, which are the same conversion)
    units.write(6, lw(ibcnt, natommol, nsum))
    _write_block(natommol, nsum, alabelt, bond)
    units.close(11)
    return


def _write_block(natommol, nsum, alabelt, bond):
    """The block of bondsclus.f that writes the averages of one molecule type to unit 11
    (written out twice in the Fortran: inside the molecule loop and after it)."""
    units.write(11, fmt(_F104, *[v for i in range(1, natommol + 1) for v in (i, Ch(alabelt[i], 3))]))
    units.write(11, lw())
    for i in range(1, natommol - 1 + 1):
        for j in range(i + 1, natommol + 1):
            bond[i, j] = bond[i, j] / f32(nsum)
            bond[j, i] = bond[j, i] / f32(nsum) - bond[i, j] * bond[i, j]
            if bond[j, i] > 0.0:
                bond[j, i] = sqrtf(bond[j, i])
            else:
                bond[j, i] = 0.0
    for i in range(1, natommol + 1):
        units.write(11, fmt(_F105, i, Ch(alabelt[i], 3), *[bond[i, j] for j in range(1, natommol + 1)]))
    units.write(11, lw())
