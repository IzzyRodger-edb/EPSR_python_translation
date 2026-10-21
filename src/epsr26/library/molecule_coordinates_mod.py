"""Translation of the routines of fortran_src/molecule_coordinates_mod.f90
(``module molecule_coordinates_mod``): ``get_coords_from_xyz``.

The module variables are in ``common/molecule_coordinates_mod.py``.
"""

from ..common import molecule_coordinates_mod as mcm
from ..fortran import (
    FortranEOF,
    FortranIOError,
    acosf,
    atan2f,
    char_array,
    iostat_of,
    findex,
    fmt,
    fstr,
    int_array,
    lw,
    read_internal,
    real_array,
    sqrtf,
    substr,
    units,
)
from ..trace import routine
from . import auxiliary_routines, parse

_F27 = "(/'get_coords_from_xyz> Attempting to read file: ',a)"


@routine
def get_coords_from_xyz(fname_run, lenfname_run):
    """molecule_coordinates_mod.f90:get_coords_from_xyz -- read the coordinates of the
    central molecule atoms from <fname_run up to its last dot>xyz into the module arrays
    (spherical polar r, theta, phi), and build the atom sequence of the molecule.
    Returns nothing.
    """
    # real :: xt,yt,zt,cosphi,sinphi
    nworddim = 201
    nstsize = 256
    # character(len=nstsize) :: sometext
    # character(len=3) :: tempchar
    ncf = int_array(nworddim)
    ncl = int_array(nworddim)
    nwords = 0
    # real, parameter :: bohr=0.52917721092   (not used)
    units.write(6, lw(substr(fname_run, 1, lenfname_run)))
    # Attempt to read coordinates of central molecule atoms
    mcm.xyzfilename = fstr(substr(fname_run, 1, auxiliary_routines.indexlastdot(fname_run)) + 'xyz', 256)
    mcm.nread = 10
    mcm.nofatoms = 0
    mcm.lenatomsequence = 0
    units.write(6, fmt(_F27, substr(mcm.xyzfilename, 1, len(mcm.xyzfilename))))
    try:
        units.open(mcm.nread, mcm.xyzfilename, status='old')
        ierr = 0
    except FortranIOError as e:
        ierr = e.iostat
    if ierr == 0:
        try:
            mcm.nofatoms = units.list_read(mcm.nread).int(mcm.nofatoms)
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = iostat_of(e)
        units.skip_record(mcm.nread)  # Blank line.
        if mcm.nofatoms > 0:
            if mcm.chsymbolt is not None:
                mcm.chsymbolt = None
                mcm.centralmolecule = None
                mcm.rt = None
                mcm.phit = None
                mcm.thetat = None
            n = mcm.nofatoms
            mcm.chsymbolt = char_array(n)
            mcm.centralmolecule = int_array(n)
            mcm.rt = real_array(n)
            mcm.phit = real_array(n)
            mcm.thetat = real_array(n)
            for it in range(1, n + 1):
                sometext = fstr(units.read_line(mcm.nread), nstsize)
                nwords = parse.parse(sometext, nwords, ncf, ncl, 10, nstsize, nworddim)
                units.write(6, lw(nwords, substr(sometext, 1, len(sometext))))
                mcm.chsymbolt[it] = fstr(substr(sometext, ncf[1], ncl[1]), 3)  # '(a3)'
                r = read_internal(substr(sometext, ncf[2], ncl[4]))
                xt = r.real()
                yt = r.real()
                zt = r.real()
                ierr = 1
                tempchar = ''
                if nwords > 4:
                    try:
                        tempchar = fstr(read_internal(substr(sometext, ncf[5], ncl[5])).word(tempchar), 3)
                        ierr = 0
                    except (FortranEOF, FortranIOError) as e:
                                    ierr = iostat_of(e)
                if ierr != 0:
                    tempchar = mcm.chsymbolt[it]
                mcm.centralmolecule[it] = min(1, findex(sometext, '#CENTRALMOLECULE'))
                mcm.rt[it] = sqrtf(xt * xt + yt * yt + zt * zt)
                mcm.thetat[it] = 0.0
                mcm.phit[it] = 0.0
                if mcm.rt[it] > 0.0:
                    mcm.thetat[it] = acosf(zt / mcm.rt[it])
                    cosphi = xt / mcm.rt[it]
                    sinphi = yt / mcm.rt[it]
                    mcm.phit[it] = atan2f(sinphi, cosphi)
                # Build the atom sequence to identify this molecule
                lentempchar = len(tempchar)
                if it == 1:
                    mcm.atomsequence = fstr(substr(tempchar, 1, lentempchar), 256)
                    mcm.lenatomsequence = lentempchar
                elif mcm.lenatomsequence < 256 - lentempchar:
                    mcm.atomsequence = fstr(substr(mcm.atomsequence, 1, mcm.lenatomsequence)
                                            + substr(tempchar, 1, lentempchar), 256)
                    mcm.lenatomsequence = mcm.lenatomsequence + lentempchar
    units.close(mcm.nread)
    units.write(6, lw(mcm.nofatoms, ' atoms read from file ', substr(mcm.xyzfilename, 1, len(mcm.xyzfilename))))
    units.write(6, lw('Atom sequence ', substr(mcm.atomsequence, 1, mcm.lenatomsequence)))
