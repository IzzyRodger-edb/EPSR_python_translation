"""Translation of fortran_src/dockato.f: ``program standalone`` (``main``), ``dockato``,
``rotatem``, ``random_rotatem``, ``search_unsaturated``, ``find_valence``,
``test_overlap``, ``find_neighbour``, ``bond_molecule``, ``replace_references``,
``add_neighbours``, ``remove_atoms``, ``clean_atom_types``.

The data of ``module dockato_arrays`` (dockato.f lines 1-20) is in
``common/dockato_arrays.py``.  ``get_rotgrp`` (declared by an INTERFACE block in
``dockato``) is the routine of ato_manipulation_rev4.f.  See docs/translation_notes/N1.md.
"""

from ..common import clustercom
from ..common import commands
from ..common import directory_name
from ..common import dockato_arrays
from ..common import molecom
from ..common import molecom_mod
from ..common import symbols
from ..common.epsr_dimension import matom, mmol, near
from ..fortran import Ch
from ..fortran import FortranEOF
from ..fortran import FortranIOError
from ..fortran import PI4ATAN
from ..fortran import cosf
from ..fortran import f32
from ..fortran import findex
from ..fortran import fmt
from ..fortran import fstr
from ..fortran import int_array
from ..fortran import iostat_of
from ..fortran import lw
from ..fortran import powf
from ..fortran import r4
from ..fortran import read_internal
from ..fortran import real_array
from ..fortran import setsub
from ..fortran import sign
from ..fortran import sinf
from ..fortran import sqrtf
from ..fortran import substr
from ..fortran import units
from ..library import ato_manipulation_rev4
from ..library import ato_read_and_save
from ..library import auxiliary_routines
from ..library import file_and_directory_commands
from ..library import minimagec
from ..library import parse
from ..library import reallocation_routines
from ..library import molecom_mod as molecom_mod_routines
from ..trace import routine

# FORMAT statements of dockato
_F699 = "(//'dockato> Search for .ato file containing docking molecule')"
_F698 = "(//'dockato> Search for .ato file to be docked onto')"
_F700 = "(/'dockato> No .ato files found in the current directory - try another directory.'/)"
_F499 = "(/'dockato> No .ato file specified - try again!'/)"
_F36 = "(/1x,i3,' new atom types in file ',a/)"
_F33 = "(1x,'Atom type ',i2,' has label ',a3)"
# 696 format(/'dockato> Specify site to be used for docking: ',$)
_F696 = "(/'dockato> Specify site to be used for docking: ')"
# 689 format(/'dockato> Specify site to be used for rotating',' this molecule: ',$)
_F689 = "(/'dockato> Specify site to be used for rotating',' this molecule: ')"
# 687 format(/'dockato> Specify site to be used for a second',' rotation of this molecule : ',$)
_F687 = "(/'dockato> Specify site to be used for a second',' rotation of this molecule : ')"
# 697 format(/'dockato> How many docking molecules are needed? ',$)
_F697 = "(/'dockato> How many docking molecules are needed? ')"
# 685 format(/'dockato> Give maximum distance from central molecule to',' dock molecules: ',$)
_F685 = "(/'dockato> Give maximum distance from central molecule to',' dock molecules: ')"
_F694 = "(/'dockato> Atoms on the docking molecule are incorrectly specified',' - re-specify types')"
# 661 format(/'dockato> Specify site to be docked onto: ',$)
_F661 = "(/'dockato> Specify site to be docked onto: ')"
# 688 format(/'dockato> Specify site to be rotated on to: ',$)
_F688 = "(/'dockato> Specify site to be rotated on to: ')"
# 686 format(/'dockato> Specify site to be rotated on to, second rotation: ',$)
_F686 = "(/'dockato> Specify site to be rotated on to, second rotation: ')"
_F684 = "(/'dockato> Specify minimum dihedral angle for second rotation ')"
_F683 = "(/'dockato> Specify maximum dihedral angle for second rotation: ')"
_F784 = "(/'dockato> Type overlaps factor (0 or yes to allow all overlaps) ')"
_F693 = "(/'dockato> Atoms on the docked molecule are incorrectly specified',' - re-specify types')"
_F599 = "(/'Problems with .ato file ',a,' Try again.'/)"
_F601 = "(/'dockato> Tried to define more molecules than allowed by dimension'/)"
_F602 = "(/'dockato> Tried to define more atoms than allowed by dimension'/)"
_F664 = "(/'dockato> ifind = ',i5)"
_F662 = "(/'dockato> Added ',i5,' molecules')"
_F663 = "(/'dockato> (no more valid sites found for docking.)')"
# 35 format(/'dockato> Type name of file to put mixture in: ',$)
_F35 = "(/'dockato> Type name of file to put mixture in: ')"
_F101 = "(a)"


def main():
    """dockato.f:standalone -- Program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    dockato()


def _read5_word(old, n):
    """``read(5,*) c`` for a character*n variable c."""
    r = units.list_read(5)
    return fstr(r.word(old), n)


@routine
def dockato():
    """dockato.f:dockato -- dock molecules of one .ato file onto sites of another.

    A specified number of the first (docking) molecule .ato file are docked onto the
    specified site on the  second (docked) molecule .ato file.

    9-Oct-2006. Modified so that when docking occurs, the added molecule becomes physically
    bonded to the docked molecule. ...
    The bonding sites are designated by Qv atoms, which are removed once docking
    has occurred. The number of neighbouring Qv atoms is used to designate the degree of
    saturation of each site.
    """
    # use reallocation_routines
    # use molecom_mod
    # use dockato_arrays
    m = molecom
    c = clustercom
    cm = commands
    da = dockato_arrays

    # character*256 readfile
    # character*4 ext
    # character*20 nooverlapstext
    # logical nooverlaps !If true docking molecule will be checked for overlaps with docked molecule.
    # names of docking atoms
    # character*3 atomv,atomd,atomvr,atomdr,atomvrr,atomdrr
    # character*1 iaccept
    # real dminrr,dmaxrr ! minimum and maximum distance to determine if a second rotation is possible
    # Variables to save filenames and atom labels
    # character*80 dockingfile,dockedfile,outputfile
    # logical fixeddihedral     !Use to indicate whether a specific dihedral has been set.
    readfile = ""
    ext = ""
    atomv = atomd = atomvr = atomdr = atomvrr = atomdrr = ""
    iaccept = ""
    dockingfile = dockedfile = outputfile = ""
    # Locals the Fortran reads before assigning on some paths (stack garbage there);
    # 0 is used here (see docs/translation_notes/N1.md)
    rmaxcell = f32(0.0)
    cellmax = f32(0.0)
    upol = f32(0.0)
    rho = f32(0.0)
    amubar = f32(0.0)
    amubar2 = f32(0.0)
    ndist = 0
    nvalencev = 0
    nvalenced = 0
    xv = yv = zv = f32(0.0)
    xd = yd = zd = f32(0.0)
    dihedralangfixed = f32(0.0)
    imv = 0  # NB: undefined in the Fortran if no docking site is ever found
    iadrr = 0
    labref = 0

    # Keep the compiler happy

    ilendockingfile = 0
    ilendockedfile = 0
    ilenatomd = 0
    ilenatomdr = 0
    ilenatomdrr = 0
    ilenatomv = 0
    ilenatomvr = 0
    ilenatomvrr = 0
    ilendminrr = 0
    ilendmaxrr = 0
    ilenoutputfile = 0
    ndock = 0
    radmax = f32(0.0)
    dminrr = f32(0.0)
    dmaxrr = f32(0.0)
    da.xadds = f32(0)
    da.yadds = f32(0)
    da.zadds = f32(0)
    da.phixs = f32(0)
    da.phiys = f32(0)
    da.phizs = f32(0)
    da.natos = 0
    da.amassts = f32(0)
    nooverlapstext = "<undefined>"

    ntest = 1
    # get name of docking file - the only first molecule in this file will be used for docking
    if cm.nwords > ntest:
        nt = ntest + 1
        dockingfile = fstr(substr(cm.sometext, cm.ncf[nt], cm.ncl[nt]), 80)
        ilendockingfile = cm.ncl[nt] - cm.ncf[nt] + 1
    # get name of atom to be used for translating the docking molecule
    ntest = ntest + 1
    if cm.nwords > ntest:
        nt = ntest + 1
        atomd = fstr(substr(cm.sometext, cm.ncf[nt], cm.ncl[nt]), 3)
        ilenatomd = cm.ncl[nt] - cm.ncf[nt] + 1
    # get name of atom to be used for rotating the docking molecule
    ntest = ntest + 1
    if cm.nwords > ntest:
        nt = ntest + 1
        atomdr = fstr(substr(cm.sometext, cm.ncf[nt], cm.ncl[nt]), 3)
        ilenatomdr = cm.ncl[nt] - cm.ncf[nt] + 1
    # get name of second atom to be used to generate a second rotation if desired
    ntest = ntest + 1
    if cm.nwords > ntest:
        nt = ntest + 1
        atomdrr = fstr(substr(cm.sometext, cm.ncf[nt], cm.ncl[nt]), 3)
        ilenatomdrr = cm.ncl[nt] - cm.ncf[nt] + 1
    # get number of molecules to be docked
    ntest = ntest + 1
    if cm.nwords > ntest:
        nt = ntest + 1
        ndock = read_internal(substr(cm.sometext, cm.ncf[nt], cm.ncl[nt])).int(ndock)
    # get maximum distance of these molecules from first molecule
    ntest = ntest + 1
    if cm.nwords > ntest:
        nt = ntest + 1
        radmax = read_internal(substr(cm.sometext, cm.ncf[nt], cm.ncl[nt])).real(radmax)
    # get name of .ato file to be docked onto.
    ntest = ntest + 1
    if cm.nwords > ntest:
        nt = ntest + 1
        dockedfile = fstr(substr(cm.sometext, cm.ncf[nt], cm.ncl[nt]), 80)
        ilendockedfile = cm.ncl[nt] - cm.ncf[nt] + 1
    # get name of atom to be translated onto
    ntest = ntest + 1
    if cm.nwords > ntest:
        nt = ntest + 1
        atomv = fstr(substr(cm.sometext, cm.ncf[nt], cm.ncl[nt]), 3)
        ilenatomv = cm.ncl[nt] - cm.ncf[nt] + 1
    # get name of atom to be rotated onto
    ntest = ntest + 1
    if cm.nwords > ntest:
        nt = ntest + 1
        atomvr = fstr(substr(cm.sometext, cm.ncf[nt], cm.ncl[nt]), 3)
        ilenatomvr = cm.ncl[nt] - cm.ncf[nt] + 1
    # get name of atom on the system being docked onto to be used to generate
    # a second rotation if desired.
    ntest = ntest + 1
    if cm.nwords > ntest:
        nt = ntest + 1
        atomvrr = fstr(substr(cm.sometext, cm.ncf[nt], cm.ncl[nt]), 3)
        ilenatomvrr = cm.ncl[nt] - cm.ncf[nt] + 1
    # get minimum dihedral angle
    ntest = ntest + 1
    if cm.nwords > ntest:
        nt = ntest + 1
        dminrr = read_internal(substr(cm.sometext, cm.ncf[nt], cm.ncl[nt])).real(dminrr)
        ilendminrr = cm.ncl[nt] - cm.ncf[nt] + 1
    # get maximum dihedral angle
    ntest = ntest + 1
    if cm.nwords > ntest:
        nt = ntest + 1
        dmaxrr = read_internal(substr(cm.sometext, cm.ncf[nt], cm.ncl[nt])).real(dmaxrr)
        ilendmaxrr = cm.ncl[nt] - cm.ncf[nt] + 1
    # Determine if docking molecule will be checked for overlaps with docked molecule
    ntest = ntest + 1
    nooverlaps = False
    if cm.nwords > ntest:
        nt = ntest + 1
        nooverlapstext = fstr(substr(cm.sometext, cm.ncf[nt], cm.ncl[nt]), 20)
    # get name of output file
    ntest = ntest + 1
    if cm.nwords > ntest:
        nt = ntest + 1
        outputfile = fstr(substr(cm.sometext, cm.ncf[nt], cm.ncl[nt]), 80)
        ilenoutputfile = cm.ncl[nt] - cm.ncf[nt] + 1

    ierr = 1
    nmix = 2
    da.natomax = 0
    da.ngrotmax = 0
    da.nrotmax = 0
    da.ndihmax = 0
    da.natos = 0
    da.ngrots = 0

    # set the number of atom and molecule types to zero at the outset

    m.ntype = 0
    c.nmoltype = 0
    ntypefirst = 0
    nmoltypefirst = c.nmoltype
    imix = 0
    while imix < nmix:
        imix = imix + 1
        if imix == 1:
            if ilendockingfile == 0:
                cm.nwords = 0
                units.write(6, fmt(_F699))
            else:
                cm.nwords = 1
                cm.ncf[1] = 1
                cm.ncl[1] = ilendockingfile
                cm.sometext = setsub(cm.sometext, 1, ilendockingfile, dockingfile, cm.nstsize)
        else:
            if ilendockedfile == 0:
                cm.nwords = 0
                units.write(6, fmt(_F698))
            else:
                cm.nwords = 1
                cm.ncf[1] = 1
                cm.ncl[1] = ilendockedfile
                cm.sometext = setsub(cm.sometext, 1, ilendockedfile, dockedfile, cm.nstsize)

        # generate the name of each .ato file in turn

        ext = fstr(".ato", 4)
        ilenext = 4
        while cm.nwords <= 0:
            iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)

            # if no .ato files are available in this directory, abandon this attempt at dockato

            if iaccept != "y":
                units.write(6, fmt(_F700))
                return
            cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
            if cm.nwords <= 0:
                units.write(6, fmt(_F499))
        readfile = fstr(substr(cm.sometext, cm.ncf[1], cm.ncl[1]), 256)
        ilenreadfile = cm.ncl[1] - cm.ncf[1] + 1

        # strip out any extension that was supplied - extension must be .ato, otherwise
        # nothing is done.

        ilendot = findex(readfile, ".") - 1
        if ilendot > 0:
            readfile = fstr(substr(readfile, 1, ilendot), 256)
            ilenreadfile = ilendot
        #
        # create the full name of the input file by adding on the directory name if not present
        #
        symbols.fnameato, symbols.ilenfnameato = file_and_directory_commands.get_full_filename(
            readfile, ilenreadfile, symbols.fnameato, symbols.ilenfnameato)
        symbols.fnameato = fstr(substr(symbols.fnameato, 1, symbols.ilenfnameato) + ".ato", 256)
        symbols.ilenfnameato = symbols.ilenfnameato + 4

        # read .ato file

        c.nmoltype = 0
        nmoltypefirst = c.nmoltype
        m.sizefactor = f32(1.0)
        (rmaxcell, cellmax, upol, rho, amubar, amubar2, ndist, ierr) = ato_read_and_save.minitc(
            symbols.fnameato, rmaxcell, cellmax, upol, rho, amubar, amubar2, ndist, ierr)
        if rho > 0.0 and ierr == 0:
            ntypefile = m.ntype - ntypefirst
            nmoltypefirst = c.nmoltype
            units.write(6, fmt(_F36, ntypefile, substr(symbols.fnameato, 1, symbols.ilenfnameato)))
            if ntypefile > 0:
                for it in range(ntypefirst + 1, m.ntype + 1):
                    units.write(6, fmt(_F33, it, Ch(m.type[it], 3)))
            ntypefirst = m.ntype

            #
            # save the parameters of the first molecule of the docking molecule file
            #
            if imix == 1:
                da.xadds = c.xadd[1]
                da.yadds = c.yadd[1]
                da.zadds = c.zadd[1]
                da.phixs = c.phix[1]
                da.phiys = c.phiy[1]
                da.phizs = c.phiz[1]
                da.natos = int(c.nato[1])
                da.amassts = c.amasst[1]

                # Ensure there is enough room to store the coordinates of this molecule.

                if da.natos > da.natomax:
                    da.natomax = da.natos
                    da.alabels = reallocation_routines.reallocate1d_c(da.alabels, 3, da.natomax)
                    da.xas = reallocation_routines.reallocate1d_r(da.xas, da.natomax)
                    da.yas = reallocation_routines.reallocate1d_r(da.yas, da.natomax)
                    da.zas = reallocation_routines.reallocate1d_r(da.zas, da.natomax)
                    da.nbondss = reallocation_routines.reallocate1d_i(da.nbondss, da.natomax)
                    da.iabonds = reallocation_routines.reallocate2d_i(da.iabonds, da.natomax, near)
                    da.tbonds = reallocation_routines.reallocate2d_r(da.tbonds, da.natomax, near)
                for ia in range(c.naf[1], c.nal[1] + 1):
                    da.alabels[ia] = m.alabel[ia]
                    da.xas[ia] = m.xa[ia]
                    da.yas[ia] = m.ya[ia]
                    da.zas[ia] = m.za[ia]
                    da.nbondss[ia] = m.nbond[ia]
                    if m.nbond[ia] > 0:
                        for ib in range(1, m.nbond[ia] + 1):
                            da.iabonds[ia, ib] = m.iabond[ia, ib]
                            da.tbonds[ia, ib] = m.tbond[ia, ib]
                #
                # rotational groups
                #
                da.ngrots = int(c.ngrot[1])
                if da.ngrots > da.ngrotmax:
                    da.ngrotmax = da.ngrots
                    da.rotgrps = reallocation_routines.reallocate1d_c(da.rotgrps, 3, da.ngrotmax)
                    da.labaxis = reallocation_routines.reallocate2d_i(da.labaxis, da.ngrotmax, 2)
                    # NB: new elements of ndihedralangs are uninitialised in the Fortran and
                    # only set below for DIH groups; 0 here (Rust dockato.rs:331)
                    da.ndihedralangs = reallocation_routines.reallocate1d_i(da.ndihedralangs, da.ngrotmax)
                    da.dihedralangs = reallocation_routines.reallocate2d_r(da.dihedralangs, da.ndihmax, da.ngrotmax)
                    da.labrotdihs = reallocation_routines.reallocate3d_i(da.labrotdihs, 2, da.ndihmax, da.ngrotmax)
                if da.ngrots > 0:
                    for ig in range(1, da.ngrots + 1):
                        da.rotgrps[ig] = molecom_mod.rotgrp[ig]
                        da.labaxis[ig, 1] = molecom_mod.labaxi[ig, 1]
                        da.labaxis[ig, 2] = molecom_mod.labaxi[ig, 2]
                        # Check if there are dihedral angles defined
                        if molecom_mod.rotgrp[ig] == "DIH":
                            da.ndihedralangs[ig] = m.ndihedralang[ig]
                            if da.ndihedralangs[ig] > da.ndihmax:
                                da.ndihmax = int(m.ndihedralang[ig])
                                da.dihedralangs = reallocation_routines.reallocate2d_r(
                                    da.dihedralangs, da.ndihmax, da.ngrotmax)
                                da.labrotdihs = reallocation_routines.reallocate3d_i(
                                    da.labrotdihs, 2, da.ndihmax, da.ngrotmax)
                            for k in range(1, m.ndihedralang[ig] + 1):
                                da.dihedralangs[k, ig] = m.dihedralang[k, ig]
                                for j in range(1, 3):
                                    da.labrotdihs[j, k, ig] = m.labrotdih[j, k, ig]
                units.write(6, lw("dockato> natos,ngrots: ", da.natos, da.ngrots))

            ifind = 0
            ifindr = 0
            ifindrr = 0
            while ifind == 0 and ifindr == 0 and ifindrr == 0:
                if imix == 1:

                    # Specify the site on this molecule to be used for docking

                    if ilenatomd == 0:
                        units.write(6, fmt(_F696), advance=False)
                        atomd = _read5_word(atomd, 3)

                    # Specify the site on this molecule to be used for rotating

                    if ilenatomdr == 0:
                        units.write(6, fmt(_F689), advance=False)
                        atomdr = _read5_word(atomdr, 3)

                    # Specify the site on this molecule to be used for a second rotation

                    if ilenatomdrr == 0:
                        units.write(6, fmt(_F687), advance=False)
                        atomdrr = _read5_word(atomdrr, 3)

                    # Determine how many of these docking molecules are needed

                    if ndock == 0:
                        units.write(6, fmt(_F697), advance=False)
                        ndock = units.list_read(5).int(ndock)

                    # Determine maximum distance from the first atom in file for docking purposes.

                    if radmax == 0.0:
                        units.write(6, fmt(_F685), advance=False)
                        radmax = units.list_read(5).real(radmax)

                    iad = 0
                    iadr = 0
                    iadrr = 0

                    # Check the current atoms that at least one site of the requested type exists

                    ia = 0
                    while ia < c.nato[1]:

                        ia = ia + 1
                        if m.alabel[ia] == atomd and ifind == 0:

                            ifind = ia

                            # Save the atom number

                            iad = ia

                            # Save the name of the molecular template file

                        if m.alabel[ia] == atomdr and ifindr == 0:

                            ifindr = ia

                            # Save the atom number

                            iadr = ia

                        if m.alabel[ia] == atomdrr and ifindrr == 0:

                            ifindrr = ia

                            # Save the atom number

                            iadrr = ia

                    # Check that one of the specified atoms on the docking molecule is a 'Qv' atom
                    # (commented out in the Fortran; format 694)

                else:

                    # Specify the site on this molecule to be used for docking on to.

                    if ilenatomv == 0:
                        units.write(6, fmt(_F661), advance=False)
                        atomv = _read5_word(atomv, 3)

                    # Specify the site on this molecule to be rotated on to

                    if ilenatomvr == 0:
                        units.write(6, fmt(_F688), advance=False)
                        atomvr = _read5_word(atomvr, 3)

                    # Specify the site on this molecule to be rotated on to

                    if ilenatomvrr == 0:
                        units.write(6, fmt(_F686), advance=False)
                        atomvrr = _read5_word(atomvrr, 3)

                    # Specify minimum dihedral angle for this atom must be for the docking molecule
                    # (0 for no second rotation)

                    if ilendminrr == 0:
                        units.write(6, fmt(_F684), advance=False)
                        dminrr = units.list_read(5).real(dminrr)

                    # Specify maximum dihedral angle this atom must be for the docking molecule
                    # (0 for no second rotation)

                    if ilendmaxrr == 0:
                        units.write(6, fmt(_F683), advance=False)
                        dmaxrr = units.list_read(5).real(dmaxrr)
                    # Check whether to allow overlaps or not
                    if findex(nooverlapstext, "<undefined>") > 0:
                        units.write(6, fmt(_F784), advance=False)
                        nooverlapstext = fstr(units.read_line(5), 20)

                    # Check that the current atoms have at least one site of the requested type
                    # and check that this atom is charged.

                    ia = 0
                    while ia < m.natom:

                        ia = ia + 1
                        if m.alabel[ia] == atomv and ifind == 0:

                            ifind = ia

                        if m.alabel[ia] == atomvr and ifindr == 0:

                            ifindr = ia

                        if m.alabel[ia] == atomvrr and ifindrr == 0:

                            ifindrr = ia

                    # Check that one of the specified atoms on the docking molecule is a 'Qv' atom
                    # (commented out in the Fortran; format 693)

        else:

            m.ntype = ntypefirst
            c.nmoltype = nmoltypefirst
            units.write(6, fmt(_F599, Ch(symbols.fnameato, 256)))
            if imix == 1:
                ilendockingfile = 0
            if imix == 2:
                ilendockedfile = 0
            imix = imix - 1

    # Maximum distance of docking molecule from origin

    # Now proceed to dock molecules, 1 molecule at a time

    it = 0
    imoladd = 0

    # Counter for the search for unsaturated sites.

    ivcount = 0

    # Check whether overlaps are allowed or not

    # The overlap factor is set to the value contained in nooverlapstext, unless this value is a
    # character string, in which case overlapfactor is set to 0.0
    if findex(nooverlapstext, "yes") > 0:
        nooverlapstext = "0.0"
    elif findex(nooverlapstext, "nooverlaps") > 0:
        nooverlapstext = "1.0"
    try:
        da.overlapfactor = read_internal(nooverlapstext).real(da.overlapfactor)
        ierr = 0
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)
    if ierr != 0:
        da.overlapfactor = f32(0.0)
    nooverlaps = da.overlapfactor > 0.0
    units.write(6, fmt('("dockato> Overlap factor: ",f10.5)', da.overlapfactor))

    while it < ndock and ivcount < m.natom:

        it = it + 1

        # Save the current numbers of atoms, molecules and rotational groups.

        natomsave = m.natom
        natomrsave = m.natomr
        da.nadd = m.natom
        units.write(6, lw("dockato> nadd", da.nadd))
        nmolsave = c.nmol
        ngrrotsave = molecom_mod.ngrrot
        # nrottotsave=nrottot (set but never used)
        dmaxsave = c.dmax

        # Add the new molecule at the end of the current list of molecules.

        ia = int(c.nal[nmolsave])
        im = nmolsave
        ig = ngrrotsave
        if im < mmol:
            im = im + 1
        else:
            units.write(6, fmt(_F601))
            return

        # Centre of mass coordinates

        c.xadd[im] = da.xadds
        c.yadd[im] = da.yadds
        c.zadd[im] = da.zadds
        c.phix[im] = da.phixs
        c.phiy[im] = da.phiys
        c.phiz[im] = da.phizs
        c.nato[im] = da.natos
        c.amasst[im] = da.amassts
        m.natom = natomsave + da.natos

        # Atomic coordinates relative to centre of mass

        c.naf[im] = ia + 1
        units.write(6, lw("dockato> natos,ngrots: ", da.natos, da.ngrots))
        for is_ in range(1, da.natos + 1):
            if ia < matom:
                ia = ia + 1
            else:
                units.write(6, fmt(_F602))
                return
            m.mref[ia] = im
            m.alabel[ia] = da.alabels[is_]
            ato_manipulation_rev4.ftype(ia)
            #
            # do not count atoms which begin with 'Q' or 'q' - these are assumed
            # invisible to the scattering radiation and experiment (although they
            # can be used to calculate other structural quantities if necessary
            #
            if m.qtype[ia] != 1:
                m.natomr = m.natomr + 1
            m.xa[ia] = da.xas[is_]
            m.ya[ia] = da.yas[is_]
            m.za[ia] = da.zas[is_]
            m.nbond[ia] = da.nbondss[is_]
            if m.nbond[ia] > 0:
                for ib in range(1, m.nbond[ia] + 1):
                    m.iabond[ia, ib] = da.iabonds[is_, ib] + da.nadd
                    m.tbond[ia, ib] = da.tbonds[is_, ib]
        c.nal[im] = ia
        #
        # rotational groups
        #
        c.ngrot[im] = da.ngrots
        if c.ngrot[im] > 0:
            molecom_mod.ngrrot = molecom_mod.ngrrot + int(c.ngrot[im])
            molecom_mod_routines.testsize_rotgrps()
            for ih in range(1, c.ngrot[im] + 1):
                ig = ig + 1
                # save molecule in which this rotational group occurs
                molecom_mod.mrot[ig] = im
                units.write(6, lw("dockato> ig, nadd", ig, da.nadd))
                molecom_mod.rotgrp[ig] = da.rotgrps[ih]
                # for rotational groups read atom numbers of two atoms to form axis of rotation
                for j in range(1, 3):
                    molecom_mod.labaxi[ig, j] = da.labaxis[ih, j] + da.nadd
                m.ndihedralang[ig] = da.ndihedralangs[ih]
                for k in range(1, m.ndihedralang[ig] + 1):
                    m.dihedralang[k, ig] = da.dihedralangs[k, ih]
                    for j in range(1, 3):
                        m.labrotdih[j, k, ig] = da.labrotdihs[j, k, ih] + da.nadd
        m.natom = ia
        c.nmol = im
        # ngrrot=ig
        units.write(6, lw("dockato > ", m.natom, m.natomr, c.nmol, molecom_mod.ngrrot,
                          *[c.ngrot[i] for i in range(1, c.nmol + 1)]))

        # New value of dmax

        c.dmax = dmaxsave * powf(f32(m.natomr) / f32(natomrsave), r4(1.0) / r4(3.0))

        # Of the available atoms search for a suitable site to be docked onto and a site
        # suitable for docking.

        ifind = 0
        iav = 0
        iavr = natomsave
        iavrr = natomsave
        iad = m.natom
        iadr = m.natom
        iadrr = m.natom

        # Identify the next unsaturated site of the specified type within the atoms to be
        # docked onto.

        while iav < natomsave and ifind == 0:

            if iad >= m.natom:

                iad = natomsave
                iav, nvalencev, xv, yv, zv = search_unsaturated(iav, natomsave, nvalencev, atomv, atomvr, xv, yv, zv)
                units.write(6, lw("dockato> 1 ", iav, nvalencev, xv, yv, zv))

            if nvalencev > 0:

                # O.K. so found an unsaturated site.
                # Now identify an unsaturated site on the docking molecule of specified type

                # Molecule number for this site

                imv = int(m.mref[iav])
                iad, nvalenced, xd, yd, zd = search_unsaturated(iad, m.natom, nvalenced, atomd, atomdr, xd, yd, zd)

                units.write(6, lw("dockato> 2 ", iad, nvalenced, xd, yd, zd))

                if nvalenced > 0:

                    # OK so this atom pair apparently meets the requirements for a valid docking.
                    # (see the Fortran for the full comment)
                    # So the next thing to do is to perform a translation of this molecule as
                    # required, and then check for overlaps.

                    imd = int(m.mref[iad])
                    units.write(6, lw(iav, iad, xv, yv, zv, xd, yd, zd))
                    c.xadd[imd] = da.xadds + xv - xd
                    c.yadd[imd] = da.yadds + yv - yd
                    c.zadd[imd] = da.zadds + zv - zd

                    # Now perform rotations as required, if possible

                    # Search for a site on the docking molecule to be used for defining the rotation.

                    iadr = 0
                    tbondrot = f32(1000.0)
                    ib = 0
                    while ib < m.nbond[iad]:

                        # We choose the atom of the correct type which is the shortest distance
                        # from iad.

                        ib = ib + 1
                        tbondtest = abs(m.tbond[iad, ib])
                        if m.alabel[m.iabond[iad, ib]] == atomdr and tbondtest < tbondrot:
                            iadr = int(m.iabond[iad, ib])
                            tbondrot = tbondtest
                    units.write(6, lw(iad, iadr, tbondrot))

                    # Repeat this search on the molecule to be docked onto. In this case we are
                    # looking for an atom of type atomvr which is a neighbour of iav at a
                    # distance of tbondrot.

                    iavr = 0
                    ib = 0
                    while ib < m.nbond[iav] and iavr == 0:

                        ib = ib + 1
                        ivref = int(m.iabond[iav, ib])
                        if m.alabel[ivref] == atomvr and abs(m.tbond[iav, ib]) == tbondrot:

                            iavr = ivref

                    units.write(6, lw(iav, iavr, tbondrot))
                    if iadr > 0 and iavr > 0:

                        # If iadr and iavr are defined we perform the first rotation

                        ierr = rotatem(imd, iad, iadr, imv, iav, iavr, ierr)

                        # Now attempt to perform 2nd rotation if required and possible

                        # Now perform a second rotation if possible about the axis defined by
                        # iad to iadr

                        iadrr, iavrr, dihedralangfixed = random_rotatem(
                            atomdrr, iad, iadr, iadrr, atomvrr, iav, iavr, iavrr, dminrr, dmaxrr, dihedralangfixed)

                        units.write(6, lw(iav, iavr, iavrr, iad, iadr, iadrr, tbondrot))

                        # Test for overlaps if required

                        ioverlap = 0
                        if nooverlaps:
                            ioverlap = test_overlap(ioverlap, imd, 0, natomsave, iav, iavr, iavrr, iad, iadr,
                                                    iadrr, da.overlapfactor)

                        if ioverlap == 0:

                            ifind = 1

                        else:

                            ifind = 0

                            # If this rotation is not allowed, restore the coordinates of
                            # the docking molecule

                            c.xadd[imd] = da.xadds
                            c.yadd[imd] = da.yadds
                            c.zadd[imd] = da.zadds

                            ia = int(c.naf[imd]) - 1
                            for is_ in range(1, c.nato[imd] + 1):
                                ia = ia + 1
                                m.xa[ia] = da.xas[is_]
                                m.ya[ia] = da.yas[is_]
                                m.za[ia] = da.zas[is_]

        # If ifind.eq.0 at this point, we have searched all possible pairs of sites,
        # and not one is allowed.
        # Therefore we have to abandon this docking attempt.

        # If no site was found, ignore this added molecule.

        units.write(6, fmt(_F664, ifind))
        if ifind == 0:

            m.natom = natomsave
            m.natomr = natomrsave
            c.nmol = nmolsave
            molecom_mod.ngrrot = ngrrotsave
            c.dmax = dmaxsave

        else:

            # Molecule docking successful. Therefore bond the new molecule to the old. This means
            # eliminating atoms iad and iadr, and giving there respective neighbours to the atoms
            # they are bonded to.

            fixeddihedral = dminrr == dmaxrr
            iavrr, iadrr = bond_molecule(iav, iavr, iavrr, iad, iadr, iadrr, fixeddihedral, dihedralangfixed)
            units.write(6, lw("dockato> ", c.nmol, m.natom, m.natomr, c.nato[imv], c.nal[imv], c.nato[imd]))

            imoladd = imoladd + 1
            m.dmaxmol = m.sizefactor * c.dmax

            # If search was successful reset docked site counter to ensure it rescans the current
            # site, in case another molecule can be docked here.

            iav = iav - 1

    # We now need to go through all rotational groups and setup the lists of atoms to be
    # rotated within each group.

    # First need to get the directly bonded atoms before we can do this,
    # since get_rotgrp works only on the directly bonded atoms

    ato_manipulation_rev4.get_direct_bonds()

    molecom_mod.nrottot = 0
    molecom_mod.nrottotmax = 0
    for ig in range(1, molecom_mod.ngrrot + 1):
        ir1 = int(molecom_mod.labaxi[ig, 1])
        ir2 = int(molecom_mod.labaxi[ig, 2])
        da.labrott1, da.nrott1 = ato_manipulation_rev4.get_rotgrp(ir1, ir2, da.labrott1, da.nrott1)
        da.labrott2, da.nrott2 = ato_manipulation_rev4.get_rotgrp(ir2, ir1, da.labrott2, da.nrott2)
        # Use the smallest of the two groups
        if da.nrott1 < da.nrott2:
            molecom_mod.labaxi[ig, 1] = ir1
            molecom_mod.labaxi[ig, 2] = ir2
            molecom_mod.nrot[ig] = da.nrott1
            molecom_mod.rotf[ig] = molecom_mod.nrottot + 1
            molecom_mod.nrottot = molecom_mod.nrottot + int(molecom_mod.nrot[ig])
            molecom_mod.rotl[ig] = molecom_mod.nrottot
            molecom_mod_routines.testsize_labrot()
            for i in range(1, molecom_mod.nrot[ig] + 1):
                labref = int(molecom_mod.rotf[ig]) + i - 1
                units.write(6, lw("dockato> get labrot> ", ig, molecom_mod.rotf[ig], molecom_mod.rotl[ig],
                                  molecom_mod.nrottot, molecom_mod.nrottotmax, labref))
                molecom_mod.labrot[labref] = da.labrott1[i]
        else:
            molecom_mod.labaxi[ig, 1] = ir2
            molecom_mod.labaxi[ig, 2] = ir1
            molecom_mod.nrot[ig] = da.nrott2
            molecom_mod.rotf[ig] = molecom_mod.nrottot + 1
            molecom_mod.nrottot = molecom_mod.nrottot + int(molecom_mod.nrot[ig])
            molecom_mod.rotl[ig] = molecom_mod.nrottot
            molecom_mod_routines.testsize_labrot()
            for i in range(1, molecom_mod.nrot[ig] + 1):
                labref = int(molecom_mod.rotf[ig]) + i - 1
                units.write(6, lw("dockato> get labrot> ", ig, molecom_mod.rotf[ig], molecom_mod.rotl[ig],
                                  molecom_mod.nrottot, molecom_mod.nrottotmax, labref))
                molecom_mod.labrot[labref] = da.labrott2[i]

    units.write(6, fmt(_F662, imoladd))
    if imoladd < ndock:
        units.write(6, fmt(_F663))

    #
    # get the new filename to put the mixture in
    #
    if ilenoutputfile > 0:
        cm.nwords = 1
        cm.ncf[1] = 1
        cm.ncl[1] = ilenoutputfile
        cm.sometext = setsub(cm.sometext, 1, ilenoutputfile, outputfile, cm.nstsize)
    else:
        cm.nwords = 0
    while cm.nwords <= 0:
        units.write(6, fmt(_F35), advance=False)
        try:
            cm.sometext = fstr(units.read_line(5), cm.nstsize)
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = iostat_of(e)
        cm.nwords = parse.parse(cm.sometext, cm.nwords, cm.ncf, cm.ncl, 5, cm.nstsize, cm.nworddim)
        if cm.nwords <= 0:
            units.write(6, fmt(_F499))
    readfile = fstr(substr(cm.sometext, cm.ncf[1], cm.ncl[1]), 256)
    ilenreadfile = cm.ncl[1] - cm.ncf[1] + 1
    #
    # strip out any extension that was supplied - extension must be .ato, otherwise
    # nothing is done.
    #
    ilendot = findex(readfile, ".") - 1
    if ilendot > 0:
        readfile = fstr(substr(readfile, 1, ilendot), 256)
        ilenreadfile = ilendot
    symbols.fnameato = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                            + substr(readfile, 1, ilenreadfile) + ".ato", 256)
    symbols.ilenfnameato = directory_name.nchardir + ilenreadfile + 4

    # Reset the molecule template file for this molecule to the current filename, so that
    # any changes to this file will be reflected in any .ato files which are derived from
    # this file

    c.moltemplate[c.mtype[imv]] = fstr(substr(readfile, 1, ilenreadfile), 256)

    #
    # write the new .ato file
    #
    ato_read_and_save.moloutc(symbols.fnameato)
    return


_F701 = "('rotatem> Molecules do not overlap - cannot rotate')"
_F702 = "('rotatem> No rotation is requested')"


@routine
def rotatem(imd, iad, iadr, imv, iav, iavr, ierr):
    """dockato.f:rotatem -- rotate molecule imd about atom iad onto atom iavr of imv.

    Rotates molecule imd about atom iad so that atom iadr on imd overlaps atom iavr
    on molecule imv. It is assumed the two molecules have been previously shifted so that
    atom iad on imd overlaps atom iav on imv. If not the routine returns without doing
    anything. Minimum image convention is not used for this check so this method is only
    to be used in conjunction with dockato.
    Returns ierr.
    """
    m = molecom
    c = clustercom
    scapro = auxiliary_routines.scapro
    vecpro = auxiliary_routines.vecpro
    xa = m.xa
    ya = m.ya
    za = m.za

    dsq = real_array(2)

    # vvector defines the vector on the existing molecule on to which the docking molecule
    # is to be rotated
    # rvector defines the vector on the docking molecule which is to be rotated to coincide
    # with vvector.
    # pvector is the vector about which the rotation occurs, namely the vector normal to the
    # plane defined by vvector and rvector

    vvector = real_array(3)
    rvector = real_array(3)
    pvector = real_array(3)
    # xvector is the vector from centre of rotation to an individual atom
    # xpvector is the component of this vector along the axis of rotation
    # xrvector is the component of this vector perpendicular to the axis of rotation
    # xrpvector is the vector perpendicular to both axis of rotation and xrvector

    xvector = real_array(3)
    xpvector = real_array(3)
    xrvector = real_array(3)
    xrpvector = real_array(3)

    # First check that atom iad does indeed overlap atom iav

    units.write(6, lw(imv, iav, iavr, imd, iad, iadr))
    units.write(6, lw(m.itype[iav], m.itype[iavr], m.itype[iad], m.itype[iadr]))

    ierr = 1
    xadiff = c.xadd[imv] - c.xadd[imd] + xa[iav] - xa[iad]
    yadiff = c.yadd[imv] - c.yadd[imd] + ya[iav] - ya[iad]
    zadiff = c.zadd[imv] - c.zadd[imd] + za[iav] - za[iad]
    dsq[1] = xadiff * xadiff + yadiff * yadiff + zadiff * zadiff
    qradiussum = f32(6.0) * (m.qradius[m.itype[iad]] + m.qradius[m.itype[iav]])
    units.write(6, lw(xadiff, yadiff, zadiff, sqrtf(dsq[1]), qradiussum))
    if sqrtf(dsq[1]) > qradiussum:
        units.write(6, fmt(_F701))
        return ierr
    elif iavr == 0 or iadr == 0:
        units.write(6, fmt(_F702))
        return ierr

    # generate unit vectors for vvector and rvector, measured from the atom
    # being docked onto

    vvector[1] = c.xadd[imv] + xa[iavr] - c.xadd[imd] - xa[iad]
    vvector[2] = c.yadd[imv] + ya[iavr] - c.yadd[imd] - ya[iad]
    vvector[3] = c.zadd[imv] + za[iavr] - c.zadd[imd] - za[iad]
    rvector[1] = xa[iadr] - xa[iad]
    rvector[2] = ya[iadr] - ya[iad]
    rvector[3] = za[iadr] - za[iad]
    vnorm = sqrtf(scapro(vvector, vvector))
    rnorm = sqrtf(scapro(rvector, rvector))
    for i in range(1, 4):
        vvector[i] = vvector[i] / vnorm
        rvector[i] = rvector[i] / rnorm

    # angle between these vectors

    cphi = scapro(vvector, rvector)
    if abs(cphi) > 1.0:
        cphi = sign(f32(1.0), cphi)
    sphi = sqrtf(f32(1.0) - cphi * cphi)

    # No need to rotate if the atoms already overlap

    units.write(6, lw(cphi, sphi))
    if sphi > 0.001:

        # define the vector normal to the plane of vvector and rvector about which the molecule
        # will be rotated

        vecpro(rvector, vvector, pvector)

        # Define the unit vector in this direction

        rab = sqrtf(scapro(pvector, pvector))
        units.write(6, lw(rab, *[pvector[i] for i in range(1, 4)]))

        pvector[1] = pvector[1] / rab
        pvector[2] = pvector[2] / rab
        pvector[3] = pvector[3] / rab
        units.write(6, lw(*[pvector[i] for i in range(1, 4)]))
        # step through atoms in group
        for ia in range(c.naf[imd], c.nal[imd] + 1):
            if ia != iad:
                # Vector from the centre of rotation to this atom
                xvector[1] = xa[ia] - xa[iad]
                xvector[2] = ya[ia] - ya[iad]
                xvector[3] = za[ia] - za[iad]
                # Component of this vector parallel to axis of rotation
                frac = scapro(pvector, xvector)
                for i in range(1, 4):
                    xpvector[i] = frac * pvector[i]
                units.write(6, lw(ia, *[v for i in range(1, 4) for v in (xvector[i], xpvector[i])]))

                # Component of vector perpendicular to axis of rotation
                for i in range(1, 4):
                    xrvector[i] = xvector[i] - xpvector[i]
                # Vector perpendicular to plane defined by axis of rotation and xrvector
                vecpro(pvector, xrvector, xrpvector)
                # Generate rotated coordinates. The vector component parallel to the axis of
                # rotation is not rotated
                xa[ia] = xa[iad] + xpvector[1]
                ya[ia] = ya[iad] + xpvector[2]
                za[ia] = za[iad] + xpvector[3]
                xa[ia] = xa[ia] + xrvector[1] * cphi + xrpvector[1] * sphi
                ya[ia] = ya[ia] + xrvector[2] * cphi + xrpvector[2] * sphi
                za[ia] = za[ia] + xrvector[3] * cphi + xrpvector[3] * sphi

    elif cphi < 1.0:

        # Only can get here if sphi=0 and cphi=-1. In that case we simply invert the coordinates
        # about the centre of rotation. (Until I think of something better to do!)

        for ia in range(c.naf[imd], c.nal[imd] + 1):

            xa[ia] = xa[iad] - (xa[ia] - xa[iad])
            ya[ia] = ya[iad] - (ya[ia] - ya[iad])
            za[ia] = za[iad] - (za[ia] - za[iad])

    ierr = 0
    return ierr


_F100_random_rotatem = "(/'random_rotatem> Beginning...',2(1x,i5,1x,i5,1x,a3)/)"
_F101_random_rotatem = "(/'random_rotatem> Atoms selected: ',6(1x,i5))"
_F102_random_rotatem = "(/'random_rotatem> Angle rotated (cosine and sine) ',1x,e13.6,1x,e13.6)"
_F103_random_rotatem = "(/'random_rotatem> Angle old (cosine and sine) ',1x,e13.6,1x,e13.6)"
_F104_random_rotatem = "(/'random_rotatem> Angle check (cosine and sine) ',1x,e13.6,1x,e13.6/)"


@routine
def random_rotatem(atomdrr, iad, iadr, iadrr, atomvrr, iav, iavr, iavrr, dminrr, dmaxrr, dihedralangfixed):
    """dockato.f:random_rotatem -- random dihedral rotation about the axis iad-iadr.

    Generate a random dihedral rotation of the docking molecule about the axis defined
    by iad,iadr.  dminrr and dmaxrr give the minimum and maximum dihedral angles requested.
    The dihedral angle is measured using the existing axis of rotation and, plus the atom
    of type atomvrr on the docked molecule and the atom of type atomdrr on the docking
    molecule.
    Returns (iadrr, iavrr, dihedralangfixed).
    """
    m = molecom
    c = clustercom
    scapro = auxiliary_routines.scapro
    vecpro = auxiliary_routines.vecpro
    xa = m.xa
    ya = m.ya
    za = m.za

    # vectors to be used for rotating molecules

    rvector = real_array(3)
    vvector = real_array(3)
    rvpvector = real_array(3)
    dvector = real_array(3)
    rdpvector = real_array(3)
    drot = real_array(3)
    # NB: cphiold, sphiold, sphinew are not set on every path in the Fortran (they are then
    # only written to standard output); 0 here
    cphiold = f32(0.0)
    sphiold = f32(0.0)
    sphinew = f32(0.0)

    units.write(6, fmt(_F100_random_rotatem, iad, iadr, Ch(atomdrr, 3), iav, iavr, Ch(atomvrr, 3)))

    # Coordinates of axis of rotation

    imd = int(m.mref[iad])

    # Vector defining axis of rotation

    rvector[1] = xa[iadr] - xa[iad]
    rvector[2] = ya[iadr] - ya[iad]
    rvector[3] = za[iadr] - za[iad]

    # Convert to a unit vector in the same direction

    rnorm = sqrtf(scapro(rvector, rvector))
    for i in range(1, 4):
        rvector[i] = rvector[i] / rnorm

    # Angle of rotation, selected at random between dminrr and dmaxrr

    pi = PI4ATAN
    arg = f32(2.0) * pi * (dminrr + auxiliary_routines.ran1() * (dmaxrr - dminrr)) / f32(360.0)
    dihedralangfixed = arg * f32(180.0) / pi
    sphistep = sinf(arg)
    cphistep = cosf(arg)

    # Search for an appropriate atom on the docking molecule. We choose the one that is closest
    # to iad

    imd = int(m.mref[iad])
    iat = int(c.naf[imd]) - 1
    iadrr = 0
    dsqmin = f32(1.0e30)

    while iat < c.nal[imd] and atomdrr != "0":

        iat = iat + 1

        # The atom must be of the correct type and cannot be either the rotation centre (iad)
        # or the atom already rotated (iadr).

        if iat != iad and iat != iadr and m.alabel[iat] == atomdrr:

            # Vector defining the position of this atom relative to the centre of rotation

            dvector[1] = xa[iat] - xa[iad]
            dvector[2] = ya[iat] - ya[iad]
            dvector[3] = za[iat] - za[iad]

            # Square of distance to iad

            dsq = scapro(dvector, dvector)
            if dsq < dsqmin:

                # Cross product of rvector onto dvector gives a vector perpendicular to both

                vecpro(rvector, dvector, rdpvector)

                # Normalise the perpendicular vector

                rdpnorm = sqrtf(scapro(rdpvector, rdpvector))
                for i in range(1, 4):
                    rdpvector[i] = rdpvector[i] / rdpnorm

                # Dot product with rvector should not be >1 or <-1

                cphi = scapro(rdpvector, rvector)
                if cphi > -1.0 and cphi < 1.0:

                    dsqmin = dsq
                    iadrr = iat

    # Search for an appropriate atom on the docked molecule. We choose the atom
    # which is closest to iavr

    imv = int(m.mref[iav])
    iat = int(c.naf[imv]) - 1
    iavrr = 0
    dsqmin = f32(1.0e30)

    while iat < c.nal[imv] and atomvrr != "0":

        iat = iat + 1

        # The atom must be of the correct type and cannot be either the rotation centre (iav)
        # or the atom already rotated onto (iavr).

        if iat != iav and iat != iavr and m.alabel[iat] == atomvrr:

            # Vector defining the position of this atom relative to the axis of rotation

            vvector[1] = xa[iat] - xa[iav]
            vvector[2] = ya[iat] - ya[iav]
            vvector[3] = za[iat] - za[iav]

            # Square of distance to iad

            dsq = scapro(vvector, vvector)
            if dsq < dsqmin:

                # Cross product of rvector onto vvector gives a vector perpendicular to both

                vecpro(rvector, vvector, rvpvector)

                # Normalise the perpendicular vector

                rvpnorm = sqrtf(scapro(rvpvector, rvpvector))
                for i in range(1, 4):
                    rvpvector[i] = rvpvector[i] / rvpnorm

                # Dot product with rvector should not be 1 or -1, otherwise vvector is parallel
                # with rvector and so cannot be used to rotate.

                cphi = scapro(rvpvector, rvector)
                units.write(6, fmt(_F100_random_rotatem))
                units.write(6, lw(iat, cphi))
                if cphi > -1.0 and cphi < 1.0:

                    dsqmin = dsq
                    iavrr = iat

    # Now determine angle of rotation

    if iadrr > 0 and iavrr > 0:

        # In this case we attempt to rotate molecule imd so that the angle between rdpvector
        # and rvpvector = arg.

        # Get the current angle between rdpvector and rvpvector

        # Cosine of the angle is given by the dot product.

        cphiold = scapro(rdpvector, rvpvector)
        if abs(cphiold) >= 1.0:
            cphiold = sign(f32(1.0), cphiold)

        # Need to determine the sine. Form a vector perpendicular to both
        # rvpvector and rdpvector - this will be parallel or antiparallel to rvector

        vecpro(rdpvector, rvpvector, drot)

        # If we have to rotate backwards, this will be antiparallel to rvector, otherwise it
        # will be parallel. Since rvector is a unit vector, the dot product wih rvector will
        # give the sine of the angle between rdpvector and rvpvector

        sphiold = scapro(drot, rvector)

        # Define the new cosine and sine required

        cphinew = cphiold * cphistep + sphiold * sphistep
        sphinew = sphiold * cphistep - cphiold * sphistep

    elif iadrr > 0:

        # In this case, if iavrr is not found, simply perform a rotation of the docking
        # molecule by the amount phinew

        sphinew = sphistep
        cphinew = cphistep

    else:

        # Otherwise do nothing

        cphinew = f32(2.0)

    # Now perform rotation as requested

    if cphinew >= -1.0 and cphinew <= 1.0:

        iat = int(c.naf[imd]) - 1
        while iat < c.nal[imd]:

            iat = iat + 1
            if iat != iad and iat != iadr:

                # Vector defining the position of this atom relative to the centre of rotation

                dvector[1] = xa[iat] - xa[iad]
                dvector[2] = ya[iat] - ya[iad]
                dvector[3] = za[iat] - za[iad]

                # Component of dvector along the axis of rotation

                frac = scapro(dvector, rvector)

                # Cross product of rvector onto dvector gives a vector perpendicular to both

                vecpro(rvector, dvector, rdpvector)

                for i in range(1, 4):

                    # Vector parallel to axis of rotation

                    dpar = frac * rvector[i]

                    # Vector perpendicular to axis of rotation

                    dperp = dvector[i] - dpar

                    # Create the rotated vector

                    drot[i] = dpar + dperp * cphinew + rdpvector[i] * sphinew

                # Generate rotated coordinates

                xa[iat] = xa[iad] + drot[1]
                ya[iat] = ya[iad] + drot[2]
                za[iat] = za[iad] + drot[3]

    # As a test, if iavrr.gt.0, check that the new angle between iavrr and iadrr is
    # as requested.

    cphitest = f32(99.0)
    sphitest = f32(99.0)
    if iavrr > 0:

        # Vector defining the position of this atom relative to the centre of rotation

        dvector[1] = xa[iadrr] - xa[iad]
        dvector[2] = ya[iadrr] - ya[iad]
        dvector[3] = za[iadrr] - za[iad]

        # Cross product of rvector onto dvector gives a vector perpendicular to both

        vecpro(rvector, dvector, rdpvector)

        # Normalise the perpendicular vector

        rdpnorm = sqrtf(scapro(rdpvector, rdpvector))
        for i in range(1, 4):
            rdpvector[i] = rdpvector[i] / rdpnorm

        # Cosine of the angle is given by the dot product.

        cphitest = scapro(rdpvector, rvpvector)
        if abs(cphitest) >= 1.0:
            cphitest = sign(f32(1.0), cphitest)

        # Need to determine the the sine. Form a vector perpendicular to both
        # rvpvector and rdpvector - this will be parallel or antiparallel to rvector

        vecpro(rdpvector, rvpvector, drot)

        # If we have to rotate backwards, this will be antiparallel to rvector, otherwise it
        # will be parallel

        sphitest = scapro(drot, rvector)

    units.write(6, fmt(_F101_random_rotatem, iad, iadr, iadrr, iav, iavr, iavrr))
    units.write(6, fmt(_F102_random_rotatem, cphinew, sphinew))
    if iavrr > 0:
        units.write(6, fmt(_F103_random_rotatem, cphiold, sphiold))
        units.write(6, fmt(_F104_random_rotatem, cphitest, sphitest))

    return iadrr, iavrr, dihedralangfixed


@routine
def search_unsaturated(iav, iavlast, nvalence, atomv, atomvr, xv, yv, zv):
    """dockato.f:search_unsaturated -- first unsaturated site of the specified type.

    On entry iav is the number of the atom to search from. The first search should have
    iav=0.  On exit iav is the atom number of the atom found suitable for docking.
    atomv is the label of the site to be docked onto; xv,yv,zv are the coordinates of
    this docking site.
    Returns (iav, nvalence, xv, yv, zv).
    """
    m = molecom
    c = clustercom
    # NB: imv is written before it is assigned when the loop does not run (stack garbage in
    # the Fortran, standard output only); 0 here
    imv = 0
    bondfound = f32(0.0)

    # We search the neighbours of iav to see whether there any at the existing value of iav.
    # This is so that if a previous search found this
    # atom, but discarded because of overlaps, we do not search it again.

    nvalence = 0
    nvalencefound = 0
    bondsearch = f32(1000.0)
    ia = iav

    while ia < iavlast:

        # Search for an atom of the specified type

        ia = ia + 1
        imv = int(m.mref[ia])

        units.write(6, lw("search_unsaturated> ", ia, iavlast, " ", Ch(m.alabel[ia], 3), " ", Ch(atomv, 3),
                          " ", Ch(atomvr, 3)))

        if m.alabel[ia] == atomv:

            # Check the neighbours of ia and sum the number of ligands as specified by atomvr
            # atoms in the neighbourhood of iav

            nvalencefound, bondfound = find_valence(ia, nvalencefound, atomvr, bondfound)

            # If a non-zero valence value was identified, select the shortest bond distance found

            if nvalencefound > 0 and bondfound < bondsearch:
                iav = ia
                bondsearch = bondfound
                nvalence = nvalencefound

    # If zero valence if found, we save the last value of ia, so that it is not searched again

    if nvalence == 0:

        iav = ia

    units.write(6, lw("search_unsaturated> ", iav, imv, nvalence, " ", Ch(m.alabel[iav], 3), " ",
                      Ch(atomvr, 3)))

    # Coordinates of this atom

    imv = int(m.mref[iav])
    xv = c.xadd[imv] + m.xa[iav]
    yv = c.yadd[imv] + m.ya[iav]
    zv = c.zadd[imv] + m.za[iav]

    return iav, nvalence, xv, yv, zv


# 100 format(/'find_valence> ',$)
_F100_find_valence = "(/'find_valence> ')"


@routine
def find_valence(iav, nvalence, atomvr, bondsearch):
    """dockato.f:find_valence -- number of atomvr neighbours of iav at the shortest distance.

    Search through the neighbours of iav and calculate the number of overlaps,
    mean position of the charge, and total charge.
    Returns (nvalence, bondsearch).
    """
    m = molecom

    # Look at all the bonded neighbours of iav, if any, and find the shortest distance
    # to atoms of type atomvr. This distance is assumed to be the distance of local valence
    # atoms.

    ib = 0
    ifind = 0
    bondsearch = f32(1000.0)

    while ib < m.nbond[iav]:

        ib = ib + 1
        if m.alabel[m.iabond[iav, ib]] == atomvr:

            ifind = ifind + 1
            if abs(m.tbond[iav, ib]) < bondsearch:
                bondsearch = abs(m.tbond[iav, ib])

    units.write(6, fmt(_F100_find_valence), advance=False)
    units.write(6, lw(bondsearch, Ch(atomvr, 3)))

    # Now search for neighbours of type atomvr which are at this minimum distance.

    ib = 0
    nvalence = 0

    while ib < m.nbond[iav]:

        ib = ib + 1

        if m.alabel[m.iabond[iav, ib]] == atomvr and abs(m.tbond[iav, ib]) == bondsearch:

            nvalence = nvalence + 1

    units.write(6, fmt(_F100_find_valence), advance=False)
    units.write(6, lw(nvalence))

    return nvalence, bondsearch


_F100_test_overlap = "(/'test_overlap> ')"


@routine
def test_overlap(ioverlap, imd, natomfirst, natomsave, iav, iavr, iavrr, iad, iadr, iadrr, overlapfactor):
    """dockato.f:test_overlap -- check for overlaps, according to the Lennard-Jones parameters.

    Returns ioverlap.
    """
    m = molecom
    c = clustercom

    ioverlap = 0

    # Step through the atoms of the new molecule

    inew = int(c.naf[imd]) - 1

    while inew < c.nal[imd] and ioverlap == 0:

        inew = inew + 1
        ityped = int(m.itype[inew])

        # Coordinates of this atom

        xt = c.xadd[imd] + m.xa[inew]
        yt = c.yadd[imd] + m.ya[inew]
        zt = c.zadd[imd] + m.za[inew]

        # Test whether this atom is a neighbour of iad, iadr or iadrr

        iadneighb = find_neighbour(iad, inew, 0)
        iadrneighb = find_neighbour(iadr, inew, 0)
        iadrrneighb = find_neighbour(iadrr, inew, 0)  # noqa: F841 (set but not used, as in the Fortran)

        # Step through the other atoms and make sure none of them come within the prescribed
        # distance

        iat = natomfirst
        while iat < natomsave and ioverlap == 0:

            iat = iat + 1

            # Test whether this atom is a neighbour of iad, iadr or iadrr

            iavneighb = find_neighbour(iav, iat, 0)
            iavrneighb = find_neighbour(iavr, iat, 0)
            iavrrneighb = find_neighbour(iavrr, iat, 0)  # noqa: F841 (set but not used, as in the Fortran)

            imt = int(m.mref[iat])
            itt = int(m.itype[iat])
            epsilont = sqrtf(m.epscor[itt] * m.epscor[ityped])

            # Only do the test if epsilont is non-zero, and iat is not a neighbour of iavr, and
            # iat is not a neighbour of iav.

            if epsilont > 0.0 and ((iavneighb == 0 or iadrneighb == 0)
                                   or (iavrneighb == 0 or iadneighb == 0)):

                sigmat = f32(0.5) * overlapfactor * (m.sigcor[itt] + m.sigcor[ityped]) / c.dmax
                xx = (c.xadd[imt] + m.xa[iat] - xt) / c.dmax
                yy = (c.yadd[imt] + m.ya[iat] - yt) / c.dmax
                zz = (c.zadd[imt] + m.za[iat] - zt) / c.dmax
                xx, yy, zz = minimagec.minimagec(xx, yy, zz)
                dtest = sqrtf(xx * xx + yy * yy + zz * zz)
                if dtest < sigmat:
                    ioverlap = 1
                    units.write(6, fmt(_F100_test_overlap), advance=False)
                    units.write(6, lw(iat, imt, itt, inew, imd, ityped, dtest, sigmat))

    return ioverlap


@routine
def find_neighbour(iad, inew, iadneighb):
    """dockato.f:find_neighbour -- 1 if inew is iad or one of its bonded neighbours.

    Returns iadneighb.
    """
    m = molecom

    iadneighb = 0

    # We include the case when the atom specified (inew) is equal to the base atom specified
    # (iad) as a neighbour, since this atom should not be included in overlap calculation.

    if iad == inew:
        iadneighb = 1
    if iad > 0:

        ib = 0
        while ib < m.nbond[iad] and iadneighb == 0:

            ib = ib + 1
            if m.iabond[iad, ib] == inew:
                iadneighb = 1

    return iadneighb


@routine
def bond_molecule(iav, iavr, iavrr, iad, iadr, iadrr, fixeddihedral, dihedralangfixed):
    """dockato.f:bond_molecule -- bond the docked molecule imd to molecule imv.

    Adds the list of atoms in the molecule given by iad,iadr, iadrr, to the molecule given
    by iav, iavr, iavrr.
    Returns (iavrr, iadrr).
    """
    # use reallocation_routines
    # use molecom_mod
    # use dockato_arrays
    m = molecom
    c = clustercom
    mm = molecom_mod
    da = dockato_arrays

    # List of atoms to be kept and erradicated

    ikeep = int_array(near)
    ilose = int_array(near)

    # NB: ig is undefined in the Fortran when neither of the two blocks below that set it
    # is executed (no rotational groups on the docking molecule and nato(imv) <= 2); 0 here
    ig = 0

    # First set all the coordinates of the atoms in the new molecule to be consistent with the
    # centre of mass coordinates of the old molecule.

    imd = int(m.mref[iad])
    imv = int(m.mref[iav])
    units.write(6, lw("bond_molecule> 1 ", c.nmol, *[c.nato[j] for j in range(1, c.nmol + 1)], m.natom))
    units.write(6, lw("bond_molecule> 2 ", iav, iavr, iavrr, iad, iadr, iadrr))
    units.write(6, lw("bond_molecule> 3 ", imv, imd, c.ngrot[imv], da.natos, da.ngrots, mm.ngrrot))

    # Transfer any existing rotational groups on imd to imv

    if da.ngrots > 0:

        # Find out the last group corresponding to imv

        iglast = 0
        im = 0
        while im < imv:

            im = im + 1
            iglast = iglast + int(c.ngrot[im])

        # Shift all the existing groups above this up by ngrots.

        ig = mm.ngrrot
        ig1 = ig - da.ngrots
        while ig1 > iglast:

            # NB: the Fortran copies mref (the atom -> molecule map) here, not mrot
            m.mref[ig] = m.mref[ig1]
            mm.rotgrp[ig] = mm.rotgrp[ig1]
            for j in range(1, 3):
                mm.labaxi[ig, j] = mm.labaxi[ig1, j]
            m.ndihedralang[ig] = m.ndihedralang[ig1]
            for k in range(1, m.ndihedralang[ig] + 1):
                m.dihedralang[k, ig] = m.dihedralang[k, ig1]
                for j in range(1, 3):
                    m.labrotdih[j, k, ig] = m.labrotdih[j, k, ig1]
            ig = ig - 1
            ig1 = ig1 - 1

        # Finally insert the saved rotational groups into the space created

        c.ngrot[imv] = c.ngrot[imv] + da.ngrots
        ig = iglast
        ih = 0

        while ih < da.ngrots:

            ig = ig + 1
            ih = ih + 1
            m.mref[ig] = imv
            mm.rotgrp[ig] = da.rotgrps[ih]
            for j in range(1, 3):
                mm.labaxi[ig, j] = da.labaxis[ih, j] + da.nadd
            m.ndihedralang[ig] = da.ndihedralangs[ih]
            for k in range(1, m.ndihedralang[ig] + 1):
                m.dihedralang[k, ig] = da.dihedralangs[k, ih]
                for j in range(1, 3):
                    m.labrotdih[j, k, ig] = da.labrotdihs[j, k, ih] + da.nadd
        units.write(6, lw("bond_molecule> 4 ", imv, imd, c.ngrot[imv], da.ngrots, mm.ngrrot))

    # If iavrr or iadrr are zero it means no dihedral angle was defined. ...
    # From 5 Sept 2012 we will always generate a random rotational group about a bond unless
    # the dihedral angle has been fixed
    # From 20/01/2014, a random rotational group is generated if the dihedral angle has not
    # been fixed, Otherwise a dihedral group will be generated.

    if c.nato[imv] > 2:
        c.ngrot[imv] = c.ngrot[imv] + 1
        mm.ngrrot = mm.ngrrot + 1
        molecom_mod_routines.testsize_rotgrps()

        # First find out the last group corresponding to imv

        iglast = 0
        im = 0
        while im < imv:

            im = im + 1
            iglast = iglast + int(c.ngrot[im])

        # Now shift all the existing groups above this up by 1

        ig = mm.ngrrot
        while ig > iglast:

            ig1 = ig - 1
            mm.rotgrp[ig] = mm.rotgrp[ig1]
            for j in range(1, 3):
                mm.labaxi[ig, j] = mm.labaxi[ig1, j]
            m.ndihedralang[ig] = m.ndihedralang[ig1]
            for k in range(1, m.ndihedralang[ig] + 1):
                m.dihedralang[k, ig] = m.dihedralang[k, ig1]
                for j in range(1, 3):
                    m.labrotdih[j, k, ig] = m.labrotdih[j, k, ig1]
            ig = ig1

        # Now insert the new rotational group

        mm.rotgrp[ig] = "ROT"
        mm.labaxi[ig, 1] = iadr
        mm.labaxi[ig, 2] = iav
        m.ndihedralang[ig] = 0

        # Ensure a finite step size is given

        if m.stepri <= 0.0:
            m.stepri = f32(1.0)

    # If both iadrr and iavrr are .gt. 0 then we add in a bond between these two atoms.
    # However after rotation it is possible these and other atoms now overlap, so we need to
    # eliminate duplicate atoms in this case. If iadrr and iavrr are not included in the list
    # of overlaps, then 20/01/2014 we set up a dihedral group between them. In all cases we
    # always lose iad and iavr and keep iav and iadr

    noverlap = 2
    ikeep[1] = iav
    ilose[1] = iad
    ikeep[2] = iadr
    ilose[2] = iavr

    ioverlaprr = 0
    if iadrr > 0 and iavrr > 0.0:

        # First determine if any atoms on imd overlaps any atoms on molecule imv - if so the
        # atom on imd will be amalgamated with the corresponding atom on imv

        iadt = int(c.naf[imd]) - 1
        while iadt < c.nal[imd]:

            iadt = iadt + 1

            # Do not include atoms iad and iadr in the test

            if iadt != iad and iadt != iadr:

                x1 = c.xadd[imd] + m.xa[iadt]
                y1 = c.yadd[imd] + m.ya[iadt]
                z1 = c.zadd[imd] + m.za[iadt]

                iavt = int(c.naf[imv]) - 1
                while iavt < c.nal[imv]:

                    iavt = iavt + 1
                    if iavt != iav and iavt != iavr:

                        x2 = c.xadd[imv] + m.xa[iavt]
                        y2 = c.yadd[imv] + m.ya[iavt]
                        z2 = c.zadd[imv] + m.za[iavt]
                        xdiff = (x2 - x1) / c.dmax
                        ydiff = (y2 - y1) / c.dmax
                        zdiff = (z2 - z1) / c.dmax
                        xdiff, ydiff, zdiff = minimagec.minimagec(xdiff, ydiff, zdiff)
                        dist = sqrtf(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff) * c.dmax
                        dmin = m.qradius[m.itype[iadt]] + m.qradius[m.itype[iavt]]
                        if dist <= dmin:

                            noverlap = noverlap + 1
                            ikeep[noverlap] = iavt
                            ilose[noverlap] = iadt

                            # If this pair includes iadrr and iavrr, then record the fact

                            if iadt == iadrr and iavt == iavrr:
                                ioverlaprr = 1

    # Shift all the coordinates on imd so that they are with reference to imv as the centre.

    for ia in range(c.naf[imd], c.nal[imd] + 1):

        m.xa[ia] = m.xa[ia] + c.xadd[imd] - c.xadd[imv]
        m.ya[ia] = m.ya[ia] + c.yadd[imd] - c.yadd[imv]
        m.za[ia] = m.za[ia] + c.zadd[imd] - c.zadd[imv]
        m.mref[ia] = imv

    # Change references to atoms on imd as required by overlaps

    for ioverlap in range(1, noverlap + 1):

        replace_references(int(ilose[ioverlap]), int(ikeep[ioverlap]))
        if iavrr == ilose[ioverlap]:
            iavrr = int(ikeep[ioverlap])
        if iadrr == ilose[ioverlap]:
            iadrr = int(ikeep[ioverlap])

    # Transfer neighbours from imd to imv

    for ioverlap in range(1, noverlap + 1):

        add_neighbours(int(ilose[ioverlap]), int(ikeep[ioverlap]))

    # If atoms iavrr and iadrr do not overlap determine the distance between these two atoms,
    # and set up a new bond between them, provided this pair does not already exist.
    # From 5 Sept 2012 we only create this distance if a particular
    # dihedral angle is specified.
    # From 20/01/2014 we create a dihedral group instead of a distance.

    if ioverlaprr == 0 and fixeddihedral:

        # First check that this pair does not already appear in each others neighbour lists

        ib = 0
        ifound = 0
        while ib < m.nbond[iavrr] and ifound == 0:

            ib = ib + 1
            if m.iabond[iavrr, ib] == iadrr:
                ifound = iadrr

        ib = 0
        while ib < m.nbond[iadrr] and ifound == 0:

            ib = ib + 1
            if m.iabond[iadrr, ib] == iavrr:
                ifound = iavrr

        if ifound == 0:

            if fixeddihedral:
                mm.rotgrp[ig] = "DIH"
                m.ndihedralang[ig] = 1
                m.dihedralang[1, ig] = dihedralangfixed
                m.labrotdih[1, 1, ig] = iadrr
                m.labrotdih[2, 1, ig] = iavrr

            # (the old code that created a bond between iavrr and iadrr is commented out in
            # the Fortran)

    # Set up the new number of atoms and mass for molecule imv

    c.nato[imv] = c.nato[imv] + c.nato[imd]
    c.nal[imv] = c.nal[imv] + c.nato[imd]
    c.amasst[imv] = c.amasst[imv] + c.amasst[imd]

    # This process should have removed all references to the overlapping atoms so we can
    # remove them from the list of atoms.

    for ioverlap in range(1, noverlap + 1):
        remove_atoms(int(ilose[ioverlap]), imv)
        # Make sure we update the list references to atoms lost and kept
        io = 0
        while io < noverlap:
            io = io + 1
            if ikeep[io] > ilose[ioverlap]:
                ikeep[io] = ikeep[io] - 1
            if ilose[io] > ilose[ioverlap]:
                ilose[io] = ilose[io] - 1

    # Need to step through the rotation groups and make sure the rotation axes are not
    # duplicated

    igb = 0
    while igb < mm.ngrrot:
        igb = igb + 1
        ig = igb
        while ig < mm.ngrrot:
            ig = ig + 1
            if ((mm.labaxi[ig, 1] == mm.labaxi[igb, 1] and mm.labaxi[ig, 2] == mm.labaxi[igb, 2])
                    or (mm.labaxi[ig, 2] == mm.labaxi[igb, 1] and mm.labaxi[ig, 1] == mm.labaxi[igb, 2])):
                # Found two axes the same - eliminate the second.
                mm.ngrrot = mm.ngrrot - 1
                igt = ig - 1
                while igt < mm.ngrrot:

                    igt = igt + 1
                    igt1 = igt + 1
                    mm.rotgrp[igt] = mm.rotgrp[igt1]
                    for j in range(1, 3):
                        mm.labaxi[igt, j] = mm.labaxi[igt1, j]
                    m.ndihedralang[igt] = m.ndihedralang[igt1]
                    if m.ndihedralang[igt] > 0:
                        for k in range(1, m.ndihedralang[igt] + 1):
                            m.dihedralang[k, igt] = m.dihedralang[k, igt1]
                            for j in range(1, 3):
                                m.labrotdih[j, k, igt] = m.labrotdih[j, k, igt1]
                # Decrement ig by 1, so that it looks at the new rotational group in this place
                ig = ig - 1

    # Will need to recalculate the number of rotational groups for each molecule

    for im in range(1, c.nmol + 1):
        c.ngrot[im] = 0
    ig = 0
    while ig < mm.ngrrot:
        ig = ig + 1
        im1 = int(m.mref[mm.labaxi[ig, 1]])
        im2 = int(m.mref[mm.labaxi[ig, 2]])
        units.write(6, lw(ig, mm.labaxi[ig, 1], im1, mm.labaxi[ig, 2], im2))
        if im1 == im2:
            c.ngrot[im1] = c.ngrot[im1] + 1
            mm.mrot[ig] = im1

    # Check the existing types to see if any of previous types do not now exist in the .ato
    # file.

    itb = 0

    while itb < m.ntype:

        itb = itb + 1
        ia = 0
        itfound = 0
        while ia < m.natom and itfound == 0:
            ia = ia + 1
            if m.itype[ia] == itb:
                itfound = itb

        # Remove this type if it is not found

        if itfound == 0:
            m.ntype = m.ntype - 1
            it = itb - 1
            while it < m.ntype:
                it = it + 1
                it1 = it + 1
                m.type[it] = m.type[it1]
                m.atsymbol[it] = m.atsymbol[it1]
                m.iexchange[it] = m.iexchange[it1]
                m.epscor[it] = m.epscor[it1]
                m.sigcor[it] = m.sigcor[it1]
                m.amassa[it] = m.amassa[it1]
                m.charge[it] = m.charge[it1]
                m.qradius[it] = m.qradius[it1]

            # Now shift the type numbers of the atoms if necessary

            ia = 0
            while ia < m.natom:
                ia = ia + 1
                if m.itype[ia] > itb:
                    m.itype[ia] = m.itype[ia] - 1

    # Shift all internal atomic coordinates to be relative to centre of mass

    ato_manipulation_rev4.centreofmass(imv)

    # Finally drop the number of molecules by 1.

    c.nmol = c.nmol - 1

    units.write(6, lw("bond_molecule> 5 ", c.nmol, m.natom, m.natomr, m.ntype, mm.ngrrot))

    return iavrr, iadrr


@routine
def replace_references(iad, iav):
    """dockato.f:replace_references -- replaces references to iad with references to iav."""
    # use reallocation_routines
    # use molecom_mod
    m = molecom
    mm = molecom_mod

    imd = int(m.mref[iad])
    imv = int(m.mref[iav])

    # Need to change the type number as appropriate

    m.itype[iad] = m.itype[iav]

    # Substitute iav for iad in the bond reference numbers

    for ia in range(1, m.natom + 1):
        ib = 0
        while ib < m.nbond[ia]:
            ib = ib + 1
            if m.iabond[ia, ib] == iad:
                m.iabond[ia, ib] = iav

    # Now the rotational groups

    ig = 0
    while ig < mm.ngrrot:
        ig = ig + 1
        if mm.mrot[ig] == imd:
            mm.mrot[ig] = imv
        for j in range(1, 3):
            if mm.labaxi[ig, j] == iad:
                mm.labaxi[ig, j] = iav
        if m.ndihedralang[ig] > 0:
            for k in range(1, m.ndihedralang[ig] + 1):
                for j in range(1, 3):
                    if m.labrotdih[j, k, ig] == iad:
                        m.labrotdih[j, k, ig] = iav

    return


@routine
def add_neighbours(iad, iav):
    """dockato.f:add_neighbours -- adds the neighbours of iad to iav.

    Adds the neighbours of iad to iav, but only if this neighbour does not already appear in
    the list of neighbours for this atom.
    """
    # use reallocation_routines
    # use molecom_mod
    m = molecom

    # The neighbours of iad will become the new neighbours of iav

    ib = 0
    while ib < m.nbond[iad]:
        ib = ib + 1
        iaref = int(m.iabond[iad, ib])

        # First check that this atom is not already listed as a neighbour to iav: only
        # add it to the list if it is not already present.

        ibc = 0
        ifound = 0
        while ibc < m.nbond[iav] and ifound == 0:
            ibc = ibc + 1
            if m.iabond[iav, ibc] == iaref:
                ifound = 1
        if ifound == 0:
            m.nbond[iav] = m.nbond[iav] + 1
            m.iabond[iav, m.nbond[iav]] = iaref
            m.tbond[iav, m.nbond[iav]] = m.tbond[iad, ib]
        ibc = 0
        ifound = 0
        while ibc < m.nbond[iaref] and ifound == 0:
            ibc = ibc + 1
            if m.iabond[iaref, ibc] == iav:
                ifound = 1
        if ifound == 0:
            m.nbond[iaref] = m.nbond[iaref] + 1
            m.iabond[iaref, m.nbond[iaref]] = iav
            m.tbond[iaref, m.nbond[iaref]] = m.tbond[iad, ib]

    return


@routine
def remove_atoms(iad, imv):
    """dockato.f:remove_atoms -- removes atom iad which is no longer in use.

    imv is the molecule number on which iad occurs
    """
    # use reallocation_routines
    # use molecom_mod
    m = molecom
    c = clustercom
    mm = molecom_mod

    # integer :: typesave
    # logical :: ifoundtype

    # New numbers of atoms in molecule imv

    c.nato[imv] = c.nato[imv] - 1
    c.nal[imv] = c.nal[imv] - 1
    if m.qtype[iad] != 1:
        m.natomr = m.natomr - 1

    # Shift all the atoms number references above iad down by 1

    ia = iad
    while ia < m.natom:
        ia1 = ia + 1
        m.xa[ia] = m.xa[ia1]
        m.ya[ia] = m.ya[ia1]
        m.za[ia] = m.za[ia1]
        m.alabel[ia] = m.alabel[ia1]
        m.itype[ia] = m.itype[ia1]
        m.qtype[ia] = m.qtype[ia1]
        m.nbond[ia] = m.nbond[ia1]
        m.mref[ia] = m.mref[ia1]
        ib = 0
        while ib < m.nbond[ia]:
            ib = ib + 1
            m.iabond[ia, ib] = m.iabond[ia1, ib]
            m.tbond[ia, ib] = m.tbond[ia1, ib]
        ia = ia1

    # Now shift all bond reference numbers down by 1 if greater than iad

    for ia in range(1, m.natom + 1):
        ib = 0
        while ib < m.nbond[ia]:
            ib = ib + 1
            if m.iabond[ia, ib] > iad:
                m.iabond[ia, ib] = m.iabond[ia, ib] - 1

    # Repeat the last step for the rotational groups.

    ig = 0
    while ig < mm.ngrrot:
        ig = ig + 1
        for j in range(1, 3):
            if mm.labaxi[ig, j] > iad:
                mm.labaxi[ig, j] = mm.labaxi[ig, j] - 1
        if m.ndihedralang[ig] > 0:
            for k in range(1, m.ndihedralang[ig] + 1):
                for j in range(1, 3):
                    if m.labrotdih[j, k, ig] > iad:
                        m.labrotdih[j, k, ig] = m.labrotdih[j, k, ig] - 1

    # New numbers of atom in the box

    m.natom = m.natom - 1

    clean_atom_types()  # Clean up the list of atom types to correspond with what is actually in the box

    return


@routine
def clean_atom_types():
    """dockato.f:clean_atom_types -- remove atom types no longer present in the box.

    Removes atom type iatomtype which is no longer in use. It is assumed this has already
    been checked.
    """
    # use reallocation_routines
    # use molecom_mod
    m = molecom
    mm = molecom_mod

    # integer iatomtype,it,it1
    # logical ifoundatomtype

    iatomtype = 0
    while iatomtype < m.ntype:
        iatomtype = iatomtype + 1
        ia = 1
        ifoundatomtype = m.type[iatomtype] == m.alabel[ia]
        while ia < m.natom and not ifoundatomtype:
            ia = ia + 1
            ifoundatomtype = m.type[iatomtype] == m.alabel[ia]
        if not ifoundatomtype:
            units.write(6, lw("clean_atom_type> ", iatomtype, m.ntype,
                              *[v for it in range(1, m.ntype + 1) for v in (" ", Ch(m.type[it], 3))]))
            if iatomtype < m.ntype:
                for it in range(iatomtype, m.ntype):
                    it1 = it + 1
                    m.type[it] = m.type[it1]
                    m.qtypea[it] = m.qtypea[it1]
                    m.iatomicnotype[it] = m.iatomicnotype[it1]
                    m.atsymbol[it] = m.atsymbol[it1]
                    m.iexchange[it] = m.iexchange[it1]
                    m.epscor[it] = m.epscor[it1]
                    m.sigcor[it] = m.sigcor[it1]
                    m.amassa[it] = m.amassa[it1]
                    m.charge[it] = m.charge[it1]
                    m.qradius[it] = m.qradius[it1]
                    mm.natomsphere[it] = mm.natomsphere[it1]
                    mm.rhosphere[it] = mm.rhosphere[it1]
                    mm.radminsphere[it] = mm.radminsphere[it1]
                    mm.radmaxsphere[it] = mm.radmaxsphere[it1]
                for ia in range(1, m.natom + 1):
                    if m.itype[ia] > iatomtype:
                        m.itype[ia] = m.itype[ia] - 1
                m.ntype = m.ntype - 1
            else:
                m.ntype = m.ntype - 1
    return
