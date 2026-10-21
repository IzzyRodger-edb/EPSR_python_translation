"""Translation of fortran_src/makeato.f: ``program standalone`` and ``subroutine makeato``
-- makes a .ato file containing one molecule of one atom.

NB: the Fortran never calls ``ctype``, so the coarse-graining arrays ``natomsphere`` and
``pseudofacmol`` (module molecom_mod) are unallocated when ``moloutc`` reads them: the
reference program dies with SIGSEGV inside ``moloutc`` and, the output being buffered,
leaves an empty .ato file (docs/rust_port/TOOLS.md, EPSR-test-ex30).  See ``_sigsegv``.
"""

import sys

from ..common import clustercom, commands, directory_name, epsr_dimension, lattice, molecom, molecom_mod, symbols
from ..fortran import Ch, FortranStop, f32, findex, fmt, fstr, int_array, powf, r4, sqrtf, substr, units
from ..library import ato_read_and_save
from ..library import auxiliary_routines
from ..library import file_and_directory_commands
from ..library import parse
from ..trace import routine

# FORMAT statements ('$' dropped: written with advance=False)
_F100 = "(/'makeato> Type name of .ato file to create: ')"
_F101 = "(a)"
_F499 = "(/'makeato> No .ato file specified! - routine will exit'/)"
_F102 = "(/'makeato> MAKEATO will create a .ato file containing 1 molecule',/' which itself contains just 1 atom'/)"
_F103 = "(/'makeato> Give number of different atom types in molecule: ')"
_F1041 = "(/'makeato> Atom type ',i2,' has label: ',a3)"
_F104 = "(/'makeato> For atom type ',i2,' give it a label (up to 3 chars.): ')"
_F121 = "(a)"
_F106 = "(/'makeato> For atom type ',a3,':-')"
_F107 = "(/' Give its epsilon [kJ/mol]: ')"
_F108 = "(/' Give its sigma [A]: ')"
_F109 = "(/' Give its atomic mass [amu]: ')"
_F110 = "(/' Give its electronic charge [e] ')"
_F111 = "(/' Give its chemical symbol: ')"
_F112 = ("(/'makeato> The following atoms have been defined:-'//' At. no., type no., label, chemical symbol,'"
         ",' (x,y,z) coordinates'/)")
_F127 = "(2(1x,i3),2(1x,a3),3(1x,f10.5))"
_F115 = "(/'makeato> Give two atom numbers',' and the distance between them [A].'//'(Type 0 0 0 when finished): ')"
_F117 = ("(/'makeato> Type three atom numbers and angle between them'//' (Middle atom is taken as centre).'"
         "/' (Type 0 0 0 0 when finished): ')")
_F118 = "(/'makeato> Not enough bonds defined for atom ',i3)"
_F119 = "(/'makeato> Suitable bonds between',3(1x,i3),' not found')"
_F124 = "(/'makeato> Type the temperature [K]: ')"
_F125 = "(/'makeato> and atomic number density [per A**3]: ')"


def main():
    """makeato.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    makeato()


@routine
def makeato():
    """makeato.f:makeato -- makes a molecule or headgroup."""
    # use molecom_mod
    mcomp = epsr_dimension.mcomp
    ncf = commands.ncf
    ncl = commands.ncl
    # character*256 readfile ; character*4 ext ; integer*4 ntypes(mcomp),ilenext
    ntypes = int_array(mcomp)
    # external cosd,sind  -> EPSR's own (auxiliary_routines)
    #
    # generate the name of the .ato file
    #
    ext = ".ato"  # noqa: F841
    ilenext = 4  # noqa: F841
    if commands.nwords < 2:
        commands.nwords = 0
        while commands.nwords <= 0:
            units.write(6, fmt(_F100), advance=False)
            commands.sometext = fstr(units.read_line(5), commands.nstsize)
            commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 5, commands.nstsize,
                                          commands.nworddim)
            if commands.nwords <= 0:
                units.write(6, fmt(_F499))
                return
        readfile = fstr(substr(commands.sometext, ncf[1], ncl[1]), 256)
        ilenreadfile = ncl[1] - ncf[1] + 1
    else:
        readfile = fstr(substr(commands.sometext, ncf[2], ncl[2]), 256)
        ilenreadfile = ncl[2] - ncf[2] + 1
    #
    # strip out any extension that was supplied - extension must be .ato, otherwise
    # nothing is done.
    #
    ilenuptodot = findex(readfile, ".") - 1
    if ilenuptodot > 0:
        readfile = fstr(substr(readfile, 1, ilenuptodot), 256)
        ilenreadfile = ilenuptodot
    symbols.fnameato = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                            + substr(readfile, 1, ilenreadfile) + ".ato", 256)
    symbols.ilenfnameato = directory_name.nchardir + ilenreadfile + 4
    pi = f32(3.141592653)  # noqa: F841
    units.write(6, fmt(_F102))
    molecom.natom = 0
    clustercom.nmol = 1
    molecom.ntype = 0
    natom1 = 1  # noqa: F841
    ntype1 = 1  # noqa: F841
    # !      write(6,103)
    # !      read(5,*) ntype
    molecom.ntype = 1
    molecom.type[1] = fstr(substr(readfile, 1, 3), 3)
    for i in range(1, molecom.ntype + 1):
        if i == 1:
            units.write(6, fmt(_F1041, i, Ch(molecom.type[i], 3)))
        else:
            units.write(6, fmt(_F104, i), advance=False)
            molecom.type[i] = fstr(units.read_fmt(5, _F121)[0], 3)

        # !            write(6,105)
        # !            read(5,*) ntypes(i)
        ntypes[i] = 1
        units.write(6, fmt(_F106, Ch(molecom.type[i], 3)))
        units.write(6, fmt(_F107), advance=False)
        molecom.epscor[i] = units.list_read(5).real(molecom.epscor[i])
        units.write(6, fmt(_F108), advance=False)
        molecom.sigcor[i] = units.list_read(5).real(molecom.sigcor[i])
        units.write(6, fmt(_F109), advance=False)
        molecom.amassa[i] = units.list_read(5).real(molecom.amassa[i])
        units.write(6, fmt(_F110), advance=False)
        molecom.charge[i] = units.list_read(5).real(molecom.charge[i])
        units.write(6, fmt(_F111), advance=False)
        molecom.atsymbol[i] = fstr(units.list_read(5).word(molecom.atsymbol[i]), 3)
        molecom.qradius[i] = 0.0
        molecom.iexchange[i] = 0
        for j in range(1, ntypes[i] + 1):
            molecom.natom = molecom.natom + 1
            natom = molecom.natom
            molecom.alabel[natom] = molecom.type[i]
            molecom.xa[natom] = 0.0
            molecom.ya[natom] = 0.0
            molecom.za[natom] = 0.0
            molecom.alabel[natom] = molecom.type[i]
            molecom.itype[natom] = i
            molecom.nbond[natom] = 0
    units.write(6, fmt(_F112))
    for i in range(1, molecom.natom + 1):
        units.write(6, fmt(_F127, i, molecom.itype[i], Ch(molecom.alabel[i], 3), Ch(molecom.atsymbol[i], 3),
                           molecom.xa[i], molecom.ya[i], molecom.za[i]))
    #
    # set up bond distances for pairs of atom types. If no bond distance is defined
    # for a particular pair of atoms, then no constraint will be defined for
    # this bond
    #
    # !      write(6,113)
    nbond = molecom.nbond
    iabond = molecom.iabond
    tbond = molecom.tbond
    bangle = f32(0.0)  # NB: uninitialised in the Fortran; only read inside the (never entered) angle loop
    while True:
        # 502 continue
        # !      write(6,114)
        # !      read(5,*) ia1,ia2,bond
        ia1 = 0
        ia2 = 0
        bond = f32(0)
        while ia1 > 0 and ia2 > 0:
            nbond[ia1] = nbond[ia1] + 1
            nbond[ia2] = nbond[ia2] + 1
            iabond[ia1, nbond[ia1]] = ia2
            iabond[ia2, nbond[ia2]] = ia1
            tbond[ia1, nbond[ia1]] = bond
            tbond[ia2, nbond[ia2]] = bond
            units.write(6, fmt(_F115), advance=False)
            r = units.list_read(5)
            ia1 = r.int(ia1)
            ia2 = r.int(ia2)
            bond = r.real(bond)
        #
        # write out the bonds defined for each atom
        #
        # !     do ia=1,natom
        # !            write(6,116) alabel(ia),nbond(ia)
        # !     1,(alabel(iabond(ia,i)),iabond(ia,i),tbond(ia,i),i=1,nbond(ia))
        # !      end do
        #
        # set up bond angles for triplets of atom types. If no bond angle is defined
        # for a particular triplet of atom types, then no constraint will be defined for
        # this triplet.
        #
        # !      write(6,117)
        # !      read(5,*) ia1,ia2,ia3,bangle
        ia1 = 0
        ia2 = 0
        ia3 = 0
        goto502 = False
        while ia1 > 0 and ia2 > 0 and ia3 > 0:
            #
            # first check that atom ia2 has at least 2 bonds defined
            #
            if nbond[ia2] < 2:
                units.write(6, fmt(_F118, ia2))
                goto502 = True
                break
            #
            # Now check that both ia1 and ia3 have been defined as bonded to ia2
            #
            nfound = 0
            bond21 = f32(0.0)
            bond23 = f32(0.0)
            for ia in range(1, nbond[ia2] + 1):
                if iabond[ia2, ia] == ia1:
                    nfound = nfound + 1
                    bond21 = tbond[ia2, ia]
                elif iabond[ia2, ia] == ia3:
                    nfound = nfound + 1
                    bond23 = tbond[ia2, ia]
            #
            # check that 2 bonds were found (if more than 2 were found it will take
            # the last ones found
            #
            if nfound < 2:
                units.write(6, fmt(_F119, ia1, ia2, ia3))
                goto502 = True
                break
            #
            # O.K. so set up distance between ia1 and ia3
            #
            b13 = bond21 * bond21 + bond23 * bond23 - f32(2) * bond21 * bond23 * auxiliary_routines.cosd(bangle)
            b13 = sqrtf(b13)
            nbond[ia1] = nbond[ia1] + 1
            nbond[ia3] = nbond[ia3] + 1
            iabond[ia1, nbond[ia1]] = ia3
            iabond[ia3, nbond[ia3]] = ia1
            tbond[ia1, nbond[ia1]] = b13
            tbond[ia3, nbond[ia3]] = b13
            units.write(6, fmt(_F117), advance=False)
            r = units.list_read(5)
            ia1 = r.int(ia1)
            ia2 = r.int(ia2)
            ia3 = r.int(ia3)
            bangle = r.real(bangle)
        if not goto502:
            break
    #
    # write out the bonds defined for each atom
    #
    # !      (commented out in the Fortran)
    #
    # determine if there are any rotational groups
    #
    ig = 0
    molecom_mod.ngrrot = 0
    molecom_mod.ngrrotmax = 0
    molecom_mod.nrottot = 0
    molecom_mod.nrottotmax = 0
    # ! NOTE: If it is ever reinstated, the following code will need to be adjusted to take account
    # ! of the revised method of indexing the rotational groups
    # ! (the rotational-group input of the Fortran is commented out)
    clustercom.ngrot[1] = ig
    #
    # set some parameters ready for writing the .ato file
    #
    clustercom.nmol = 1
    clustercom.xadd[1] = 0.0
    clustercom.yadd[1] = 0.0
    clustercom.zadd[1] = 0.0
    clustercom.phix[1] = 0.0
    clustercom.phiy[1] = 0.0
    clustercom.phiz[1] = 0.0
    clustercom.nato[1] = molecom.natom
    clustercom.naf[1] = 1
    clustercom.nal[1] = molecom.natom
    molecom.stepmi = f32(0.1)
    if clustercom.ngrot[1] > 0:
        molecom.stepri = f32(1.0)
    else:
        molecom.stepri = f32(0.0)
    molecom.stepr = f32(1.0)
    molecom.stepm = f32(2.0)
    molecom.tempb = f32(65)
    units.write(6, fmt(_F124), advance=False)
    molecom.simtk = units.list_read(5).real(molecom.simtk)
    units.write(6, fmt(_F125), advance=False)
    # rho is an implicitly typed real local here (epsr16.inc is not included)
    rho = units.list_read(5).real(f32(0.0))
    #     write(6,128,advance='no') / read(5,*) ecore     (commented out)
    #     write(6,129,advance='no') / read(5,*) dcore     (commented out)
    # !      write(6,130,advance='no') / read(5,*) tempb   (commented out)
    molecom.tempb = f32(65)
    molecom.ecore = f32(0.0)
    molecom.dcore = f32(0.0)
    #
    # convert density to box size
    #
    clustercom.volume = f32(molecom.natom) / rho
    clustercom.dmax = powf(clustercom.volume, r4(1.0) / r4(3.0))
    lattice.cubic = True
    #
    # setup molecular size factor
    #
    molecom.dmaxmol = clustercom.dmax
    # !Set up template file
    clustercom.nmoltype = 1
    clustercom.moltemplate[clustercom.nmoltype] = fstr(substr(readfile, 1, ilenreadfile), 256)
    clustercom.wholerot[clustercom.nmoltype] = 0.0
    clustercom.wholetrans[clustercom.nmoltype] = 0.01
    try:
        ato_read_and_save.moloutc(symbols.fnameato)
    except TypeError:
        if molecom_mod.natomsphere is not None and molecom_mod.pseudofacmol is not None:
            raise
        _sigsegv()
    return


def _sigsegv():
    """The reference program dies here (SIGSEGV in moloutc: natomsphere/pseudofacmol are
    unallocated, makeato never calls ctype).  Reproduce the outcome: the buffered output
    of the open units is lost (the .ato file stays empty, as moloutc's open left it) and
    the program ends with the signal's exit status."""
    units.reset()
    sys.stderr.write("\nProgram received signal SIGSEGV: Segmentation fault - invalid memory reference.\n")
    raise FortranStop(139, "makeato: moloutc read the unallocated natomsphere/pseudofacmol (Fortran SIGSEGV)")
