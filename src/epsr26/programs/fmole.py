"""Translation of fortran_src/fmole.f: ``program standalone``, ``subroutine fmole``,
``initpt_fmole`` and ``update_ato`` -- a molecular "shake" of all the molecules in a
cluster (intramolecular Monte Carlo only), repeated the specified number of times.

Unlike the Rust port, the ``makemole`` call fmole makes for every molecule template
(regenerating ``<template>.ato`` from ``<template>.mol`` when that file is present) is
translated: it calls library/ato_manipulation_rev4.makemole exactly as the Fortran does.
"""

from ..common import clustercom, commands, epsr16, epsr_dimension, molecom, molecom_mod, potential_par, symbols
from ..common import coordconstraints_mod
from ..fortran import (
    Ch,
    FortranEOF,
    FortranIOError,
    char_array,
    expf,
    f32,
    findex,
    findex_padded,
    fmt,
    fstr,
    iostat_of,
    lw,
    read_internal,
    substr,
    units,
)
from ..library import ato_manipulation_rev4
from ..library import ato_read_and_save
from ..library import file_and_directory_commands
from ..library import minimagec
from ..library import molecom_mod as molecom_mod_lib
from ..library import molfit
from ..library import parse
from ..library import reallocation_routines
from ..trace import routine

# FORMAT statements of fmole ('$' dropped: written with advance=False)
_F499 = "(/'fcluster> No .ato file specified! - routine will exit'/)"
_F497 = "(/'fmole> Type the number of times to perform the shake: ')"
_F599 = "(/'fmole> Problems with .ato file ',a/'Routine will exit now'/)"
_F595 = "(/'fmole> ',a/)"
_F598 = "(/'fmole> Step size is set to zero!'/)"
_F596 = "(/'fmole> Starting intramolecular energy: ',(1x,e12.5)/)"
_F597 = ("(/'fmole> Iteration ',i0,' Intramolecular energy: ',(1x,e12.5)"
         "/1x,'No. of moves tried: ',i0,' No. of moves rejected: ',i0)")
# FORMAT statements of update_ato
_UA_F597 = "(/'update_ato> ',a)"
_UA_F107 = "(a)"
_UA_F100 = "(1x,a3,1x,a3)"
_UA_F596 = "(/'update_ato> Could not open template file :',a)"


def main():
    """fmole.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    fmole()


@routine
def fmole():
    """fmole.f:fmole -- perform a molecular "shake" on all the molecules in a cluster.
    This is repeated the specified number of times.

    2006-11-12 Made a small change to initpt_fmole routine, so that truncation function
    applied to potential, not gradient of potential.
    """
    # use molecom_mod ; use cellneighbourlist ; use coordconstraints_mod
    from ..library import cellneighbourlist

    mcomp = epsr_dimension.mcomp
    ncf = commands.ncf
    ncl = commands.ncl
    # character*256 readfile ; character*4 ext ; character*1 iaccept
    iaccept = ""  # NB: not initialised in the Fortran; search_for_files sets it
    # integer*4 nmoltypesave,ilenreadfile,num_threads,rotfreq
    # character*256 moltemplatesave(mcomp)
    moltemplatesave = char_array(mcomp)
    # real ftethered ; fcoordconstraints, ftest, ftestb (implicit real): set by fmolec1
    ftethered = f32(0.0)
    fcoordconstraints = f32(0.0)
    #
    # generate the name of the .ato file
    #
    ext = ".ato"
    ilenext = 4
    if commands.nwords < 2:
        commands.nwords = 0
        while commands.nwords <= 0:
            iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)
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
    # determine the number of times to do the shake
    #
    ntimes = 0
    ierr = 1
    if commands.nwords > 2:
        try:
            ntimes = read_internal(substr(commands.sometext, ncf[3], ncl[3])).int(ntimes)
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = iostat_of(e)
    while ierr != 0:
        units.write(6, fmt(_F497), advance=False)
        try:
            ntimes = units.list_read(5).int(ntimes)
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = iostat_of(e)
    #
    # determine the frequency for neighbour list update
    #
    nrenew = 0
    # !      ierr=1
    # !      if(nwords.gt.3) then
    # !         read(sometext(ncf(4):ncl(4)),*,iostat=ierr) nrenew
    # !      endif
    # !      do while(ierr.ne.0)
    # !         write(6,496)
    # !496   format(/'fmole> Type the frequency to update the neighbour list (0): ',$)
    # !         read(5,*,iostat=ierr) nrenew
    # !      end do
    if nrenew < 1 or nrenew > ntimes:
        nrenew = ntimes
    #
    # strip out any extension that was supplied for the .ato file
    # - extension must be .ato, otherwise nothing is done.
    #
    ilendot = findex(readfile, ".") - 1
    if ilendot > 0:
        readfile = fstr(substr(readfile, 1, ilendot), 256)
        ilenreadfile = ilendot
    units.write(6, lw(substr(readfile, 1, ilenreadfile)))
    #
    # create the full name of the input file by adding on the directory name if not present
    #
    symbols.fnameato, symbols.ilenfnameato = file_and_directory_commands.get_full_filename(
        readfile, ilenreadfile, symbols.fnameato, symbols.ilenfnameato)
    symbols.fnameato = fstr(substr(symbols.fnameato, 1, symbols.ilenfnameato) + ".ato", 256)
    symbols.ilenfnameato = symbols.ilenfnameato + 4
    units.write(6, lw(substr(readfile, 1, ilenreadfile)))
    units.write(6, lw(substr(symbols.fnameato, 1, symbols.ilenfnameato)))

    # Save this for the second read

    readfile = symbols.fnameato
    ilenreadfile = symbols.ilenfnameato

    # initialise number of atom and molecular types to zero

    molecom.ntype = 0
    clustercom.nmoltype = 0
    coordconstraints_mod.ncoordconstraints = 0
    #
    # read .ato file
    #
    molecom.sizefactor = f32(1.0)
    (epsr16.rmaxcell, epsr16.cellmax, epsr16.upol, epsr16.rho, epsr16.amubar, epsr16.amubar2, epsr16.ndist,
     ierr) = ato_read_and_save.minitc(symbols.fnameato, epsr16.rmaxcell, epsr16.cellmax, epsr16.upol,
                                      epsr16.rho, epsr16.amubar, epsr16.amubar2, epsr16.ndist, ierr)
    if epsr16.rho == 0.0 or ierr != 0:
        units.write(6, fmt(_F599, Ch(symbols.fnameato, 256)))
        return

    # Save the names of the template files.

    nmoltypesave = clustercom.nmoltype
    for it in range(1, clustercom.nmoltype + 1):
        moltemplatesave[it] = clustercom.moltemplate[it]

    # For each template, run makemole to ensure the corresponding .ato has the latest
    # values.

    for it in range(1, nmoltypesave + 1):

        ilenmoltemplate = findex_padded(moltemplatesave[it], " ", 256) - 1
        if ilenmoltemplate > 0:

            commands.sometext = fstr("makemole " + substr(moltemplatesave[it], 1, ilenmoltemplate),
                                     commands.nstsize)
            commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 5, commands.nstsize,
                                          commands.nworddim)
            if commands.nwords > 1:
                units.write(6, fmt(_F595, substr(commands.sometext, ncf[1], ncl[2])))
                ato_manipulation_rev4.makemole()

    # Makemole may have overwritten some values of the original ato file so now need
    # to read the original ato file again.

    symbols.fnameato = readfile
    symbols.ilenfnameato = ilenreadfile
    molecom.sizefactor = f32(1.0)
    molecom.ntype = 0
    clustercom.nmoltype = 0
    (epsr16.rmaxcell, epsr16.cellmax, epsr16.upol, epsr16.rho, epsr16.amubar, epsr16.amubar2, epsr16.ndist,
     ierr) = ato_read_and_save.minitc(symbols.fnameato, epsr16.rmaxcell, epsr16.cellmax, epsr16.upol,
                                      epsr16.rho, epsr16.amubar, epsr16.amubar2, epsr16.ndist, ierr)

    # Update this ato file from the template files

    update_ato(ierr)

    ato_read_and_save.normalise_coordinates()
    #
    # initialize a repulsive exponential potential based on ecore and dcore
    # - this is to push non-bonded atoms as far apart as possible
    #
    rminpt = f32(4.0)
    rmaxpt = f32(6.0)
    initpt_fmole(rminpt, rmaxpt)
    #
    # if stepmi is set to zero, nothing will happen
    #
    if molecom.stepmi <= 0.0 and molecom.stepri <= 0.0:
        units.write(6, fmt(_F598))
        return
    #
    # save the step sizes and set the others to zero
    #
    stepms = molecom.stepm
    stepris = molecom.stepri
    steprs = molecom.stepr
    tempdihsave = molecom.tempdih
    # !      tempdih=0.0 !For intramolecular refinement we ignore any dihedral angles  - these need to be refined
    # ! including the inter-molecular interactions.
    molecom.stepm = f32(0.0)
    molecom.stepr = f32(0.0)
    if clustercom.nmol > 1:
        molecom.stepri = f32(0.0)  # Only if there is more than 1 molecule, since intermolecular forces are not defined in fmole.
    #
    # calculate the initial intramolecular energy
    #
    ftest = f32(0.0)
    ftestb = f32(0.0)

    # return in an error occurred in fmolec1

    #     if(ierr.ne.0) return
    #
    # force the molecular "shake"
    #
    intra = 1
    rotfreq = 1
    molecom.refine_intra = False
    molecom.coul_intra = False
    molecom.refine_intra_fac = f32(0.0)
    molecom.coul_intra_fac = f32(0.0)
    # !Set the requested rejection rate to 0.8
    epsr16.rejrate = f32(0.8)
    nisum = intra
    ireset = 0
    in_ = nrenew
    it = 0

    # ! set the maximum separation of atoms to build neighbour list. In this case
    # ! no intermolecular interactions are included

    sepmax = f32(0.0)
    nmolcellin = 10

    # ! Initialise cell neighbour list

    cellneighbourlist.get_date_and_time("Start initialise_cellneighbourlist#", 0)
    cellneighbourlist.initialise_cellneighbourlist(nmolcellin, clustercom.xadd, clustercom.yadd, clustercom.zadd,
                                                   clustercom.nmol, sepmax)
    cellneighbourlist.get_date_and_time("Finish initialise_cellneighbourlist#", 1)

    while it < ntimes:
        it = it + 1
        if in_ < nrenew:
            in_ = in_ + 1
        else:
            in_ = 1
            #
            # At beginning and every nrenew iterations we renew energy of distribution
            # and reset the neighbour list
            #
            num_threads = 0
            ato_manipulation_rev4.adjust_intra_wts()  # Adjust and set intramolecular energy weights
            ierr, num_threads, ftest, ftestb, ftethered, fcoordconstraints = ato_manipulation_rev4.fmolec1(
                ierr, num_threads, ftest, ftestb, ftethered, fcoordconstraints)
            if ierr != 0:
                return
        if it == 1:
            units.write(6, fmt(_F596, ftestb))
        #         write(6,596) ftestb

        # Ensure molecules stay within box.

        for im in range(1, clustercom.nmol + 1):
            clustercom.xadd[im], clustercom.yadd[im], clustercom.zadd[im] = minimagec.minimagec(
                clustercom.xadd[im], clustercom.yadd[im], clustercom.zadd[im])
            if clustercom.tethered[im]:
                clustercom.phix[im], clustercom.phiy[im], clustercom.phiz[im] = minimagec.minimagec(
                    clustercom.phix[im], clustercom.phiy[im], clustercom.phiz[im])
        # !         write(6,*) 'fmole> ',it
        ftest, ftestb, nisum = molfit.molfit(ftest, ftestb, nisum, intra, rotfreq, ireset)
        if 100 * (it // 100) == it or it == ntimes:
            for j in range(1, 3):
                clustercom.nmoveit[j] = 0
                for imoltype in range(1, clustercom.nmoltype + 1):
                    clustercom.nmoveit[j] = (clustercom.nmoveit[j] + clustercom.nmovei[j, imoltype]
                                             + clustercom.ninter[j, imoltype])
            units.write(6, fmt(_F597, it, ftestb, clustercom.nmoveit[1], clustercom.nmoveit[2]))
            ato_manipulation_rev4.adjust_intra_wts()  # Adjust and set intramolecular energy weights
    #
    # restore step sizes
    #
    molecom.stepm = stepms
    molecom.stepri = stepris
    molecom.stepr = steprs
    molecom.tempdih = tempdihsave
    #
    # write out the new molecule coordinates
    #
    ato_read_and_save.unnormalise_coordinates()
    ato_read_and_save.moloutc(symbols.fnameato)

    units.write(6, lw("Done fmole"))

    return


@routine
def initpt_fmole(rminpt, rmaxpt):
    """fmole.f:initpt_fmole -- r grid, truncation functions, zeroed empirical potentials
    and the repulsive molecular potential (generate_mp) based on ecore and dcore."""
    # use molecom_mod
    pp = potential_par
    pi = f32(3.141592653)  # noqa: F841
    #
    # maximum value of a sum - this it to prevent overflows
    #
    summax = expf(f32(25.0))  # noqa: F841
    #
    # generate x scae and renormalized x scale and potential damping function at low r
    #
    epsr16.npt = 1000
    epsr16.cellst = f32(0.5) * clustercom.dmax / f32(epsr16.npt)
    epsr16.npt2 = epsr16.npt - 1
    epsr16.nptmin = 1
    epsr16.nptmax = epsr16.npt2
    for j in range(1, epsr16.npt + 1):
        epsr16.x[j] = f32(j - 1) * epsr16.cellst
        if epsr16.x[j] <= rminpt:
            epsr16.nptmin = j
        if epsr16.x[j] <= rmaxpt:
            epsr16.nptmax = j
        epsr16.xstar[j] = epsr16.x[j] / clustercom.dmax
    for j in range(1, epsr16.npt2 + 1):
        xx = epsr16.x[j]
        pp.wfac[j] = ato_manipulation_rev4.truncate(xx, rminpt, rmaxpt)
        pp.wfacc[j] = ato_manipulation_rev4.truncatec(xx, rmaxpt)
    pp.wfac[epsr16.npt] = 0.0
    #
    # set the empirical potentials to zero
    #
    for ip in range(1, molecom.npairs + 1):
        for j in range(1, epsr16.npt + 1):
            pp.xav[j, ip] = 0.0
            pp.xavs[j, ip] = 0.0
            pp.xav1[j, ip] = 0.0
            pp.xavs1[j, ip] = 0.0
            pp.epot[j, ip] = 0.0
            pp.epot1[j, ip] = 0.0
            pp.ecoul[j, ip] = 0.0
            pp.ecoul1[j, ip] = 0.0
        pp.irefpcof[ip] = ip
        pp.rminex[ip] = 0.0
        pp.pcoefex[ip] = 0.0
        pp.rmaxex[ip] = 0.0
        pp.pcoefexmax[ip] = 0.0

    # ! Generate molecular potential (based on ecore and dcore)

    molecom.power = f32(12.0)
    ato_manipulation_rev4.generate_mp()

    return


def _read_items(text, targets):
    """List-directed internal read of len(targets) items (('int'|'real', array, index))
    with iostat=: items are stored one by one until an error or the end of the record,
    as gfortran does.  Returns the iostat value."""
    r = read_internal(text)
    try:
        for kind, arr, idx in targets:
            if kind == "int":
                arr[idx] = r.int(arr[idx])
            else:
                arr[idx] = r.real(arr[idx])
    except (FortranEOF, FortranIOError) as e:
        return iostat_of(e)
    return 0


@routine
def update_ato(ierrout):
    """fmole.f:update_ato -- updates the .ato file, using the latest template files.

    Note that if an attempt is made to change any of the atom types or the number of
    different types, nothing will be changed from the point at which the error was
    detected onwards for that molecule.
    """
    # use reallocation_routines ; use molecom_mod
    near = epsr_dimension.near
    ncf = commands.ncf
    ncl = commands.ncl
    # integer*4 ierrout ; character*20 ext ; integer*4 ilenext
    # integer*4 icurrenttemplate,imoltype,ntypetemplate,natomtemplate,ilenmoltemplate
    # character(len=3), allocatable :: alabeltemplate,typetemplate,rotgrptemplate
    alabeltemplate = typetemplate = rotgrptemplate = None
    # character(len=3) :: labeltest,symboltest
    labeltest = ""  # NB: uninitialised in the Fortran; set by the read below
    symboltest = ""
    # integer, allocatable :: nbondtemplate,nrottemplate ; iabondtemplate,labrottemplate ; labaxitemplate
    nbondtemplate = nrottemplate = None
    iabondtemplate = labrottemplate = None
    labaxitemplate = None
    # real, allocatable :: tbondtemplate,dihedralangs
    tbondtemplate = dihedralangs = None
    # integer, allocatable :: ndihedralangs ; labrotdihs
    ndihedralangs = None
    labrotdihs = None
    # NB: natomtemplate and nrotgrptemplate are uninitialised before the first template
    # is read (the first molecule always reads one: icurrenttemplate starts at 0)
    natomtemplate = 0
    nrotgrptemplate = 0
    ierr = 0

    # Label of the current template file

    icurrenttemplate = 0
    # !      if(.not.allocated(alabeltemplate)) allocate(alabeltemplate(1),typetemplate(1),rotgrptemplate(1))

    # Step through the molecules. Identify its type and read in the correct molecular template
    # if it is different from the previous molecule

    # Number the rotational groups

    molecom_mod.ngrrot = 0
    molecom_mod.nrottot = 0
    ig = 0

    # Extension of template file

    ext = ".ato"
    ilenext = 4

    for im in range(1, clustercom.nmol + 1):

        imoltype = int(clustercom.mtype[im])

        # Only try to read a template file if one is specified

        if imoltype != icurrenttemplate:

            # Set the initial number of atoms in the template file to zero.
            # This is to ensure nothing is done in the event that a template file does not exist.

            natomtemplate = 0

            # Assign the extension to the template file

            icurrenttemplate = imoltype
            ilenmoltemplate = len(clustercom.moltemplate[imoltype])

            # Check that a template file is specified

            if ilenmoltemplate > 0:

                # O.K so open it and get the molecular parameters

                symbols.fname_template = fstr(substr(clustercom.moltemplate[imoltype], 1, ilenmoltemplate)
                                              + substr(ext, 1, ilenext), 256)
                symbols.ilenfname_template = ilenmoltemplate + ilenext
                units.write(6, fmt(_UA_F597, substr(symbols.fname_template, 1, symbols.ilenfname_template)))

                # Convert this to the full filename

                symbols.fnametemplate, symbols.ilenfnametemplate = file_and_directory_commands.get_full_filename(
                    symbols.fname_template, symbols.ilenfname_template, symbols.fnametemplate,
                    symbols.ilenfnametemplate)

                # Open it and read the relevant atom types, bond lengths and rotational groups.

                try:
                    units.open(10, symbols.fnametemplate, status="old")
                    ierr = 0
                except (FortranEOF, FortranIOError) as e:
                    ierr = iostat_of(e)

                if ierr == 0:

                    # Decide whether this is a cubic or a non-cubic simulation box

                    try:
                        commands.sometext = fstr(units.read_line(10), commands.nstsize)
                        ierr = 0
                    except (FortranEOF, FortranIOError) as e:
                        ierr = iostat_of(e)
                    commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 5, commands.nstsize,
                                                  commands.nworddim)
                    if commands.nwords != 3:

                        # Not cubic, read 3 more lines

                        for _k in range(3):
                            try:
                                units.skip_record(10)
                                ierr = 0
                            except (FortranEOF, FortranIOError) as e:
                                ierr = iostat_of(e)

                    else:

                        # ! Otherwise read 1 more lines

                        try:
                            units.skip_record(10)
                            ierr = 0
                        except (FortranEOF, FortranIOError) as e:
                            ierr = iostat_of(e)

                    ntypetemplate = 0
                    nrotmax = 0
                    ndihmax = 0

                    # Number of atoms in the template molecule

                    natomtemplate = units.list_read(10).int(natomtemplate)
                    alabeltemplate = reallocation_routines.reallocate1d_c(alabeltemplate, 3, natomtemplate)
                    typetemplate = reallocation_routines.reallocate1d_c(typetemplate, 3, natomtemplate)
                    nbondtemplate = reallocation_routines.reallocate1d_i(nbondtemplate, natomtemplate)
                    iabondtemplate = reallocation_routines.reallocate2d_i(iabondtemplate, natomtemplate, near)
                    tbondtemplate = reallocation_routines.reallocate2d_r(tbondtemplate, natomtemplate, near)

                    # Run through the atoms in the template files and pick up the atom types and bond
                    # length parameters

                    for ia in range(1, natomtemplate + 1):

                        try:
                            alabeltemplate[ia] = units.read_fmt(10, _UA_F100)[0]
                            ierr = 0
                        except (FortranEOF, FortranIOError) as e:
                            ierr = iostat_of(e)
                        it = 0
                        ic = 0

                        # Save the number of atom types in this molecule. Reminder: it is not possible to change
                        # the atom types when updating the molecule template.

                        while ic < ntypetemplate and it == 0:
                            ic = ic + 1
                            if alabeltemplate[ia] == typetemplate[ic]:
                                it = ic

                        # If it.eq.0, then no atom type of this kind has been found so far in the template

                        if it == 0:
                            ntypetemplate = ntypetemplate + 1
                            typetemplate[ntypetemplate] = alabeltemplate[ia]
                        units.skip_record(10)
                        r = units.list_read(10)
                        nbondtemplate[ia] = r.int(nbondtemplate[ia])
                        for ja in range(1, nbondtemplate[ia] + 1):
                            iabondtemplate[ia, ja] = r.int(iabondtemplate[ia, ja])
                            tbondtemplate[ia, ja] = r.real(tbondtemplate[ia, ja])

                    # Read rotational groups

                    nrotgrptemplate = units.list_read(10).int(nrotgrptemplate)
                    if nrotgrptemplate > 0:
                        nrotmax = 0
                        ndihmax = 0
                        rotgrptemplate = reallocation_routines.reallocate1d_c(rotgrptemplate, 3, nrotgrptemplate)
                        labaxitemplate = reallocation_routines.reallocate2d_i(labaxitemplate, nrotgrptemplate, 2)
                        nrottemplate = reallocation_routines.reallocate1d_i(nrottemplate, nrotgrptemplate)
                        ndihedralangs = reallocation_routines.reallocate1d_i(ndihedralangs, nrotgrptemplate)
                        for ir in range(1, nrotgrptemplate + 1):
                            try:
                                rotgrptemplate[ir] = units.read_fmt(10, _UA_F100)[0]
                                ierr = 0
                            except (FortranEOF, FortranIOError) as e:
                                ierr = iostat_of(e)
                            r = units.list_read(10)
                            labaxitemplate[ir, 1] = r.int(labaxitemplate[ir, 1])
                            labaxitemplate[ir, 2] = r.int(labaxitemplate[ir, 2])
                            if rotgrptemplate[ir] == "DIH":
                                ndihedralangs[ir] = units.list_read(10).int(ndihedralangs[ir])
                                if ndihedralangs[ir] > ndihmax:
                                    ndihmax = int(ndihedralangs[ir])
                                    dihedralangs = reallocation_routines.reallocate2d_r(dihedralangs, ndihmax,
                                                                                        nrotgrptemplate)
                                    labrotdihs = reallocation_routines.reallocate3d_i(labrotdihs, 2, ndihmax,
                                                                                      nrotgrptemplate)
                                for k in range(1, ndihedralangs[ir] + 1):
                                    r = units.list_read(10)
                                    labrotdihs[1, k, ir] = r.int(labrotdihs[1, k, ir])
                                    labrotdihs[2, k, ir] = r.int(labrotdihs[2, k, ir])
                                    dihedralangs[k, ir] = r.real(dihedralangs[k, ir])
                            else:
                                ndihedralangs[ir] = 0
                            try:
                                commands.sometext = fstr(units.read_line(10), commands.nstsize)
                                ierr = 0
                            except (FortranEOF, FortranIOError) as e:
                                ierr = iostat_of(e)
                            commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 5,
                                                          commands.nstsize, commands.nworddim)
                            nrottemplate[ir] = read_internal(substr(commands.sometext, ncf[1], ncl[1])).int(
                                nrottemplate[ir])
                            if nrottemplate[ir] > nrotmax:
                                nrotmax = int(nrottemplate[ir])
                                labrottemplate = reallocation_routines.reallocate2d_i(labrottemplate, nrotmax,
                                                                                      nrotgrptemplate)
                            r = read_internal(substr(commands.sometext, ncf[2], ncl[commands.nwords]))
                            for jr in range(1, commands.nwords - 1 + 1):
                                labrottemplate[jr, ir] = r.int(labrottemplate[jr, ir])
                            # !If there are more values to read, then read them directly
                            if nrottemplate[ir] > commands.nwords - 1:
                                r = units.list_read(10)
                                for jr in range(commands.nwords, nrottemplate[ir] + 1):
                                    labrottemplate[jr, ir] = r.int(labrottemplate[jr, ir])

                    # Read the potential parameters in case they have changed. It will be
                    # assumed these occur in the same order as they are found in the molecule

                    for it in range(1, ntypetemplate + 1):

                        try:
                            vals = units.read_fmt(10, _UA_F100)
                            labeltest, symboltest = vals[0], vals[1]
                            ierr = 0
                        except (FortranEOF, FortranIOError) as e:
                            ierr = iostat_of(e)

                        # Find the type number corresponding to this label

                        ic = 0
                        is_ = 0
                        while ic < molecom.ntype and is_ == 0:

                            ic = ic + 1
                            if labeltest == molecom.type[ic]:
                                is_ = ic

                        # Read the parameters if a matching type has been found.

                        if is_ > 0:

                            molecom.atsymbol[is_] = symboltest
                            # !                        read(10,*) epscor(is),sigcor(is)
                            # !     *,amassa(is),charge(is),qradius(is),natomsphere(is),radminsphere(is),radmaxsphere(is),pseudofacmol(is)
                            try:
                                commands.sometext = fstr(units.read_line(10), commands.nstsize)
                                ierr = 0
                            except (FortranEOF, FortranIOError) as e:
                                ierr = iostat_of(e)
                            commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 11,
                                                          commands.nstsize, commands.nworddim)
                            ierr = _read_items(substr(commands.sometext, ncf[1], ncl[5]),
                                               [("real", molecom.epscor, is_), ("real", molecom.sigcor, is_),
                                                ("real", molecom.amassa, is_), ("real", molecom.charge, is_),
                                                ("real", molecom.qradius, is_)])
                            # !Get the coarse graining parameters if present
                            if commands.nwords > 5:
                                # !If this is negative it represents a density, not the number of atoms in the sphere
                                ierr = _read_items(substr(commands.sometext, ncf[6], ncl[6]),
                                                   [("real", molecom_mod.natomsphere, is_)])
                            else:
                                molecom_mod.natomsphere[is_] = 1.0
                            if commands.nwords > 7:
                                ierr = _read_items(substr(commands.sometext, ncf[7], ncl[8]),
                                                   [("real", molecom_mod.radminsphere, is_),
                                                    ("real", molecom_mod.radmaxsphere, is_)])
                                if molecom_mod.radmaxsphere[is_] <= 0.0:
                                    # In case default value is read in, set this to radius of atom
                                    molecom_mod.radmaxsphere[is_] = f32(0.5) * molecom.sigcor[is_]
                            else:
                                molecom_mod.radminsphere[is_] = 0.0
                                molecom_mod.radmaxsphere[is_] = f32(0.5) * molecom.sigcor[is_]
                            if commands.nwords > 8:
                                ierr = _read_items(substr(commands.sometext, ncf[9], ncl[9]),
                                                   [("real", molecom_mod.pseudofacmol, is_)])
                            else:
                                molecom_mod.pseudofacmol[is_] = 1.0

                        else:

                            # If not ignore the next line

                            units.skip_record(10)

                else:

                    # Signal that the specified template file was not found

                    units.write(6, fmt(_UA_F596, substr(symbols.fnametemplate, 1, symbols.ilenfnametemplate)))

                units.close(10)

        # O.K. so the template file has been read. Now try to apply it to the current molecule.

        # First check that the number of atoms in the template file is equal to the number
        # of atoms in the current molecule

        if clustercom.nato[im] == natomtemplate:

            # Now step through the atoms and if the labels correspond, set the neighbour numbers
            # and bond lengths

            nadd = clustercom.naf[im] - 1

            for ia in range(clustercom.naf[im], clustercom.nal[im] + 1):

                itref = ia - nadd

                # Only change the bond parameters if the labels match

                if molecom.alabel[ia] == alabeltemplate[itref]:

                    molecom.nbond[ia] = nbondtemplate[itref]
                    ja = 0
                    while ja < molecom.nbond[ia]:
                        ja = ja + 1
                        molecom.iabond[ia, ja] = iabondtemplate[itref, ja] + nadd
                        molecom.tbond[ia, ja] = tbondtemplate[itref, ja]

            # Rotational groups

            clustercom.ngrot[im] = nrotgrptemplate
            molecom_mod.ngrrot = molecom_mod.ngrrot + int(clustercom.ngrot[im])
            molecom_mod_lib.testsize_rotgrps()
            ir = 0
            while ir < clustercom.ngrot[im]:

                ir = ir + 1
                ig = ig + 1
                molecom_mod.rotgrp[ig] = rotgrptemplate[ir]
                molecom_mod.labaxi[ig, 1] = labaxitemplate[ir, 1] + nadd
                molecom_mod.labaxi[ig, 2] = labaxitemplate[ir, 2] + nadd
                molecom_mod.nrot[ig] = nrottemplate[ir]
                molecom_mod.rotf[ig] = molecom_mod.nrottot + 1
                molecom_mod.nrottot = molecom_mod.nrottot + int(molecom_mod.nrot[ig])
                molecom_mod.rotl[ig] = molecom_mod.nrottot
                molecom_mod_lib.testsize_labrot()
                for j in range(1, molecom_mod.nrot[ig] + 1):
                    labref = molecom_mod.rotf[ig] + j - 1
                    molecom_mod.labrot[labref] = labrottemplate[j, ir] + nadd
                if molecom_mod.rotgrp[ig] == "DIH":
                    molecom.ndihedralang[ig] = ndihedralangs[ir]
                    for k in range(1, molecom.ndihedralang[ig] + 1):
                        molecom.dihedralang[k, ig] = dihedralangs[k, ir]
                        for j in range(1, 3):
                            molecom.labrotdih[j, k, ig] = labrotdihs[j, k, ir] + nadd
                else:
                    molecom.ndihedralang[ig] = 0

        else:

            # Remember to update the rotational groups reference integers, in case they are
            # not accessed by the template file

            molecom_mod.ngrrot = molecom_mod.ngrrot + int(clustercom.ngrot[im])
            ig = molecom_mod.ngrrot

    return
