"""Translation of fortran_src/epsr_standalone_rev3.f, part 1 (lines 1-2522):
``program epsr_standalone`` and ``subroutine runepsr17`` -- the EPSR driver.

The other routines of the Fortran file (``interp``, ``inita``, ``initg``, ``sumgr``,
``ptgnsq``, ``get_wts``, ``gensq``, ``lsparabola1``, ...) are translated in
``epsr_standalone_rev3_part2.py`` / ``_part3.py`` and merged into this module by the
integrator; ``runepsr17`` calls them by bare name.

The library modules (and the Fortran-module data modules) are imported inside the
routines, because several of them are being translated at the same time as this file;
the module objects are used exactly as with a top-level import (``molfit.molfit(...)``).

Translation of fortran_src/epsr_standalone_rev3.f lines 2523-4736 (package J).

interp, inita, j1qr, j1qr_mod, gaussianqr2, lorentzianqr2, pofr, initg, entrpa, entrpb,
entrpc, sumgr, sumsq, sumpe, initpt, ptgnsq, repulsivept, norma, openf, get_wts, wtype,
wxinit, wxtype, wxcreate.

This file is merged by the integrator into programs/epsr_standalone_rev3.py (runepsr17,
package I, before it; gensq ... multiply_pcoefex, package K, after it).  Routines of the
same Fortran file are called by bare name: uniform, generate_rp, generate_ep,
calc_gaussq, calc_pwrexpq, fit_expansion and invert_cwt are defined in _part3.

The Fortran dummy argument ``in`` is a Python keyword and is spelled ``in_``.  The compat layer is
imported by name (not ``*``) so that ruff reports undefined names; see docs/translation_notes/J.md.

Translation of fortran_src/epsr_standalone_rev3.f, lines 4737-6584 (package K).

gensq, fit_expansion, expand, calc_gaussq, calc_gaussr, calc_gaussr_grad, calc_pwrexpq,
calc_pwrexpr, calc_pwrexpr_grad, dnext, ftback, uniform, invert_cwt, generate_rp,
generate_ep, lsparabola, lsparabola1, multiply_pot_coeffs, multiply_pcoefex.

To be merged by the integrator into programs/epsr_standalone_rev3.py (runepsr17, package
I; lines 2523-4736 are package J's _part2.py).  Routines of the same Fortran file are
called by bare name.  auxiliary_routines and ato_manipulation_rev4 are imported inside
the routines that use them because those modules did not exist yet when this part was
written; the integrator may hoist the imports to module level.
See docs/translation_notes/K.md.
"""

import math
from .. import trace
from ..common import clustercom
from ..common import commands
from ..common import directory_name
from ..common import epsr16
from ..common import epsrinput_par
from ..common import iresetnumbers
from ..common import lattice
from ..common import molecom
from ..common import potential_par
from ..common import sumvalues
from ..common import symbols
from ..common import version
from ..common import weights
from ..common import x_rayf0
from ..common import molecom_mod
from ..common import containerforces_mod as containerforces_mod_data
from ..common import lngamma_mod as lngamma_mod_data
from ..common.epsr_dimension import mfile
from ..common.epsr_dimension import mpairs
from ..common.epsr_dimension import mpts
from ..common.epsr_dimension import mpz
from ..common.epsr_dimension import mqd
from ..common.epsr_dimension import mpairs as _mpairs
from ..fortran import *  # noqa: F403
from ..fortran import PI4ATAN
from ..fortran import Ch
from ..fortran import FortranEOF
from ..fortran import FortranIOError
from ..fortran import char_array
from ..fortran import expf
from ..fortran import f32
from ..fortran import findex
from ..fortran import findex_padded
from ..fortran import fmt
from ..fortran import fstr
from ..fortran import idiv
from ..fortran import ifix
from ..fortran import int_array
from ..fortran import lw
from ..fortran import nint
from ..fortran import pow3
from ..fortran import powf
from ..fortran import r4
from ..fortran import read_internal
from ..fortran import real_array
from ..fortran import runtime
from ..fortran import sign
from ..fortran import sqrtf
from ..fortran import substr
from ..fortran import units
from ..fortran import FortranStop
from ..fortran import cosf
from ..fortran import dble_array
from ..fortran import dsqrt
from ..fortran import logf
from ..fortran import pad
from ..fortran import sinf
from ..library import ato_manipulation_rev4
from ..library import auxiliary_routines
from ..library import file_and_directory_commands
from ..library import parse
from ..library import reallocation_routines
from ..library import containerforces_mod
from ..library import lngamma_mod
from ..trace import routine


# =============================================================================
# from epsr_standalone_rev3.py
# =============================================================================

# FORMAT statements of runepsr17
_F2224 = "(i7,3g14.7)"
_F352 = "(' Read ',i5,' points from file ',a)"
_F958 = "(1x,6(1x,i5))"
_F912 = "((1x,a))"
_F1999 = "(2(1x,e13.6),2(1x,i16),1x,e13.6)"
_F331 = "(a2)"
_F992 = "('(3x,a,1h-,a,',i2.2,'x,a,1h-,a,',i2.2,'x)')"
_F594 = "(10(1x,e13.6))"
_F146 = "(/2x,'Current date: ',a2,':',a2,':',a2,' on ',a2,'-',a2,'-',a4/2x,'EPSR version ',a,/)"
_F147 = (
    "(1x,' Trial no. ',i5,5x,' Iteration no. ',i5"
    ",5x,' No. of configurations in sum ',i5"
    ",/1x,' Chi-square ',e13.6,' R-space coeff.: ',e13.6,' Abs. empirical energy: ',1x,e13.6"
    ",/1x,' Quality factor ',e13.6,' Control level, baseline and threshold: ',3(e13.6,1x)"
    ",/1x,' Pair energies: repulsive ',e13.6,' harmonic ',e13.6,' per molecule'"
    ",/1x,' Tethered energy: 'e13.6,' per molecule'"
    ",/1x,' Coordination constraint energy: 'e13.6,' per molecule'"
    ",/1x,' Three-body energy: 'e13.6,' per molecule'"
    ",/1x,' Three-body pressure: 'e13.6,' kbar'"
    ",/1x,' No. of .erg lines: ',i7,' Average gradient of quality factor: ',5(e13.6,1x)"
    ")"
)
_F145 = "(1x,' Sizefactor ',f10.5,' Number density ',e12.5)"
_F148 = "((1x,' R-factor for data file ',i2,':',1x,e11.4))"
_F149 = "((1x,' R-factor for partial ',i2,':',1x,e11.4))"
_F219 = (
    "(1x,' Noise in average g(r) ',e13.6"
    "/1x,' Intermolecular energy and U_pol ',e13.6,1x,e13.6,' kJ/mole'"
    "/1x,' Intermolecular pressure ',e13.6,' kbar'"
    "/1x,' Average dipole moment and std. ',e13.6,1x,e13.6,' D'"
    # c     1/1x,' Exponential energy ',e13.6,' kJ/mole'
    # c     1/1x,' Exponential pressure ',e13.6,' kbar'
    "/1x,' Av. intermolecular energy ',e13.6,'+-',e13.6,' kJ/mole'"
    "/1x,' Av. intermolecular pressure ',e13.6,'+-',e13.6,' kbar')"
)
_F220 = "(1x,' Number of rejects in MOLFIT:-'/((1x,4(i7,1x))))"
_F278 = "(i5/(6(1x,e12.5)))"


def _ios(e):
    """The iostat value of an I/O exception (FortranEOF -> -1)."""
    return -1 if isinstance(e, FortranEOF) else e.iostat


def _stop_requested():
    """True once the user pressed Ctrl+C in the shell during this run (shell/signals.py)."""
    from ..shell import signals

    return signals.stop_requested()


def _iteration_finished(iter, niter):
    """Hook called right after each EPSR iteration of runepsr17 completes (after the
    ``3000 continue``).  Returns True when the loop must end (graceful Ctrl+C)."""
    return iter >= niter


def main():
    """epsr_standalone_rev3.f:epsr_standalone -- program to perform the EPSR command.

    2009/3/30 Attempt to make this a standalone program
    """
    from ..library import file_and_directory_commands

    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    runepsr17()


@routine
def runepsr17():
    """epsr_standalone_rev3.f:runepsr17 -- EPSR17_3: Empirical Potential Structure
    Refinement, fitting in Q space (setup, the ``do 3000 iter=1,niter`` loop, results).

    Standalone version of EPSR using only ASCII input and output files.  See the long
    history comment at the top of the Fortran routine.
    """
    # Library routines and Fortran-module data (imported here: see the module docstring)
    from ..common import diffuse_mod, molecom_mod, threebody_mod, timer
    from ..library import (
        ato_manipulation_rev4,
        ato_read_and_save,
        auxiliary_routines,
        cellneighbourlist,
        file_and_directory_commands,
        minimagec,
        molfit,
        parse,
        read_and_save_files,
        setup_epsr_rev1,
    )

    # c EPSR17_3.FOR
    # c Version 17_3 of EPSR (29/01/2004) for fitting in Q space.
    # (The full history comment of the Fortran is not repeated here; see
    #  fortran_src/epsr_standalone_rev3.f lines 19-153.)

    # use timer
    # use cellneighbourlist !Stores and processes the cell neighbour list
    # use diffuse_mod
    # use threebody_mod
    # use molecom_mod

    # ---------------------------------------------------------------- local arrays
    z = real_array(mpz)
    q = real_array(mpts)
    datai = real_array(mpts, mfile)
    errs = real_array(mpts, mfile)
    ufile = real_array(mpts, mfile)
    ugdif = real_array(mpz, mfile)
    ugdifsmo = real_array(mpz, mfile)  # noqa: F841
    qwind = real_array(2, mpts, mpairs)
    rfac = real_array(mpairs)
    datat = real_array(mpts, mfile)
    backsc = real_array(mpts, mfile)
    wav = real_array(2, mpz, mpairs)
    errsav = real_array(mpts, mpairs)
    soqism = real_array(mpts, mfile)
    rshmin = real_array(mfile)
    szeros = real_array(mfile)
    szerostemp = real_array(mfile)
    rfacd = real_array(mfile)
    potnoi = real_array(mfile)
    energf = real_array(mfile)
    pressf = real_array(mfile)
    efac = real_array(mpairs)  # noqa: F841
    tweak = real_array(mfile)
    datap = real_array(mpts, mpairs)
    datagr = real_array(mpz, mpairs)
    datagr2 = real_array(mpz, mpairs)
    dataigr = real_array(mpz, mfile)
    grlim = real_array(mpairs)
    grzerolim = real_array(mfile)
    efilereq = real_array(mfile)

    qd = real_array(mqd)
    data = real_array(mqd)
    err = real_array(mqd)
    gofzism = real_array(mpz, mfile)
    gz2ism = real_array(mpz, mfile)

    # c Temporary store of simulated intermolecular partial structure factors
    soqintertemp = real_array(mpts, mpairs)
    dts = real_array(mpts)
    dti = real_array(mpts)

    # c nrtype is the type of dataset being read - see below.
    # c nqminp and nqmaxp are used in the Fourier transform of the estimated partial structure
    # c factors.
    nrtype = int_array(mfile)
    nqmind = int_array(mfile)
    nqmaxd = int_array(mfile)
    nqminp = int_array(mpairs)
    nqmaxp = int_array(mpairs)

    infili = char_array(mfile)  # character*256
    fnamew = char_array(mfile)  # character*256
    fnamedata = char_array(mfile)  # character*256
    fnamewts = char_array(mfile)  # character*256
    titlepsf = char_array(mpairs)  # character*29
    titledat = char_array(mfile)  # character*29

    # ! Variables for Bragg peak calculation if needed
    qbragg = None  # real, allocatable, dimension(:)
    intenbragg = None  # real, allocatable, dimension(:,:)
    intenconv = None  # real, allocatable, dimension(:,:)

    # ! Variables for calculating the Porod scattering
    porodamp = real_array(mfile)
    porodadd = real_array(mfile)

    # ---------------------------------------------------------------- scalars
    # NB: the Fortran does not initialise the following locals before (possibly) reading
    # or writing them; 0 reproduces the reference (Rust epsr.rs: Control::default(),
    # Summary::default(), ctrllevel/ctrlbaseline = 0.0, parab = None -> 0.0).
    ierr = 0
    navmaxparab = 0
    nparabstartfactor = 0
    rmaxgrsave = f32(0.0)
    nparab = 0
    aparab = f32(0.0)
    bparab = f32(0.0)
    cparab = f32(0.0)
    aav = f32(0.0)
    bav = f32(0.0)
    lastqf = f32(0.0)
    nextqf = f32(0.0)
    epcup = f32(0.0)
    epccoeff = f32(0.0)
    ctrllevel = f32(0.0)
    ctrlbaseline = f32(0.0)
    fcurrent = f32(0.0)
    eminmaxdir = f32(0.0)  # NB: misspelt eminmaxdif in a write(6,*); never set in the Fortran
    rfactot = f32(0.0)
    chisq = f32(0.0)
    oiseb = f32(0.0)
    sentb = f32(0.0)
    energd = f32(0.0)
    pressd = f32(0.0)
    energabs = f32(0.0)
    energthreeb = f32(0.0)
    pressthreeb = f32(0.0)
    ftethered = f32(0.0)
    fcoordconstraints = f32(0.0)
    foldsv = f32(0.0)
    foldsb = f32(0.0)
    hashsearch = ""  # character*10
    newphib = f32(0.0)
    newthetac = f32(0.0)
    newphic = f32(0.0)
    alphaold = f32(0.0)
    betaold = f32(0.0)
    gamold = f32(0.0)
    # NB: overallfac is passed to ptgnsq but never set; ptgnsq does not use it.
    overallfac = f32(0.0)
    thresholde = f32(0.0)
    porodpwr = f32(0.0)

    # c
    # c EPSR begins....
    # c

    pi = PI4ATAN
    degconv = pi / r4(180.0)  # noqa: F841

    elimitfac = r4(0.02)  # noqa: F841

    # c
    # c set up the EPSR symbol tables
    # c
    setup_epsr_rev1.setup_symbols_epsr()

    # c
    # c set up the EPSR input file
    # c
    nin = 11
    nout = 6
    ierrout = 0
    nout, ierrout = file_and_directory_commands.setup_input_file(nin, nout, ierrout)
    # c
    # c if an error occurred on setup, exit now
    # c
    if ierrout != 0:
        return

    # c
    # c read the input file
    # c
    nout, ierrout = read_and_save_files.read_input_file(nin, nout, ierrout)
    # c
    # c if an error was encountered on read, exit now
    # c
    if ierrout != 0:
        return

    # c
    # c now read the symbol_values into their appropriate values in the
    # c program.
    # c
    icount = symbols.nsymboluptodata
    ifilecount = 1
    iref = file_and_directory_commands.ireff(icount, ifilecount)
    # c NFILE = number of data files to fit (30 maximum)
    epsr16.nfile = read_internal(symbols.symbol_value[iref]).int(epsr16.nfile)
    efilereqsum = f32(0.0)
    for in_ in range(1, epsr16.nfile + 1):
        # c
        # c INFILE = name of data file for S(Q) data
        # c
        iref = iref + 1
        infili[in_] = fstr(symbols.symbol_value[iref], 256)
        # c
        # c add on the directory name
        # c
        ilen = findex_padded(infili[in_], " ", 256) - 1
        fnamedata[in_] = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                              + substr(infili[in_], 1, ilen), 256)
        # c
        # c FNAMEW = name of file containing weighting factors for this dataset
        # c
        iref = iref + 1
        fnamew[in_] = fstr(symbols.symbol_value[iref], 256)
        # c
        # c add on the directory name
        # c
        ilen = findex_padded(fnamew[in_], " ", 256) - 1
        fnamewts[in_] = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                             + substr(fnamew[in_], 1, ilen), 256)
        # c
        # c type of input file (1,2, 3 or 4, 5)
        # c     NRTYPE=1, Q, S(Q) ASCII values, 1 set per line
        # c     NRTYPE=2, Q, S(Q), E(Q) ASCII values, 1 set per line
        # c     NRTYPE=3,4 N, Q, S(Q), E(Q) ASCII values, 1 set per line, 4 blank lines
        # c(GENIE SHOW DATA format, NRTYPE=3 corresponds to histogram format, 4 = point
        # c mode)
        # c     NRTYPE=5, Q, S(Q), E(Q) ASCII values, 1 set per line, using GUDRUN .DCS format (histogram)
        # c RSHMIN = minimum radius for this dataset - used in removing the background
        # c SZEROS is the Q=0 limit for this file. If non-zero it is used to extrapolate the data
        # c from the smallest Q data point to Q=0. If zero then the last data point is used as the
        # c q=0 value.
        # c TWEAK is a simple factor to multiply the data by prior to any processing. Normally 1.0.
        # c EFILEREQ is the relative importance of this dataset in calculating the inversion matrix
        # c - increasing it should improve the fit to this dataset in principle, at the expense of
        # c the fit to the other datasets. Set it to zero if you want to calculate the simulated
        # c pattern for this dataset but not use the data to calculate the empirical potential. EFILEREQ can take any value but it will be normalised internally to
        # c prevent the inversion matrix becoming indeterminate.
        # c
        iref = iref + 1
        nrtype[in_] = read_internal(symbols.symbol_value[iref]).int(nrtype[in_])
        iref = iref + 1
        rshmin[in_] = read_internal(symbols.symbol_value[iref]).real(rshmin[in_])
        iref = iref + 1
        szeros[in_] = read_internal(symbols.symbol_value[iref]).real(szeros[in_])
        iref = iref + 1
        tweak[in_] = read_internal(symbols.symbol_value[iref]).real(tweak[in_])
        iref = iref + 1
        efilereq[in_] = read_internal(symbols.symbol_value[iref]).real(efilereq[in_])
        efilereqsum = efilereqsum + efilereq[in_]
        iref = iref + 1
        porodadd[in_] = read_internal(symbols.symbol_value[iref]).real(porodadd[in_])

    # c NB efilereq cannot be zero for every dataset, otherwise the inversion matrix will crash if feedback > 0.0

    if efilereqsum <= 0.0:
        efilereq[1] = 1.0

    icount = 1
    ifilecount = 1
    iref = file_and_directory_commands.ireff(icount, ifilecount)
    # c
    # c FEEDBACK = confidence in data factor. Must be 0.0 <= FEEDBACK <= 1.0. Should only be set
    # c to 1.0 if the matrix of weight coefficients is well conditioned.
    # c Can be increased or decreased to emphasize or de-emphasize the data.
    # c POTFAC decides whether potential will be refined and how much of the difference
    # c between data and fit are fed back into the EP at each iteration. Typically POTFAC = 0.1
    # c or 0.2. If oscillations develop in total energy of simulation, probably means POTFAC is
    # c too large, so reduce its amplitude. However if POTFAC is too small then simulation will
    # c take a long time to find a fit. POTFAC=0.0 means no EP refinement, (but simulation will
    # c proceed as usual).
    # c PRESSREQ is the requested absolute energy (kJ/mole) of empirical potential.
    # c This can be used to control degree of fit and energy/pressure for all datasets together.
    # c SIZEFACTOR (should be >=1) is used to shrink molecules during initial refinement when
    # c molecular entanglement may occur. Advisable to set potfac to 0 if sizefactor > 1
    # c
    potential_par.feedback = read_internal(symbols.symbol_value[iref]).real(potential_par.feedback)
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    potential_par.potfacamp = read_internal(
        substr(commands.sometext, commands.ncf[1], commands.ncl[1])).real(potential_par.potfacamp)
    if commands.nwords > 1:
        potential_par.potfac = read_internal(
            substr(commands.sometext, commands.ncf[2], commands.ncl[2])).real(potential_par.potfac)
    else:
        potential_par.potfac = potential_par.potfacamp
    if potential_par.potfacamp == 0.0:
        potential_par.potfac = f32(0.0)
    if potential_par.potfac > 1.0:
        potential_par.potfac = f32(1.0)
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    pressreq = read_internal(substr(commands.sometext, commands.ncf[1], commands.ncl[1])).real(f32(0.0))
    potential_par.absolute_energy = pressreq
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    ereqmin = read_internal(substr(commands.sometext, commands.ncf[1], commands.ncl[1])).real(f32(0.0))
    if ereqmin < 0.0:
        ereqmin = f32(0.0)
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    ereqmax = read_internal(substr(commands.sometext, commands.ncf[1], commands.ncl[1])).real(f32(0.0))
    if ereqmax < 0.0:
        ereqmax = f32(0.0)
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    ereqstepread = read_internal(substr(commands.sometext, commands.ncf[1], commands.ncl[1])).real(f32(0.0))
    ereqstepread = max(ereqstepread, f32(0.0))
    # !      ereqstepread=min(ereqstepread,1.0)
    potential_par.ereqstep = ereqstepread
    potential_par.erequnit = potential_par.ereqstep
    iref = iref + 1
    epcthresh = read_internal(symbols.symbol_value[iref]).real(f32(0.0))
    epcthreshuse = epcthresh  # noqa: F841
    iref = iref + 1
    ereqbias = read_internal(symbols.symbol_value[iref]).real(f32(0.0))
    iref = iref + 1
    r = read_internal(symbols.symbol_value[iref])
    navmaxparab = r.int(navmaxparab)
    nparabstartfactor = r.int(nparabstartfactor)
    if navmaxparab <= 0 or nparabstartfactor < 1:  # If a problem with either value, set the defaults.
        navmaxparab = 50
        nparabstartfactor = 1
    nparabstart = navmaxparab * nparabstartfactor
    potential_par.erequnit = potential_par.ereqstep
    iref = iref + 1
    epcrmin = read_internal(symbols.symbol_value[iref]).real(f32(0.0))
    iref = iref + 1
    epcfrac = read_internal(symbols.symbol_value[iref]).real(f32(0.0))
    iref = iref + 1
    num_threads = read_internal(symbols.symbol_value[iref]).int(0)
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    porodq = read_internal(substr(commands.sometext, commands.ncf[1], commands.ncl[1])).real(f32(0.0))
    if commands.nwords < 2 or porodq <= 0.0:
        porodamp[:] = 0.0
        porodpwr = f32(0.0)
    else:
        porodpwr = read_internal(substr(commands.sometext, commands.ncf[2], commands.ncl[2])).real(porodpwr)
    iref = iref + 1
    nmolcellin = read_internal(symbols.symbol_value[iref]).int(0)
    # !Determine if pressure control is required
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    control_pressure = False
    pacontrol = False
    pbcontrol = False
    pccontrol = False
    ppcontrol = False
    control_p_value = f32(0.0)
    control_p_step = f32(0.0)
    if commands.nwords > 0:
        word1 = substr(commands.sometext, commands.ncf[1], commands.ncl[1])
        # FIX (FORTRAN_BUGS.md 1.5): the Fortran accepts 'ca' here but its axis branch
        # below tests 'ac', so 'ca' selected no axis and 'ac' was rejected; as in the Rust
        # (pressure.rs) both 'ca' and 'ac' are accepted and select the c and a axes.
        control_pressure = commands.nwords > 2 and (
            word1 == "a"
            or word1 == "b"
            or word1 == "c"
            or word1 == "ab"
            or word1 == "bc"
            or word1 == "ca"
            or word1 == "ac"
            or word1 == "abc"
            or word1 == "all"
            or word1 == "p")
        if control_pressure:
            try:
                r = read_internal(substr(commands.sometext, commands.ncf[2], commands.ncl[3]))
                control_p_value = r.real(control_p_value)
                control_p_step = r.real(control_p_step)
                ierr = 0
            except (FortranEOF, FortranIOError) as e:
                ierr = _ios(e)
            control_pressure = ierr == 0 and control_p_step > 0.0
        if control_pressure:
            if word1 == "a":
                pacontrol = True
                pbcontrol = False
                pccontrol = False
                ppcontrol = False
            elif word1 == "b":
                pacontrol = False
                pbcontrol = True
                pccontrol = False
                ppcontrol = False
            elif word1 == "c":
                pacontrol = False
                pbcontrol = False
                pccontrol = True
                ppcontrol = False
            elif word1 == "ab":
                pacontrol = True
                pbcontrol = True
                pccontrol = False
                ppcontrol = False
            elif word1 == "bc":
                pacontrol = False
                pbcontrol = True
                pccontrol = True
                ppcontrol = False
            elif word1 == "ac" or word1 == "ca":  # FIX (FORTRAN_BUGS.md 1.5): Fortran tests 'ac' only
                pacontrol = True
                pbcontrol = False
                pccontrol = True
                ppcontrol = False
            elif word1 == "abc" or word1 == "all":
                pacontrol = True
                pbcontrol = True
                pccontrol = True
                ppcontrol = False
            elif word1 == "p":
                # !This represents pressure control using the repulsive potential
                pacontrol = False
                pbcontrol = False
                pccontrol = False
                ppcontrol = True
    npcontrols = 0
    if pacontrol:
        npcontrols = npcontrols + 1
    if pbcontrol:
        npcontrols = npcontrols + 1
    if pccontrol:
        npcontrols = npcontrols + 1
    if ppcontrol:
        npcontrols = npcontrols + 1
    # c Determine if to refine the intra-molecular structure
    iref = iref + 1
    molecom.refine_intra = False
    molecom.coul_intra = False
    molecom.refine_intra_fac = f32(0.0)
    molecom.coul_intra_fac = f32(0.0)
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    if commands.nwords > 0:
        try:
            molecom.refine_intra_fac = read_internal(
                substr(commands.sometext, commands.ncf[1], commands.ncl[1])).real(molecom.refine_intra_fac)
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = _ios(e)
    if molecom.refine_intra_fac > 0.0:
        molecom.refine_intra = True
    if commands.nwords > 1:
        try:
            molecom.coul_intra_fac = read_internal(
                substr(commands.sometext, commands.ncf[2], commands.ncl[2])).real(molecom.coul_intra_fac)
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = _ios(e)
    if molecom.coul_intra_fac > 0.0:
        molecom.coul_intra = True
    symbols.symbol_value[iref] = fstr(fmt("(f8.3,1x,f8.3)", molecom.refine_intra_fac, molecom.coul_intra_fac), 256)
    ierr = 0
    iref = iref + 1
    irefsizefactor = iref
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    molecom.sizefactor = read_internal(
        substr(commands.sometext, commands.ncf[1], commands.ncl[1])).real(molecom.sizefactor)
    if commands.nwords > 1:
        molecom.decline_rate = read_internal(
            substr(commands.sometext, commands.ncf[2], commands.ncl[2])).real(molecom.decline_rate)
        # ! Make sure crazy values are ignored
        if molecom.decline_rate > 1.0 or molecom.decline_rate < 0.1:
            molecom.decline_rate = r4(0.9)
    else:
        molecom.decline_rate = r4(0.9)
    # ! Set the threshold energy below which the size factor is reduced.
    if commands.nwords > 2:
        thresholde = read_internal(substr(commands.sometext, commands.ncf[3], commands.ncl[3])).real(thresholde)
    else:
        thresholde = f32(0.0)
    # !Save the values of stepmi and stepri
    stepmisave = molecom.stepmi
    steprisave = molecom.stepri
    # !If sizefactor > 1.0, then do not intramolecular rotations and refinement.
    # ! Also do not update repulsive or empirical potentials.
    if molecom.sizefactor > 1.0:
        potential_par.potfac = f32(0.0)
        molecom.stepmi = f32(0.0)
        # !         stepri=0.0
    elif potential_par.potfacamp > 0.0 and potential_par.potfac <= 0.0:
        potential_par.potfac = potential_par.potfacamp
    # !      if(sizefactor.gt.2.0) then
    # !         stepri=0.0
    # !      end if
    units.write(6, fmt('("Steps: ",4(1x,e12.5))', molecom.stepmi, molecom.stepri, stepmisave, steprisave))
    if potential_par.feedback > 1.0:
        potential_par.feedback = r4(0.9)
    if potential_par.feedback < 0.0:
        potential_par.feedback = f32(0.0)
    # c do not attempt to refine potential if there is no weight given to the data
    if potential_par.feedback <= 0.0:
        potential_par.potfacamp = f32(0.0)
        potential_par.potfac = f32(0.0)
    # c POLAR = polarizablity (not currently used)
    polar = f32(0.0)

    # c Set the data importance factors to the requested value. These have to be
    # c normalised to the number of non-zero values to avoid overemphasizing one dataset
    # c which could make the inversion matrix indeterminate.

    nonzero = 0
    sumnonzero = f32(0.0)
    for in_ in range(1, epsr16.nfile + 1):
        if efilereq[in_] > 0.0:
            nonzero = nonzero + 1
            sumnonzero = sumnonzero + efilereq[in_]

    # c Create the normalised dataset importance factors, EFACM

    if nonzero > 0 and sumnonzero > 0.0:
        sumnonzero = sumnonzero / f32(nonzero)
    else:
        sumnonzero = f32(1.0)
    for in_ in range(1, epsr16.nfile + 1):
        potential_par.efacm[in_] = efilereq[in_] / sumnonzero

    # c NQ = No. of Q values (Maximum = 1000)
    # c QSTEP = requested step in  Q
    iref = iref + 1
    nq = read_internal(symbols.symbol_value[iref]).int(0)
    if nq > mpts:
        nq = mpts
    iref = iref + 1
    qstep = read_internal(symbols.symbol_value[iref]).real(f32(0.0))
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    molecom_mod.qscaleinputfactor = read_internal(
        substr(commands.sometext, commands.ncf[1], commands.ncl[1])).real(molecom_mod.qscaleinputfactor)
    if commands.nwords > 1:
        molecom_mod.selfbroadpwr = read_internal(
            substr(commands.sometext, commands.ncf[2], commands.ncl[2])).real(molecom_mod.selfbroadpwr)
    else:
        molecom_mod.selfbroadpwr = f32(0.0)
    molecom_mod.coarsegraining = molecom_mod.qscaleinputfactor > 1.0
    # c NEGAT not currently operational - hold at zero
    negat = 0  # noqa: F841
    # c IRESET = 1 to set EP to zero and reinitialize certain accummulators
    # c     else IRESET = 0.
    iref = iref + 1
    ireset = read_internal(symbols.symbol_value[iref]).int(0)
    if ireset == 1:
        pressreq = f32(0.0)
        potential_par.ereqstep = f32(0.0)
        potential_par.erequnit = f32(0.0)
    elif ireset == 2:
        # !         pressreq=0.0
        if potential_par.ereqstep > 0.0:
            pressreq = potential_par.erequnit
    potential_par.absolute_energy = pressreq
    iref = iref + 1
    iinit = read_internal(symbols.symbol_value[iref]).int(0)
    # c NTIMES = no. of trials per iteration (typically 1, but could be more)
    # c NITER = no. of iterations (typically 1)
    # c NSUMS is the number of previous iterations that have already been accumulated
    # c     If NSUMS is -1, then a running average is used to
    # c     calculate g(r) and S(Q), as per average energy and average pressure. (See comments in
    # c subroutine INITA.)
    iref = iref + 1
    ntimes = read_internal(symbols.symbol_value[iref]).int(0)
    iref = iref + 1
    niter = read_internal(symbols.symbol_value[iref]).int(0)
    iref = iref + 1
    epsrinput_par.nsums = read_internal(symbols.symbol_value[iref]).int(epsrinput_par.nsums)
    # c INTRA is the number of intermolecular moves that are made before an INTRA-
    # c molecular move (molecular vibration) is made - typically INTRA=50 or 100.
    # ! rotfreq is the number of iterations between performing internal molecular rotations
    # c INTER is the number of iterations the potential is averaged over
    iref = iref + 1
    intra = read_internal(symbols.symbol_value[iref]).int(0)
    iref = iref + 1
    rotfreq = read_internal(symbols.symbol_value[iref]).int(0)
    iref = iref + 1
    inter = read_internal(symbols.symbol_value[iref]).int(0)
    # c INTER cannot be less than 1
    if inter < 1:
        inter = 1
    # c number of iterations to occur after intramolecular disordering before
    # c accumulating energy and pressure averages.
    nisumm = ifix(r4(0.3) * f32(intra))
    # c RHO no. density in scattering units per A**3 (this is determined by
    # c the .ATO file - this value is written here for information only)
    iref = iref + 1
    epsr16.rho = read_internal(symbols.symbol_value[iref]).real(epsr16.rho)
    # c NDIST = number of r-space distributions (automatically set to npairs)
    # c CELLST = step size in linear region of r - controls number of r values
    # c
    iref = iref + 1
    epsr16.cellst = read_internal(symbols.symbol_value[iref]).real(epsr16.cellst)
    # ! Increase in rmax compared to default (0.5*cell)
    iref = iref + 1
    try:
        rmaxgrsave = read_internal(symbols.symbol_value[iref]).real(rmaxgrsave)
        ierr = 0
    except (FortranEOF, FortranIOError) as e:
        ierr = _ios(e)
    epsrinput_par.rmaxgr = rmaxgrsave
    iref = iref + 1
    ngrsamples = read_internal(symbols.symbol_value[iref]).int(0)
    # c The following factors should not be used in this version of the program.
    efacs = f32(0.0)  # noqa: F841
    pfaclj = f32(0.0)  # noqa: F841
    gfacs = f32(0.0)
    gfacmn = f32(0.0)
    gfacmx = f32(0.0)
    if gfacs > 0.0 and gfacs < gfacmn:
        gfacs = gfacmn
    if gfacs > 0.0 and gfacs > gfacmx:
        gfacs = gfacmx
    gfac = gfacs
    # c resolution function of data (Gaussian) FWHM(Q) = FWHM + FWHMQ*Q
    iref = iref + 1
    fwhm = read_internal(symbols.symbol_value[iref]).real(f32(0.0))
    iref = iref + 1
    fwhmq = read_internal(symbols.symbol_value[iref]).real(f32(0.0))
    # c WTENTP is not used
    # c NSMOOp = defines whether a background based on Rmin will be subtracted
    # c from the data prior to fitting the potential
    # c (0 for no background subtraction)
    wtentp = f32(0.0)  # noqa: F841
    iref = iref + 1
    nsmoop = read_internal(symbols.symbol_value[iref]).int(0)
    # c read broadening of Q space data using revised Lorch function
    iref = iref + 3
    epsrinput_par.lorchqwidth = f32(0.0)
    epsrinput_par.lorchqwidthq = f32(0.0)
    try:
        r = read_internal(symbols.symbol_value[iref])
        epsrinput_par.lorchqwidth = r.real(epsrinput_par.lorchqwidth)
        epsrinput_par.lorchqwidthq = r.real(epsrinput_par.lorchqwidthq)
        ierr = 0
    except (FortranEOF, FortranIOError) as e:
        ierr = _ios(e)
    # !Read broadening and maximum Q for Bragg peak calculation, and Debye-Waller factor.
    # !Use 0 for qbraggmax if you do not want to do this calculation
    iref = iref + 1
    r = read_internal(symbols.symbol_value[iref])
    qbraggwidth = r.real(f32(0.0))
    qbraggmax = r.real(f32(0.0))
    # !Read multiplicities for a, b and c axes. This so that when calculating Bragg peak intensities
    # !proper account is taken of the number of unit cells used to build the EPSR box
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    r = read_internal(substr(commands.sometext, commands.ncf[1], commands.ncl[3]))
    multa = r.int(0)  # multipliticities along a, b and c
    multb = r.int(0)
    multc = r.int(0)
    if commands.nwords > 5:
        # use this axis if non-zero, otherwise no Bragg values along this axis
        r = read_internal(substr(commands.sometext, commands.ncf[4], commands.ncl[6]))
        usea = r.int(0)
        useb = r.int(0)
        usec = r.int(0)
        if usea != 0:
            usea = 1  # Must be either 0 or 1
        if useb != 0:
            useb = 1  # Must be either 0 or 1
        if usec != 0:
            usec = 1  # Must be either 0 or 1
    else:
        usea = 1
        useb = 1
        usec = 1
    # c read minimum Q value to calculate Bragg intensities and radius value to change
    # c from local structure to long range structure
    iref = iref + 1
    r = read_internal(symbols.symbol_value[iref])
    hklqmin = r.real(f32(0.0))
    epsrinput_par.rminbragg = r.real(epsrinput_par.rminbragg)
    debyewaller = r.real(f32(0.0))
    # !Ensure the Bragg intensities are not calculated if rminbragg < 0.0
    if epsrinput_par.rminbragg <= 0.0:
        qbraggmax = f32(0.0)
    # !If qbraggmax .gt. 0, ensure rminbragg is large enough to not affect amplitude at qbraggmax
    # !      if(qbraggmax.gt.0.0) rminbragg=max(rminbragg,30.0/qbraggmax)
    # !Allocate required arrays
    nbragg = 0
    nbraggup = 0
    if qbraggmax > 0.0:
        nbragg = nint(qbraggmax / qstep) + 1
        # !For the convolution we need to include more Bragg peaks so there is no
        # !discontinuity at qbraggmax
        nbraggup = nint((qbraggmax + r4(3.0) / epsrinput_par.rminbragg) / qstep) + 1
        if nbraggup > nq:
            nbraggup = nq
    # !Specify whether diffuse scattering will be calculated. If so we need the number of unit cells along each axis, and h,k,l pattern
    # ! e.g. 1 2 0 would mean any h or k for which l=0, 1 1 2 wouild mean any h,k,l for which h=k, etc.
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    diffuse_mod.diffusedefined = commands.nwords > 5
    if diffuse_mod.diffusedefined:
        r = read_internal(symbols.symbol_value[iref])
        diffuse_mod.multdiffa = r.int(diffuse_mod.multdiffa)
        diffuse_mod.multdiffb = r.int(diffuse_mod.multdiffb)
        diffuse_mod.multdiffc = r.int(diffuse_mod.multdiffc)
        diffuse_mod.ldiffmax = r.int(diffuse_mod.ldiffmax)
        diffuse_mod.ldiffstep = r.int(diffuse_mod.ldiffstep)
        diffuse_mod.mdiffstep = r.int(diffuse_mod.mdiffstep)
        if diffuse_mod.ldiffmax > 40:
            # Practical limit based on current accuracy of dlms (generalised spherical harmonics)
            diffuse_mod.ldiffmax = 40
    else:
        diffuse_mod.ndiffuse = 0
    # ! Variables for lattice refinement if needed
    iref = iref + 1
    nlatrefine = read_internal(symbols.symbol_value[iref]).int(0)
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    alatstep = f32(0)
    blatstep = f32(0)
    clatstep = f32(0)
    alphastep = f32(0)
    betastep = f32(0)
    gamstep = f32(0)
    if commands.nwords > 5 and nlatrefine > 0:
        try:
            r = read_internal(substr(commands.sometext, commands.ncf[1], commands.ncl[6]))
            alatstep = r.real(alatstep)
            blatstep = r.real(blatstep)
            clatstep = r.real(clatstep)
            alphastep = r.real(alphastep)
            betastep = r.real(betastep)
            gamstep = r.real(gamstep)
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = _ios(e)
        if ierr != 0:
            alatstep = f32(0)
            blatstep = f32(0)
            clatstep = f32(0)
            alphastep = f32(0)
            betastep = f32(0)
            gamstep = f32(0)
            nlatrefine = 0
    refinelattice = nlatrefine > 0
    # c read requested Monte Carlo rejection ratio
    iref = iref + 1
    epsr16.rejrate = read_internal(symbols.symbol_value[iref]).real(epsr16.rejrate)
    # c read minimum Q values for Fourier transforms and potential fit, qminm and qminp. If qminp is zero, the minimum
    # ! Q for fitting the potential is qminm.
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                  commands.nstsize, commands.nworddim)
    if commands.nwords > 0:
        try:
            epsr16.qminm = read_internal(substr(commands.sometext, commands.ncf[1], commands.ncl[1])).real(epsr16.qminm)
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = _ios(e)
        if commands.nwords > 1:
            try:
                epsr16.qminp = read_internal(
                    substr(commands.sometext, commands.ncf[2], commands.ncl[2])).real(epsr16.qminp)
                ierr = 0
            except (FortranEOF, FortranIOError) as e:
                ierr = _ios(e)
        else:
            epsr16.qminp = f32(0.0)
    else:
        epsr16.qminm = f32(0.0)
        epsr16.qminp = f32(0.0)
    symbols.symbol_value[iref] = fstr(fmt("(2(1x,f10.5))", epsr16.qminm, epsr16.qminp), 256)
    facpress = r4(0.001)  # noqa: F841
    facpress3 = r4(0.01)  # noqa: F841
    # c
    # c define pressure and energy factors not currently used, but needed.
    # c
    eljfac = f32(1.0)  # noqa: F841
    pljfac = f32(1.0)
    pifac = r4(1.0) / sqrtf(r4(2.0) * pi)  # noqa: F841
    epsr16.cell = epsr16.cellmax
    datafc = r4(4.0) * pi * epsr16.rho  # noqa: F841
    # c
    # c write the repulsive potential coefficients
    # c
    # c     do ip=1,npairs
    # c           write(6,*) n1extr(ip),n2extr(ip),dcor(ip)
    # c     end do
    # c
    # c initialize distributions
    # c
    units.write(6, lw(epsr16.rmaxcell, epsr16.cellmax))
    nz = 0
    iinit, nz = inita(ireset, iinit, nq, q, qstep, nz, z, fwhm, fwhmq, intra)  # noqa: F821
    units.write(6, lw(epsrinput_par.rminpt, epsrinput_par.rmaxpt))

    # !Find integer corresponding to minimum Q
    nqminm = 1
    for i in range(1, nq + 1):
        if q[i] <= epsr16.qminm:
            nqminm = i
    # c
    # c Set up the weights matrix and its inverse
    # c
    get_wts(q, nq, fnamewts, iinit)  # noqa: F821
    # c
    # c set the value of nesumt to inter
    # c
    if iinit != 0:
        inter = epsrinput_par.nesumt
    else:
        epsrinput_par.nesumt = inter

    # c
    # c normalise atom and molecular coordinates to box size
    # c
    ato_read_and_save.normalise_coordinates()
    # ! Get the maximum dimension of any molecule
    ato_manipulation_rev4.get_maxmolsize()
    units.write(6, lw("epsr> Maximum molecular dimension (A): ", molecom.maxmolsize * clustercom.dmax))
    # c
    # c Ignore energy from polarisation unless needed
    # c
    epsr16.upol = epsr16.upol * epsrinput_par.rcharge
    # c
    # c initial g(r) sums
    # c
    units.write(6, lw(epsrinput_par.nsum, epsrinput_par.nesum, epsrinput_par.nsums))
    if epsrinput_par.nsums >= 0:
        epsrinput_par.nsum = initg(nq, nz, symbols.fnamerun, symbols.ilenfnamerun,  # noqa: F821
                                   epsrinput_par.nsum, ireset)
    else:
        epsrinput_par.nesum = initg(nq, nz, symbols.fnamerun, symbols.ilenfnamerun,  # noqa: F821
                                    epsrinput_par.nesum, ireset)

    # c initialize reference and empirical potentials

    initpt(ireset, epsrinput_par.rminpt, epsrinput_par.rmaxpt, epsrinput_par.expecf, polar,  # noqa: F821
           szeros, rshmin, epsrinput_par.rbroad)

    # ! Determine position of data for Porod scattering, if applicable. Assumes q increases monotonically with iq.
    nporodq = 0
    if porodq > 0.0:
        foundporodq = False
        porodqdif = porodq
        porodqmin = porodqdif
        iq = 0
        while iq < nq and not foundporodq:
            iq = iq + 1
            porodqdif = abs(porodq - q[iq])
            foundporodq = porodqdif > porodqmin
            if not foundporodq:
                porodqmin = porodqdif
                nporodq = iq

    # C read in data to be fitted

    sumerr = f32(0.0)  # noqa: F841
    for in_ in range(1, epsr16.nfile + 1):
        infile = fnamedata[in_]  # character*256

        # c get the length of this name including spaces

        ilendot = findex(infile, ".")
        ilentot = ilendot + findex(substr(infile, ilendot + 1, 256), " ") - 1
        units.write(6, lw(nrtype[in_], " ", substr(infile, 1, ilentot)))
        units.open(8, infile, status="old")

        # c Read a line to see if there are any problems

        try:
            units.skip_record(8)
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = _ios(e)

        # c for Genie style files skip the next 4 lines

        if nrtype[in_] > 2 and nrtype[in_] < 5:
            for i in range(1, 4 + 1):
                try:
                    units.skip_record(8)
                    ierr = 0
                except (FortranEOF, FortranIOError) as e:
                    ierr = _ios(e)

        # c for GUDRUN files, ignore any lines that have # in them

        if nrtype[in_] > 4:
            try:
                r = units.list_read(8)
                hashsearch = fstr(r.word(hashsearch), 10)
                ierr = 0
            except (FortranEOF, FortranIOError) as e:
                ierr = _ios(e)
            while findex(hashsearch, "#") > 0 and ierr == 0:
                try:
                    r = units.list_read(8)
                    hashsearch = fstr(r.word(hashsearch), 10)
                    ierr = 0
                except (FortranEOF, FortranIOError) as e:
                    ierr = _ios(e)

            # c If we got here there was no hash in the line read, so treat the line as good.

        # c Backspace the data file 1 line to the first good line

        units.backspace(8)

        ilen = findex(infile, ".")

        # c read in data to be fitted

        nqd = 0
        while nqd < mqd and ierr == 0:
            nqd1 = nqd + 1
            if nrtype[in_] == 1:
                try:
                    r = units.list_read(8)
                    qd[nqd1] = r.real(qd[nqd1])
                    data[nqd1] = r.real(data[nqd1])
                    ierr = 0
                except (FortranEOF, FortranIOError) as e:
                    ierr = _ios(e)
                # c     write(6,*) nqd,nqd1,qd(nqd1),data(nqd1)
                err[nqd1] = 0.0
            elif nrtype[in_] == 2 or nrtype[in_] >= 5:
                try:
                    r = units.list_read(8)
                    qd[nqd1] = r.real(qd[nqd1])
                    data[nqd1] = r.real(data[nqd1])
                    err[nqd1] = r.real(err[nqd1])
                    ierr = 0
                except (FortranEOF, FortranIOError) as e:
                    ierr = _ios(e)
            else:
                try:
                    n, qd[nqd1], data[nqd1], err[nqd1] = units.read_fmt(8, _F2224)
                    ierr = 0
                except (FortranEOF, FortranIOError) as e:
                    ierr = _ios(e)
            if ierr == 0:
                nqd = nqd1

        # c Ignore the last point which may be rubbish

        if nrtype[in_] > 2:
            nqd = nqd - 1

        # c for histogram files, calculate the correct q value

        if nrtype[in_] == 3 or nrtype[in_] == 5:
            q1 = qd[1]
            for i in range(1, nqd + 1):
                q2 = qd[i + 1]
                qd[i] = r4(0.5) * (q1 + q2)
                q1 = q2

        # c multiply data by density factor so that all data are per atom, and input Q values by qscale input factor

        for i in range(1, nqd + 1):
            # !                  qd(i)=qd(i)*qscaleinputfactor
            data[i] = data[i] * epsr16.rhofac[in_] * tweak[in_]

        units.write(6, fmt(_F352, nqd, substr(infile, 1, ilentot)))

        # c interpolate input data onto same scale as calculation

        nqmind[in_], nqmaxd[in_] = interp(nqd, nq, qd, q, data, datai, in_, nqmind[in_], nqmaxd[in_],  # noqa: F821
                                          mqd, mpts, mfile)
        nqmind[in_], nqmaxd[in_] = interp(nqd, nq, qd, q, err, errs, in_, nqmind[in_], nqmaxd[in_],  # noqa: F821
                                          mqd, mpts, mfile)

        # c prevent nqmax from being .ge. nq (This is for smoothing purposes in PTGNSQ)

        epsr16.nqmin[in_] = nqmind[in_]
        epsr16.nqmax[in_] = nqmaxd[in_]
        if epsr16.nqmax[in_] >= nq:
            epsr16.nqmax[in_] = nq - 1
        if epsr16.nqmin[in_] < nqminm:
            epsr16.nqmin[in_] = nqminm
        epsr16.nqpts[in_] = epsr16.nqmax[in_] - epsr16.nqmin[in_] + 1
        units.write(6, fmt(_F958, nqmind[in_], nqmaxd[in_], epsr16.nqmin[in_], epsr16.nqmax[in_],
                           epsr16.nqpts[in_], nq))

        # ! Set the Porod amplitude for this dataset

        if nporodq >= epsr16.nqmin[in_] and nporodq <= epsr16.nqmax[in_]:
            porodamp[in_] = datai[nporodq, in_] + porodadd[in_]
            units.write(6, fmt('("Porod amp. set to ",f10.5," for dataset ",i5)', porodamp[in_], in_))
        else:
            porodamp[in_] = 0.0
    # ! initialise the r-factor to write
    rfacwrite = f32(0)
    # c
    # c open file to save energy, pressure and R-factor values
    # c
    fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + ".erg", 256)
    units.write(6, fmt(_F912, substr(fname, 1, symbols.ilenfnamerun + 4)))
    openf(9, fname)  # noqa: F821
    # !      open(35,file='dump.erg',status='unknown')
    # c
    # c skip to the end of the file unless iinit.eq.1
    # c
    if iinit != 1:

        # c calculate the coefficients of the least squares parabola which goes through the R-factor
        # c and absolute energies

        fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + ".qdr", 256)
        (nparab, aparab, bparab, cparab, rfacwrite, aav, bav, lastqf, nextqf) = lsparabola1(  # noqa: F821
            9, fname, nparab, navmaxparab, epcfrac, aparab, bparab, cparab, rfacwrite, aav, bav, lastqf, nextqf)
        ierr = 0
        # !Skip to end of file dump file
        # !         do while(ierr.eq.0)
        # !            read(35,*,iostat=ierr)
        # !         end do
        # !         backspace(35)

    # c
    # c begin outer loop
    # c
    stdevo = f32(0.0)  # noqa: F841
    stde1o = f32(0.0)
    qalfco = f32(0.0)  # noqa: F841
    # c
    # c initialize control on size of fluctuations
    # c
    trat = f32(0.0)  # noqa: F841
    faco = f32(1.0)  # noqa: F841
    senmin = r4(1.0e20)  # noqa: F841

    # do 3000 iter=1,niter
    for iter in range(1, niter + 1):
        trace.context["iter"] = iter
        wtensm = f32(0.0)  # noqa: F841
        sumrf = f32(0.0)  # noqa: F841
        sument = f32(0.0)  # noqa: F841
        topsm = f32(0.0)  # noqa: F841
        botsm = f32(0.0)  # noqa: F841
        boasum = f32(0.0)  # noqa: F841
        qalsum = f32(0.0)  # noqa: F841
        qalsm2 = f32(0.0)  # noqa: F841
        nqals = 0  # noqa: F841
        # c
        # c reset s(q) accumulators
        # c
        if ireset == 1 or iter == 1:
            for j in range(1, 2 + 1):
                for imoltype in range(1, clustercom.nmoltype + 1):
                    clustercom.nmovei[j, imoltype] = 0
                    clustercom.ninter[j, imoltype] = 0
                    clustercom.nwhole[j, imoltype] = 0
                    clustercom.nmovem[j, imoltype] = 0

        # ! set the maximum separation of molecules. This corrected for the maximum size of any molecule

        sepmax = epsr16.xstar[epsr16.nptmax] + molecom.maxmolsize + r4(3.0) / clustercom.dmax
        units.write(6, fmt("(a,f10.5)", "epsr> Maximum separation for neighbourlist = ", sepmax * clustercom.dmax))

        # ! Initialise cell neighbour list

        cellneighbourlist.get_date_and_time("Start initialise_cellneighbourlist#", 0)
        cellneighbourlist.initialise_cellneighbourlist(nmolcellin, clustercom.xadd, clustercom.yadd, clustercom.zadd,
                                                       clustercom.nmol, sepmax)
        cellneighbourlist.get_date_and_time("Finish initialise_cellneighbourlist#", 1)
        # ! Ensure molecules are within box.
        for im in range(1, clustercom.nmol + 1):
            clustercom.xadd[im], clustercom.yadd[im], clustercom.zadd[im] = minimagec.minimagec(
                clustercom.xadd[im], clustercom.yadd[im], clustercom.zadd[im])
            if clustercom.tethered[im]:
                clustercom.phix[im], clustercom.phiy[im], clustercom.phiz[im] = minimagec.minimagec(
                    clustercom.phix[im], clustercom.phiy[im], clustercom.phiz[im])
        # ! Initialise various energies
        foldsv = f32(0.0)
        foldsb = f32(0.0)
        ftethered = f32(0.0)
        fcoordconstraints = f32(0.0)
        ierr, num_threads, foldsv, foldsb, ftethered, fcoordconstraints = ato_manipulation_rev4.fmolec1(
            ierr, num_threads, foldsv, foldsb, ftethered, fcoordconstraints)
        units.write(6, lw(ierr, num_threads, foldsv, foldsb, ftethered, fcoordconstraints,
                          *[threebody_mod.sume_triplet_dist[id] / f32(clustercom.nmol)
                            for id in range(1, threebody_mod.nthreeb + 1)]))
        units.write(35, fmt("(a,i4,a,2(1x,e13.6),a,4(1x,e13.6))", "#Refinement: ", iter, " Initial energies: ",
                            foldsv, foldsb, " Sizefactor et al.: ", molecom.sizefactor, molecom.stepmi,
                            molecom.stepri, pressreq))

        idone = 0  # noqa: F841
        # c
        # c set up theoretical R-factor and calculate renormalization on error bars
        # c
        theorf = f32(1.0)  # noqa: F841

        # ! Prevent intramolecular rotations if sizefactor > 1 and foldsv > abs(foldsb).
        # ! This to stop molecules adopting unwarranted configurations when intermolecular energy is very large
        # ! (e.g. after a change in sizefactor)
        # !      if(sizefactor.gt.1.0.and.foldsv.gt.abs(foldsb)) then
        # !         stepri=0.0
        # !      else
        # !         stepri=steprisave
        # !      end if
        units.write(6, fmt('("Steps: ",4(1x,e12.5))', molecom.stepmi, molecom.stepri, stepmisave, steprisave))
        # !      return

        # c exit if an error occurred in fmolec1

        if ierr != 0:
            return
        # c     write(6,*) nmol,foldsv*nmol,foldsb*nmol

        cellneighbourlist.get_date_and_time("Start molfit iterations#", 0)
        # do 2000 itimes=1,ntimes
        for itimes in range(1, ntimes + 1):
            # c     write(6,*) 'Got to molfit'
            # ! Initialise the intramolecular energy weights
            ato_manipulation_rev4.initialise_timer()
            foldsv, foldsb, epsrinput_par.nisum = molfit.molfit(foldsv, foldsb, epsrinput_par.nisum, intra,
                                                                rotfreq, ireset)
            # !         sumetriplets=0.0
            # !         do im=1,nmol
            # !            entest=ergold_triplet(im)
            # !            entest=en_mole_triplets(im,xadd,yadd,zadd,phix,phiy,phiz,xa,ya,za,mref,tethered)
            # !            if(im.eq.1) write(6,*) im,ergold_triplet(im),ergnew_triplet(im),entest
            # !            sumetriplets=sumetriplets+entest
            # !         end do
            # !         sumetriplets=sumetriplets/real(nmol)
            # !If nisum.eq.1, the molecular shake has occurred, so recalculate the energies
            if epsrinput_par.nisum == 1:
                ato_manipulation_rev4.adjust_intra_wts()
                units.write(6, fmt(_F1999, foldsv, foldsb, timer.sumtime1, timer.sumtime2))
                # c Ensure molecules stay within box.
                for im in range(1, clustercom.nmol + 1):
                    clustercom.xadd[im], clustercom.yadd[im], clustercom.zadd[im] = minimagec.minimagec(
                        clustercom.xadd[im], clustercom.yadd[im], clustercom.zadd[im])
                    if clustercom.tethered[im]:
                        clustercom.phix[im], clustercom.phiy[im], clustercom.phiz[im] = minimagec.minimagec(
                            clustercom.phix[im], clustercom.phiy[im], clustercom.phiz[im])
                foldsv = f32(0.0)
                foldsb = f32(0.0)
                ftethered = f32(0.0)
                fcoordconstraints = f32(0.0)
                ierr, num_threads, foldsv, foldsb, ftethered, fcoordconstraints = ato_manipulation_rev4.fmolec1(
                    ierr, num_threads, foldsv, foldsb, ftethered, fcoordconstraints)
            units.write(6, fmt(_F1999, foldsv, foldsb, timer.sumtime1, timer.sumtime2))
            # c     write(6,*) nmol,foldsv*nmol,foldsb*nmol
            if epsrinput_par.nisum == 1 and molecom.stepmi > 0.0:
                units.write(35, fmt("(a)", "#Intramolecular shake> "))
            if (epsrinput_par.nisum != 1
                    and (epsrinput_par.nisum - 1) == rotfreq * idiv(epsrinput_par.nisum - 1, rotfreq)
                    and molecom.stepri > 0.0):
                units.write(35, fmt("(a)", "#Intramolecular rotation> "))
            units.write(35, fmt("(i6,2(1x,e13.6),2(1x,i16))", epsrinput_par.nisum, foldsv, foldsb,
                                timer.sumtime1, timer.sumtime2))
        # 2000 continue
        cellneighbourlist.get_date_and_time("Finish molfit iterations#", 1)
        # ! Ensure molecules are within box.
        # !      do im=1,nmol
        # !         call minimagec(xadd(im),yadd(im),zadd(im))
        # !         if(tethered(im)) call minimagec(phix(im),phiy(im),phiz(im))
        # !      end do
        # ! Check final energies are consistent
        # ! Initialise various energies
        # !      foldsv=0.0
        # !      foldsb=0.0
        # !      ftethered=0.0
        # !      call fmolec1(ierr,num_threads,foldsv,foldsb,ftethered)
        # !      write(35,'(a,i4,a,2(1x,e13.6),a,3(1x,e13.6))') '#Refinement: ',iter,' Final energies: ',foldsv,foldsb

        # ! set the maximum separation of atoms

        # !      sepmax=xstar(npt1)

        # ! Initialise cell neighbour list

        # !      call get_date_and_time('Start initialise_cellneighbourlist#',0)
        # !      call initialise_cellneighbourlist(nmolcell,xadd,yadd,zadd,nmol,sepmax)
        # !      call get_date_and_time('Finish initialise_cellneighbourlist#',1)

        # ! *** The following section is needed if we are to refine lattice parameters

        oldfit = f32(0.0)
        alatold = lattice.alat
        blatold = lattice.blat
        clatold = lattice.clat
        alphaold, betaold, gamold = auxiliary_routines.convertphibthetacphic(
            lattice.phib, lattice.thetac, lattice.phic, alphaold, betaold, gamold)
        alatnew = alatold
        blatnew = blatold
        clatnew = clatold
        alphanew = alphaold
        betanew = betaold
        gamnew = gamold
        ilatrefine = 0
        while ilatrefine < nlatrefine + 1:
            cellneighbourlist.get_date_and_time("epsr> Start atogen#", 0)
            ngrsamples = ato_manipulation_rev4.atogen(wav, ngrsamples)
            cellneighbourlist.get_date_and_time("epsr> Finished atogen#", 1)
            units.write(6, lw("Got to gensq"))
            gensq(nq, q, qwind, wav)  # noqa: F821
            # !Get the Bragg intensities if required
            if nbragg > 0:
                units.write(6, lw(qstep, qbraggmax, qbraggwidth, multa, multb, multc, nbragg, nbraggup))
                if qbragg is not None:
                    qbragg = None
                if intenbragg is not None:
                    intenbragg = None
                if intenconv is not None:
                    intenconv = None
                qbragg = real_array(nbraggup)
                intenbragg = real_array(molecom.npairs, nbraggup)
                intenconv = real_array(molecom.npairs, nbragg)
                for ib in range(1, nbraggup + 1):
                    qbragg[ib] = qstep * f32(ib - 1)
                lattice.twopi = r4(2.0) * pi
                # !Output file name
                outputfile = fstr(substr(symbols.fname_run, 1, symbols.ilenfname_run), 256)
                ncharoutput = symbols.ilenfname_run
                fname = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                             + substr(outputfile, 1, ncharoutput) + ".inten", 256)
                units.open(25, fname, status="unknown")
                # !            if(rminbragg.le.0.0.or.rminbragg.gt.x(npt1)) rminbragg=x(npt1)
                ato_read_and_save.unnormalise_coordinates()
                ato_manipulation_rev4.bragginten(qstep, qbraggmax, qbraggwidth, fwhmq, multa, multb, multc,
                                                 usea, useb, usec, hklqmin, nbraggup, qbragg, intenbragg)
                ato_read_and_save.normalise_coordinates()
                # !Subtract the single molecule scattering and normalise to the number of unit cells
                units.write(25, lw("# ", molecom.npairs, nbragg, nbraggup))
                multfactor = f32(multa * multb * multc)
                for ib in range(1, nbraggup + 1):
                    for ip in range(1, molecom.npairs + 1):
                        intenbragg[ip, ib] = intenbragg[ip, ib] / multfactor - qwind[1, ib, ip]
                for ib in range(1, nbragg + 1):
                    for ip in range(1, molecom.npairs + 1):
                        intenconv[ip, ib] = 0.0
                # !Apply Debye-Waller factor to Bragg intensities
                if debyewaller > 0.0:
                    for ib in range(1, nbraggup + 1):
                        qbraggfactor = qbragg[ib] * debyewaller
                        qbraggfactor = r4(0.5) * qbraggfactor * qbraggfactor
                        if qbraggfactor < 40.0:
                            qbraggfactor = expf(-qbraggfactor)
                        else:
                            qbraggfactor = f32(0.0)
                        for ip in range(1, molecom.npairs + 1):
                            intenbragg[ip, ib] = qbraggfactor * intenbragg[ip, ib]
                # !Step through Q values and add broadened intensities to convolution function
                gfactor = abs(epsrinput_par.rminbragg) / (r4(2.0) * (powf(lattice.twopi, r4(1.5))))
                gfactor0 = (pow3(abs(epsrinput_par.rminbragg))) / ((powf(lattice.twopi, r4(1.5))))
                for ib in range(1, nbraggup + 1):
                    sumfactor = r4(4.0) * pi * qbragg[ib] * qstep
                    for it in range(1, nbragg + 1):
                        sum = f32(0.0)
                        # !Gaussian broadening function
                        if qbragg[ib] > 0.0 and qbragg[it] > 0.0:
                            diff = (qbragg[ib] - qbragg[it]) * epsrinput_par.rminbragg
                            epsr16.add = (qbragg[ib] + qbragg[it]) * epsrinput_par.rminbragg
                            if abs(diff) < 8.0:
                                sum = expf(-(r4(0.5) * diff * diff))
                            if abs(epsr16.add) < 8.0:
                                sum = sum - expf(-(r4(0.5) * epsr16.add * epsr16.add))
                            sum = sum * gfactor / (qbragg[it])
                        else:
                            epsr16.add = qbragg[ib] * epsrinput_par.rminbragg
                            diff = qbragg[it] * epsrinput_par.rminbragg
                            sum = r4(0.5) * (diff * diff + epsr16.add * epsr16.add)
                            if sum < 32.0:
                                sum = gfactor0 * expf(-sum)
                            else:
                                sum = f32(0.0)
                            sum = sum * qbragg[ib]
                        sum = sum * sumfactor
                        for ip in range(1, molecom.npairs + 1):
                            intenconv[ip, it] = intenconv[ip, it] + (qwind[2, ib, ip] - intenbragg[ip, ib]) * sum
                # !Subtract the convolution function
                for ip in range(1, molecom.npairs + 1):
                    for ib in range(1, nbragg + 1):
                        # !               intenbragg(ip,ib)=intenbragg(ip,ib)-intenconv(ip,ib)
                        # !               qwind(2,ib,ip)=qwind(2,ib,ip)+intenbragg(ip,ib)
                        qwind[2, ib, ip] = intenbragg[ip, ib] + intenconv[ip, ib]
                        # !               qwind(2,ib,ip)=intenbragg(ip,ib)-qwind(1,ib,ip)
                for ib in range(1, nbragg + 1):
                    vals = []
                    for ip in range(1, molecom.npairs + 1):
                        vals += [qwind[2, ib, ip], intenbragg[ip, ib], intenconv[ip, ib]]
                    units.write(25, lw(qbragg[ib], *vals))
                units.close(25)
                # ! Output the diffuse scattering if requested
                epsrinput_par.nsums = ato_manipulation_rev4.sharm_diffuse(qstep, nq, epsrinput_par.nsums,
                                                                          epsr16.nfile)
            if refinelattice and ilatrefine <= nlatrefine:
                # ! a) Now calculate fit to data if refining lattice parameters
                rfacraw = f32(0.0)
                for in_ in range(1, epsr16.nfile + 1):
                    rfacdraw = f32(0.0)
                    # ! Determine whether this is an X-ray dataset
                    ix = int(weights.xfile[in_])
                    for i in range(1, nq + 1):
                        soqism[i, in_] = 0.0
                        backsc[i, in_] = 0.0
                        sum = f32(0.0)
                        sumraw = f32(0.0)
                        if ix > 0:
                            for ic in range(1, molecom.npairs + 1):
                                sumraw = f32(sumraw + x_rayf0.weightx[i, ic, ix] * (qwind[1, i, ic] + qwind[2, i, ic]))
                            # ! Add on the X-ray background due to mismatched form factors.
                            sum = f32(sum + x_rayf0.backx[i, ix])
                            sumraw = f32(sumraw + x_rayf0.backx[i, ix])
                        else:
                            for ic in range(1, molecom.npairs + 1):
                                sumraw = f32(sumraw + weights.ciwt[in_, ic] * qwind[1, i, ic]
                                             + weights.cwt[in_, ic] * qwind[2, i, ic])
                        if i >= nqmind[in_] and i <= nqmaxd[in_]:
                            diffraw = datai[i, in_] - sumraw
                        else:
                            diffraw = f32(0.0)
                        if i >= epsr16.nqmin[in_] and i <= epsr16.nqmax[in_]:
                            rfacdraw = rfacdraw + diffraw * diffraw
                    if epsr16.nqpts[in_] > 0:
                        rfacdraw = rfacdraw / f32(epsr16.nqpts[in_])
                    rfacraw = rfacraw + rfacdraw
                rfacraw = rfacraw / f32(epsr16.nfile)
                newfit = rfacraw
                if ilatrefine == 0:
                    units.write(6, fmt("(a)", "refinelattice> Starting values:-"))
                    units.write(6, fmt("(a,6(1x,e12.6))", "refinelattice> Initial lattice pars.: ",
                                       alatold, blatold, clatold, alphaold, betaold, gamold))
                    units.write(6, fmt("(a,2(1x,e12.6))", "refinelattice> Initial fit: ", newfit))
                else:
                    units.write(6, fmt("(a,i5)", "refinelattice> End of iteration: ", ilatrefine))
                    units.write(6, fmt("(a,6(1x,e12.6))", "refinelattice> Old lattice pars.: ",
                                       alatold, blatold, clatold, alphaold, betaold, gamold))
                    units.write(6, fmt("(a,6(1x,e12.6))", "refinelattice> New lattice pars.: ",
                                       alatnew, blatnew, clatnew, alphanew, betanew, gamnew))
                    units.write(6, fmt("(a,2(1x,e12.6))", "refinelattice> New fit and old fit: ", newfit, oldfit))
                # ! b) If not the first iteration, compare newfit with oldfit and if better, or oldfit = 0, accept new parameters. Otherwise restore old parameters.
                if ilatrefine == 0:
                    oldfit = newfit
                elif newfit <= oldfit:  # Fit improved so save lattice parameters
                    oldfit = newfit
                    alatold = alatnew
                    blatold = blatnew
                    clatold = clatnew
                    alphaold = alphanew
                    betaold = betanew
                    gamold = gamnew
                    if ilatrefine == nlatrefine:
                        # If this is the last iteration, no need to go back and recalculate S(Q).
                        ilatrefine = ilatrefine + 1
                else:  # Restore original lattice parameters
                    alatfac = alatold / alatnew
                    blatfac = blatold / blatnew
                    clatfac = clatold / clatnew
                    newphib, newthetac, newphic = auxiliary_routines.convertalphabetagamma(
                        alphaold, betaold, gamold, newphib, newthetac, newphic)
                    ato_read_and_save.unnormalise_coordinates()
                    ato_manipulation_rev4.change_boxsize(alatfac, blatfac, clatfac, newphib, newthetac, newphic)
                    ato_read_and_save.normalise_coordinates()
                # ! c) Choose new lattice parameters and try a new box shape/size - skip this step if the last iteration
                if ilatrefine < nlatrefine:
                    alatnew = alatold + alatstep * (r4(2.0) * auxiliary_routines.ran1() - r4(1.0)) / clustercom.dmax
                    alatfac = alatnew / alatold
                    if lattice.cubic:
                        blatnew = alatnew
                        clatnew = alatnew
                        blatfac = alatfac
                        clatfac = alatfac
                        alphanew = r4(90.0)
                        betanew = r4(90.0)
                        gamnew = r4(90.0)
                    else:
                        blatnew = blatold + blatstep * (r4(2.0) * auxiliary_routines.ran1() - r4(1.0)) / clustercom.dmax
                        blatfac = blatnew / blatold
                        clatnew = clatold + clatstep * (r4(2.0) * auxiliary_routines.ran1() - r4(1.0)) / clustercom.dmax
                        clatfac = clatnew / clatold
                        alphanew = alphaold + alphastep * (r4(2.0) * auxiliary_routines.ran1() - r4(1.0))
                        betanew = betaold + betastep * (r4(2.0) * auxiliary_routines.ran1() - r4(1.0))
                        gamnew = gamold + gamstep * (r4(2.0) * auxiliary_routines.ran1() - r4(1.0))
                    newphib, newthetac, newphic = auxiliary_routines.convertalphabetagamma(
                        alphanew, betanew, gamnew, newphib, newthetac, newphic)
                    # ! unnormalise the coordinates before changing box size and shape
                    ato_read_and_save.unnormalise_coordinates()
                    # ! Need to ensure density does not change
                    rhoold = f32(molecom.natomr) / (lattice.alatx * lattice.blaty * lattice.clatz)
                    ato_manipulation_rev4.change_boxsize(alatfac, blatfac, clatfac, newphib, newthetac, newphic)
                    rhonew = f32(molecom.natomr) / (lattice.alatx * lattice.blaty * lattice.clatz)
                    vol = rhonew / rhoold
                    dmaxfac = powf(vol, r4(1.0) / r4(3.0))
                    alatfac = dmaxfac
                    blatfac = dmaxfac
                    clatfac = dmaxfac
                    newphib = lattice.phib
                    newthetac = lattice.thetac
                    newphic = lattice.phic
                    ato_manipulation_rev4.change_boxsize(alatfac, blatfac, clatfac, newphib, newthetac, newphic)
                    ato_read_and_save.normalise_coordinates()
                    alatnew = alatnew * alatfac
                    blatnew = blatnew * blatfac
                    clatnew = clatnew * blatfac  # NB: blatfac (= clatfac here), as in the Fortran
            ilatrefine = ilatrefine + 1
            # ! e) Go back to *** (atogen) and recalculate fit to data

        # ! **** End of lattice parameter refinement.

        units.write(6, lw("Got to entropc"))
        oiseb = entrpc(oiseb, wav)  # noqa: F821
        nb = 12  # noqa: F841
        nb2 = 15  # noqa: F841
        # c
        # c sum distributions
        # c
        epsrinput_par.nsum = epsrinput_par.nsum + 1
        if epsrinput_par.nsums >= 0:
            sumvalues.fc1 = float(r4(1.0) / f32(epsrinput_par.nsum))
            sumvalues.fc2 = float(f32(epsrinput_par.nsum - 1)) * sumvalues.fc1
        else:
            sumvalues.fc1 = float(r4(1.0) / f32(epsrinput_par.nesum + 1))
            sumvalues.fc2 = float(f32(epsrinput_par.nesum)) * sumvalues.fc1

        # ! Update uniform atom distribution if needed
        # !      if(nisum.ge.intra) then
        # !            ycoh(1)=0.0
        # !            atry=5.0/27.0
        # !            call wholecoord(alatx,dmax)
        # !            call wholecoord(blaty,dmax)
        # !            call wholecoord(clatz,dmax)
        # !            call uniform(ycoh,cellst,rho,cell,mpz)
        # !            call fractcoord(alatx,dmax)
        # !            call fractcoord(blaty,dmax)
        # !            call fractcoord(clatz,dmax)
        # !            nisum=0
        # !      end if
        # c
        # c sum Q-space distributions
        # c
        units.write(6, lw("Got to sumsq"))
        sumsq(nq, qwind)  # noqa: F821
        # c
        # c sum r-space distributions
        # c
        units.write(6, lw("Got to sumgr"))
        sumgr(nz, wav)  # noqa: F821
        # c
        # c calculate entropy of average distribution
        # c
        units.write(6, lw("Got to entrpb"))
        sentb = entrpb(sentb, epsrinput_par.nsum)  # noqa: F821
        # c     endif
        stde1o = oiseb / epsr16.rho  # noqa: F841
        # c
        # c calculate intermolecular energy and pressure
        # c
        # !      if(nthreeb.gt.0) write(6,*) ' Standard deviation (rad^2) for each triplet: ',(sume_triplet_dist(id)/real(nmol),id=1,nthreeb)
        cellneighbourlist.get_date_and_time("epsr> Start press#", 0)
        energd, pressd, energabs, energthreeb, pressthreeb = ato_manipulation_rev4.press(
            energd, pressd, energf, pressf, energabs, energthreeb, pressthreeb)
        # !      absolute_energy=energabs
        cellneighbourlist.get_date_and_time("epsr> Finished press#", 1)
        # !      if(nthreeb.gt.0) write(6,*) ' Standard deviation (rad^2) for each triplet: ',(sume_triplet_dist(id)/real(nmol),id=1,nthreeb)
        # c
        # c correct energy for static polarisation
        # c
        energd = energd + epsr16.upol
        # c
        # c reset accumulators if necessary
        # c
        if epsrinput_par.nesum == 0:
            epsrinput_par.ensum = f32(0.0)
            epsrinput_par.ensum2 = f32(0.0)
            epsrinput_par.prsum = f32(0.0)
            epsrinput_par.prsum2 = f32(0.0)
        # c
        # c only sum energies if nisum .ge. nisumm  - this is to prevent incorrect values which
        # c occur after an intramolecular disordering move from biassing the overall
        # c energy
        # c
        if epsrinput_par.nisum >= nisumm:
            epsrinput_par.nesum = epsrinput_par.nesum + 1
            if epsrinput_par.nesum > epsrinput_par.nesumt:
                epsrinput_par.nesum = epsrinput_par.nesumt
            nesum1 = epsrinput_par.nesum - 1
            epsrinput_par.ensum = (f32(nesum1) * epsrinput_par.ensum + energd) / f32(epsrinput_par.nesum)
            ensumc = energd - epsrinput_par.ensum
            epsrinput_par.ensum2 = (f32(nesum1) * epsrinput_par.ensum2 + ensumc * ensumc) / f32(epsrinput_par.nesum)
            epsrinput_par.prsum = (f32(nesum1) * epsrinput_par.prsum + pressd) / f32(epsrinput_par.nesum)
            prsumc = pressd - epsrinput_par.prsum
            epsrinput_par.prsum2 = (f32(nesum1) * epsrinput_par.prsum2 + prsumc * prsumc) / f32(epsrinput_par.nesum)
        enav = epsrinput_par.ensum
        enav2 = epsrinput_par.ensum2
        if enav2 > 0.0:
            enav2 = sqrtf(enav2)
        else:
            enav2 = f32(0.0)
        prav = epsrinput_par.prsum
        prav2 = epsrinput_par.prsum2
        if prav2 > 0.0:
            prav2 = sqrtf(prav2)
        else:
            prav2 = f32(0.0)
        # c
        # c sum the contributions from the empirical potential
        # c
        epsum = f32(0)
        epsum1 = f32(0.0)
        for in_ in range(1, epsr16.nfile + 1):
            epsum = epsum + abs(energf[in_])
            epsum1 = epsum1 + energf[in_]
        # c
        # c set up intermolecular data
        # c
        for in_ in range(1, epsr16.nfile + 1):

            ix = int(weights.xfile[in_])
            if ix > 0:

                # c If X-ray data, need to apply the Q-dependent weights

                for iq in range(1, nq + 1):
                    # c
                    # c sum the intramolecular terms and subtract them from the diffraction data
                    # c
                    if iq >= nqmind[in_] and iq <= nqmaxd[in_]:
                        sum = f32(0.0)
                        for ip in range(1, molecom.npairs + 1):
                            sum = f32(sum + x_rayf0.weightx[iq, ip, ix] * epsr16.soqintra[iq, ip])
                        datat[iq, in_] = datai[iq, in_] - sum
                    else:
                        datat[iq, in_] = 0.0

                # !Remember to subtract the intramolecular term from the zero limit as well!
                szerostemp[in_] = 0.0
                if szeros[in_] != 0.0:
                    sum = f32(0.0)
                    for ip in range(1, molecom.npairs + 1):
                        sum = f32(sum + x_rayf0.weightx[1, ip, ix] * epsr16.soqintra[1, ip])
                    szerostemp[in_] = szeros[in_] - sum

            else:

                for iq in range(1, nq + 1):

                    if iq >= nqmind[in_] and iq <= nqmaxd[in_]:
                        # c
                        # c sum the intramolecular terms and subtract them from the diffraction data
                        # c
                        sum = f32(0.0)
                        for ip in range(1, molecom.npairs + 1):
                            sum = f32(sum + weights.ciwt[in_, ip] * epsr16.soqintra[iq, ip])
                        datat[iq, in_] = datai[iq, in_] - sum
                    else:
                        datat[iq, in_] = 0.0

                # !Remember to subtract the intramolecular term from the zero limit as well!
                szerostemp[in_] = 0.0
                if szeros[in_] != 0.0:
                    sum = f32(0.0)
                    for ip in range(1, molecom.npairs + 1):
                        sum = f32(sum + weights.ciwt[in_, ip] * epsr16.soqintra[1, ip])
                    szerostemp[in_] = szeros[in_] - sum

        # c Set up difference data - simulation and subtract a background from data if needed

        for id in range(1, epsr16.nfile + 1):

            # c Form the simulated fit for this dataset

            # c If X-ray data, need to apply the Q-dependent weights

            ix = int(weights.xfile[id])

            # c Set up the gradient for extrapolating to Q=0

            if szeros[id] != 0.0 and epsr16.nqmin[id] > 1:

                grads0 = (datat[epsr16.nqmin[id], id] - szerostemp[id]) / (q[epsr16.nqmin[id]] - q[1])
                conss0 = szerostemp[id]

            else:

                grads0 = f32(0.0)
                conss0 = datat[epsr16.nqmin[id], id]

            for i in range(1, nq + 1):
                # c
                # c Sum the intermolecular partials for this dataset at this q value
                # c
                sum = f32(0.0)
                if ix > 0:
                    for ic in range(1, molecom.npairs + 1):
                        sum = f32(sum + x_rayf0.weightx[i, ic, ix] * epsr16.soqinter[i, ic])
                else:
                    for ic in range(1, molecom.npairs + 1):
                        sum = f32(sum + weights.cwt[id, ic] * epsr16.soqinter[i, ic])

                # c Calculate the difference between data and fit. For Q-values below qmin we take the first
                # c valid data point as the data below Qmin, unless szeros is non-zero.

                if i < epsr16.nqmin[id]:
                    if szeros[id] != 0.0:
                        dts[i] = conss0 + grads0 * (q[i] - q[1]) - sum
                    else:
                        dts[i] = 0.0
                elif i >= epsr16.nqmin[id] and i <= epsr16.nqmax[id]:
                    dts[i] = datat[i, id] - sum
                    # !Add in the Porod scattering if required
                    if porodamp[id] > 0.0:
                        dts[i] = dts[i] + porodamp[id] * powf(porodq / q[i], porodpwr)
                else:
                    dts[i] = 0.0

            # c
            # c form a smooth background to data if requested and subtract it
            # c
            if nsmoop > 0:
                ftback(id, rshmin[id], q, dts, dti)  # noqa: F821
                for i in range(1, epsr16.nqmax[id] + 1):
                    if i >= epsr16.nqmin[id]:
                        datat[i, id] = datat[i, id] - dts[i]

        # c Set up the current simulated partial structure factors

        for ip in range(1, molecom.npairs + 1):

            for iq in range(1, nq + 1):

                soqintertemp[iq, ip] = qwind[2, iq, ip]
                # c           soqintertemp(iq,ip)=soqinter(iq,ip)

        # c
        # c update potential
        # c
        wtentt = epsr16.rho
        # c     subroutine ptgnsq(nsum,q,wav,data,qwind
        # c     *,rminpt,rmaxpt,szerostemp,overallfac,gfac
        # c     *,rminp,potnoi,wtentp
        # c     *,rbroad,expecf,nsmoop)
        units.write(6, lw("Got to ptgnsq"))
        # ! Only calculate the new potential if the quality factor got worse at the previous iteration, or automatic adjustment of ereq is off, or
        # ! we are initial phase of automatic adjustment of ereq.
        # !      if(pressreq.gt.0.0.and.(ereqstep.eq.0.0.or.nextqf.gt.lastqf.or.nparab.lt.nparabstart)) call ptgnsq(nesum,nq,q,wav,datat,soqintertemp
        if pressreq > 0.0:
            ptgnsq(epsrinput_par.nesum, nq, q, wav, datat, soqintertemp,  # noqa: F821
                   epsrinput_par.rminpt, epsrinput_par.rmaxpt, szerostemp, overallfac, gfac,
                   rshmin, potnoi, wtentt,
                   epsrinput_par.rbroad, epsrinput_par.expecf, nsmoop)
        # !Regenerate the repulsive potential if potfac > 0
        repulsivept(potential_par.potfac)  # noqa: F821
        sumefa = f32(0.0)  # noqa: F841
        nsig = 0  # noqa: F841
        sum1 = f32(0.0)
        sum2 = f32(0.0)
        for imoltype in range(1, clustercom.nmoltype + 1):
            sum1 = sum1 + f32(clustercom.nwhole[1, imoltype]) + f32(clustercom.nmovem[1, imoltype])
            sum2 = sum2 + f32(clustercom.nwhole[2, imoltype]) + f32(clustercom.nmovem[2, imoltype])
        if sum1 > 0.0:
            rfrac = sum2 / sum1  # noqa: F841
        else:
            rfrac = f32(0.0)  # noqa: F841

        rhoout = epsr16.rho * pljfac * pljfac * pljfac  # noqa: F841
        # c
        # c for last attempt values are written to enable continuation
        # c of present attempt
        # c calculate chisq at same time
        # c
        eminmax = f32(0.0)
        feminmax = f32(0.0)  # noqa: F841
        # Python shell: Ctrl+C during epsr makes this iteration the last one, so the results
        # below are written and every file is saved consistently (as if niter had been iter).
        if niter > iter and _stop_requested():
            niter = iter
            units.write(6, lw("runepsr17> Ctrl+C: this is the last iteration; writing results and saving files."))
        if iter == niter:
            # c
            # c calculate the potential energies based on the mean g(r) and
            # c current potentials
            # c
            sumpe(nz, z, symbols.fnamerun, symbols.ilenfnamerun)  # noqa: F821

            # c
            # c invert the data according to the inversion matrix
            # c
            for ic in range(1, molecom.npairs + 1):
                # c
                # c set the column title for this pair of atoms
                # c
                atom1 = fstr(potential_par.tlab1[potential_par.irefpcof[ic]], 4)
                atom2 = fstr(potential_par.tlab2[potential_par.irefpcof[ic]], 4)
                ilent1 = findex_padded(atom1, " ", 4) - 1
                ilent2 = findex_padded(atom2, " ", 4) - 1
                titlepsf[ic] = fstr("   " + substr(atom1, 1, ilent1) + "-" + substr(atom2, 1, ilent2), 29)
                # c     write(6,*) ilent1,ilent2,tlab1(ic)(1:3),tlab2(ic)(1:3),titlepsf(ic)(1:20)

                # c set the assymptotic limit for g(r)

                grlim[ic] = 1.0

                # c Minimum and maximum Q values for transform of partial structure factor data

                nqminp[ic] = 1
                nqmaxp[ic] = nq

                iqinv = 2
                for iq in range(1, nq + 1):
                    datt = f32(0.0)
                    if weights.nxfile > 0:

                        # c Ensure the current q values of the inversion matrix correspond to this Q value

                        while weights.qinv[iqinv] < iq and iqinv < weights.nqinv:
                            iqinv = iqinv + 1

                        # c Interpolation ratio

                        ratinter = (f32(iq - int(weights.qinv[iqinv - 1]))
                                    / f32(int(weights.qinv[iqinv]) - int(weights.qinv[iqinv - 1])))
                        ratinter1 = r4(1.0) - ratinter
                        for in_ in range(1, epsr16.nfile + 1):

                            # c Interpolate (linearly) inversion matrix to this Q value

                            cwttemp = f32(float(ratinter1) * weights.cwttq[iqinv - 1, ic, in_]
                                          + float(ratinter) * weights.cwttq[iqinv, ic, in_])

                            if iq < epsr16.nqmin[in_]:
                                # c
                                # c fill in the data below nqmin
                                # c
                                term = datat[epsr16.nqmin[in_], in_] * cwttemp
                            else:
                                term = datat[iq, in_] * cwttemp
                            if potential_par.wsum[in_] <= 0.0:
                                term = f32(0.0)
                            datt = datt + term
                        # c
                        # c add on the extra terms from the simulation
                        # c
                        if potential_par.feedback < 1.0:
                            for in_ in range(1, molecom.npairs + 1):
                                inref = in_ + epsr16.nfile
                                cwttemp = f32(float(ratinter1) * weights.cwttq[iqinv - 1, ic, inref]
                                              + float(ratinter) * weights.cwttq[iqinv, ic, inref])
                                datt = datt + epsr16.soqinter[iq, in_] * cwttemp
                    else:
                        ratinter1 = f32(1.0)
                        for in_ in range(1, epsr16.nfile + 1):

                            # c Interpolate (linearly) inversion matrix to this Q value

                            cwttemp = f32(float(ratinter1) * weights.cwttq[1, ic, in_])
                            if iq < epsr16.nqmin[in_]:
                                # c
                                # c fill in the data below nqmin
                                # c
                                term = datat[epsr16.nqmin[in_], in_] * cwttemp
                            else:
                                term = datat[iq, in_] * cwttemp
                            if potential_par.wsum[in_] <= 0.0:
                                term = f32(0.0)
                            datt = datt + term
                        # c
                        # c add on the extra terms from the simulation
                        # c
                        if potential_par.feedback < 1.0:
                            for in_ in range(1, molecom.npairs + 1):
                                inref = in_ + epsr16.nfile
                                cwttemp = f32(float(ratinter1) * weights.cwttq[1, ic, inref])
                                datt = datt + epsr16.soqinter[iq, in_] * cwttemp
                    datap[iq, ic] = datt
            # c
            # c set up titles for data output files
            # c
            for in_ in range(1, epsr16.nfile + 1):
                ilent2 = findex(infili[in_], ".") - 1
                ilent1 = ilent2 - 28
                if ilent1 < 1:
                    ilent1 = 1
                titledat[in_] = fstr(substr(infili[in_], ilent1, ilent2), 29)
            # c
            # c axis titles
            # c
            titler = "   r"  # character*14
            titleq = "   Q"  # character*14
            # c
            # c perform direct Fourier transform on partial structure factors and inverted data ready for output
            # c
            ato_manipulation_rev4.stog(molecom_mod.rhocg, q, datap, z, datagr, nqminp, nqmaxp, nz, molecom.npairs,
                                       mpts, mpz, mpairs, grlim)
            chisq = f32(0.0)
            for ic in range(1, molecom.npairs + 1):
                rfac[ic] = 0.0
                for i in range(1, nq + 1):
                    a = epsr16.soqinter[i, ic]  # noqa: F841
                    b = epsr16.sq2inter[i, ic]  # noqa: F841
                    if datap[i, ic] != 0.0:
                        dif = datap[i, ic] - epsr16.soqinter[i, ic]
                        errsav[i, ic] = dif
                        rfac[ic] = rfac[ic] + dif * dif
                    else:
                        errsav[i, ic] = 0.0
                rfac[ic] = rfac[ic] / f32(nq)
                chisq = chisq + rfac[ic]
            infile = infili[epsr16.nfile]
            ilen = findex(infile, ".")
            # c
            # c write fit to partial structure factors
            # c
            ext1 = fstr(fmt(_F331, ".f"), 2)
            fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(ext1, 1, 2), 256)
            fname = auxiliary_routines.w_int(nq, fname, molecom.npairs, q, epsr16.soqinter, epsr16.sq2inter,
                                             titleq, titlepsf, mpts, mpairs)

            ext1 = fstr(fmt(_F331, ".s"), 2)
            fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(ext1, 1, 2), 256)
            fname = auxiliary_routines.w_int(nq, fname, molecom.npairs, q, epsr16.soqintra, epsr16.sq2intra,
                                             titleq, titlepsf, mpts, mpairs)
            # c
            # c write difference partial structure factors minus fit
            # c
            ext1 = fstr(fmt(_F331, ".d"), 2)
            fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(ext1, 1, 2), 256)
            # c
            # c Output will be written between Qmin and Qmax of input datafiles
            # c
            lpt1 = nq
            lpt2 = 1
            for in_ in range(1, epsr16.nfile + 1):
                lpt1 = min(lpt1, int(nqmind[in_]))
                lpt2 = max(lpt2, int(nqmaxd[in_]))
            fname = auxiliary_routines.w_intr(lpt1, lpt2, fname, molecom.npairs, q, errsav, epsr16.sq2inter,
                                              titleq, titlepsf, mpts, mpairs)
            # c
            # c write estimated partial structure factors
            # c
            ext1 = fstr(fmt(_F331, ".q"), 2)
            fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(ext1, 1, 2), 256)
            fname = auxiliary_routines.w_intr(lpt1, lpt2, fname, molecom.npairs, q, datap, epsr16.sq2inter,
                                              titleq, titlepsf, mpts, mpairs)
            # c
            # c write out the input data after interpolation by the program - the so-called
            # c "total" data and after background subtraction
            # c
            # c
            # c set up intermolecular data
            # c
            for in_ in range(1, epsr16.nfile + 1):

                # c Identify any X-ray data sets

                ix = int(weights.xfile[in_])

                for i in range(1, nq + 1):
                    # c
                    # c sum the intramolecular terms and subtract them from the diffraction data
                    # c
                    sum = f32(0.0)
                    if i >= nqmind[in_] and i <= nqmaxd[in_]:
                        if ix > 0:
                            for ic in range(1, molecom.npairs + 1):
                                sum = f32(sum + x_rayf0.weightx[i, ic, ix] * epsr16.soqintra[i, ic])
                        else:
                            for ic in range(1, molecom.npairs + 1):
                                sum = f32(sum + weights.ciwt[in_, ic] * epsr16.soqintra[i, ic])
                    datat[i, in_] = datat[i, in_] + sum
            ext1 = fstr(fmt(_F331, ".t"), 2)
            fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(ext1, 1, 2), 256)
            fname = auxiliary_routines.w_intr(lpt1, lpt2, fname, epsr16.nfile, q, datat, errs,
                                              titleq, titledat, mpts, mfile)
            # c
            # c calculate the fit to the input data
            # c
            rfactot = f32(0.0)
            rfacraw = f32(0.0)
            for in_ in range(1, epsr16.nfile + 1):
                rfacd[in_] = 0.0
                rfacdraw = f32(0.0)

                # c Determine whether this is an X-ray dataset

                ix = int(weights.xfile[in_])

                for i in range(1, nq + 1):
                    soqism[i, in_] = 0.0
                    backsc[i, in_] = 0.0
                    sum = f32(0.0)
                    sumraw = f32(0.0)

                    if ix > 0:
                        for ic in range(1, molecom.npairs + 1):
                            sum = f32(sum + x_rayf0.weightx[i, ic, ix]
                                      * (epsr16.soqintra[i, ic] + epsr16.soqinter[i, ic]))
                            sumraw = f32(sumraw + x_rayf0.weightx[i, ic, ix]
                                         * (qwind[1, i, ic] + qwind[2, i, ic]))

                        # c Add on the X-ray background due to mismatched form factors.

                        sum = f32(sum + x_rayf0.backx[i, ix])
                        sumraw = f32(sumraw + x_rayf0.backx[i, ix])
                    else:
                        for ic in range(1, molecom.npairs + 1):
                            sum = f32(sum
                                      + weights.ciwt[in_, ic] * epsr16.soqintra[i, ic]
                                      + weights.cwt[in_, ic] * epsr16.soqinter[i, ic])
                            sumraw = f32(sumraw
                                         + weights.ciwt[in_, ic] * qwind[1, i, ic]
                                         + weights.cwt[in_, ic] * qwind[2, i, ic])
                    # !Add in the Porod scattering if required
                    if q[i] > 0.0 and porodamp[in_] > 0.0:
                        porod = porodamp[in_] * powf(porodq / q[i], porodpwr)
                        sum = sum + porod
                        sumraw = sumraw + porod

                    soqism[i, in_] = sum
                    ufile[i, in_] = sumraw
                    if i >= nqmind[in_] and i <= nqmaxd[in_]:
                        backsc[i, in_] = datai[i, in_] - soqism[i, in_]
                        diffraw = datai[i, in_] - sumraw
                    else:
                        backsc[i, in_] = 0.0
                        diffraw = f32(0.0)
                    if i >= epsr16.nqmin[in_] and i <= epsr16.nqmax[in_]:
                        rfacd[in_] = rfacd[in_] + backsc[i, in_] * backsc[i, in_]
                        rfacdraw = rfacdraw + diffraw * diffraw
                if epsr16.nqpts[in_] > 0:
                    rfacd[in_] = rfacd[in_] / f32(epsr16.nqpts[in_])
                    rfacdraw = rfacdraw / f32(epsr16.nqpts[in_])
                else:
                    rfacd[in_] = 0.0
                    rfacdraw = f32(0.0)
                rfactot = rfactot + rfacd[in_]
                rfacraw = rfacraw + rfacdraw
            rfactot = rfactot / f32(epsr16.nfile)
            rfacraw = rfacraw / f32(epsr16.nfile)
            # c
            # c write out fit to total data
            # c
            ext1 = fstr(fmt(_F331, ".u"), 2)
            fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(ext1, 1, 2), 256)
            fname = auxiliary_routines.w_intr(1, nq, fname, epsr16.nfile, q, soqism, errs,
                                              titleq, titledat, mpts, mfile)
            # c
            # c write out difference between fit to total data and total data
            # c
            ext1 = fstr(fmt(_F331, ".v"), 2)
            fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(ext1, 1, 2), 256)
            fname = auxiliary_routines.w_intr(lpt1, lpt2, fname, epsr16.nfile, q, backsc, errs,
                                              titleq, titledat, mpts, mfile)
            # c
            # c write out site-site partial g(r)s as calculated by the program
            # c
            ext1 = fstr(fmt(_F331, ".g"), 2)
            fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(ext1, 1, 2), 256)
            fname = auxiliary_routines.w_int(nz, fname, molecom.npairs, z, epsr16.gorinter, epsr16.gr2inter,
                                             titler, titlepsf, mpz, mpairs)
            ext1 = fstr(fmt(_F331, ".l"), 2)
            fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(ext1, 1, 2), 256)
            fname = auxiliary_routines.w_int(nz, fname, molecom.npairs, z, epsr16.entsig, epsr16.entsig,
                                             titler, titlepsf, mpz, mpairs)
            ext1 = fstr(fmt(_F331, ".y"), 2)
            fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(ext1, 1, 2), 256)
            fname = auxiliary_routines.w_int(nz, fname, molecom.npairs, z, epsr16.gorintra, epsr16.gr2intra,
                                             titler, titlepsf, mpz, mpairs)
            # c
            # c write out the Fourier transform of the estimated partial structure factors
            # c
            ext1 = fstr(fmt(_F331, ".r"), 2)
            fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(ext1, 1, 2), 256)
            fname = auxiliary_routines.w_int(nz, fname, molecom.npairs, z, datagr, epsr16.gr2inter,
                                             titler, titlepsf, mpz, mpairs)

            # c combine partial g(r)'s according to weights matrix

            for in_ in range(1, epsr16.nfile + 1):
                for i in range(1, nz + 1):
                    sum1 = f32(0.0)
                    sum2 = f32(0.0)
                    for ic in range(1, molecom.npairs + 1):
                        sum = f32(weights.ciwt[in_, ic] * epsr16.gorintra[i, ic])
                        sum1 = sum1 + sum
                        sum2 = f32(sum2 + sum + weights.cwt[in_, ic] * epsr16.gorinter[i, ic])
                    gz2ism[i, in_] = sum1

                    # c subtract the large r limit so result oscillates about zero.

                    # c                 gofzism(i,in)=sum2-xwt(in)

                # c set the large r limit for Fourier transform of data

                grzerolim[in_] = 0.0

            # c Form the Fourier transform of the fits to the data over the same Q range as the
            # c data.

            # !            write(6,*) 'epsr> 2013: ',rhocg,q(10),nqmin,nqmax,nz,nfile,grzerolim
            ato_manipulation_rev4.stog(epsr16.rho, q, soqism, z, gofzism, epsr16.nqmin, epsr16.nqmax, nz,
                                       epsr16.nfile, mpts, mpz, mfile, grzerolim)
            # c
            # c write out the combined partial g(r)'s
            # c
            ext1 = fstr(fmt(_F331, ".x"), 2)
            fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(ext1, 1, 2), 256)
            fname = auxiliary_routines.w_int(nz, fname, epsr16.nfile, z, gofzism, gz2ism,
                                             titler, titledat, mpz, mfile)
            # c
            # c perform direct Fourier transform on input data ready for output
            # c
            ato_manipulation_rev4.stog(epsr16.rho, q, datai, z, dataigr, epsr16.nqmin, epsr16.nqmax, nz,
                                       epsr16.nfile, mpts, mpz, mfile, grzerolim)
            # c
            # c write out direct Fourier transform of input data.
            # c
            ext1 = fstr(fmt(_F331, ".w"), 2)
            fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(ext1, 1, 2), 256)
            fname = auxiliary_routines.w_int(nz, fname, epsr16.nfile, z, dataigr, gz2ism,
                                             titler, titledat, mpz, mfile)
            # c
            # c write out the reference potentials
            # c
            ext1 = fstr(fmt(_F331, ".o"), 2)
            fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(ext1, 1, 2), 256)
            fname = auxiliary_routines.w_int(nz, fname, molecom.npairs, z, potential_par.ecoul, potential_par.epot,
                                             titler, titlepsf, mpz, mpairs)
            # c
            # c write out the empirical potentials
            # c
            ext1 = fstr(fmt(_F331, ".p"), 2)
            fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(ext1, 1, 2), 256)
            fname = auxiliary_routines.w_int(nz, fname, molecom.npairs, z, potential_par.xavs, potential_par.xav,
                                             titler, titlepsf, mpz, mpairs)
            # ! Calculate a correlation coefficient between empirical potential and difference between simulated
            # ! and estimated total g(r)s
            epccoeff = f32(0.0)
            epcup = f32(0.0)
            epcdn = f32(0.0)  # noqa: F841
            eminmax = f32(0.0)  # noqa: F841
            # !            if(energabs.gt.0.0) then
            # !Transform raw simulated data in Q space to r space, using same Q ranges as data
            # (datagr2(mpz,mpairs) is passed to the dummy gofr(mpz,mfile): same leading dimension)
            ato_manipulation_rev4.stog(epsr16.rho, q, ufile, z, datagr2, epsr16.nqmin, epsr16.nqmax, nz,
                                       epsr16.nfile, mpts, mpz, mfile, grzerolim)
            nugdif = 0
            for in_ in range(1, epsr16.nfile + 1):
                ix = int(weights.xfile[in_])
                for i in range(1, nz + 1):
                    ugdif[i, in_] = 0.0
                    # ! Get the weighted sum of gofrs at this z value. Only calculate the r-space coefficient if the sum is .ne. 0.0
                    gortemp = f32(0.0)
                    for ip in range(1, molecom.npairs + 1):
                        if ix > 0:
                            gortemp = f32(gortemp + x_rayf0.weightx[1, ip, ix] * wav[2, i, ip])
                        else:
                            gortemp = f32(gortemp + weights.cwt[in_, ip] * wav[2, i, ip])
                    gortemp = gortemp / epsr16.ycoh[i]
                    if gortemp != 0.0 and z[i] > epcrmin:
                        prod = datagr2[i, in_] - dataigr[i, in_]
                        ugdif[i, in_] = prod
                # !                  epcfac=1.5
                # !                  delnz=nint(epcfac*pi/(q(nqmax(in))*cellst))+1 !Form a smoothed version of the data averaging over half the FT truncation width.
                # !                  nterms=2*delnz+1
                # !                  ifind=0
                for i in range(1, nz + 1):
                    # !                     ugdifsmo(i,in)=0.0
                    # !                     if(ugdif(i,in).ne.0.0.and.ifind.eq.0) ifind=i !Ignore any low r region where ugdif is zero
                    # !                     if(ifind.ne.0) then
                    # !                        sum1=0.0 !Zero the smoothed ugdif
                    # !                        do j=-delnz,delnz !step over the range of z values to be summed, reflecting at each end to keep total area constant
                    # !                           iref=i+j
                    # !                           if(iref.lt.ifind) then
                    # !                              iref=2*ifind-iref
                    # !                           else if(iref.gt.nz) then
                    # !                              iref=2*nz-iref
                    # !                           end if
                    # !                           sum1=sum1+ugdif(iref,in)
                    # !                        end do
                    # !                        ugdifsmo(i,in)=sum1/real(nterms)!Form top-hat smoothed version of data
                    # !                     end if
                    if z[i] > epcrmin and ugdif[i, in_] != 0.0:
                        nugdif = nugdif + 1
                        prod = ugdif[i, in_]
                        epcup = epcup + prod * prod
                        # !                        prod=ugdifsmo(i,in)
                        # !                        epcdn=epcdn+prod*prod
            epcup = epcup / f32(nugdif)
            # !               fractemp=rfacraw/(rfacraw+epcup)
            # !               epcdn=epcdn/real(nugdif) !Ostensibly this represents unfittable systematic effects - not currently used
            # !               epccoeff=epcfrac*fractemp*epcup+(1.0-epcfrac)*(1.0-fractemp)*rfacraw
            epccoeff = epcfrac * epcup + (r4(1.0) - epcfrac) * rfacraw
            lastqf = nextqf
            nextqf = epccoeff
            # !               if(epccoeff.gt.0.0) then
            # !                  bias=(epcup+epcdn)/epccoeff
            # !               else
            # !                  bias=-2.0
            # !               end if
            if potential_par.ereqstep > 0.0 and potential_par.potfac > 0.0:
                factor = f32(0.0)
                grad = f32(0.0)
                if nparab > nparabstart and nparabstart > 0:
                    # ! Use the parabola to define the best value of ereq (pressreq) if the parabola has a definite minimum,
                    # ! the minimum is below the current value of energabs, and energabs.gt.1.0.
                    if cparab != 0.0:
                        # Fitted correlation at current absolute energy
                        fcurrent = aparab + energabs * (bparab + energabs * cparab)
                        epcthreshuse = epcthresh
                        # !                        eminmax=-bparab/(2.0*cparab) !Position of minimum or maximum ugdif correlation.
                        # !                        feminmax=aparab+eminmax*(bparab+eminmax*cparab) !Fitted correlation at min/max energy
                        grad = bparab + r4(2.0) * cparab * energabs  # Gradient of fitted correlation at current energy
                        # ! The following code will tend to reduce ereq if the parabola is convex: this implies we are in a flatish region so should therefore be reducing ereq.
                        # !                     if(cparab.lt.0.0.and.abs(1.0-eminmax/energabs).lt.0.2) then !Only do this if the maximum energy is within 20% of current energy.
                        # !                        epcthreshuse=max(feminmax,epcthresh)
                        # !                     end if
                        if aparab < epcthreshuse:
                            grad = epcthreshuse - aparab  # Force backwards step in this case.
                        # Strength of the bias is determined by the size of the difference between current energy and min/max energy
                        eminmaxdif = sign(ereqbias, grad)
                        if abs(eminmaxdif) != 0.0:
                            # Choose from a distribution of step sizes, biassed by size of eminmaxdif.
                            term = r4(1.0) + eminmaxdif * (eminmaxdif - r4(4.0) * auxiliary_routines.ran1() + r4(2.0))
                            factor = r4(2.0) * potential_par.ereqstep * (r4(1.0) - sqrtf(term)) / eminmaxdif
                        else:
                            # Choose a step at random in the range +-ereqstep
                            factor = r4(2.0) * potential_par.ereqstep * (r4(2.0) * auxiliary_routines.ran1() - r4(1.0))
                    elif bparab != 0.0:
                        fcurrent = aparab + energabs * bparab  # Fitted correlation at current absolute energy
                        ctrllevel = aparab
                        ctrlbaseline = fcurrent
                        # !                        ctrlbaseline=fcurrent+epcthresh*bav*energabs
                        # !                        gradav1=epcthreshuse*gradav
                        # !                        grad=bparab+gradav1 !Gradient of fitted correlation at current energy
                        grad = ctrlbaseline - ctrllevel
                        # Strength of the bias is determined by the size of the difference between current energy and min/max energy
                        eminmaxdif = sign(ereqbias, grad)
                        # !                        if(fcurrent.lt.epcthreshuse) then
                        # !                           grad=epcthreshuse-fcurrent !Force backwards step in this case.
                        # !                           eminmaxdif=sign(ereqbias,grad) !Strength of the bias is determined by the size of the difference between current energy and min/max energy
                        # !                        else
                        # !                        if(gradav1.gt.0.0) then
                        # !                           eminmaxdif=ereqbias*grad/gradav1
                        # !                        end if
                        if abs(eminmaxdif) != 0.0:
                            # !                           term=1.0+eminmaxdif*(eminmaxdif-4.0*ran1()+2.0) !Choose from a distribution of step sizes, biassed by size of eminmaxdif.
                            # !                           factor=ereqstep*(1.0-sqrt(term))/eminmaxdif
                            factor = -potential_par.ereqstep * sign(r4(1.0), grad)
                        else:
                            # Choose a step at random in the range +-ereqstep
                            factor = potential_par.ereqstep * (r4(2.0) * auxiliary_routines.ran1() - r4(1.0))
                        units.write(6, lw(fcurrent, ctrllevel, ctrlbaseline, grad, eminmaxdir, factor))
                    else:
                        factor = potential_par.ereqstep
                    # !                     if(nextqf.gt.lastqf.and.grad.lt.0.0) then !Only change the EP amplitude if the fit got better. Always change the amplitude if the gradient is +ve.
                    if nextqf > lastqf:  # Only change the EP amplitude if the fit got better.
                        factor = f32(0.0)
                else:
                    # Assign initial ereq factor - this is such as to prevent epccoeff getting smaller than epcthresh
                    factor = potential_par.ereqstep
                # !                  factor=max(factor,-0.05)
                # !                  factor=min(factor,+0.05)
                if factor > 0.0:
                    factor = max(potential_par.erequnit, factor)  # Ensure the positive step is always at least 1 unit.
                if energabs > r4(2.0) * potential_par.erequnit and factor < 0.0:
                    factor = min(-potential_par.erequnit, factor)  # Ensure the negative step is always at least -1 unit.
                # !                  if(pressreq.lt.100.0) then
                # Keep ereq integer multiples of erequnit
                pressreq = potential_par.erequnit * f32(nint((energabs + factor) / potential_par.erequnit))
                # !                  else
                # !                     pressreq=erequnit*real(nint((energabs*(1.0+0.01*factor))/erequnit))
                # !                  end if
                # !                  if(pressreq.le.erequnit.and.nparab.gt.navmaxparab.and.nsums.lt.0) iinit=-1! Signal that you want to reset the distributions
                # !                  if(nparab.gt.nparabstart.and.nparabstart.gt.0) then
                if pressreq <= max(ereqmin, potential_par.erequnit) and factor < 0.0:
                    pressreq = ereqmin + f32(10) * potential_par.erequnit  # To stop getting stuck at minimum ereq.
                if ereqmax > 0.0:
                    if pressreq >= ereqmax:
                        pressreq = ereqmin + r4(0.75) * (ereqmax - ereqmin)  # To stop getting stuck at maximum ereq.
                    pressreq = min(pressreq, ereqmax)
                # !                  end if
                units.write(6, lw("epsr-automatic> ", nparab, nparabstart, energabs, pressreq))
                # !               end if
                # !               write(ext1,331) '.c'
                # !               fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
                # !               call w_int(nz,fname,nfile,z,ugdif,ugdifsmo,titler,titledat,mpz,mfile)
            # c Calculate and write out the intermolecule running coordination numbers. These are
            # c multiplied by the atomic fraction of the second atom of each pair. The first column
            # c corresponds to b about a, while the second column corresponds to a about b.

            # c We need to change the column titles appropriately

            for ip in range(1, molecom.npairs + 1):
                # c
                # c set the column title for this pair of atoms
                # c
                atom1 = fstr(potential_par.tlab1[potential_par.irefpcof[ip]], 4)
                atom2 = fstr(potential_par.tlab2[potential_par.irefpcof[ip]], 4)
                ilent1 = findex_padded(atom1, " ", 4) - 1
                ilent2 = findex_padded(atom2, " ", 4) - 1

                # c Calculate the number of spaces between titles

                icharttitle = ilent1 + ilent2 + 1
                iremain = 14 - icharttitle - 3
                ibetween = iremain + 4
                formattitle = fstr(fmt(_F992, ibetween, iremain), 256)
                ilenformat = 12 + 2 + 10 + 2 + 2  # noqa: F841
                titlepsf[ip] = fstr(fmt(formattitle, substr(atom1, 1, ilent1), substr(atom2, 1, ilent2),
                                        substr(atom2, 1, ilent2), substr(atom1, 1, ilent1)), 29)

                datagr[1, ip] = 0
                datagr2[1, ip] = 0
                jt1 = int(molecom.n2extr[ip])
                jt2 = int(molecom.n1extr[ip])
                for j in range(2, nz + 1):
                    cprod1 = molecom.cfrac[jt1] * epsr16.ycoh[j]
                    cprod2 = molecom.cfrac[jt2] * epsr16.ycoh[j]
                    datagr[j, ip] = datagr[j - 1, ip] + epsr16.gorinter[j, ip] * cprod1
                    datagr2[j, ip] = datagr2[j - 1, ip] + epsr16.gorinter[j, ip] * cprod2
            ext1 = fstr(fmt(_F331, ".z"), 2)
            fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(ext1, 1, 2), 256)
            fname = auxiliary_routines.w_int(nz, fname, molecom.npairs, z, datagr, datagr2,
                                             titler, titlepsf, mpz, mpairs)

            # C WRITE CURRENT DATA TO INPUT FILE IN A FORM SUITABLE TO BE READ AS NEXT
            # C INPUT FILE

            norma(ireset, niter, niter)  # noqa: F821
            # c     write(6,*) 'Done Norma'
            # c
            # c reset the values of ireset and iinit
            # c
            if ireset != 0:
                ireset = 0
            if iinit > 0:
                iinit = 0
            if iinit < 0:
                iinit = 1
            rfacwrite = rfactot
        else:
            # FORTRAN_BUGS.md section 2: on non-final iterations the .erg column 3 is the
            # rfacwrite of the .erg history (lsparabola1), or 0 (Rust epsr.rs rfacraw_last).
            rfacraw = rfacwrite
        # c
        # c output summary of results
        # c
        # c
        # c write out current energy and pressure and Rfactor and density
        # c
        units.write(9, fmt(_F594, energd, pressd, rfacraw, energabs, epsr16.rho, epcup, epccoeff,
                           aparab, bparab, cparab))
        # c     write(6,594) energd,pressd,rfactot,rhoout
        extfile = ".out"  # character*10
        fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(extfile, 1, 4), 256)
        openf(10, fname)  # noqa: F821
        wtent = f32(0.0)  # noqa: F841
        currentdate, currenttime, currentzone, _values = runtime.date_and_time()
        currentdate = fstr(currentdate, 8)
        currenttime = fstr(currenttime, 10)
        currentzone = fstr(currentzone, 5)
        units.write(10, fmt(_F146, substr(currenttime, 1, 2), substr(currenttime, 3, 4), substr(currenttime, 5, 6),
                            substr(currentdate, 7, 8), substr(currentdate, 5, 6), substr(currentdate, 1, 4),
                            Ch(version.version, 20)))
        units.write(10, fmt(_F147, ntimes, iter, epsrinput_par.nsum, rfacwrite, epcup, energabs,
                            epccoeff, ctrllevel, ctrlbaseline, epcthresh,
                            foldsv, foldsb, ftethered, fcoordconstraints, energthreeb, pressthreeb,
                            nparab, bav))
        if threebody_mod.nthreeb > 0:
            units.write(10, lw(" Energy for each triplet: ",
                               *[threebody_mod.sume_triplet_dist[id] / f32(clustercom.nmol)
                                 for id in range(1, threebody_mod.nthreeb + 1)]))
        units.write(10, fmt(_F145, molecom.sizefactor, epsr16.rho))
        units.write(10, fmt('("  Net charge = ",f10.5)', molecom.sumcharge))
        vals = []
        for j in range(1, epsr16.nfile + 1):
            vals += [j, rfacd[j]]
        units.write(10, fmt(_F148, *vals))
        # !      write(10,149) (j,rfac(j),j=1,npairs)
        units.write(10, fmt(_F219, sentb / epsr16.rho, energd, epsr16.upol, pressd, epsr16.amubar, epsr16.amubar2,
                            enav, enav2, prav, prav2))
        for j in range(1, 2 + 1):
            clustercom.nmoveit[j] = 0
            clustercom.nintert[j] = 0
            clustercom.nwholet[j] = 0
            clustercom.nmovemt[j] = 0
            for imoltype in range(1, clustercom.nmoltype + 1):
                clustercom.nmoveit[j] = clustercom.nmoveit[j] + clustercom.nmovei[j, imoltype]
                clustercom.nintert[j] = clustercom.nintert[j] + clustercom.ninter[j, imoltype]
                clustercom.nwholet[j] = clustercom.nwholet[j] + clustercom.nwhole[j, imoltype]
                clustercom.nmovemt[j] = clustercom.nmovemt[j] + clustercom.nmovem[j, imoltype]
        vals = []
        for j in range(1, 2 + 1):
            vals += [clustercom.nmoveit[j], clustercom.nintert[j], clustercom.nwholet[j], clustercom.nmovemt[j]]
        units.write(10, fmt(_F220, *vals))
        units.close(10)
        # c
        # c write input file
        # c
        # c
        # c first save the internal values to their symbol counterparts
        # c

        text = fstr("potfac", 20)
        icount = 0
        icount = file_and_directory_commands.search_for_symbol(text, icount)
        iref = file_and_directory_commands.ireff(icount, 1)
        symbols.symbol_value[iref] = fstr(fmt("(2(1x,f10.6))", potential_par.potfacamp, potential_par.potfac), 256)
        iref = iref + 1
        if pressreq > 100000.0:
            pressreq = r4(100000.0)
        symbols.symbol_value[iref] = fstr(fmt("(e12.5)", pressreq), 256)
        iref = iref + 1
        symbols.symbol_value[iref] = fstr(fmt("(f10.2)", ereqmin), 256)
        iref = iref + 1
        symbols.symbol_value[iref] = fstr(fmt("(f10.2)", ereqmax), 256)
        iref = iref + 1
        symbols.symbol_value[iref] = fstr(fmt("(f6.3,1x,f6.3)", ereqstepread, potential_par.erequnit), 256)
        iref = iref + 3
        symbols.symbol_value[iref] = fstr(fmt("(i5,1x,i5)", navmaxparab, nparabstartfactor), 256)
        if epsrinput_par.nsums >= 0:
            epsrinput_par.nsums = epsrinput_par.nsum
        text = fstr("nsumt", 20)
        icount = 0
        icount = file_and_directory_commands.search_for_symbol(text, icount)
        iref = file_and_directory_commands.ireff(icount, 1)
        symbols.symbol_value[iref] = fstr(fmt("(i5)", epsrinput_par.nsums), 256)
        text = fstr("inter", 20)
        icount = 0
        icount = file_and_directory_commands.search_for_symbol(text, icount)
        iref = file_and_directory_commands.ireff(icount, 1)
        symbols.symbol_value[iref] = fstr(fmt("(i5)", inter), 256)
        iref = file_and_directory_commands.ireff(iresetnumbers.iresetref, 1)
        symbols.symbol_value[iref] = fstr(fmt("(i2)", ireset), 256)
        iref = file_and_directory_commands.ireff(iresetnumbers.iinitref, 1)
        symbols.symbol_value[iref] = fstr(fmt("(i2)", iinit), 256)
        # c
        # c unnormalise the coordinates before writing the .ato file
        # c
        ato_read_and_save.unnormalise_coordinates()
        # c     write(6,*) 'Done unnormalise'

        # c Generate new sizefactor if needed

        sizefactorold = molecom.sizefactor  # noqa: F841
        if molecom.sizefactor > 1.0 and foldsv < thresholde:
            molecom.sizefactor = molecom.sizefactor * molecom.decline_rate
            if molecom.sizefactor < 1.0:
                molecom.sizefactor = f32(1.0)
        # !If a larger positive energy excursion occurs, increase the effective box size.
        # !      if(sizefactor.gt.1.0.and.foldsv.gt.abs(foldsb)) then
        # !         sizefactor=sizefactor/decline_rate
        # !      end if
        symbols.symbol_value[irefsizefactor] = fstr(
            fmt("(f10.5,1x,f10.5,1x,e12.5)", molecom.sizefactor, molecom.decline_rate, thresholde), 256)
        # !      text='qstep'
        # !      call search_for_symbol(text,icount)
        # !      iref=ireff(icount,1)
        # !      write(symbol_value(iref),'(f10.5)') qstep
        text = fstr("rmaxgr", 20)
        # NB: icount is not reset to 0 here in the Fortran (it keeps the 'inter' value)
        icount = file_and_directory_commands.search_for_symbol(text, icount)
        iref = file_and_directory_commands.ireff(icount, 1)
        symbols.symbol_value[iref] = fstr(fmt("(e13.6)", rmaxgrsave), 256)
        units.write(6, fmt('("Steps: ",4(1x,e12.5))', molecom.stepmi, molecom.stepri, stepmisave, steprisave))
        # !Restore intramolecular step sizes
        molecom.stepmi = stepmisave
        molecom.stepri = steprisave

        # !Adjust boxsize to control pressure if requested

        if control_pressure and npcontrols > 0:
            # !Set up the smallest value that can be added to 1.0 to give a new value
            # !         smallest=tiny(1.0)
            # !         oneplusrat=1.0+smallest
            # !         do while (oneplusrat.eq.1.0)
            # !            smallest=smallest*10.0
            # !            oneplusrat=1.0+smallest
            # !         end do
            rat = pressd - control_p_value
            if abs(rat) > 1.0:
                rat = sign(r4(1.0), rat)
            rat = control_p_step * rat
            # !         if(abs(rat).lt.smallest) rat=sign(smallest,rat)
            # !         if(abs(rat).gt.ratplimit) rat=sign(ratplimit,rat)
            dmaxfac = powf(r4(1.0) + rat, r4(1.0) / f32(npcontrols))
            afac = f32(1.0)
            bfac = f32(1.0)
            cfac = f32(1.0)
            pfac = f32(1.0)
            if pacontrol:
                afac = dmaxfac
            if pbcontrol:
                bfac = dmaxfac
            if pccontrol:
                cfac = dmaxfac
            if ppcontrol:
                pfac = r4(1.0) / dmaxfac
            units.write(6, lw(afac, bfac, cfac, pfac))
            if ppcontrol:
                # !Only do this if potfac is greater than zero, i.e. potential is being modified
                if potential_par.potfac > 0.0:
                    multiply_pcoefex(pfac)  # noqa: F821
            else:
                newphib = lattice.phib
                newthetac = lattice.thetac
                newphic = lattice.phic
                units.write(6, lw(lattice.alat, lattice.blat, lattice.clat, lattice.phib, lattice.thetac, lattice.phic))
                ato_manipulation_rev4.change_boxsize(afac, bfac, cfac, newphib, newthetac, newphic)
                units.write(6, lw(lattice.alat, lattice.blat, lattice.clat, lattice.phib, lattice.thetac, lattice.phic))
                uniform(epsr16.ycoh, epsr16.cellst, epsr16.rho, epsrinput_par.rmaxgr, mpz)  # noqa: F821
        # !Make sure any changes to the repulsive and Gaussian parameters are recorded
        # !(Need to implement this - currently only new values of pcoefex will be saved.)

        read_and_save_files.save_files_others(14)

        # c     write(6,*) 'Done save_files_others'

        # c
        # c Normalise coordinates ready for the next iteration
        # c
        # ! Apply current sizefactor in case it has changed
        molecom.dmaxmol = clustercom.dmax * molecom.sizefactor
        ato_read_and_save.normalise_coordinates()
        # ! Recalculate the potential in case dmax has changed
        # !      call initpt(ireset,rminpt,rmaxpt,expecf,polar,szeros,rshmin,rbroad)

        # c     write(6,*) 'Done normalise'

        # ------------------------------------------------------------ 3000 continue
        # End of EPSR iteration `iter` (the Fortran label 3000 closes the iteration loop).
        if _iteration_finished(iter, niter):
            break
    units.close(9)
    # !      close(35)
    # c     write(6,*) 'Done close(9)'
    return


# =============================================================================
# from epsr_standalone_rev3_part2.py
# =============================================================================

# FORMAT statements
_F_INITA_959 = "(5(1x,i5))"
_F_SUMPE_200 = "(3(1x,e13.6))"
_F_SUMPE_201 = "(2(1x,a3))"
_F_PTGNSQ_332 = "(a,i2.2,a)"
_F_PTGNSQ_59 = "(1x,'Minimum radius integer of potential:',(1x,i5))"
_F_OPENF_100 = "(1x,' Open error in OPENF',1x,a50,1x,a10)"
_F_GETWTS_108 = "(1x,a3,1x,a4)"
_F_GETWTS_109 = "(/'get_wts> Exceeded limit of ',i2,' X-ray files.')"
_F_GETWTS_107 = "(1x,a3,1x,a3)"
_F_WTYPE_100 = "(1x,'Mismatched types between atomic coordinates and','neutron weights for file ',i5)"
_F_WXTYPE_100 = "(1x,'Mismatched types between atomic coordinates and','X-ray weights for file ',i5)"
_F_WXCREATE_100 = "(13(1x,e13.6))"


@routine
def interp(nqd, nzpts, qd, qi, data, alambi, in_, nqmin, nqmax, nqdim, nqoutd, mfile):
    """epsr_standalone_rev3.f:interp -- interpolates data on grid qd onto grid qi.

    Output in alambi(:, in).
    Returns (nqmin, nqmax).
    """
    #
    # INTERPOLATES DATA ON GRID QD ONTO GRID QI. OUTPUT IN ALAMBI
    #
    q1 = qd[1]
    q2 = qd[2]
    d1 = data[1]
    d2 = data[2]
    n = 2
    i = 1
    nqmin = 1
    while True:
        # 5 CONTINUE
        qt = qi[i]
        if q1 > qt:
            # GO TO 15
            # 15 CONTINUE
            i = i + 1
            nqmin = i
            continue  # GO TO 5
        # 10 CONTINUE
        goto30 = False
        while q2 < qt:
            # GO TO 20
            # 20 CONTINUE
            n = n + 1
            if n > nqd:
                goto30 = True
                break
            d1 = d2
            q1 = q2
            q2 = qd[n]
            d2 = data[n]
            # GO TO 10
        if goto30:
            break
        grad = (d2 - d1) / (q2 - q1)
        alambi[i, in_] = d1 + (qt - q1) * grad
        i = i + 1
        if i > nzpts:
            break  # GO TO 30
        # GO TO 5
    # 30 CONTINUE
    nqmax = i - 1
    return nqmin, nqmax


#
# initialize distribution A
#
@routine
def inita(ireset, iinit, nq, q, qstep, nz, z, fwhm, fwhmq, intra):
    """epsr_standalone_rev3.f:inita -- initialize distribution A.

    Reads the EP parameters from the symbol table, the accumulators from the .inpa file,
    sets up the r and Q grids, the Fourier transform matrix foqols and (when coarse
    graining) the spherical form factors.
    Returns (iinit, nz); q and z are filled in place.
    """
    # real nitterms,niitermstotal,niitermsself ; integer ncoeffpread,nterms
    pi = PI4ATAN
    epsr16.entwt[1] = f32(0.0)
    #     write(6,966) nfile,ndist
    #
    # Now read the empirical potential parameters from the relevant symbol_values
    #
    # ROVERLAP = minimum allowed separation between atoms on different molecules
    # RMINFAC = factor on sigmas used to set up minimum separation for indivdidual atom pairs
    # RMINPT = point at which the truncation of the potential beings, e.g. 8 A
    # RMAXPT = limiting radius to which potential is calculated. Typically this might be 10A
    # RBROAD = Decay applied to potential perturbation (Normally small or zero)
    # EXPECF = hardness of Exponential potential for controlling g(r) below rmin (normally 0.1)
    # expecfpwr = power to increase expecf as a function of minimum distance: expecf=expecf*d**(expecfpwr)
    # RCHARGE = 1 to calculate energy due to polarisation of molecule. Currently this only works
    # correctly for SPC/E water. Normally set to 0.
    # NB: ifilecount is an undeclared, uninitialised local of inita.  ireff only uses it when
    # iqminref lies in the per-dataset symbol range, which it does not (iqminref <=
    # nsymboluptodata); 1 is passed, as the Rust does (structure.rs inita: ireff(iqminref, 1)).
    ifilecount = 1
    sv = symbols.symbol_value
    iref = file_and_directory_commands.ireff(iresetnumbers.iqminref, ifilecount)
    iref = iref + 1
    epsrinput_par.roverlap = read_internal(sv[iref]).real(epsrinput_par.roverlap)
    iref = iref + 1
    epsrinput_par.rminfac = read_internal(sv[iref]).real(epsrinput_par.rminfac)
    iref = iref + 1
    epsrinput_par.rminpt = read_internal(sv[iref]).real(epsrinput_par.rminpt)
    epsrinput_par.rminpt = max(r4(1.0), epsrinput_par.rminpt)  # Make sure silly values can't be specified
    iref = iref + 1
    epsrinput_par.rmaxpt = read_internal(sv[iref]).real(epsrinput_par.rmaxpt)
    epsrinput_par.rmaxpt = max(r4(2.0), epsrinput_par.rmaxpt)  # Make sure silly values can't be specified
    iref = iref + 1
    epsrinput_par.rbroad = read_internal(sv[iref]).real(epsrinput_par.rbroad)
    iref = iref + 1
    epsrinput_par.expecf = r4(0.3)
    epsrinput_par.expecfpwr = f32(0.0)
    commands.sometext = fstr(sv[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 11,
                                  commands.nstsize, commands.nworddim)
    if commands.nwords > 0:
        epsrinput_par.expecf = read_internal(
            substr(commands.sometext, commands.ncf[1], commands.ncl[1])).real(epsrinput_par.expecf)
    if commands.nwords > 1:
        epsrinput_par.expecfpwr = read_internal(
            substr(commands.sometext, commands.ncf[2], commands.ncl[2])).real(epsrinput_par.expecfpwr)
    # Decide whether or not to use a harmonic repulsive force
    iref = iref + 1
    potential_par.harmonic = False
    if sv[iref] == "harmonic":
        potential_par.harmonic = True
    else:
        sv[iref] = "exponential"
    iref = iref + 1
    # This decides whether to use Gaussian or 'modmorse' for additional potentials
    potential_par.addpotential = fstr(read_internal(sv[iref]).word(potential_par.addpotential), 20)
    iref = iref + 1
    epsrinput_par.rcharge = read_internal(sv[iref]).real(epsrinput_par.rcharge)
    iref = iref + 1
    # POWER = power for repulsive term: (Normally 12, must be > 6, otherwise defaults to 12)
    molecom.power = read_internal(sv[iref]).real(molecom.power)
    if molecom.power <= 6:
        molecom.power = f32(12)
    #
    # Read number of coefficients, no. of iterations to refine potential, separation in r
    # between coefficients, maximum radius to calculate coefficients,
    # acceptance factor, and width of functions
    #
    iref = iref + 1
    potential_par.pdstep = r4(0.1)
    potential_par.pdmin = r4(1.0) + potential_par.pdstep / epsrinput_par.rminpt
    # pdmax is set by programme to the truncation in the potential
    potential_par.pdmax = epsrinput_par.rmaxpt
    # ncoeffp is set to 10 times rmaxpt or mcoeff, whichever is smaller
    #      ncoeffp=min(int(10.0*rmaxpt+0.0001),mcoeff)
    potential_par.ncoeffp = min(
        (ifix(epsrinput_par.rminpt / potential_par.pdstep)
         + ifix(logf(epsrinput_par.rmaxpt / epsrinput_par.rminpt) / logf(potential_par.pdmin))),
        potential_par.mcoeff)
    units.write(6, lw("inita> ncoeffp = ", potential_par.ncoeffp))
    ncoeffpread = 0
    ncoeffpread = read_internal(sv[iref]).int(ncoeffpread)  # noqa: F841
    iref = iref + 1
    pdmaxread = f32(0.0)
    pdmaxread = read_internal(sv[iref]).real(pdmaxread)  # noqa: F841
    iref = iref + 1
    pdstepread = f32(0.0)
    pdstepread = read_internal(sv[iref]).real(pdstepread)  # noqa: F841
    iref = iref + 1
    potential_par.npitss = read_internal(sv[iref]).int(potential_par.npitss)
    iref = iref + 1
    potential_par.paccept = read_internal(sv[iref]).real(potential_par.paccept)
    iref = iref + 1
    # Select whether to use Gaussian or Poisson representation of EP.
    potential_par.gaussian = False
    if sv[iref] == "T" or sv[iref] == "t" or sv[iref] == "Gaussian":
        potential_par.gaussian = True
        sv[iref] = "Gaussian"
    else:
        sv[iref] = "Poisson"
    iref = iref + 1
    potential_par.psigma2read = read_internal(sv[iref]).real(potential_par.psigma2read)
    iref = iref + 1
    potential_par.refpotfac = read_internal(sv[iref]).real(potential_par.refpotfac)
    #      pdstep=pdmax/ncoeffp
    potential_par.psigma2 = potential_par.psigma2read
    # To prevent overflows, there is a limit on psigma2
    if potential_par.psigma2 < r4(0.004):
        potential_par.psigma2 = r4(0.004)
        potential_par.psigma2read = potential_par.psigma2

    if ireset == 0:
        #
        # only need to open the .inpa file if iinit.eq.0
        #
        if iinit == 0:
            extfile = fstr(".inpa", 10)
            fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(extfile, 1, 5), 256)
            ierr = 0
            try:
                units.open(10, fname, status="old")
            except FortranIOError as e:
                ierr = e.iostat
            #
            # if there are problems with the inpa file, then set iinit to 1 and do not attempt to read any
            # more from this file
            #
            if ierr != 0:
                iinit = 1
        #
        # if iinit is still zero, then carry on reading .inpa file
        #
        if iinit == 0:
            #
            # OTHER parameters relate to calculating average energy and pressure
            #
            r = units.list_read(10)
            epsrinput_par.nisum = r.int(epsrinput_par.nisum)
            epsrinput_par.nesum = r.int(epsrinput_par.nesum)
            epsrinput_par.nesumt = r.int(epsrinput_par.nesumt)
            epsrinput_par.ensum = r.real(epsrinput_par.ensum)
            epsrinput_par.ensum2 = r.real(epsrinput_par.ensum2)
            epsrinput_par.prsum = r.real(epsrinput_par.prsum)
            epsrinput_par.prsum2 = r.real(epsrinput_par.prsum2)
            #
            # read in existing site-site coefficients from .PCOF file
            #
            for id in range(1, epsr16.nfile + 1):
                units.skip_record(10)
                r = units.list_read(10)
                potential_par.npcoeff[id] = r.int(potential_par.npcoeff[id])
                potential_par.pstepd[id] = r.real(potential_par.pstepd[id])
                epsr16.sigma2[id] = r.real(epsr16.sigma2[id])
                #
                # read the difference coefficients - these represent the difference between the data
                # of the present simulation and the fit.
                #
                r = units.list_read(10)
                for ic in range(1, int(potential_par.npcoeff[id]) + 1):
                    potential_par.pcoefd[ic, id] = r.real(potential_par.pcoefd[ic, id])
            #
            # read the inverted weights matrix for this potential
            #
            r = units.list_read(10)
            potential_par.feedbackold = r.real(potential_par.feedbackold)
            for in_ in range(1, epsr16.nfile + 1):
                potential_par.efacold[in_] = r.real(potential_par.efacold[in_])
            weights.nqinv = r.int(weights.nqinv)
            nnew = epsr16.nfile + molecom.npairs
            for iq in range(1, weights.nqinv + 1):
                for ip in range(1, molecom.npairs + 1):
                    r = units.list_read(10)
                    weights.qinv[iq] = r.int(weights.qinv[iq])
                    r = units.list_read(10)
                    for id in range(1, nnew + 1):
                        weights.cwttq[iq, ip, id] = r.dble(weights.cwttq[iq, ip, id])
        else:
            epsrinput_par.nsum = 0
            epsrinput_par.nesum = 0
            epsrinput_par.nesumt = 5
            epsrinput_par.nsums = -1
            epsrinput_par.nisum = intra
            #
            # initialise the number of coefficients for each dataset and
            # set the potential coefficients to zero
            #
            for id in range(1, epsr16.nfile + 1):
                potential_par.npcoeff[id] = potential_par.ncoeffp
                epsr16.sigma2[id] = potential_par.psigma2
                for ic in range(1, int(potential_par.npcoeff[id]) + 1):
                    potential_par.pcoefd[ic, id] = f32(0.0)
                potential_par.pstepd[id] = r4(0.002)
        units.close(10)
        #
        # initialise certain accumulators
        #
        if epsrinput_par.nsums >= 0:
            epsrinput_par.nsum = epsrinput_par.nsums
        else:
            epsrinput_par.nsum = 0
        potential_par.npits = potential_par.npitss
    else:
        epsrinput_par.nsum = 0
        epsrinput_par.nesum = 0
        epsrinput_par.nesumt = 5
        epsrinput_par.nsums = -1
        epsrinput_par.nisum = intra
        #
        # initialise the number of coefficients for each dataset and
        # set the potential coefficients to zero
        #
        for id in range(1, epsr16.nfile + 1):
            potential_par.npcoeff[id] = potential_par.ncoeffp
            epsr16.sigma2[id] = potential_par.psigma2
            for ic in range(1, int(potential_par.npcoeff[id]) + 1):
                potential_par.pcoefd[ic, id] = f32(0.0)
            potential_par.pstepd[id] = r4(0.002)
        potential_par.npits = potential_par.npitss
    #
    # rmaxgr should be no bigger than half the box size
    #
    if epsrinput_par.rmaxgr > r4(0.5) * epsr16.cellmax or epsrinput_par.rmaxgr <= 0.0:
        epsrinput_par.rmaxgr = r4(0.5) * epsr16.cellmax
    units.write(6, lw("epsr-inita> rmaxgr,rmaxcell,cell,cellmax = ", epsrinput_par.rmaxgr, epsr16.rmaxcell,
                      epsr16.cell, epsr16.cellmax))
    #
    # set up x values (=r values). For x>rlog the step size increases exponentially
    #
    x = epsr16.x
    epsr16.npt = 0
    x3dold = f32(0.0)  # noqa: F841
    celstt = epsr16.cellst
    stpfac = r4(1.0) + celstt / epsr16.rlog  # noqa: F841
    celmin = epsr16.cell  # noqa: F841
    x1dold = r4(-0.5) * celstt
    cellhf = r4(0.25) * celstt
    rminsv = x1dold  # noqa: F841
    x[1] = f32(0.0)
    j = 2
    epsr16.npt1 = 2
    epsr16.npt2 = 2
    epsr16.nptmin = 2
    epsr16.nptmax = 2
    x1d = x1dold + celstt
    rmaxgr = epsrinput_par.rmaxgr
    rminpt = epsrinput_par.rminpt
    rmaxpt = epsrinput_par.rmaxpt
    while j < mpz and x1d < rmaxgr:
        x3d = x1d * x1d * x1d / f32(3)
        #        x3d=x1d
        x[j] = x1d - cellhf
        epsr16.entwt[j] = f32(1) / x[j]
        x1dold = x1d
        x3dold = x3d  # noqa: F841
        #        if(x1d.gt.rlog) celstt=celstt*stpfac
        cellhf = r4(0.5) * celstt
        #
        # save the integer corresponding to the range of the potential
        #
        if x1d < rminpt:
            epsr16.nptmin = j + 1
        if x1d < rmaxpt:
            epsr16.nptmax = j + 1
        #
        # save the index for x values below the specified maximum value for
        # valid Fourier transforms
        #
        if x1d < rmaxgr:
            epsr16.npt1 = j
        #     write(6,*) j,x(j),nptmin,nptmax,npt1
        j = j + 1
        x1d = x1dold + celstt
    #
    # set the last value of x
    #
    x[j] = x1d - cellhf
    #     write(6,*) j,x(j),nptmin,nptmax,npt1
    epsr16.npt = j
    epsr16.npt2 = j - 2
    if epsr16.npt1 > epsr16.npt2:
        epsr16.npt1 = epsr16.npt2
    units.write(6, fmt(_F_INITA_959, epsr16.nptmin, epsr16.nptmax, epsr16.npt1, epsr16.npt2, epsr16.npt))
    decay = x[epsr16.nptmax] - x[epsr16.nptmin]  # noqa: F841
    if epsrinput_par.nesumt < 1:
        epsrinput_par.nesumt = 1
    if epsrinput_par.nesum > epsrinput_par.nesumt - 1:
        epsrinput_par.nesum = epsrinput_par.nesumt - 1
    if epsrinput_par.nisum > intra:
        epsrinput_par.nisum = intra
    #
    # calculate uniform distribution for the box
    # To do this try to ensure that at least one element in the nearest shell (27
    # boxes within 1 step from origin, including one at origin) is
    # sampled on average. We repeat this calculation every intra MC moves so that
    # it does not unduly bias the summed g(r).
    #
    # NB: volume is the clustercom COMMON variable, not a local.
    clustercom.volume = lattice.alatx * lattice.blaty * lattice.clatz
    rhotemp = f32(molecom.natomr) / clustercom.volume
    uniform(epsr16.ycoh, epsr16.cellst, rhotemp, epsrinput_par.rmaxgr, mpz)  # noqa: F821
    # (the commented-out code that read the uniform distribution from a .uni file is
    # omitted here; see the Fortran lines 2854-2880)

    #
    # define the q scale on which calculation is made
    #
    for i in range(1, nq + 1):
        q[i] = qstep * f32(i - 1)
    #
    # sigma for instrumental resolution
    #
    # sqrt(2.*alog(2.0)) is folded by gfortran in single precision
    sqrt2ln2 = sqrtf(r4(2.0) * logf(r4(2.0)))
    sigres = r4(0.5) * fwhm / sqrt2ln2
    sigrsq = r4(0.5) * fwhmq / sqrt2ln2
    #
    # calculate transform matrices
    #
    pi2 = f32(2) * pi  # noqa: F841
    trfac = r4(4.0) * pi * rhotemp
    lorchlim = r4(4.4934) / epsrinput_par.rmaxgr
    ycoh = epsr16.ycoh
    foqols = epsr16.foqols
    lorchqwidth = epsrinput_par.lorchqwidth
    lorchqwidthq = epsrinput_par.lorchqwidthq
    npt1 = epsr16.npt1
    trfac_d = float(trfac)
    for i in range(1, nq + 1):
        qq = q[i]
        r1 = f32(0)
        j1qr1 = j1qr(qq, r1)
        rstep = x[1] - r1
        for j in range(1, npt1 + 1):
            rr = x[j]
            rstep = rr - r1
            r2 = rr + rstep
            j1qr2 = j1qr(qq, r2)
            if ycoh[j] > 0.0:
                foqols[j, i] = trfac_d * (j1qr2 - j1qr1) / float(ycoh[j])
            else:
                foqols[j, i] = f32(0.0)
            r1 = r2
            j1qr1 = j1qr2
            qsig = sigres + qq * sigrsq
            qrs = qsig * rr
            qr2 = r4(0.5) * qrs * qrs
            sums = f32(0.0)
            if qr2 < 30.0:
                sums = expf(-qr2)
                if lorchqwidth > 0.0 or lorchqwidthq > 0.0:
                    bwidth = rr * lorchlim * (lorchqwidth + lorchqwidthq * qq)
                    sums = sums * pofr(bwidth)
            foqols[j, i] = foqols[j, i] * sums
    # If coarse graining, calculate the spherical form factors for each atom type
    if molecom_mod.coarsegraining:
        nmol = clustercom.nmol
        natom = molecom.natom
        ntype = molecom.ntype
        npairs = molecom.npairs
        itype = molecom.itype
        mref = molecom.mref
        naf = clustercom.naf
        nal = clustercom.nal
        m = molecom_mod
        # Calculate modified scale factors based on number of atoms in each molecule
        m.scalefacmoluse = reallocation_routines.reallocate1d_r(m.scalefacmoluse, nmol)
        for im in range(1, nmol + 1):
            # Get the total number of atoms in this molecule
            totalatoms = f32(0.0)
            for ia in range(int(naf[im]), int(nal[im]) + 1):
                totalatoms = totalatoms + m.natomsphere[itype[ia]] * m.pseudofacmol[itype[ia]]
            # selffac is used because intra-molecular atomic separations can be significant, especially for larger molecules
            m.scalefacmoluse[im] = powf(totalatoms, m.selfbroadpwr)  # May need to use cube root? Or extra smearing factor?
        qfactor = r4(4.0) * pi / r4(3.0)  # Factor to get average number of atoms per spherical shell
        #         qfactor=4.0*pi*(qscaleinputfactor**3) !Factor to get average number of atoms per spherical shell
        #         qfactor=(qscaleinputfactor**3) !Factor to get average number of atoms per spherical shell
        # In this first part we determine how many spherical form factor bins there will be for each atom type (nscalefacmol), and the bin no. for each molecule
        # (All the atoms in a molecule will belong to this same bin, since the scale factor applies to the molecule, not the individual atoms.)
        # Determine the minimum and maximum scale factors for each molecular type
        if m.scalefactypemin is not None:
            m.scalefactypemin = None
            m.scalefactypemax = None
            m.nsffbintype = None
        m.scalefactypemin = real_array(ntype)
        m.scalefactypemax = real_array(ntype)
        m.nsffbintype = int_array(ntype)
        m.scalefactypemin[:] = r4(1.0)
        m.scalefactypemax[:] = r4(1.0)
        for ia in range(1, natom + 1):
            im = mref[ia]
            it = itype[ia]
            m.scalefactypemin[it] = min(m.scalefacmol[im], m.scalefactypemin[it])
            m.scalefactypemax[it] = max(m.scalefacmol[im], m.scalefactypemax[it])
        m.maxnscalefactype = 0
        m.nsffbintype[:] = 0
        # Get the number of form factor bins for each atom type
        for it in range(1, ntype + 1):
            # No. of scale factor bins for each atom type
            m.nsffbintype[it] = ifix((m.scalefactypemax[it] - m.scalefactypemin[it]) / m.scalefacratio) + 1
            if m.nsffbintype[it] > m.maxnscalefactype:
                m.maxnscalefactype = int(m.nsffbintype[it])
            #            write(6,*) 'inita> Minimum and maximum scale factors for atom type: ',it,scalefactypemin(it),scalefactypemax(it),nsffbintype(it)
            # Now search through the atom types, and assign the
        # Assign the sphere form factor bin for each molecule, and create a store for the no. terms in each form factor bin for each atom type
        if m.assignatomscalefac is not None:
            m.assignatomscalefac = None
        m.assignatomscalefac = int_array(natom)
        for ia in range(1, natom + 1):
            im = mref[ia]
            it = itype[ia]
            ibin = ifix((m.scalefacmol[im] - m.scalefactypemin[it]) / m.scalefacratio) + 1
            m.assignatomscalefac[ia] = ibin
            #            write(6,*) 'inita> ',ia,ibin
        if m.sphereformfactor is not None:
            m.sphereformfactor = None
            m.intraformfactor = None
            m.intraformfactorself = None
            m.ntermssff = None
            m.ntermsintra = None
            m.ntermsintraself = None
        # NB: fracntermssff is allocated below but missing from the deallocate list above, so a
        # second call of inita with coarse graining fails in the Fortran with a runtime error.
        if m.fracntermssff is not None:
            raise FortranStop(2, "Fortran runtime error: Attempting to allocate already allocated variable "
                                 "'fracntermssff'")
        maxn = m.maxnscalefactype
        m.sphereformfactor = real_array(nq, maxn, ntype)
        m.intraformfactor = real_array(nq, npairs)
        m.intraformfactorself = real_array(nq, ntype)
        m.ntermssff = int_array(maxn, ntype)
        m.fracntermssff = real_array(maxn, ntype)
        m.ntermsintra = int_array(npairs)
        m.ntermsintraself = int_array(ntype)
        # Zero all the type form factors: will form an average over the different radii for each type
        m.ntermssff[:] = 0
        m.fracntermssff[:] = 0
        m.ntermsintra[:] = 0
        m.ntermsintraself[:] = 0
        m.sphereformfactor[:] = 0.0
        m.intraformfactor[:] = 0.0
        m.intraformfactorself[:] = 0.0
        sff = m.sphereformfactor
        for ia in range(1, natom + 1):
            im = mref[ia]
            it = itype[ia]
            #            qfactorit=qfactor*cfraccg(it)
            qfactorit = qfactor * m.rhosphere[it]
            radlocalmin = m.radminsphere[it]
            radlocalmax = m.radmaxsphere[it]
            radlocalmin3 = radlocalmin * radlocalmin * radlocalmin  # noqa: F841
            radlocalmax3 = radlocalmax * radlocalmax * radlocalmax  # noqa: F841
            #            radlocalmin=0.0
            #            radlocalmax=0.5*sigcor(it)
            ibin = m.assignatomscalefac[ia]
            m.ntermssff[ibin, it] = m.ntermssff[ibin, it] + 1
            scalefac = m.scalefacmol[im]
            selffac = m.scalefacmoluse[im]
            #            selffac=1.0
            radusemin = radlocalmin * scalefac
            radusemax = radlocalmax * scalefac
            radusemin3 = qfactorit * (radusemin * radusemin * radusemin)
            radusemax3 = qfactorit * (radusemax * radusemax * radusemax)
            radusemin = radlocalmin * selffac
            radusemax = radlocalmax * selffac
            #            radusemin=radlocalmin
            #            radusemax=radlocalmax
            for i in range(1, nq + 1):
                qq = q[i]
                # Weight the sphere form factors by the volume of each sphere
                sff[i, ibin, it] = sff[i, ibin, it] + (radusemax3 * gaussianqr2(qq, radusemax)
                                                       - radusemin3 * gaussianqr2(qq, radusemin))
                # Weight the sphere form factors by the average volume of this sphere
                #               sphereformfactor(i,it)=sphereformfactor(i,it)+qfactor*rhosphere(it)*(radlocalmax3*gaussianqr2(qrmax)-radlocalmin3*gaussianqr2(qrmin))
        sumatoms = f32(0)  # Will need to normalise everything per atom
        for it in range(1, ntype + 1):
            nterms = 0
            sumatomsbin = f32(0.0)
            for ib in range(1, int(m.nsffbintype[it]) + 1):
                nterms = nterms + int(m.ntermssff[ib, it])
                # Q = 0 value corresponds to number of atoms in this type of spherical shell
                sumatomsbin = sumatomsbin + sff[1, ib, it]
                if m.ntermssff[ib, it] > 1:
                    for i in range(1, nq + 1):
                        # Mean sphere form factor for this form factor bin
                        sff[i, ib, it] = sff[i, ib, it] / f32(m.ntermssff[ib, it])
            #            write(6,*) 'inita> sumatomsbin: ',it,cfrac(it),nsffbintype(it),sumatomsbin,nterms
            # Get fractional ntermssff for each atom type
            if nterms > 0:
                # Normalise to number of terms. Don't need atomic fractions here because nterms is already proportional to atom fraction
                sumatoms = sumatoms + molecom.cfrac[it] * sumatomsbin / f32(nterms)
                for ib in range(1, int(m.nsffbintype[it]) + 1):
                    m.fracntermssff[ib, it] = f32(m.ntermssff[ib, it]) / f32(nterms)
            else:
                for ib in range(1, int(m.nsffbintype[it]) + 1):
                    m.fracntermssff[ib, it] = f32(0.0)
            #            write(6,*) 'inita> :',it,(fracntermssff(ib,it),ib=1,nsffbintype(it))
        sumatoms = sumatoms * f32(molecom.natomr) / f32(nmol)  # convert to average number of atoms per molecule
        units.write(6, lw("inita> Average number of atoms per molecule: ", sumatoms))
        if sumatoms > 0.0:
            sfactor = sqrtf(r4(1.0) / sumatoms)  # Eventually everything is normalised to number of atoms in box
            # sphereformfactor=sfactor*sphereformfactor (element by element, float32)
            sff[...] = sfactor * sff
        units.open(10, "formfactors.txt", status="unknown")
        for i in range(1, nq + 1):
            units.write(10, lw(q[i], *[sff[i, 1, it] for it in range(1, ntype + 1)]))
        units.close(10)
        # Now do the intra-molecular terms
        iff = m.intraformfactor
        iffs = m.intraformfactorself
        for im in range(1, nmol + 1):
            scalefac = m.scalefacmol[im]
            # selffac is used because intra-molecular atomic separations can be significant, especially for larger molecules
            #            selffac=scalefacmoluse(im) !May need to use cube root? Or extra smearing factor?
            selffac = m.scalefacmol[im]  # May need to use cube root? Or extra smearing factor?
            for ia in range(int(naf[im]), int(nal[im]) + 1):
                it = itype[ia]
                qfactorit = qfactor * m.rhosphere[it]
                raditmin = m.radminsphere[it] * scalefac
                raditmax = m.radmaxsphere[it] * scalefac
                raditmin3 = qfactorit * (raditmin * raditmin * raditmin)
                raditmax3 = qfactorit * (raditmax * raditmax * raditmax)
                #               raditminself=raditmin
                #               raditmaxself=raditmax
                raditmin = raditmin * selffac
                raditmax = raditmax * selffac
                for ja in range(ia, int(nal[im]) + 1):
                    jt = itype[ja]
                    ip = molecom.npairref[it, jt]
                    m.ntermsintra[ip] = m.ntermsintra[ip] + 1
                    if it == jt:
                        # Total number atoms of this type in pseudo molecule
                        nitterms = m.natomsphere[it] * m.pseudofacmol[it]
                        if nitterms > 0:  # If there are no self terms for this atom, ignore it
                            # Total number of pairs of type it,it in pseudomolecule
                            niitermstotal = nitterms * (nitterms - f32(1))
                            # Number true "self" pairs in the pseudo molecule
                            niitermsself = m.pseudofacmol[it] * m.natomsphere[it] * (m.natomsphere[it] - f32(1))
                            # Fractional contribution from true self terms
                            selfcorrself = niitermsself / (nitterms * nitterms)
                            #                        selfcorr=0.5*(niitermstotal-niitermsself)/(nitterms*nitterms) !Fractional contribution from remainder
                            selfcorr = (niitermstotal - niitermsself) / (nitterms * nitterms)  # Fractional contribution from remainder
                            #                     write(6,*) 'inita> ',it,nitterms,niitermstotal,niitermsself,selfcorrself,selfcorr
                            m.ntermsintraself[it] = m.ntermsintraself[it] + 1
                            for i in range(1, nq + 1):
                                qq = q[i]
                                formfactorit = (raditmax3 * gaussianqr2(qq, raditmax)
                                                - raditmin3 * gaussianqr2(qq, raditmin))
                                iffs[i, it] = iffs[i, it] + selfcorrself * formfactorit * formfactorit
                                iff[i, ip] = iff[i, ip] + selfcorr * formfactorit * formfactorit
                    else:
                        qfactorjt = qfactor * m.rhosphere[jt]
                        radjtmin = m.radminsphere[jt] * scalefac
                        radjtmax = m.radmaxsphere[jt] * scalefac
                        radjtmin3 = qfactorjt * (radjtmin * radjtmin * radjtmin)
                        radjtmax3 = qfactorjt * (radjtmax * radjtmax * radjtmax)
                        radjtmin = radjtmin * selffac
                        radjtmax = radjtmax * selffac
                        for i in range(1, nq + 1):
                            qq = q[i]
                            formfactorit = (raditmax3 * gaussianqr2(qq, raditmax)
                                            - raditmin3 * gaussianqr2(qq, raditmin))
                            formfactorjt = (radjtmax3 * gaussianqr2(qq, radjtmax)
                                            - radjtmin3 * gaussianqr2(qq, radjtmin))
                            iff[i, ip] = iff[i, ip] + formfactorit * formfactorjt
        #         write(6,*) 'inita> ',(ntermsintra(ip),ip=1,npairs)
        #         write(6,*) 'inita> ',(ntermsintraself(ip),ip=1,ntype)
        # Normalise form factors to number of atoms
        for ip in range(1, npairs + 1):
            if m.ntermsintra[ip] > 0 and sumatoms > 0.0:
                sfactor = r4(1.0) / (f32(m.ntermsintra[ip]) * sumatoms)
                #               sfactor=1.0/(ntermsintra(ip))
                for i in range(1, nq + 1):
                    iff[i, ip] = sfactor * iff[i, ip]
        for it in range(1, ntype + 1):
            if m.ntermsintraself[it] > 0 and sumatoms > 0.0:
                ip = molecom.npairref[it, it]  # noqa: F841
                sfactor = r4(1.0) / (f32(m.ntermsintraself[it]) * sumatoms)
                #               sfactor=1.0/(ntermsintra(ip)*sumatoms)
                for i in range(1, nq + 1):
                    iffs[i, it] = sfactor * iffs[i, it]
        units.open(10, "intraformfactors.txt", status="unknown")
        for i in range(1, nq + 1):
            units.write(10, lw(q[i], *[iff[i, ip] for ip in range(1, npairs + 1)]))
        # NB: unit 10 is reopened on another file without a close; the OPEN closes
        # intraformfactors.txt first (units.open does the same).
        units.open(10, "intraformfactorsself.txt", status="unknown")
        for i in range(1, nq + 1):
            units.write(10, lw(q[i], *[iffs[i, it] for it in range(1, ntype + 1)]))
        units.close(10)
    #
    # set z scale identical to x scale - inherited from earlier version
    #
    nz = epsr16.npt1
    for i in range(1, nz + 1):
        z[i] = x[i]
    #      stop
    return iinit, nz


@routine
def j1qr(q, r):
    """epsr_standalone_rev3.f:j1qr -- integral of r sin(Qr)/Q dr (real*8 result).

    Returns j1qr.
    """
    if q > 0.0:
        qr = float(q * r)
        j1qr = (math.sin(qr) - qr * math.cos(qr)) / float(q * q * q)
    else:
        j1qr = float(r * r * r / r4(3.0))
    return j1qr


@routine
def j1qr_mod(q, r):
    """epsr_standalone_rev3.f:j1qr_mod -- j1qr divided by (qr)**3 (real*8 result).

    Returns j1qr_mod.
    """
    qr = float(q * r)
    j1qr_mod = 0.0
    if qr > 0.0:
        j1qr_mod = (math.sin(qr) - qr * math.cos(qr)) / (qr * qr * qr)
    else:
        # 1.0/3.0 is folded in single precision, then converted to real*8
        j1qr_mod = float(r4(1.0) / r4(3.0))
    return j1qr_mod


@routine
def gaussianqr2(q, r):
    """epsr_standalone_rev3.f:gaussianqr2 -- Gaussian form factor for Gaussian density
    distribution across sphere.

    Returns gaussianqr2.
    """
    gaussianqr2 = r4(1.0)
    if r <= 0.0:
        return gaussianqr2  # If sigma is zero, uniform scattering
    qr = q * r
    #      arg=0.64309*qr !ratio of sigma to radius so that r=0 density is same as uniform sphere (sqrt(2/pi)/3)**(1/3)
    arg = r4(0.51042225) * qr  # Gives volume of Gaussian half that of sphere, and about 73% of mass inside specified sphere radius
    arg = r4(0.5) * arg * arg
    if arg < 30.0:
        # The 3 is because the number of atoms in the sphere is 4*pi (already included) / 3 (not included) * R**3.
        gaussianqr2 = expf(-arg)
    else:
        gaussianqr2 = f32(0.0)
    return gaussianqr2


@routine
def lorentzianqr2(q, r):
    """epsr_standalone_rev3.f:lorentzianqr2 -- (Lorentzian) form factor across sphere.

    Returns lorentzianqr2.
    """
    lorentzianqr2 = f32(0.0)
    if r <= 0.0:
        return lorentzianqr2  # If sigma is zero, there can be no scattering
    qr = q * r
    #      arg=0.64309*qr !ratio of sigma to radius so that r=0 density is same as uniform sphere (sqrt(2/pi)/3)**(1/3)
    arg = r4(0.59647) * qr  # Gives volume of Gaussian half that of sphere, and about 73% of mass inside specified sphere radius
    if arg < 30.0:
        # The 3 is because the number of atoms in the sphere is 4*pi (already included) / 3 (not included) * R**3.
        lorentzianqr2 = expf(-arg) / f32(3)
    return lorentzianqr2


@routine
def pofr(delr):
    """epsr_standalone_rev3.f:pofr -- Fourier transform of a uniform sphere (real*4).

    Returns pofr.
    """
    if delr <= 0.0:
        pofr = r4(1.0)
    else:
        delr3 = r4(3.0) / (delr * delr * delr)
        pofr = delr3 * (sinf(delr) - delr * cosf(delr))
    return pofr


#
# initialize g(r) sums
#
@routine
def initg(nq, nz, inpfile, ileninpfile, nsum, ireset):
    """epsr_standalone_rev3.f:initg -- initialize g(r) sums.

    Reads the accumulated g(r) and S(Q) (.g .y .f .s files) when nsum > 0, otherwise
    zeroes them.
    Returns nsum.
    """
    qarray = real_array(mpts)
    rarray = real_array(mpz)
    ierr = 0
    # NB: lptin is an uninitialised local; r_int only assigns it.
    lptin = 0
    #     write(6,*) nsum
    if nsum > 0:
        #
        # read in existing fits in r and Q space
        #
        ext = fstr(fmt("('.g')"), 2)
        fname = fstr(substr(inpfile, 1, ileninpfile) + substr(ext, 1, 2), 256)
        lptin, fname, ierr = auxiliary_routines.r_int(lptin, fname, molecom.npairs, rarray, epsr16.gorinter,
                                                      epsr16.gr2inter, mpz, mpairs, ierr)

        ext = fstr(fmt("('.y')"), 2)
        fname = fstr(substr(inpfile, 1, ileninpfile) + substr(ext, 1, 2), 256)
        lptin, fname, ierr = auxiliary_routines.r_int(lptin, fname, molecom.npairs, rarray, epsr16.gorintra,
                                                      epsr16.gr2intra, mpz, mpairs, ierr)

        ext = fstr(fmt("('.f')"), 2)
        fname = fstr(substr(inpfile, 1, ileninpfile) + substr(ext, 1, 2), 256)
        lptin, fname, ierr = auxiliary_routines.r_int(lptin, fname, molecom.npairs, qarray, epsr16.soqinter,
                                                      epsr16.sq2inter, mpts, mpairs, ierr)

        ext = fstr(fmt("('.s')"), 2)
        fname = fstr(substr(inpfile, 1, ileninpfile) + substr(ext, 1, 2), 256)
        lptin, fname, ierr = auxiliary_routines.r_int(lptin, fname, molecom.npairs, qarray, epsr16.soqintra,
                                                      epsr16.sq2intra, mpts, mpairs, ierr)

    else:

        #
        # zero gorinter
        #
        nsum = 0
        for id in range(1, molecom.npairs + 1):
            for j in range(1, nz + 1):
                epsr16.gorintra[j, id] = f32(0.0)
                epsr16.gr2intra[j, id] = f32(0.0)
                epsr16.gorinter[j, id] = f32(0.0)
                epsr16.gr2inter[j, id] = f32(0.0)
        for in_ in range(1, molecom.npairs + 1):
            for i in range(1, nq + 1):
                epsr16.soqintra[i, in_] = f32(0.0)
                epsr16.sq2intra[i, in_] = f32(0.0)
                epsr16.soqinter[i, in_] = f32(0.0)
                epsr16.sq2inter[i, in_] = f32(0.0)
    return nsum


#
# calculate entropy of distribution
#
@routine
def entrpa(n1, n2, dt, entrop, mpts):
    """epsr_standalone_rev3.f:entrpa -- calculate entropy of distribution dt(n1:n2).

    Returns entrop.
    """
    entrop = f32(0.0)
    yav1 = dt[n1 - 1]
    yav2 = dt[n1]
    ysum = f32(0.0)
    for j in range(n1, n2 + 1):
        yav3 = dt[j + 1]
        ysum = ysum + r4(1.0)
        rj = r4(0.5) * abs(yav3 - yav1)
        pj = r4(0.25) * (r4(2.0) * yav2 + yav1 + yav3)
        dif = yav2 - pj
        difa = abs(dif)
        if difa < r4(0.5) * rj:
            sent = dif * dif / rj
        else:
            sent = difa - r4(0.25) * rj
        entrop = entrop + sent
        yav1 = yav2
        yav2 = yav3
    if ysum > 0.0:
        entrop = entrop / ysum
    return entrop


#
# calculate entropy of average distribution
#
@routine
def entrpb(entrop, nsum):
    """epsr_standalone_rev3.f:entrpb -- calculate entropy of average distribution.

    NB (FORTRAN_BUGS.md section 2, reproduced): the loop runs over the first nfile pairs
    (id=1,nfile), not over all npairs.
    Returns entrop.
    """
    if nsum <= 0:
        entrop = f32(0.0)
    else:
        n1 = 2
        n2 = epsr16.npt2
        entrop = f32(0.0)
        gorinter = epsr16.gorinter
        ycoh = epsr16.ycoh
        entwt = epsr16.entwt
        for id in range(1, epsr16.nfile + 1):
            yav1 = f32(0.0)
            yav2 = gorinter[n1, id] * ycoh[n1]
            for j in range(n1, n2 + 1):
                yav3 = gorinter[j + 1, id] * ycoh[j + 1]
                rj = r4(0.5) * abs(yav3 - yav1)
                pj = r4(0.25) * (r4(2.0) * yav2 + yav1 + yav3)
                dif = yav2 - pj
                difa = abs(dif)
                if difa < r4(0.5) * rj:
                    sent = dif * dif / rj
                else:
                    sent = difa - r4(0.25) * rj
                entrop = entrop + sent * entwt[j]
                yav1 = yav2
                yav2 = yav3
        ysuma = r4(1.0)
        entrop = entrop / f32(nsum) / ysuma
    return entrop


#
# calculate entropy of molecular distribution
#
@routine
def entrpc(entrop, wav):
    """epsr_standalone_rev3.f:entrpc -- calculate entropy of molecular distribution.

    Also sets wav(2,npt2+1,id) and wav(2,npt2+2,id) to ycoh.
    Returns entrop.
    """
    n1 = 2
    n2 = epsr16.npt2
    npt2 = epsr16.npt2
    ycoh = epsr16.ycoh
    entwt = epsr16.entwt
    entrop = f32(0.0)
    for id in range(1, molecom.npairs + 1):
        wav[2, npt2 + 1, id] = ycoh[npt2 + 1]
        wav[2, npt2 + 2, id] = ycoh[npt2 + 2]
        yav1 = f32(0.0)
        yav2 = wav[2, n1, id]
        for j in range(n1, n2 + 1):
            yav3 = wav[2, j + 1, id]
            rj = r4(0.5) * abs(yav3 - yav1)
            pj = r4(0.25) * (r4(2.0) * yav2 + yav1 + yav3)
            dif = yav2 - pj
            difa = abs(dif)
            if difa < r4(0.5) * rj:
                sent = dif * dif / rj
            else:
                sent = difa - r4(0.25) * rj
            entrop = entrop + sent * entwt[j]
            yav1 = yav2
            yav2 = yav3
    return entrop


#
# sum g(r) and error bars
#
@routine
def sumgr(nz, wav):
    """epsr_standalone_rev3.f:sumgr -- sum g(r) and error bars.

    Uses fc1, fc2 and the work arrays of common/sumvalues/.
    """
    zerolimit = dble_array(mpairs)
    sv = sumvalues
    add1 = sv.add1
    add2 = sv.add2
    sumtot1 = sv.sumtot1
    sumtot2 = sv.sumtot2
    errort = sv.errort
    ycoh = epsr16.ycoh
    gorintra = epsr16.gorintra
    gorinter = epsr16.gorinter
    gr2intra = epsr16.gr2intra
    gr2inter = epsr16.gr2inter
    for id in range(1, molecom.npairs + 1):
        atom1 = fstr(potential_par.tlab1[potential_par.irefpcof[id]], 4)
        atom2 = fstr(potential_par.tlab2[potential_par.irefpcof[id]], 4)
        #      write(6,*) 'sumgr> fc1 and fc2: ',id,nz,npt1,npt2,' ',atom1,' ',atom2
        ilent1 = findex_padded(atom1, " ", 4) - 1
        ilent2 = findex_padded(atom2, " ", 4) - 1
        titlepsf = fstr(substr(atom1, 1, ilent1) + "-" + substr(atom2, 1, ilent2), 9)  # noqa: F841
        ilent = ilent1 + ilent2 + 1  # noqa: F841
        epsr16.entsig[1, id] = f32(0.0)
        zerolimit[id] = 0.0
        for j in range(2, nz + 1):
            if ycoh[j] > 0.0:
                add1[1] = wav[1, j, id] / ycoh[j]
                add1[2] = wav[2, j, id] / ycoh[j]
            else:
                add1[1] = 0.0
                add1[2] = 0.0
            sumtot1[1] = gorintra[j, id]
            sumtot1[2] = gorinter[j, id]
            sumtot2[1] = gr2intra[j, id]
            sumtot2[2] = gr2inter[j, id]
            #      write(6,*) 'sumgr> ',id,j,nz,sumtot1(1),sumtot1(2),sumtot2(1),sumtot2(2)
            fc1 = float(sv.fc1)
            fc2 = float(sv.fc2)
            for k in range(1, 3):
                add2[k] = add1[k] * add1[k]
                sumtot2[k] = sumtot2[k] * sumtot2[k] + sumtot1[k] * sumtot1[k]
                sumtot1[k] = fc2 * sumtot1[k] + fc1 * add1[k]
                sumtot2[k] = fc2 * sumtot2[k] + fc1 * add2[k]
                errort[k] = sumtot2[k] - sumtot1[k] * sumtot1[k]
                if errort[k] > 0.0:
                    errort[k] = dsqrt(errort[k])
                else:
                    errort[k] = 0.0
            gorintra[j, id] = sumtot1[1]
            gorinter[j, id] = sumtot1[2]
            gr2intra[j, id] = errort[1]
            gr2inter[j, id] = errort[2]
            zerolimit[id] = zerolimit[id] + float(ycoh[j] * (gorinter[j, id] - r4(1.0)))
            epsr16.entsig[j, id] = zerolimit[id]
        #            write(6,'(a,i4,1x,a,1x,a,e12.5)') 'sumgr> Atom pair ',id,titlepsf(1:ilent),' zero limit: ',zerolimit(id)


#
# sum S(Q) and error bars
#
@routine
def sumsq(nq, qwind):
    """epsr_standalone_rev3.f:sumsq -- sum S(Q) and error bars.

    Uses fc1, fc2 and the work arrays of common/sumvalues/.
    """
    sv = sumvalues
    add1 = sv.add1
    add2 = sv.add2
    sumtot1 = sv.sumtot1
    sumtot2 = sv.sumtot2
    errort = sv.errort
    soqintra = epsr16.soqintra
    soqinter = epsr16.soqinter
    sq2intra = epsr16.sq2intra
    sq2inter = epsr16.sq2inter
    for id in range(1, molecom.npairs + 1):
        for j in range(1, nq + 1):
            add1[1] = qwind[1, j, id]
            add1[2] = qwind[2, j, id]
            sumtot1[1] = soqintra[j, id]
            sumtot1[2] = soqinter[j, id]
            sumtot2[1] = sq2intra[j, id]
            sumtot2[2] = sq2inter[j, id]
            fc1 = float(sv.fc1)
            fc2 = float(sv.fc2)
            for k in range(1, 3):
                add2[k] = add1[k] * add1[k]
                sumtot2[k] = sumtot2[k] * sumtot2[k] + sumtot1[k] * sumtot1[k]
                sumtot1[k] = fc2 * sumtot1[k] + fc1 * add1[k]
                sumtot2[k] = fc2 * sumtot2[k] + fc1 * add2[k]
                errort[k] = sumtot2[k] - sumtot1[k] * sumtot1[k]
                if errort[k] > 0.0:
                    errort[k] = dsqrt(errort[k])
                else:
                    errort[k] = 0.0
            soqintra[j, id] = sumtot1[1]
            soqinter[j, id] = sumtot1[2]
            sq2intra[j, id] = errort[1]
            sq2inter[j, id] = errort[2]
            #                  write(25,'(a,2(1x,i7),6(1x,e13.6))') 'sumsq-1> ',id,j,add1(1),add1(2),sumtot1(1),sumtot1(2)


#
# sum the potential energies and write a file
#
@routine
def sumpe(nz, z, inpfile, ileninpfile):
    """epsr_standalone_rev3.f:sumpe -- sum the potential energies and write the .terg file."""
    toteref = f32(0.0)
    toteemp = f32(0.0)
    totetot = f32(0.0)
    rmaxpot = r4(0.5) * clustercom.dmax
    for id in range(1, molecom.npairs + 1):
        epsr16.energyref[id] = f32(0.0)
        epsr16.energyemp[id] = f32(0.0)
        epsr16.energytot[id] = f32(0.0)
        it = molecom.n1extr[id]
        jt = molecom.n2extr[id]
        cprod = molecom.cfrac[it] * molecom.cfrac[jt]
        if it != jt:
            cprod = cprod * f32(2)
        j = 2
        #
        # only integrate up to half the box dimension since g(r) may not be
        # reliable beyond that distance
        #
        while j <= nz and z[j] <= rmaxpot:
            wt = cprod * epsr16.ycoh[j] * epsr16.gorinter[j, id]
            epsr16.energyref[id] = epsr16.energyref[id] + wt * potential_par.ecoul[j, id]
            epsr16.energyemp[id] = epsr16.energyemp[id] + wt * potential_par.xavs[j, id]
            j = j + 1
        epsr16.energytot[id] = epsr16.energyref[id] + epsr16.energytot[id]
        toteref = toteref + epsr16.energyref[id]
        toteemp = toteemp + epsr16.energyemp[id]
    totetot = toteref + toteemp
    #
    # write these results to a file
    #
    ext = fstr(fmt("('.terg')"), 5)
    fname = fstr(substr(inpfile, 1, ileninpfile) + substr(ext, 1, 5), 256)
    units.open(20, fname, status="unknown")
    units.write(20, fmt(_F_SUMPE_200, toteref, toteemp, totetot))
    for id in range(1, molecom.npairs + 1):
        units.write(20, fmt(_F_SUMPE_201, Ch(molecom.type[molecom.n1extr[id]], 3),
                            Ch(molecom.type[molecom.n2extr[id]], 3)))
        units.write(20, fmt(_F_SUMPE_200, epsr16.energyref[id], epsr16.energyemp[id], epsr16.energytot[id]))
    units.close(20)


@routine
def initpt(ireset, rminpt, rmaxpt, expecf, polar, smoos, rminp, rbroad):
    """epsr_standalone_rev3.f:initpt -- truncation functions, zero the EPs, then generate
    the reference potential and the empirical potential."""
    units.write(6, lw("initpt> rminpt,rmaxpt: ", rminpt, rmaxpt))
    pi = r4(3.141592653)  # noqa: F841
    if molecom.power > 0.0:
        ratt = powf(r4(1.0e10), r4(1.0) / molecom.power)  # noqa: F841
    #
    # maximum value of a sum - this it to prevent overflows
    #
    summax = expf(r4(25.0))  # noqa: F841
    #
    # generate renormalized x scale and potential damping function at low r
    #
    x = epsr16.x
    npt = epsr16.npt
    for j in range(1, npt + 1):
        epsr16.xstar[j] = x[j] / clustercom.dmax
    widthptc = max(r4(3.5), r4(0.2) * rmaxpt)  # noqa: F841
    for j in range(1, npt + 1):
        xx = x[j]

        # Lennard-Jones and EP truncation function

        potential_par.wfac[j] = ato_manipulation_rev4.truncate(xx, rminpt, rmaxpt)

        # Coulomb truncation function

        #         wfacc(j)=truncatecweeks(xx,widthptc)
        potential_par.wfacc[j] = ato_manipulation_rev4.truncatec(xx, rmaxpt)

    potential_par.wfac[npt] = f32(0.0)
    potential_par.wfacc[npt] = 0
    #
    # set the empirical potentials to zero, and test for the maximum distance required for the potential
    #
    for ip in range(1, molecom.npairs + 1):
        irefp = potential_par.irefpcof[ip]
        for j in range(1, npt + 1):
            if x[j] <= abs(potential_par.rminex[irefp] + r4(4.0) * expecf) and x[j] > x[epsr16.nptmax]:
                epsr16.nptmax = j
            potential_par.xav[j, ip] = f32(0.0)
            potential_par.xavs[j, ip] = f32(0.0)
            potential_par.xav1[j, ip] = f32(0.0)
            potential_par.xavs1[j, ip] = f32(0.0)
    #
    # Generate reference potential from parameters
    #
    generate_rp()  # noqa: F821
    #
    # generate the EP from the site-site coefficients
    #
    generate_ep(False)  # noqa: F821

    return


@routine
def ptgnsq(nsum, nq, q, wav, data, qwind, rminpt, rmaxpt, szeros, overallfac, gfac, rminp, potnoi, wtentp,
           rbroad, sigma0, nsmoop):
    """epsr_standalone_rev3.f:ptgnsq -- refine the empirical potential from the
    differences between the data and the simulated intermolecular partials.

    qwind is qwind(mpts,mpairs) here (the caller passes soqintertemp).  Writes delfit.dat.
    """
    # NB: qtemp is given one extra column: the delfit.dat write tests qtemp(i,nfile+1)
    # (the outer loop variable id after the loop); the local array is zero there.
    qtemp = real_array(mpts, mfile + 1)
    fittemp = real_array(mpts, mfile)
    dt = real_array(mpz)
    dts = real_array(mpz)
    rfactorep = real_array(mfile)
    nmin = int_array(mfile)
    nmini = int_array(mfile)

    summax = expf(r4(25.0))  # noqa: F841
    pi = PI4ATAN
    # Set the minimum Q for fitting the potential
    if epsr16.qminp <= 0.0:
        qminptemp = epsr16.qminm
    else:
        qminptemp = epsr16.qminp
    stepcl = epsr16.cellst  # noqa: F841
    xavmax = f32(0.0)  # noqa: F841
    qstep = q[3] - q[2]
    qrat = r4(2.0) * pi / qstep  # noqa: F841
    expmax = r4(25.0)  # noqa: F841
    potadd = potential_par.potfac  # Limit the amplitude of the perturbation
    npairs = molecom.npairs
    for id in range(1, epsr16.nfile + 1):
        nmin[id] = 2
        nmini[id] = 2

        # Form the simulated fit for this dataset

        # If X-ray data, need to apply the Q-dependent weights

        ix = weights.xfile[id]

        # Set up the gradient for extrapolating to Q=0

        nqmn = epsr16.nqmin[id]
        if szeros[id] != 0.0 and nqmn > 1:

            grads0 = (data[nqmn, id] - szeros[id]) / (q[nqmn] - q[1])
            conss0 = szeros[id]

        else:

            grads0 = f32(0.0)
            conss0 = f32(0.0)

        for i in range(1, nq + 1):
            #
            # Sum the intermolecular partials for this dataset at this q value
            #
            sum = f32(0.0)
            if ix > 0:
                for ic in range(1, npairs + 1):
                    sum = f32(float(sum) + x_rayf0.weightx[i, ic, ix] * float(qwind[i, ic]))
            else:
                for ic in range(1, npairs + 1):
                    sum = f32(float(sum) + weights.cwt[id, ic] * float(qwind[i, ic]))

            # Calculate the difference between data and fit. For Q-values below qmin we take the first
            # valid data point as the data below Qmin, unless szeros is non-zero.

            if q[i] < max(qminptemp, q[nqmn]):
                if szeros[id] != 0.0:
                    qtemp[i, id] = conss0 + grads0 * (q[i] - q[1]) - sum
                else:
                    qtemp[i, id] = f32(0.0)
            elif i >= nqmn and i <= epsr16.nqmax[id]:
                qtemp[i, id] = data[i, id] - sum
            else:
                qtemp[i, id] = f32(0.0)

    # Format for diagnostic files

    formatout = fstr(fmt(_F_PTGNSQ_332, "(", (2 * mfile + 1), "(1x,e12.5))"), 14)

    # If X-ray data are present we need to reconstruct the differences using the Q=0
    # weights. This is to avoid having to deconvolute the X-ray form factors in r-space.
    # This is possibly a controversial practice since it involves using the inversion
    # matrix and therefore eliminates any aspects of the data that are incompatible between
    # datasets. For compatibility the same is done for datasets even when there are no
    # X-ray data present.
    # (the commented-out reconstruction code, Fortran lines 3747-3866, is omitted)

    if potential_par.feedback > 0.0 and potadd > 0.0:
        for id in range(1, epsr16.nfile + 1):
            #
            # determine the minimum bin for the intermolecular potential
            #
            for j in range(1, epsr16.npt2 + 1):
                if epsr16.x[j] < rminp[id]:
                    nmini[id] = j + 1
            units.write(6, fmt(_F_PTGNSQ_59, nmini[id]))
            #
            # parameters for generating the potential
            #
            # nfreq represents the number of times during the fit the coefficients are
            # smoothed
            #
            nfreq = 10
            nterms = int(potential_par.npcoeff[id])
            stepp = potential_par.pstepd[id]
            #
            # decay=sigma2(id) controls the width of the Gaussians used to fit the data.
            # Typical value should be of order 2*pi/qmax, or a bit less. Obviously reducing
            # sigma2 will produce more detail in the fitted s(q)/g(r), but may lead to
            # greater possibility of truncation artifacts being reproduced.
            #
            decayp = epsr16.sigma2[id]
            acceptr = potential_par.paccept
            #
            # set up the transform matrix from power-exponentials to Q-space
            #
            for ic in range(1, nterms + 1):
                potential_par.pwtold[ic] = potential_par.pcoefd[ic, id]
                potential_par.pwt[ic] = potential_par.pwtold[ic]
            # NB: ncfirst, rejrat and chisqnew are uninitialised locals; the callees only
            # assign them.
            ncfirst = 0
            rejrat = f32(0.0)
            chisqnew = f32(0.0)
            if potential_par.gaussian:
                ncfirst = calc_gaussq(epsr16.nqmax[id], nterms, potential_par.pdstep, potential_par.pdmax,  # noqa: F821
                                      decayp, epsr16.rho, q, potential_par.expsq, f32(0.0), ncfirst)
            else:
                ncfirst = calc_pwrexpq(epsr16.nqmax[id], nterms, potential_par.pdmin, potential_par.pdstep,  # noqa: F821
                                       potential_par.pdmax, decayp, epsr16.rho, q, potential_par.expsq, f32(0.0),
                                       ncfirst)
            #
            # calculate perturbation to potential
            #
            # form the difference between data and fit.
            #
            for i in range(1, epsr16.nqmax[id] + 1):
                dt[i] = -(qtemp[i, id])
            #
            # fit a series of power-exponentials in q-space to this difference
            #
            if szeros[id] == 0.0:
                stepp, acceptr, rejrat, chisqnew = fit_expansion(  # noqa: F821
                    dt, potential_par.expsq, dts, potential_par.pwt, potential_par.pwtold,
                    stepp, acceptr, rejrat, chisqnew, 1, epsr16.nqmax[id], 1, nterms,
                    potential_par.npits, nfreq)
            else:
                stepp, acceptr, rejrat, chisqnew = fit_expansion(  # noqa: F821
                    dt, potential_par.expsq, dts, potential_par.pwt, potential_par.pwtold,
                    stepp, acceptr, rejrat, chisqnew, 1, epsr16.nqmax[id], 1, nterms,
                    potential_par.npits, nfreq)
            # Save fit

            for i in range(1, epsr16.nqmax[id] + 1):
                fittemp[i, id] = dts[i]

            #                 write(6,*) decayp,chisqnew
            #
            # save revised step and chisq values
            #
            potential_par.pstepd[id] = stepp
            potnoi[id] = chisqnew
            #
            # save the new fluctuation coefficients and add them to the empirical
            # potential coefficients
            #
            rate = molecom.simt  # noqa: F841
            pfact = potadd
            for ic in range(2, nterms + 1):
                #
                # save the new fluctuation coefficients
                #
                potential_par.pcoefd[ic, id] = potential_par.pwtold[ic]
                #
                # increment the site-site potential coefficients
                #
                for ip in range(1, npairs + 1):
                    irefp = potential_par.irefpcof[ip]
                    potential_par.pcoefp[ic, irefp] = (
                        float(potential_par.pcoefp[ic, irefp])
                        + float(pfact * potential_par.efacp[irefp]) * weights.cwtpot[ip, id]
                        * float(potential_par.pwtold[ic]))
            rfactorep[id] = f32(0.0)
            if szeros[id] == 0.0:
                for i in range(epsr16.nqmin[id], nq + 1):
                    if abs(qtemp[i, id]) > 0.0:
                        rfactorep[id] = rfactorep[id] + abs(fittemp[i, id] + qtemp[i, id]) / abs(qtemp[i, id])
                    else:
                        rfactorep[id] = rfactorep[id] + abs(fittemp[i, id] + qtemp[i, id])
                rfactorep[id] = rfactorep[id] / f32(nq - epsr16.nqmin[id] + 1)
            else:
                for i in range(1, nq + 1):
                    if abs(qtemp[i, id]) > 0.0:
                        rfactorep[id] = rfactorep[id] + abs(fittemp[i, id] + qtemp[i, id]) / abs(qtemp[i, id])
                    else:
                        rfactorep[id] = rfactorep[id] + abs(fittemp[i, id] + qtemp[i, id])
                rfactorep[id] = rfactorep[id] / f32(nq)
            units.write(6, lw("ptgnsq> ", id, potadd, pfact, rfactorep[id]))
        # after the do loop the Fortran id is nfile+1
        id = epsr16.nfile + 1
        units.open(10, "delfit.dat", status="unknown")
        nfile = epsr16.nfile
        for i in range(1, nq + 1):
            # NB: tests qtemp(i,id) with id = nfile+1 (zero), so the second branch is taken
            if abs(qtemp[i, id]) > 0.0:
                vals = []
                for jd in range(1, nfile + 1):
                    vals += [fittemp[i, jd], -qtemp[i, jd],
                             abs(fittemp[i, jd] + qtemp[i, jd]) / abs(qtemp[i, jd])]
                units.write(10, fmt(formatout, q[i], *vals))
            else:
                vals = []
                for jd in range(1, nfile + 1):
                    vals += [fittemp[i, jd], -qtemp[i, jd], abs(fittemp[i, jd] + qtemp[i, jd])]
                units.write(10, fmt(formatout, q[i], *vals))
        units.close(10)

        # Regenerate EP

        #           call generate_ep(.true.)
        generate_ep(False)  # noqa: F821

    #
    # if no. of iterations to calculate potential is large (ireset=1) reset
    # this value to input value for all subsequent iterations
    #
    potential_par.npits = potential_par.npitss
    return


@routine
def repulsivept(potadd):
    """epsr_standalone_rev3.f:repulsivept -- regenerates the reference potential including
    new values of the exponential repulsive terms."""
    #
    # modify the repulsive terms as required
    #
    if potadd > 0.0:
        pfact = potadd  # noqa: F841
        simt = molecom.simt
        for ip in range(1, molecom.npairs + 1):
            irefp = potential_par.irefpcof[ip]
            units.write(6, lw("repulsivept potex> ", potential_par.potex[ip], potential_par.coordex[irefp]))
            if potential_par.rminex[irefp] > 0.0:
                pchange = potential_par.potex[ip] - potential_par.coordex[irefp]
                # Limit the value of pcoefex becoming absurdly small
                if pchange == 0.0 and abs(potential_par.pcoefex[irefp]) > r4(0.0001):
                    #
                    # if potex(ic) is zero, then reduce the current pcoefex(ic) by 10% to
                    # ensure this term does not dominate the simulation
                    #
                    potential_par.potex[ip] = r4(-0.1)
                    potential_par.pcoefex[irefp] = potential_par.pcoefex[irefp] * (r4(1.0) + potential_par.potex[ip])
                else:
                    # Restrict changes to less than 0.5% of simt
                    if potential_par.coordex[irefp] > 1.0:
                        # For larger coordination numbers, the change is expressed as a fraction of the requested coordination number
                        pchange = pchange / potential_par.coordex[irefp]
                    pchange = sign(min(r4(0.005) * simt, r4(0.025) * simt * abs(pchange)), pchange)
                    #                       pchange=sign(0.001*simt,pchange)
                    potential_par.pcoefex[irefp] = potential_par.pcoefex[irefp] + pchange
            # Repeat for pcoefexmax
            potential_par.pcoefexmax[irefp] = potential_par.potexmax[ip]
            # (commented-out pcoefexmax adjustment, Fortran lines 4046-4058, omitted)
        #
        # regenerate the reference potential
        #
        generate_rp()  # noqa: F821

    return


#
# write a new input file
#
@routine
def norma(ireset, itimes, ntimes):
    """epsr_standalone_rev3.f:norma -- write current data to the .inpa file (input for the
    next run), store the EP parameters back in the symbol table, write the .uni file."""
    #
    # WRITE CURRENT DATA TO INPUT FILE IN A FORM SUITABLE TO BE READ AS NEXT
    # INPUT FILE
    #
    if itimes == ntimes:
        extfile = fstr(".inpa", 10)
        fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(extfile, 1, 5), 256)
        openf(14, fname)
        units.write(14, lw(epsrinput_par.nisum, epsrinput_par.nesum, epsrinput_par.nesumt, epsrinput_par.ensum,
                           epsrinput_par.ensum2, epsrinput_par.prsum, epsrinput_par.prsum2))
        #
        # Write the data potential parameters
        #
        for id in range(1, epsr16.nfile + 1):
            units.write(14, lw(id))
            units.write(14, lw(potential_par.npcoeff[id], potential_par.pstepd[id], epsr16.sigma2[id]))
            units.write(14, lw(*[potential_par.pcoefd[ic, id] for ic in range(1, int(potential_par.npcoeff[id]) + 1)]))
        nnew = epsr16.nfile + molecom.npairs
        units.write(14, lw(potential_par.feedback, *[potential_par.efacm[j] for j in range(1, epsr16.nfile + 1)],
                           weights.nqinv))
        for iq in range(1, weights.nqinv + 1):
            for ip in range(1, molecom.npairs + 1):
                units.write(14, lw(weights.qinv[iq], molecom.n1extr[ip], molecom.n2extr[ip]))
                units.write(14, lw(*[weights.cwttq[iq, ip, id] for id in range(1, nnew + 1)]))
        units.close(14)
        #
        # write the site-site potential parameters
        #
        # first save the current internal values to their symbol counterparts
        #
        sv = symbols.symbol_value
        iref = file_and_directory_commands.ireff(iresetnumbers.iqminref, 1)
        iref = iref + 1
        sv[iref] = fstr(lw(epsrinput_par.roverlap), 256)
        iref = iref + 1
        sv[iref] = fstr(lw(epsrinput_par.rminfac), 256)
        iref = iref + 1
        sv[iref] = fstr(lw(epsrinput_par.rminpt), 256)
        iref = iref + 1
        sv[iref] = fstr(lw(epsrinput_par.rmaxpt), 256)
        iref = iref + 1
        sv[iref] = fstr(lw(epsrinput_par.rbroad), 256)
        iref = iref + 1
        sv[iref] = fstr(fmt("(f7.3,2x,f7.3)", epsrinput_par.expecf, epsrinput_par.expecfpwr), 256)
        # reppottype has been assigned at input, so skip over it
        iref = iref + 1
        # addpotential has been assigned at input, so skip over it
        iref = iref + 1
        iref = iref + 1
        sv[iref] = fstr(lw(epsrinput_par.rcharge), 256)
        iref = iref + 1
        sv[iref] = fstr(lw(molecom.power), 256)
        iref = iref + 1
        sv[iref] = fstr(lw(potential_par.ncoeffp), 256)
        iref = iref + 1
        sv[iref] = fstr(lw(potential_par.pdmax), 256)
        iref = iref + 1
        sv[iref] = fstr(lw(potential_par.pdstep), 256)
        iref = iref + 1
        sv[iref] = fstr(lw(potential_par.npitss), 256)
        iref = iref + 1
        sv[iref] = fstr(lw(potential_par.paccept), 256)
        # Gaussian has been assigned at input, so skip over it
        iref = iref + 1
        iref = iref + 1
        sv[iref] = fstr(fmt("(f8.5)", potential_par.psigma2read), 256)
        iref = iref + 1
        sv[iref] = fstr(fmt("(f8.5)", potential_par.refpotfac), 256)
        #
        # Write the uniform distribution to a file called EPSR.uni
        #
        extfile = fstr(".uni", 10)
        fname = fstr(substr(symbols.fnamerun, 1, symbols.ilenfnamerun) + substr(extfile, 1, 4), 256)
        openf(10, fname)
        vals = []
        for j in range(1, epsr16.npt1 + 1):
            vals += [epsr16.x[j], epsr16.ycoh[j]]
        units.write(10, fmt("(2(1x,e13.6))", *vals))
        units.close(10)
    return


@routine
def openf(nunit, fname):
    """epsr_standalone_rev3.f:openf -- open a file with status 'unknown', else 'new',
    else stop."""
    stat = fstr("unknown", 10)
    ierr = 0
    try:
        units.open(nunit, fname, status=stat)
    except FortranIOError as e:
        ierr = e.iostat
    if ierr == 0:
        return
    stat = fstr("new", 10)
    ierr = 0
    try:
        units.open(nunit, fname, status=stat)
    except FortranIOError as e:
        ierr = e.iostat
    if ierr == 0:
        return
    units.write(6, fmt(_F_OPENF_100, Ch(fname, 256), Ch(stat, 10)))
    raise FortranStop(0)


@routine
def get_wts(q, nq, fnamew, ireset):
    """epsr_standalone_rev3.f:get_wts -- assign neutron (or X-ray) weightings to the
    individual atom pairs and (re)generate the inverted weights matrix when needed."""
    #
    # assign neutron weightings to individual atom pairs - these are initially
    # all set to zero
    #
    wread = real_array(2)
    axt = dble_array(5)
    bxt = dble_array(5)
    cxt = 0.0
    alphaxt = f32(0.0)
    deltaxt = f32(0.0)
    a1 = ""
    a2 = ""
    # NB: normx and iafformaff are uninitialised locals; both are read (or deduced)
    # before use.
    normx = 0
    iafformaff = 0
    nrho = 0
    npairs = molecom.npairs
    nfile = epsr16.nfile
    #
    # set some accumulators to zero
    #
    for ic in range(1, npairs + 1):
        potential_par.uwt[ic] = f32(0.0)

    # Set the X-ray form factors to zero

    wxinit(nq)

    # Number of X-ray files. This should not exceed mxfile, otherwise only the last one input will
    # be used.
    weights.nxfile = 0

    for in_ in range(1, nfile + 1):
        for ic in range(1, npairs + 1):
            for j in range(1, 3):
                potential_par.weight[j, in_, ic] = f32(0.0)
        units.open(10, fnamew[in_], status="old")

        # decide whether this is an X-ray or neutron wts file

        ilenxray1 = findex(fnamew[in_], "X-ray.wts")
        ilenxray2 = findex(fnamew[in_], "x-ray.wts")
        ilenxray3 = findex(fnamew[in_], "XWTS.wts")
        ilenxray = max(ilenxray1, ilenxray2, ilenxray3)

        if ilenxray > 0:

            #
            # read factor for dataset to convert it to per atom. For X-rays this factor will decide
            # whether data are normalised to the single atom scattering (1) or square of the mean
            # scattering length (2).
            #
            commands.sometext = fstr(units.read_line(10), 3020)
            commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                          commands.nstsize, commands.nworddim)

            if commands.nwords == 3:

                r = read_internal(substr(commands.sometext, commands.ncf[1], commands.ncl[commands.nwords]))
                epsr16.rhofac[in_] = r.real(epsr16.rhofac[in_])
                normx = r.int(normx)
                iafformaff = r.int(iafformaff)

            else:

                # input from old epsrwtsx

                r = read_internal(substr(commands.sometext, commands.ncf[1], commands.ncl[commands.nwords]))
                epsr16.rhofac[in_] = r.real(epsr16.rhofac[in_])

                # Deduce the values of rhofac(in), normx and iafformaff from this one value

                # nrho + means data are per atom, -ve data are per molecule

                if epsr16.rhofac[in_] < 0:

                    nrho = 2

                else:

                    nrho = 1

                # iafformaff is 1 means use IAFFs for normalisation, otherwise use MAFFs

                iafformaff = ifix((abs(epsr16.rhofac[in_]) + f32(1)) / f32(2))

                # In the old version normx = 1 means normalise to single atom scattering, otherwise normalise to sum of Fs
                # all squared. normx = 0 was not an option

                normx = ifix(abs(epsr16.rhofac[in_]) - f32(2 * (iafformaff - 1)))

                # However in the new version, normx = 1 means normalise to sum of Fs all squared, and normx = 2 means normalise to single atom scattering

                if normx == 1:
                    normx = 2
                else:
                    normx = 1
                units.write(6, lw("get_wts> ", epsr16.rhofac[in_], nrho, normx, iafformaff))

                epsr16.rhofac[in_] = f32(nrho)

            if weights.nxfile < x_rayf0.mxfile:
                weights.nxfile = weights.nxfile + 1
                nxfile = weights.nxfile

                # Save the normalisation of the X-ray data.

                x_rayf0.datanorm[nxfile] = ifix(epsr16.rhofac[in_])
                epsr16.rhofac[in_] = r4(1.0)
                x_rayf0.normxx[nxfile] = normx
                x_rayf0.iafformaffx[nxfile] = iafformaff
                units.write(6, lw("get_wts> ", nxfile, x_rayf0.datanorm[nxfile]))

                # Setup pointers to the X-ray form factor arrays

                weights.xfile[in_] = nxfile

                # And to the dataset these weights correspond to.

                x_rayf0.datapoint[nxfile] = in_

                units.write(6, lw("get_wts> ", in_, nxfile, weights.xfile[in_], x_rayf0.datapoint[nxfile]))

                # Read form factor coefficients for this atom.

                a1, a2, ierr = _read_a1a2(_F_GETWTS_108, a1, a2, 4)
                cxt, alphaxt, deltaxt, ierr = _read_xray_coeffs(cxt, axt, bxt, alphaxt, deltaxt)
                while ierr == 0:

                    # Assign these X-ray coefficients if the specified atom type matches an atom in the
                    # .ato file.  Any atoms which are not assigned form factor coefficients are
                    # automatically assumed not to contribute to this diffraction pattern.
                    #
                    wxtype(a1, axt, bxt, cxt, alphaxt, deltaxt, weights.nxfile, q, nq)
                    a1, a2, ierr = _read_a1a2(_F_GETWTS_108, a1, a2, 4)
                    cxt, alphaxt, deltaxt, ierr = _read_xray_coeffs(cxt, axt, bxt, alphaxt, deltaxt)

            else:

                units.write(6, fmt(_F_GETWTS_109, x_rayf0.mxfile))

        else:

            r = units.list_read(10)
            epsr16.rhofac[in_] = r.real(epsr16.rhofac[in_])
            weights.xfile[in_] = 0

            # Read partial structure factor assignment and neutron weighting for
            # intramolecular and inter molecular terms.
            #
            a1, a2, ierr = _read_a1a2(_F_GETWTS_107, a1, a2, 3)
            ierr = _read_wread(wread)
            while ierr == 0:
                #
                # assign this neutron weighting if it matches a pair of atoms in the molecule
                # Pairs which are not assigned a weighting by this input are left as zero.
                #
                wtype(a1, a2, wread, in_)
                a1, a2, ierr = _read_a1a2(_F_GETWTS_107, a1, a2, 3)
                ierr = _read_wread(wread)

        units.close(10)

    # Prepare the X-ray weighting factors if present.

    if weights.nxfile > 0:
        wxcreate(q, nq)

    cfraccg = molecom.cfraccg
    for in_ in range(1, nfile + 1):

        # sum up intermolecular weights and modulus of these for this file

        potential_par.wsum[in_] = f32(0.0)
        epsr16.xwt[in_] = f32(0.0)
        for ic in range(1, npairs + 1):
            it = molecom.n1extr[ic]
            jt = molecom.n2extr[ic]
            cprod = cfraccg[it] * cfraccg[jt]
            if it != jt:
                cprod = r4(2.0) * cprod
            weights.ciwt[in_, ic] = cprod * potential_par.weight[1, in_, ic]
            weights.cwt[in_, ic] = cprod * potential_par.weight[2, in_, ic]
            potential_par.wsum[in_] = float(potential_par.wsum[in_]) + abs(weights.cwt[in_, ic])
            epsr16.xwt[in_] = float(epsr16.xwt[in_]) + weights.cwt[in_, ic]

    # Repeat this for the X-ray weights if present

    if weights.nxfile > 0:

        for ix in range(1, weights.nxfile + 1):
            for ip in range(1, npairs + 1):
                it = molecom.n1extr[ip]
                jt = molecom.n2extr[ip]
                cprod = cfraccg[it] * cfraccg[jt]
                if it != jt:
                    cprod = r4(2.0) * cprod
                for iq in range(1, nq + 1):
                    x_rayf0.weightx[iq, ip, ix] = float(cprod) * x_rayf0.weightx[iq, ip, ix]

        # Set up the Q-values on which to calculate the inversion matrix.

        qinvstep = r4(2.0)
        nqinvstep = ifix(qinvstep / (q[2] - q[1]))

    else:

        # If no X-ray data, only need to calculate inversion at one Q value.

        nqinvstep = nq
    #
    # if ireset.eq.1, then form the inversion of the coefficients cwtpot
    #
    # Also calculate the inversion if any of the factors has changed
    #
    iflag = 0
    if ireset == 1 or abs(potential_par.feedback - potential_par.feedbackold) > r4(0.001):
        iflag = 1
    else:
        in_ = 1
        while iflag == 0 and in_ <= nfile:
            if abs(potential_par.efacm[in_] - potential_par.efacold[in_]) > r4(0.001):
                iflag = 1
            in_ = in_ + 1
    if iflag != 0:

        # Set the initial values of the inversion matrix
        nnew = npairs + nfile
        for iff in range(1, nnew + 1):
            for ip in range(1, npairs + 1):
                #
                # for the partials from the simulation, assume the initial
                # transform matrix is diagonal, as it would be if there were no
                # data
                #
                ifref = iff - nfile
                if ifref == ip:
                    weights.cwtt[ip, iff] = 1.0
                else:
                    weights.cwtt[ip, iff] = 0.0
        #
        # generate weights to be used for potential inversion at this Q value
        #

        weights.nqinv = 0
        for iqinv in range(1, nq + 1, nqinvstep):

            if weights.nqinv < weights.mqinv:
                weights.nqinv = weights.nqinv + 1
                weights.qinv[weights.nqinv] = iqinv

                units.write(6, lw("get_wts> ", weights.nqinv, q[weights.qinv[weights.nqinv]], iqinv, q[iqinv]))

                for in_ in range(1, nfile + 1):

                    ix = weights.xfile[in_]
                    if ix > 0:
                        for ic in range(1, npairs + 1):
                            it = molecom.n1extr[ic]
                            jt = molecom.n2extr[ic]
                            if cfraccg[it] > 0.0 and cfraccg[jt] > 0.0:
                                # Correct for coarse graining if present
                                cprodrat = molecom.cfrac[it] * molecom.cfrac[jt] / (cfraccg[it] * cfraccg[jt])
                            else:
                                cprodrat = f32(0.0)
                            weights.cwtpot[in_, ic] = float(cprodrat) * x_rayf0.weightx[iqinv, ic, ix]
                    else:
                        for ic in range(1, npairs + 1):
                            it = molecom.n1extr[ic]
                            jt = molecom.n2extr[ic]
                            if cfraccg[it] > 0.0 and cfraccg[jt] > 0.0:
                                # Correct for coarse graining if present
                                cprodrat = molecom.cfrac[it] * molecom.cfrac[jt] / (cfraccg[it] * cfraccg[jt])
                            else:
                                cprodrat = f32(0.0)
                            weights.cwtpot[in_, ic] = float(cprodrat) * weights.cwt[in_, ic]

                invert_cwt(weights.cwtpot, weights.cwtt, weights.carray, weights.carrayold,  # noqa: F821
                           potential_par.feedback, potential_par.efacm, nfile, npairs, mfile, mpairs)

                # Save the results into a big array!

                for in_ in range(1, nnew + 1):
                    for ic in range(1, npairs + 1):
                        weights.cwttq[weights.nqinv, ic, in_] = weights.cwtt[ic, in_]

    # For the potential inversion coefficients, halve the unlike terms to avoid doubling the
    # magnitude of the Empirical Potential.

    for in_ in range(1, nfile + 1):
        for ic in range(1, npairs + 1):
            weights.cwtpot[ic, in_] = weights.cwttq[1, ic, in_]
            if molecom.n1extr[ic] != molecom.n2extr[ic]:
                weights.cwtpot[ic, in_] = 0.5 * weights.cwtpot[ic, in_]
    return


def _read_a1a2(fmtstr, a1, a2, w2):
    """``read(10,fmtstr,iostat=ierr) a1,a2`` (a1, a2 character*3) for get_wts.

    ``w2`` is the width of the second A edit descriptor: with ``a4`` into a character*3
    variable Fortran keeps the rightmost 3 characters of the 4-character field.
    Returns (a1, a2, ierr); a1 and a2 keep their old values when the read fails.
    """
    ierr = 0
    try:
        v = units.read_fmt(10, fmtstr)
        a1 = fstr(v[0], 3)
        if w2 > 3:
            a2 = fstr(pad(v[1], w2)[w2 - 3:], 3)
        else:
            a2 = fstr(v[1], 3)
    except FortranEOF:
        ierr = -1
    except FortranIOError as e:
        ierr = e.iostat
    return a1, a2, ierr


def _read_wread(wread):
    """``read(10,*,iostat=ierr) wread(1),wread(2)`` for get_wts.  Returns ierr."""
    ierr = 0
    try:
        r = units.list_read(10)
        wread[1] = r.real(wread[1])
        wread[2] = r.real(wread[2])
    except FortranEOF:
        ierr = -1
    except FortranIOError as e:
        ierr = e.iostat
    return ierr


def _read_xray_coeffs(cxt, axt, bxt, alphaxt, deltaxt):
    """``read(10,*,iostat=ierr) cxt,(axt(i),bxt(i),i=1,5),alphaxt,deltaxt`` for get_wts.

    Items read before an error keep their new values, as in Fortran.
    Returns (cxt, alphaxt, deltaxt, ierr); axt, bxt are filled in place.
    """
    ierr = 0
    try:
        r = units.list_read(10)
        cxt = r.dble(cxt)
        for i in range(1, 6):
            axt[i] = r.dble(axt[i])
            bxt[i] = r.dble(bxt[i])
        alphaxt = r.real(alphaxt)
        deltaxt = r.real(deltaxt)
    except FortranEOF:
        ierr = -1
    except FortranIOError as e:
        ierr = e.iostat
    return cxt, alphaxt, deltaxt, ierr


@routine
def wtype(a1, a2, wread, in_):
    """epsr_standalone_rev3.f:wtype -- store the neutron weights wread of pair (a1,a2)
    for dataset in."""
    n1 = 0
    n2 = 0
    for i in range(1, molecom.ntype + 1):
        if a1 == molecom.type[i]:
            n1 = i
        if a2 == molecom.type[i]:
            n2 = i
    if n1 > 0 and n2 > 0:
        icref = molecom.npairref[n1, n2]
        for j in range(1, 3):
            potential_par.weight[j, in_, icref] = wread[j]
    else:
        units.write(6, fmt(_F_WTYPE_100, in_))
    return


@routine
def wxinit(nq):
    """epsr_standalone_rev3.f:wxinit -- initialise the X-ray form factors for all
    components to zero."""
    for ix in range(1, x_rayf0.mxfile + 1):
        for it in range(1, molecom.ntype + 1):
            for iq in range(1, nq + 1):
                x_rayf0.f0x[iq, it, ix] = 0.0
    return


@routine
def wxtype(a1, axt, bxt, cxt, alphaxt, deltaxt, in_, q, nq):
    """epsr_standalone_rev3.f:wxtype -- set up the X-ray form factor for the specified
    atom type."""
    n1 = 0
    one = 1.0
    pi4 = 16.0 * math.atan(one)
    for i in range(1, molecom.ntype + 1):
        if a1 == molecom.type[i]:
            n1 = i
    if n1 > 0:
        for iq in range(1, nq + 1):
            qx = float(q[iq]) / pi4
            qx2 = qx * qx
            sum = cxt
            for i in range(1, 6):
                term = -bxt[i] * qx2
                if term > -40.0 and term < +40.0:
                    sum = sum + axt[i] * math.exp(term)
            x_rayf0.f0x[iq, n1, in_] = sum

            # Multiply by the modification function according to
            # equation (5) of JCP, vol 113, p 9149 (2000)

            if alphaxt != 0.0 and deltaxt > 0.0:

                expon = q[iq] * q[iq] / (r4(2.0) * deltaxt * deltaxt)
                if expon < 20.0:
                    fac = f32(1.0 - (float(alphaxt) / x_rayf0.f0x[1, n1, in_]) * float(expf(-expon)))
                else:
                    fac = r4(1.0)
            else:
                fac = r4(1.0)
            x_rayf0.fm0x[iq, n1, in_] = float(fac) * sum
    else:
        units.write(6, fmt(_F_WXTYPE_100, in_))
    return


@routine
def wxcreate(q, nq):
    """epsr_standalone_rev3.f:wxcreate -- Q-dependent X-ray weights for each atom pair,
    the self-scattering background, and the files X-raywts.dat, X-rayf0x.dat."""
    # Get the ratio of atoms to molecules in case it is needed

    amolrat = f32(molecom.natomr) / f32(clustercom.nmol)
    cfraccg = molecom.cfraccg
    f0x = x_rayf0.f0x
    fm0x = x_rayf0.fm0x
    weightx = x_rayf0.weightx

    for ix in range(1, weights.nxfile + 1):

        # Determine the type of normalisation

        # nrho 1 means data are per atom, 2 data are per molecule

        nrho = int(x_rayf0.datanorm[ix])

        # iafformaff is 1 means use IAFFs for normalisation, otherwise use MAFFs

        iafformaff = int(x_rayf0.iafformaffx[ix])

        # normx = 0 means no normalisation, normx = 1 means normalise to <F>^2,
        # normx = 2 means normalise to <F^2>

        normx = int(x_rayf0.normxx[ix])
        units.write(6, lw("wxcreate> ", x_rayf0.datanorm[ix], nrho, normx, iafformaff))

        # Form the Q-dependent weights for each atom pair

        for iq in range(1, nq + 1):

            # Form the single atom term

            singlem = f32(0.0)
            singlei = f32(0.0)
            bbar = 0.0
            for it in range(1, molecom.ntype + 1):

                # Set up the single atom scattering for the modified and independent atom form factors

                if molecom.qtypea[it] == 0:
                    cf = float(cfraccg[it])
                    singlem = f32(float(singlem) + cf * (fm0x[iq, it, ix] * fm0x[iq, it, ix]))
                    singlei = f32(float(singlei) + cf * (f0x[iq, it, ix] * f0x[iq, it, ix]))
                    if iafformaff != 2:
                        if normx == 0:
                            bbar = bbar + cf
                        elif normx == 1:
                            bbar = bbar + cf * f0x[iq, it, ix]
                        elif normx == 2:
                            bbar = bbar + cf * (f0x[iq, it, ix] * f0x[iq, it, ix])
                    else:
                        if normx == 0:
                            bbar = bbar + cf
                        elif normx == 1:
                            bbar = bbar + cf * fm0x[iq, it, ix]
                        elif normx == 2:
                            bbar = bbar + cf * (fm0x[iq, it, ix] * fm0x[iq, it, ix])

            if normx == 1:
                bbar = bbar * bbar

            for ip in range(1, molecom.npairs + 1):

                n1 = molecom.n1extr[ip]
                n2 = molecom.n2extr[ip]
                weightx[iq, ip, ix] = fm0x[iq, n1, ix] * fm0x[iq, n2, ix] / bbar

                # Correct for per molecule data if needed

                if nrho == 2 and normx == 1:
                    weightx[iq, ip, ix] = weightx[iq, ip, ix] / float(amolrat)

            # Calculate background shift to simulated structure factor

            if iafformaff == 1:
                x_rayf0.backx[iq, ix] = float(singlem - singlei) / bbar
            else:
                x_rayf0.backx[iq, ix] = 0.0
            if nrho == 2 and normx == 1:
                x_rayf0.backx[iq, ix] = x_rayf0.backx[iq, ix] / float(amolrat)

        # Set the Q=0 weights

        for ip in range(1, molecom.npairs + 1):
            potential_par.weight[1, x_rayf0.datapoint[ix], ip] = weightx[1, ip, ix]
            potential_par.weight[2, x_rayf0.datapoint[ix], ip] = potential_par.weight[1, x_rayf0.datapoint[ix], ip]

        # Write an output file to check that all is well

        units.open(10, "X-raywts.dat", status="unknown")
        for iq in range(1, nq + 1):
            units.write(10, fmt(_F_WXCREATE_100, q[iq], *[weightx[iq, ip, ix] for ip in range(1, molecom.npairs + 1)]))
        units.close(10)

        # Write an output file to show the X-ray form factors used

        units.open(10, "X-rayf0x.dat", status="unknown")
        for iq in range(1, nq + 1):
            units.write(10, fmt(_F_WXCREATE_100, q[iq], *[fm0x[iq, it, ix] for it in range(1, molecom.ntype + 1)]))
        units.close(10)

    return


# =============================================================================
# from epsr_standalone_rev3_part3.py
# =============================================================================

# exp(25.0) is folded by gfortran (correctly rounded); equal to expf(25.0) (checked).
_EXP25 = expf(r4(25.0))

# 200 format(i5/(6(1x,e12.5)))  (invert_cwt)
_F200_INVERT_CWT = "(i5/(6(1x,e12.5)))"
# 100 format(7(1x,e13.6))  (lsparabola, lsparabola1)
_F100_LSPARABOLA = "(7(1x,e13.6))"


@routine
def gensq(nq, q, sofq, wav):
    """epsr_standalone_rev3.f:gensq -- partial structure factors of the current box.

    sofq(1,i,ic) is the intramolecular, sofq(2,i,ic) the intermolecular term, from the
    histograms wav through the sine-transform table foqols.
    """
    # real q(mpts),sofq(2,mpts,mpairs),wav(2,mpz,mpairs),sumwavself(mpairs)
    sumwavself = real_array(_mpairs)
    #      open(25,file='dump-25.txt',status='unknown')
    #      write(6,*) 'gensq> ',maxnscalefactype
    if molecom_mod.coarsegraining:
        molecom_mod.prodssff = real_array(molecom_mod.maxnscalefactype, molecom_mod.maxnscalefactype)
    # Set the fudge factors on the self terms
    ycohfactor = f32(1.0)  # noqa: F841
    sumwavself[:] = 0.0
    npt1 = epsr16.npt1
    foqols = epsr16.foqols
    # (speed only) element lists of the arrays read in the loops; values and summation
    # order are unchanged.
    ycoh = list(epsr16.ycoh)
    for ic in range(1, molecom.npairs + 1):
        it = molecom.n1extr[ic]
        jt = molecom.n2extr[ic]
        for j in range(2, npt1 + 1):
            sumwavself[ic] = sumwavself[ic] + wav[1, j, ic]
        wav1 = list(wav[1, :, ic])
        wav2 = list(wav[2, :, ic])
        for i in range(1, nq + 1):
            sum1 = f32(0.0)
            sum2 = f32(0.0)
            sum3 = f32(0.0)
            if molecom_mod.coarsegraining:
                # Calculate the product of form factors at this Q value
                prodssff = molecom_mod.prodssff
                sphereformfactor = molecom_mod.sphereformfactor
                for ibinia in range(1, molecom_mod.nsffbintype[it] + 1):
                    for ibinja in range(1, molecom_mod.nsffbintype[jt] + 1):
                        prodssff[ibinia, ibinja] = sphereformfactor[i, ibinia, it] * sphereformfactor[i, ibinja, jt]
            fq = list(foqols[:, i])
            for j in range(2, npt1 + 1):
                # First do intra-molecular terms
                add1 = wav1[j]
                if abs(add1) > 0.0:
                    sum1 = sum1 + add1 * fq[j]
            #            if(coarsegraining) then
            #               sum2=intraformfactor(i,ic) !This term is added to sum2 so that it gets the correct neutron weighting
            #               if(it.eq.jt.and.cfraccg(it).gt.0.0) then
            #                  sum1=intraformfactorself(i,it)/cfraccg(it) !Include the self term in this case
            #               end if
            #            end if
            # Now inter-molecular terms
            natomsphere_jt = molecom_mod.natomsphere[jt]
            for j in range(2, npt1 + 1):
                # (the coarse-grained wavinter/prodssff version is commented out in the Fortran)
                add2 = (wav2[j] - ycoh[j]) * natomsphere_jt
                sum2 = sum2 + add2 * fq[j]
                sum3 = sum3 + wav2[j]
            #
            # in the event that there are no intermolecular correlations for a particular atom pair
            # (e.g. only one molecule of that type in the box)
            # we set the simulated structure factor to zero
            #
            if sum3 <= 0.0:
                sum2 = f32(0.0)
            sofq[1, i, ic] = sum1
            sofq[2, i, ic] = sum2
    #      write(6,'((1x,i3,1x,a3,1x,a3,1x,f10.5))') (ic,type(n1extr(ic)),type(n2extr(ic)),sumwavself(ic),ic=1,npairs)
    #      close(25)
    if molecom_mod.coarsegraining:
        if molecom_mod.prodssff is not None:
            molecom_mod.prodssff = None
    return


@routine
def fit_expansion(yin, expr, yout, wt, wtold, stepr, acceptr, rejrat, chisqnew, nmin, npt, ncfirst, nterms, nits, nfreq):
    """epsr_standalone_rev3.f:fit_expansion -- Monte Carlo fit of the coefficients wt so
    that sum_i wt(i)*expr((j-1)*nterms+i) fits yin(j), j = nmin..npt.

    Consumes ran1.  Returns (stepr, acceptr, rejrat, chisqnew).
    """
    from ..library import auxiliary_routines

    # Internal variables
    # integer :: nstep,itry,ntry,its,npoints
    # real :: wtmax,chisqold,rejfac,chisqstart,chisqfinish,chisqtest
    #
    # maximum value of a sum - this it to prevent overflows
    #
    summax = _EXP25
    nstep = 100
    ntry = 150
    acceptr = r4(0.05)
    npoints = npt - nmin + 1
    #
    # assign the initial coefficients
    #
    wtmax = f32(0.0)
    for it in range(1, nterms + 1):
        wt[it] = wtold[it]
        wtmax = max(abs(wt[it]), wtmax)
    if stepr < r4(1.0e-6) * wtmax:
        stepr = r4(2.0e-6) * wtmax  # Prevent step size getting too small!
    if nits > 0:
        chisqold = f32(0.0)
        icount = 0
        chisqnew = f32(0.0)
        for i in range(nmin, npt + 1):
            term = yin[i] - expand(nterms, i, expr, wt, summax)
            chisqnew = chisqnew + term * term
        chisqold = chisqnew
        chisqtest = r4(-1.0)
        its = 0
        itry = 0
        while itry < 10 or (itry < ntry and chisqtest < r4(1.0e-4)):
            #         do while (itry.lt.ntry.and.abs(chisqtest).gt.1.0e-3.and.chisqtest.lt.0)
            itry = itry + 1
            accept = acceptr * chisqold / f32(npoints)
            rejrat = f32(0.0)
            chisqstart = chisqold / f32(npoints)
            for it in range(1, nstep + 1):
                its = its + 1  # noqa: F841
                #
                # select some new trial values
                #
                for ic in range(ncfirst, nterms + 1):
                    if abs(wtold[ic]) < 1.0:
                        wt[ic] = wtold[ic] + stepr * (2.0 * auxiliary_routines.ran1() - 1.0)
                    else:
                        wt[ic] = wtold[ic] * (1.0 + stepr * (2.0 * auxiliary_routines.ran1() - 1.0))
                chisqnew = f32(0.0)
                for i in range(nmin, npt + 1):
                    term = yin[i] - expand(nterms, i, expr, wt, summax)
                    chisqnew = chisqnew + term * term
                exponent = (chisqnew - chisqold) / accept
                if exponent > 0.0 and exponent <= 20.0:
                    test = expf(-exponent) - auxiliary_routines.ran1()
                elif exponent > 20.0:
                    test = r4(-1.0)
                else:
                    test = r4(1.0)
                if test >= 0.0:
                    #
                    # perform 1 smoothing on the coefficients to prevent them becoming noisy
                    # every nfreq times a set of moves are accepted
                    #
                    icount = icount + 1
                    if icount == nfreq:
                        icount = 0
                        auxiliary_routines.smoop(ncfirst, nterms, 1, wt)
                        #
                        # regenerate chisqnew for the modified coeffs.
                        #
                        chisqnew = f32(0.0)
                        for i in range(nmin, npt + 1):
                            term = yin[i] - expand(nterms, i, expr, wt, summax)
                            chisqnew = chisqnew + term * term
                    #
                    # accept new values
                    #
                    chisqold = chisqnew
                    for ic in range(1, nterms + 1):
                        wtold[ic] = wt[ic]
                else:
                    rejrat = rejrat + r4(1.0)
            #            chisqtest=0.0 ... (commented-out alternative test in the Fortran)
            rejrat = rejrat / f32(nstep)
            #          if(chisqnew.gt.0.0) chisqtest=chisqtest/chisqnew
            chisqfinish = chisqold / f32(npoints)
            if chisqstart > 0.0:
                chisqtest = (chisqfinish - chisqstart) / chisqstart
            else:
                chisqtest = chisqfinish - chisqstart
            if rejrat > 0.0:
                rejfac = r4(0.75) / rejrat
            else:
                rejfac = r4(2.0)
            stepr = stepr * rejfac
            #            write(6,*) 'fit_expansion> ',itry,chisqfinish,chisqtest,rejrat
        units.write(6, lw("fit_expansion> ", itry, chisqfinish, chisqtest, rejrat))
        for ic in range(1, nterms + 1):
            wt[ic] = wtold[ic]
    chisqnew = f32(0.0)
    for i in range(1, npt + 1):
        term = expand(nterms, i, expr, wt, summax)
        yout[i] = term
        if i >= nmin:
            d = yin[i] - term
            chisqnew = chisqnew + d * d
    return stepr, acceptr, rejrat, chisqnew


@routine
def expand(nterms, j, expr, coeff, summax):
    """epsr_standalone_rev3.f:expand -- sum_i expr((j-1)*nterms+i)*coeff(i), clamped at
    +-summax after every term.

    Returns expand.
    """
    # real*4 coeff(*),expr(*)
    sum = f32(0.0)
    jref0 = (j - 1) * nterms
    # (speed only) iterate over the row and the coefficients in the Fortran order
    for e, c in zip(expr[jref0 + 1:jref0 + nterms + 1], coeff[1:nterms + 1]):
        sum = sum + e * c
        if abs(sum) > summax:
            sum = sign(summax, sum)
    expand_ = sum
    return expand_


@routine
def calc_gaussq(lpt, nterms, dstep, dmax, sigma, rho, xin, expr, rmin, ncfirst):
    """epsr_standalone_rev3.f:calc_gaussq -- Q-space transforms of the Gaussian EP basis.

    Returns ncfirst.
    """
    # real*4 xin(*),expr(*)
    pi = r4(3.141592653)  # noqa: F841
    #
    # maximum value of exponent
    #
    expmax = r4(25.0)
    #
    # set up array of coefficients
    #
    for m in range(1, nterms + 1):
        #
        # Set the value of d
        #
        d = f32(m) * dstep
        #
        # save the integer of the first non-zero coefficient
        #
        if d < rmin:
            ncfirst = m + 1
        #
        # generate value of sigma at this d value
        #
        sigd = sigma
        #
        # calculate Q dependence of mth term
        #
        for i in range(1, lpt + 1):
            jref = (i - 1) * nterms + m
            q = xin[i]
            qsigd = q * sigd
            qs2 = 0.5 * qsigd * qsigd
            if qs2 < expmax:
                qd = q * d
                if qd > 0.0:
                    sinqd = sinf(qd) / qd
                else:
                    sinqd = r4(1.0)
                expr[jref] = sinqd * expf(-qs2)
            else:
                expr[jref] = 0.0
    return ncfirst


@routine
def calc_gaussr(lpt, nterms, dstep, dmax, sigmab, rho, xin, expr, rmin, ncfirst, fbroad, rbroad):
    """epsr_standalone_rev3.f:calc_gaussr -- r-space Gaussian EP basis functions.

    Returns ncfirst.
    """
    # real*4 xin(*),expr(*)
    pi = PI4ATAN
    piroot = sqrtf(0.5 * pi)
    gfac = piroot / (r4(4.0) * pi * pi * rho)
    #
    # maximum value of exponent
    #
    expmax = r4(25.0)
    #
    # set up array of coefficients
    #
    for m in range(1, nterms + 1):
        #
        # calculate d dependence of mth term
        #
        d = f32(m) * dstep
        #
        # save the integer of the first non-zero coefficient
        #
        if d < rmin:
            ncfirst = m + 1
        #
        # generate value of sigma at this d value
        #
        sigd = fbroad * sigmab
        gfacs = gfac / sigd
        for i in range(1, lpt + 1):
            jref = (i - 1) * nterms + m
            r = xin[i]
            rmd = (r - d) / sigd
            rpd = (r + d) / sigd
            rmd2 = 0.5 * rmd * rmd
            rpd2 = 0.5 * rpd * rpd
            sum = f32(0.0)
            if r > 0.0 and d > 0.0:
                gfact = gfacs / (r * d)
                if rmd2 < expmax:
                    sum = sum + expf(-rmd2)
                if rpd2 < expmax:
                    sum = sum - expf(-rpd2)
            else:
                gfact = 2.0 * gfacs / (sigd * sigd)
                if rpd2 < expmax:
                    sum = sum + expf(-rpd2)
            expr[jref] = gfact * sum
    return ncfirst


@routine
def calc_gaussr_grad(lpt, nterms, dstep, dmax, sigmab, rho, xin, expr, rmin, ncfirst, fbroad, rbroad):
    """epsr_standalone_rev3.f:calc_gaussr_grad -- r-space gradients of the Gaussian EP
    basis functions.

    Returns ncfirst.
    """
    # real*4 xin(*),expr(*)
    pi = PI4ATAN
    piroot = sqrtf(0.5 * pi)
    gfac = piroot / (r4(4.0) * pi * pi * rho)
    #
    # maximum value of exponent
    #
    expmax = r4(25.0)
    #
    # set up array of coefficients
    #
    for m in range(1, nterms + 1):
        #
        # calculate d dependence of mth term
        #
        d = f32(m) * dstep
        #
        # save the integer of the first non-zero coefficient
        #
        if d < rmin:
            ncfirst = m + 1
        #
        # generate value of sigma at this d value
        #
        sigd = fbroad * sigmab
        gfacs = gfac / sigd
        sigd2 = sigd * sigd
        for i in range(1, lpt + 1):
            jref = (i - 1) * nterms + m
            r = xin[i]
            rmd = (r - d) / sigd
            rpd = (r + d) / sigd
            rmd2 = 0.5 * rmd * rmd
            rpd2 = 0.5 * rpd * rpd
            sum = f32(0.0)
            if abs(r) > 0.0 and abs(d) > 0.0:
                gfact = gfacs / (r * d)
                if rmd2 < expmax:
                    sum = sum - (1.0 / r + (r - d) / sigd2) * expf(-rmd2)
                if rpd2 < expmax:
                    sum = sum + (1.0 / r + (r + d) / sigd2) * expf(-rpd2)
            else:
                gfact = r4(0.0)
            expr[jref] = gfact * sum
    return ncfirst


@routine
def calc_pwrexpq(lpt, nterms, dmin, dstep, dmax, sigmab, rho, xin, expr, rmin, ncfirst):
    """epsr_standalone_rev3.f:calc_pwrexpq -- Q-space transforms of the power-exponential
    EP basis functions.

    Returns ncfirst.
    """
    # real xin(*),expr(*),d,dmin,dstep,dmax,dnext,rootd,sigmab,sigma
    pi = r4(3.141592653)  # noqa: F841
    #
    # maximum value of exponent
    #
    expmax = r4(25.0)  # noqa: F841
    #
    # set up array of coefficients
    #
    d = dstep
    rootd = r4(1.0)
    for m in range(1, nterms + 1):
        #
        # calculate Q dependence of mth term
        #
        #            d=m*dstep
        #
        # save the integer of the first non-zero coefficient
        #
        if d < rmin:
            ncfirst = m + 1
        sigma = rootd * sigmab
        #            sigma=sigmab
        if sigma <= d and sigma > 0.0 and m > 1:
            n = ifix(d / sigma - 1.0)
        else:
            n = 0
        an2 = f32(n + 2)
        an2fac = 1.0 / an2
        for i in range(1, lpt + 1):
            q = xin[i]
            qsigma = q * sigma
            qs2 = qsigma * qsigma
            onemqs2 = 1.0 - qs2
            onepqs2 = 1.0 + qs2
            onepqs2s = onepqs2 * onepqs2  # noqa: F841
            rtonepqs2 = sqrtf(onepqs2)
            alpha = atanf(qsigma)
            alphan = f32(n) * alpha
            if qsigma > 0.0:
                term = 2.0 * cosf(alphan) + onemqs2 * sinf(alphan) / qsigma
            else:
                term = an2
            pn = term * an2fac / powi(rtonepqs2, 4 + n)
            jref = (i - 1) * nterms + m
            expr[jref] = pn
        d, rootd = dnext(d, dmin, dstep, rootd)
    return ncfirst


@routine
def calc_pwrexpr(lpt, nterms, dmin, dstep, dmax, sigmab, rho, xin, expr, rmin, ncfirst, fbroad, rbroad):
    """epsr_standalone_rev3.f:calc_pwrexpr -- r-space power-exponential EP basis functions.

    Returns ncfirst.
    """
    # use lngamma_mod
    # real*4 xin(*),expr(*),d,dmin,dstep,dmax,dnext,rootd,sigmab,sigma
    pi = r4(3.141592653)
    #
    # maximum value of exponent
    #
    expmax = r4(25.0)
    ntmax = ifix(dmax / sigmab + r4(3.0))
    lngamma_mod.fctrln(ntmax)
    gamln = lngamma_mod_data.gamln
    #
    # set up array of coefficients
    #
    d = dstep
    rootd = r4(1.0)
    for m in range(1, nterms + 1):
        #
        # calculate r dependence of mth term
        #
        #            d=m*dstep
        #
        # save the integer of the first non-zero coefficient
        #
        if d < rmin:
            ncfirst = m + 1
        #
        # generate value of sigma at this d value
        #
        sigma = rootd * fbroad * sigmab
        if sigma <= d and sigma > 0.0:
            n = ifix(d / sigma - 1.0)
        else:
            n = 0
        nref = n + 3
        gfac = r4(4.0) * pi * rho * sigma * sigma * sigma
        gfacln = logf(gfac)
        sigln = f32(gamln[nref] + float(gfacln))
        #
        # for the first term we add a factor exp(rmin/sigma)
        # This is to ensure low  r repulsive potential develops quickly
        #
        if m == 1:
            sigln = sigln - rmin / sigma
        for i in range(1, lpt + 1):
            jref = (i - 1) * nterms + m
            r = xin[i]
            fact = rbroad * r
            rsigma = r / sigma
            if r > 0.0:
                rsigmaln = logf(rsigma)
                expon = f32(n) * rsigmaln - rsigma - sigln - fact
            else:
                if n == 0:
                    expon = -rsigma - sigln - fact
                else:
                    #
                    # for n>0, the function is identically zero at r=0
                    #
                    expon = -2.0 * expmax
            if expon > -expmax:
                expr[jref] = expf(expon)
            else:
                expr[jref] = 0.0
        d, rootd = dnext(d, dmin, dstep, rootd)
    return ncfirst


@routine
def calc_pwrexpr_grad(lpt, nterms, dmin, dstep, dmax, sigmab, rho, xin, expr, rmin, ncfirst, fbroad, rbroad):
    """epsr_standalone_rev3.f:calc_pwrexpr_grad -- r-space gradients of the
    power-exponential EP basis functions.

    Returns ncfirst.
    """
    # use lngamma_mod
    # real xin(*),expr(*),dmin,dstep,dmax,dnext,rootd,sigmab,sigma
    pi = r4(3.141592653)
    #
    # maximum value of exponent
    #
    expmax = r4(25.0)
    ntmax = ifix(dmax / sigmab + r4(3.0))
    lngamma_mod.fctrln(ntmax)
    gamln = lngamma_mod_data.gamln
    #
    # set up array of coefficients
    #
    d = dstep
    rootd = r4(1.0)
    for m in range(1, nterms + 1):
        #
        # calculate r dependence of mth term
        #
        #            d=m*dstep
        #            write(6,*) "calc_pwrexp_grad> ", m,d,dmin,dstep,rootd
        #
        # save the integer of the first non-zero coefficient
        #
        if d < rmin:
            ncfirst = m + 1
        #
        # generate value of sigma at this d value
        #
        sigma = rootd * fbroad * sigmab
        if sigma <= d and sigma > 0.0 and m > 1:
            n = ifix(d / sigma - 1.0)
        else:
            n = 0
        nref = n + 3
        gfac = r4(4.0) * pi * rho * sigma * sigma * sigma
        gfacln = logf(gfac)
        sigln = f32(gamln[nref] + float(gfacln))
        #
        # for the first term we add a factor exp(rmin/sigma)
        # This is to ensure low  r repulsive potential develops quickly
        #
        if m == 1:
            sigln = sigln - rmin / sigma
        for i in range(1, lpt + 1):
            jref = (i - 1) * nterms + m
            r = xin[i]
            fact = rbroad * r
            rsigma = r / sigma
            if r > 0.0:
                rsigmaln = logf(rsigma)
                expon = f32(n) * rsigmaln - rsigma - sigln - fact
                #
                # add on the term to give the gradient of the potential
                #
                fact = f32(n) / r - 1.0 / sigma - rbroad
            else:
                if n == 0:
                    expon = -rsigma - sigln - fact
                    fact = -(1.0 / sigma + rbroad)
                elif n == 1:
                    expon = -rsigma - sigln - fact
                    fact = 1.0 / sigma
                else:
                    #
                    # for n>1, the gradient is identically zero at r=0
                    #
                    expon = -2.0 * expmax
                    fact = r4(0.0)
            if expon > -expmax:
                expr[jref] = fact * expf(expon)
            else:
                expr[jref] = 0.0
        d, rootd = dnext(d, dmin, dstep, rootd)
    return ncfirst


@routine
def dnext(d, dfac, dstep, rootd):
    """epsr_standalone_rev3.f:dnext -- next basis-function position; the spacing grows
    geometrically (factor dfac) once that exceeds dstep.

    Returns (dnext, rootd).
    """
    # real dnext,d,dstep,dfac,dtest,rootd
    dnext_ = d + dstep
    dtest = d * dfac
    if dtest > dnext_:
        dnext_ = dtest
        rootd = rootd * sqrtf(dfac)
    return dnext_, rootd


@routine
def ftback(in_, rmin, q, sofq, grt):
    """epsr_standalone_rev3.f:ftback -- the part of sofq that comes from r < rmin.

    perform a direct transform on sofq for r values, and use rmin to
    generate a background to sofq. Result is in sofq.
    (The Fortran dummy ``in`` is ``in_`` here: ``in`` is a Python keyword.)
    """
    # dimension q(*),sofq(*),grt(*)
    nqmin = epsr16.nqmin[in_]
    nqmax = epsr16.nqmax[in_]
    x = epsr16.x
    rho = epsr16.rho
    for i in range(nqmin, nqmax + 1):
        sofq[i] = sofq[i] * q[i]
    pi = r4(3.141592653)
    qstep = q[4] - q[3]
    sfac = 0.5 * qstep / pi / pi / rho
    jsave = 2
    for j in range(2, epsr16.npt1 + 1):
        rr = x[j]
        if rr >= rmin:
            break  # go to 200
        jsave = j
        sum = f32(0.0)
        for i in range(nqmin, nqmax + 1):
            qq = q[i]
            qr = qq * rr
            sum = sum + sofq[i] * sinf(qr)
        grt[j] = sfac * sum
    # 200 continue
    #
    # calculate the S(Q) from the truncated g(r)
    #
    gfac = r4(4.0) * pi * rho
    for i in range(nqmin, nqmax + 1):
        qq = q[i]
        if qq > 0.0:
            r1 = f32(0.0)
            sum = f32(0.0)
            for j in range(2, jsave + 1):
                rr = x[j]
                r2 = 2.0 * rr - r1
                rstep = r2 - r1
                qr = qq * rr
                sum = sum + grt[j] * sinf(qr) * rstep
                r1 = r2
            sofq[i] = gfac * sum / qq
        else:
            r1 = f32(0.0)
            sum = f32(0.0)
            for j in range(2, jsave + 1):
                rr = x[j]
                r2 = 2.0 * rr - r1
                rstep = r2 - r1
                sum = sum + grt[j] * rr * rstep
                r1 = r2
            sofq[i] = gfac * sum
    return


@routine
def uniform(ycoh, xst, rho, rmaxgr, ndim):
    """epsr_standalone_rev3.f:uniform -- the uniform distribution (number of atoms per
    shell for density rho) out to rmaxgr.

    calculates the uniform distribution for a unit cell of arbitrary shape
    """
    # include 'lattice.inc'   (nothing from it is used)
    # integer*8 ntries
    # dimension ycoh(ndim)

    # now calculate the region out to the edge of the box analytically

    pi = PI4ATAN
    rhofac = r4(4.0) * pi * rho / r4(3.0)
    x1 = f32(0.0)
    x13 = f32(0.0)
    xstep = xst * r4(0.5)
    x2 = x1 + xstep
    ycoh[1] = 0.0
    j = 1
    while x2 < rmaxgr and j < ndim:
        j = j + 1
        x2 = x1 + xstep
        x23 = x2 * x2 * x2
        ycoh[j] = rhofac * (x23 - x13)
        xstep = xst
        x1 = x2
        x13 = x23
    return


@routine
def invert_cwt(cwt, cwtt, carray, carrayold, feedback, efacm, nfile, npairs, mfile, mpairs):
    """epsr_standalone_rev3.f:invert_cwt -- Monte Carlo "smart transpose" cwtt of the
    weights matrix cwt (consumes ran1).

    generate the smart tranpose of cwt to take account of the different weights
    of particular partials in different datasets

    (speed only) The real*8 arrays cwt, cwtt, carray and carrayold are worked on as
    nested Python lists (``cwt_[i][j]`` is cwt(i,j); Python floats are real*8, so the
    arithmetic is unchanged) and copied back into the caller's arrays at the end.
    """
    from ..library import auxiliary_routines

    ran1 = auxiliary_routines.ran1  # bound at entry (tracing state at the call is kept)
    # real*8 cwt(mpairs,mpairs),cwtt(mpairs,mpairs+mfile)
    # 1,carray(mpairs+mfile),carrayold(mpairs+mfile,mpairs+mfile)
    # real*8 temp,sum,sumold,total,totalold,tdif,add
    # real*4 efacm(mfile)
    cwt_ = cwt.tolist()
    cwtt_ = cwtt.tolist()
    carray_ = carray.tolist()
    carrayold_ = carrayold.tolist()
    #
    # the chosen method is to select the values of the transpose matrix by Monte Carlo,
    # calculate the cross-correlation matrix, Carray, and aim towards diagonal elements
    # of this array to be unity, while off diagonal units are zero. To develop a single
    # solution when the problem is overdetermined, the sum of the squares of the tranpose
    # coefficients is kept as small as possible.
    #
    niter = 50000
    ntimes = 20  # noqa: F841
    #
    # FEEDBACK is the feedback factor used in the simulation. FEEDBACK must not be > 1.0
    #
    feedback1 = r4(1.0) - feedback
    feedback1d = float(feedback1)  # feedback1 as it enters the real*8 expressions
    #
    # we will assume there is, in an addition to the data, a distinct
    # set of partials available (in this case from the simulation itself)
    # to make up for the lack of distinct information in the data
    #
    # this may affect the values of the inversion matrix, but is needed to
    # force a unique solution. However since we use the inversion matrix only
    # on the DIFFERENCE between data a fit, the inversion matrix for these
    # pseudo datasets will have no effect on the values of the empirical potential
    # since the difference for the pseudo data is by definition zero.
    #
    nnew = npairs + nfile
    temp = float(r4(1.0) / f32(nnew))
    #
    # write out input coefficients and form average modulo(weight)
    #
    iadd = 0  # noqa: F841
    units.write(6, lw())
    for i in range(1, nfile + 1):
        for j in range(1, npairs + 1):
            cwt_[i][j] = cwt_[i][j] * float(feedback) * float(efacm[i])
        units.write(6, fmt(_F200_INVERT_CWT, i, *[cwt_[i][j] for j in range(1, npairs + 1)]))
    nfile1 = nfile + 1  # noqa: F841
    #
    # keep refining until a sensible solution is found
    #
    totaltest = r4(1.0e-15)
    #
    # sum the cross array
    #
    total = 0.0
    for i1 in range(1, nnew + 1):
        for i2 in range(1, nnew + 1):
            if i1 <= nfile:
                sum = 0.0
                for j in range(1, npairs + 1):
                    sum = sum + cwt_[i1][j] * cwtt_[j][i2]
            else:
                j = i1 - nfile
                sum = feedback1d * cwtt_[j][i2]
            carrayold_[i1][i2] = sum
            if i1 == i2:
                t = sum - 1.0
                total = total + t * t
            else:
                total = total + sum * sum
    totalold = total
    #
    # set starting step size
    #
    amp = f32(float(r4(0.01)) * total)
    #
    # new fitting temperature
    #
    temp = (totalold - float(nfile)) / float(nnew * nnew)
    itimes = 1
    iter = 0
    nreject = 0
    while (totalold - float(nfile)) >= float(totaltest) and itimes <= 1000:
        #        write(6,*) itimes,iter,nreject,totalold,temp
        for i2 in range(1, nnew + 1):
            #
            # generate new values
            #
            for j in range(1, npairs + 1):
                iter = iter + 1
                # add=(2.0*ran1()-1)*amp: real*4 expression, stored in real*8
                add = float((2.0 * ran1() - 1.0) * amp)
                #
                # calculate the cross array for this set of values
                #
                for i1 in range(1, nfile + 1):
                    carray_[i1] = carrayold_[i1][i2] + cwt_[i1][j] * add
                #
                # include the extra term for the diagonal cwts
                #
                i1 = j + nfile
                carray_[i1] = carrayold_[i1][i2] + feedback1d * add
                #
                # calculate the change in total for the changes in these specific elements of carray
                #
                tdif = 0.0
                for i1 in range(1, nfile + 1):
                    sum = carray_[i1]
                    sumold = carrayold_[i1][i2]
                    if i1 == i2:
                        t = sum - 1.0
                        told = sumold - 1.0
                        tdif = tdif + t * t - told * told
                    else:
                        tdif = tdif + sum * sum - sumold * sumold
                i1 = j + nfile
                sum = carray_[i1]
                sumold = carrayold_[i1][i2]
                if i1 == i2:
                    t = sum - 1.0
                    told = sumold - 1.0
                    tdif = tdif + t * t - told * told
                else:
                    tdif = tdif + sum * sum - sumold * sumold
                tdif = tdif / temp
                #
                # perform Monte Carlo test to accept or reject this move
                #
                iaccept = 0
                if tdif <= 0.0:
                    iaccept = 1
                elif tdif > 20.0:
                    iaccept = 0
                elif float(ran1()) <= math.exp(-tdif):
                    iaccept = 1
                else:
                    iaccept = 0
                if iaccept == 1:
                    #
                    # move accepted - save the new coefficients
                    #
                    total = total + tdif * temp
                    totalold = total
                    cwtt_[j][i2] = cwtt_[j][i2] + add
                    for i1 in range(1, nfile + 1):
                        carrayold_[i1][i2] = carray_[i1]
                    #
                    # include the extra term for the diagonal cwts
                    #
                    i1 = j + nfile
                    carrayold_[i1][i2] = carray_[i1]
                else:
                    nreject = nreject + 1
                    #           write(6,*) iter,nreject,totalold,tdif,temp
        #        write(6,*) itimes,iter,nreject,totalold,temp
        #
        # modify step size to give 50% rejection rate
        #
        if iter >= niter:
            if nreject > 0:
                fac = r4(0.5) * f32(iter - nreject) / f32(nreject)
            else:
                fac = r4(4.0)
            if fac < 0.95:
                fac = r4(0.95)
            if fac > 1.05:
                fac = r4(1.05)
            amp = amp * fac
            #           write(6,*)
            #           write(6,*) itimes,total,totalold,totalold-nfile
            #           write(6,*) nreject,iter,fac,amp,temp
            itimes = itimes + 1
            iter = 0
            nreject = 0
            #
            # sum the cross array
            #
            total = 0.0
            for i1 in range(1, nnew + 1):
                for i2 in range(1, nnew + 1):
                    if i1 <= nfile:
                        sum = 0.0
                        for j in range(1, npairs + 1):
                            sum = sum + cwt_[i1][j] * cwtt_[j][i2]
                    else:
                        j = i1 - nfile
                        sum = feedback1d * cwtt_[j][i2]
                    carrayold_[i1][i2] = sum
                    if i1 == i2:
                        t = sum - 1.0
                        total = total + t * t
                    else:
                        total = total + sum * sum
            totalold = total
            #
            # new fitting temperature
            #
            temp = (totalold - float(nfile)) / float(nnew * nnew)
    # 200 format(i5/(6(1x,e12.5)))
    units.write(6, lw())
    #     do j=1,npairs
    #        write(6,200) j,(cwtt(j,i),i=1,nnew)
    #     end do
    #
    # calculate the cross array to demonstrate that we have indeed separated
    # all the partials - this array should be the unit matrix
    #
    for j1 in range(1, npairs + 1):
        for j2 in range(1, npairs + 1):
            sum = 0.0
            #
            # first sum over the data files
            #
            for i in range(1, nfile + 1):
                sum = sum + cwtt_[j1][i] * cwt_[i][j2]
            #
            # add on the pseudo terms
            #
            i = j2 + nfile
            carray_[j2] = sum + cwtt_[j1][i] * feedback1d
        units.write(6, fmt(_F200_INVERT_CWT, j1, *[carray_[j2] for j2 in range(1, npairs + 1)]))
        #
        # modify the cwtt for i.le.nfile
        #
        for i in range(1, nfile + 1):
            cwtt_[j1][i] = cwtt_[j1][i] * float(feedback) * float(efacm[i])
        #
        # modify the cwtt(j,i) for i .gt. nfile
        #
        for i in range(nfile + 1, nnew + 1):
            cwtt_[j1][i] = cwtt_[j1][i] * feedback1d
    # (copy the work lists back into the caller's arrays)
    cwt[...] = cwt_
    cwtt[...] = cwtt_
    carray[...] = carray_
    carrayold[...] = carrayold_
    return


@routine
def generate_rp():
    """epsr_standalone_rev3.f:generate_rp -- reference potential (n-6 LJ + Coulomb, extra
    Gaussians / modified Morse, repulsive core) and its gradients for every pair.
    """
    from ..library import ato_manipulation_rev4

    # use threebody_mod
    # use containerforces_mod !These are attractive forces designed to prevent particular atom pairs from becoming too far apart.
    # NB: ``add`` is a COMMON variable (epsr16.inc), so it is assigned as epsr16.add.
    x = epsr16.x
    xstar = epsr16.xstar
    npt = epsr16.npt
    npt2 = epsr16.npt2
    pp = potential_par

    # real width,pos,gradlimit,loggradlimit,ddiff

    # Limit on gradient of repulsive force
    gradlimit = epsrinput_par.expecfpwr
    if gradlimit > 0.0:
        loggradlimit = logf(gradlimit)
    else:
        gradlimit = r4(0.0)
        loggradlimit = r4(0.0)

    # Generate molecular potential (based on ecore and dcore)

    ato_manipulation_rev4.generate_mp()

    power = molecom.power
    pi = r4(3.141592653)  # noqa: F841
    if power > 0.0:
        ratt = powf(r4(1.0e10), r4(1.0) / power)  # noqa: F841
    #
    # maximum value of a sum - this it to prevent overflows
    #
    summax = _EXP25  # noqa: F841
    nterms = pp.ncoeffp  # noqa: F841
    decayp = pp.psigma2  # noqa: F841

    #
    # set up L-J, Coulomb, Morse and Gaussian potentials between pairs.
    #
    # first define prefactor for the well depth
    #
    rat = r4(6.0) / power
    rat1 = sqrtf(rat)
    rat2 = r4(1.0) / rat1
    power1 = (power + f32(6)) / (power - f32(6))
    rat3 = powf(rat1, power1)
    prefac = (rat2 - rat1) * rat3
    prefac = r4(1.0) / prefac
    for ic in range(1, molecom.npairs + 1):
        irefp = pp.irefpcof[ic]
        units.write(6, lw("generate_rp> ", ic, irefp, pp.rminex[irefp], pp.rmaxex[irefp]))
        it = molecom.n1extr[ic]
        jt = molecom.n2extr[ic]
        #
        # Coulomb interactions are treated as between uniform spherical charge clouds
        # of diameter equal to qradius(it)
        #
        ra = molecom.qradius[it]
        rb = molecom.qradius[jt]
        ect = molecom.ecor[ic]
        qct = molecom.qcor[ic] * pp.coulfac[irefp]
        dc = molecom.dcor[ic]
        x1 = x[1]
        x2 = x[2]
        for j in range(2, npt + 1):
            x2 = x[j]
            #
            # modified L-J potential n-6 where n>6
            #
            ddif = dc / x2
            if ddif > 0.0 and ect > 0.0:
                fac = prefac * ect
                xx = logf(ddif)
                pwr1 = r4(0.5) * (power + r4(6.0)) * xx
                if abs(pwr1) > 20.0:
                    pwr1 = sign(r4(20.0), pwr1)
                pwr2 = r4(0.5) * (power - r4(6.0)) * xx
                if abs(pwr2) > 40.0:
                    pwr2 = sign(r4(40.0), pwr2)
                sumlj = r4(2.0) * fac * sinhf(pwr2) * expf(pwr1)
            else:
                sumlj = r4(0.0)
            xc = x2
            sumc = qct * ato_manipulation_rev4.chsph(ra, rb, xc)

            # Apply termination functions

            sumlj = pp.refpotfac * sumlj * pp.wfac[j]
            sumc = pp.refpotfac * sumc * pp.wfacc[j]

            # Save results

            pp.epot[j, ic] = sumlj + sumc
            pp.ecoul[j, ic] = sumlj + sumc
            #                  if (coul_intra) xav(j,ic)=xav(j,ic)+coul_intra_fac*sumc

            # Add in Gaussians or modified Morse functions (Only the first can be a modified Morse function)

            ig = 0
            while ig < pp.npgauss[irefp]:

                ig = ig + 1

                wpg = abs(pp.wpgauss[ig, irefp])
                if wpg > 0.0:
                    rpg = abs(pp.rpgauss[ig, irefp])
                    hpg = pp.hpgauss[ig, irefp]
                    epsr16.add = r4(0.0)
                    if ig == 1 and pp.addpotential == "modmorse":
                        width = wpg
                        pos = rpg
                        g6, width = ato_manipulation_rev4.gauss6pot(x[j], pos, width)
                        epsr16.add = abs(hpg) * g6
                    else:
                        rdiff = (x[j] - rpg) / wpg
                        term = logf(abs(hpg)) - r4(0.5) * rdiff * rdiff
                        if abs(term) < 25.0:
                            epsr16.add = sign(r4(1.0), hpg) * expf(term)
                    #     write(6,*) j,x(j),rpg,wpg,hpg,rdiff,term,add
                    #                        epot(j,ic)=epot(j,ic)+add*wfac(j)
                    pp.ecoul[j, ic] = pp.ecoul[j, ic] + epsr16.add * pp.wfac[j]
            # Try this! Put full reference potential into intramolecular terms
            if molecom.coul_intra:
                pp.xav[j, ic] = molecom.coul_intra_fac * pp.ecoul[j, ic]
            x1 = x2
        #
        # add on the extra repulsive term for the intermolecular potential
        # then calculate the potential gradients.
        #
        # only do something if intensity below rminex(ip) was found
        #
        if abs(pp.pcoefex[irefp]) > 0.0 and abs(pp.rminex[irefp]) > 0.0:
            expecfuse = epsrinput_par.expecf * sqrtf(abs(pp.rminex[irefp]))
            if not pp.harmonic:
                # Exponential repulsive core
                rminextest = abs(pp.rminex[irefp])
                #                  call calc_pwrexpr(npt2,1,pdstep,pdmax
                #     *,expecf,rho,x,expgr,rminextest,ncfirst,1.0,rbroad)
                for j in range(1, npt2 + 1):
                    ddiff = (rminextest - x[j]) / expecfuse
                    if gradlimit > 0.0 and ddiff >= loggradlimit:  # Use a linear term to prevent exponential divergences
                        epsr16.add = pp.pcoefex[irefp] * gradlimit * (r4(1.0) + ddiff - loggradlimit)
                    else:
                        epsr16.add = pp.pcoefex[irefp] * expf(ddiff)
                    #                        epot(j,ic)=epot(j,ic)+add*wfac(j)
                    pp.ecoul[j, ic] = pp.ecoul[j, ic] + epsr16.add * pp.wfac[j]
            else:
                # Gaussian repulsive core
                rminextest = abs(pp.rminex[irefp]) + expecfuse
                for j in range(1, npt2 + 1):
                    ddiff = (rminextest - x[j]) / expecfuse
                    if ddiff > 0.0:
                        if gradlimit > 0.0 and ddiff >= gradlimit:  # Use a linear term to prevent exponential divergences
                            epsr16.add = r4(0.5) * pp.pcoefex[irefp] * gradlimit * (f32(2) * ddiff - gradlimit)
                        else:
                            epsr16.add = r4(0.5) * pp.pcoefex[irefp] * ddiff * ddiff
                        #                        epot(j,ic)=epot(j,ic)+add*wfac(j)
                        pp.ecoul[j, ic] = pp.ecoul[j, ic] + epsr16.add * pp.wfac[j]
        pp.epot[1, ic] = pp.epot[2, ic]
        pp.ecoul[1, ic] = pp.ecoul[2, ic]
        pp.xav[1, ic] = pp.xav[2, ic]
        #
        # set up gradients of these potentials
        #
        sumlj1 = pp.epot[1, ic]
        sumc1 = pp.ecoul[1, ic]
        sumx1 = pp.xav[1, ic]
        x1 = xstar[1]
        x2 = xstar[2]
        for j in range(2, npt + 1):
            x2 = xstar[j]
            xdif = x2 - x1
            sumlj2 = pp.epot[j, ic]
            sumc2 = pp.ecoul[j, ic]
            sumx2 = pp.xav[j, ic]
            pp.epot1[j - 1, ic] = (sumlj2 - sumlj1) / xdif
            pp.ecoul1[j - 1, ic] = (sumc2 - sumc1) / xdif
            pp.xav1[j - 1, ic] = (sumx2 - sumx1) / xdif
            x1 = x2
            sumlj1 = sumlj2
            sumc1 = sumc2
            sumx1 = sumx2
        pp.epot1[npt, ic] = 0.0
        pp.ecoul1[npt, ic] = 0.0
        pp.xav[npt, ic] = 0.0  # (xav, not xav1, as in the Fortran)
        # Get a list of atom pairs for this type. If rminex(irefp).ge.0.0 then the list is empty.
        containerforces_mod.get_containerforce_list(
            clustercom.nmol, ic, pp.rmaxex[irefp], it, jt, molecom.itype, molecom.qtype, molecom.natom,
            molecom.qtypea, molecom.ntype)
    if containerforces_mod_data.numcontainerforces > 0:
        for ic in range(1, containerforces_mod_data.numcontainerforces + 1):
            units.write(6, lw("generate_rp> ", ic, containerforces_mod_data.pairrefstore[1, ic],
                              containerforces_mod_data.pairrefstore[2, ic]))
    return


@routine
def generate_ep(changepotfac):
    """epsr_standalone_rev3.f:generate_ep -- generate EP from coefficients: rebuild xavs,
    its gradient xavs1 and the total potential xavtot (+ gradient) from pcoefp.
    """
    from ..library import ato_manipulation_rev4

    pp = potential_par
    x = epsr16.x
    xstar = epsr16.xstar
    npt = epsr16.npt
    npt2 = epsr16.npt2
    # logical, intent(in) :: changepotfac
    # real energabs,pressfac,abs_energy_ep,absediff

    # First evaluate the empirical potential with current coefficients
    #
    # maximum value of a sum - this it to prevent overflows
    #
    summax = _EXP25
    #
    # sum up the EP's over the different datafiles
    #
    nterms = pp.ncoeffp
    decayp = pp.psigma2
    #     write(6,*) nterms,decayp,rbroad,rho,xstar(10),dmax
    #
    # set up the transform matrix from power-exponentials to r-space
    #
    # NB: ncfirst is an uninitialised local in the Fortran; with rmin = 0.0 the callee
    # never assigns it and nothing reads it.  0 is used.
    ncfirst = 0
    if pp.gaussian:
        ncfirst = calc_gaussr_grad(npt2, nterms, pp.pdstep, pp.pdmax, decayp, epsr16.rho, x, pp.expgr, r4(0.0),
                                   ncfirst, r4(1.0), epsrinput_par.rbroad)
    else:
        ncfirst = calc_pwrexpr_grad(npt2, nterms, pp.pdmin, pp.pdstep, pp.pdmax, decayp, epsr16.rho, x, pp.expgr,
                                    r4(0.0), ncfirst, r4(1.0), epsrinput_par.rbroad)
    #
    # step through the partials and calculate the empirical potential
    #
    for ip in range(1, molecom.npairs + 1):
        irefp = pp.irefpcof[ip]
        #
        # Set up the potential coefficients, ignoring the first one
        #
        for ic in range(2, nterms + 1):
            pp.pwtold[ic] = pp.pcoefp[ic, irefp]
        pp.pwtold[1] = 0.0

        #
        # generate the potential corresponding to these coefficients
        #
        sum1 = f32(0.0)  # noqa: F841
        sums1 = f32(0.0)
        for j in range(1, npt2 + 1):
            term = expand(nterms, j, pp.expgr, pp.pwtold, summax)
            pp.xavs1[j, ip] = term * pp.pcoefp[1, irefp]

    # Calculate the potential from the gradients.
    #
    for ip in range(1, molecom.npairs + 1):

        for j in range(1, npt2 + 1):
            #
            # gradients need to be dimensionless, and multiply by truncation function
            #
            pp.xavs1[j, ip] = pp.xavs1[j, ip] * clustercom.dmax
        #
        # from gradients, calculate empirical potentials, integrating
        # from the last point
        #
        pp.xavs[npt, ip] = 0.0
        pp.xavs1[npt, ip] = 0.0
        x1 = xstar[npt]
        xts1 = pp.xavs1[npt, ip]
        for j in range(1, npt2 + 1):
            jref = npt - j
            x2 = xstar[jref]
            xts2 = pp.xavs1[jref, ip]
            pp.xavs[jref, ip] = pp.xavs[jref + 1, ip] - r4(0.5) * (x1 - x2) * (xts1 + xts2)
            x1 = x2
            xts1 = xts2

        # Multiply empirical potential by truncation function

        for j in range(1, npt2 + 1):
            #               xavs(j,ip)=xavs(j,ip)
            pp.xavs[j, ip] = pp.xavs[j, ip] * pp.wfac[j]
        #            if(ip.eq.1) write(6,*) 'generate_ep> ',ip,irefp,(xavs(j,ip),j=1,20)
    # Determine the absolute energy of the EP
    energabs = ato_manipulation_rev4.abs_energy_ep()
    pressfac = r4(1.0)
    if abs(pp.absolute_energy) == 0.0:
        pressfac = r4(0.0)
        energabs = r4(0.0)
    elif abs(energabs) > 0.0:
        if abs(energabs) < abs(pp.erequnit):
            pressfac = abs(pp.erequnit) / abs(energabs)
        else:
            pressfac = abs(pp.absolute_energy) / abs(energabs)
            if pressfac > 1.0 and pp.ereqstep == 0.0:
                pressfac = r4(1.0)
    units.write(6, lw("generate_ep> ", pp.absolute_energy, energabs, pressfac))
    # Apply potential factor
    multiply_pot_coeffs(pressfac)
    # Now re-evaluate the derivatives with these truncated potentials, multiplying by supplied factor
    for ip in range(1, molecom.npairs + 1):
        irefp = pp.irefpcof[ip]
        # Multiply by potential factor if required and create the total potential
        if pp.npgauss[irefp] == 0:
            for j in range(1, npt2 + 1):
                pp.xavs[j, ip] = pp.xavs[j, ip] * pressfac
                pp.xavtot[j, ip] = pp.ecoul[j, ip] + pp.xavs[j, ip]
        else:
            for j in range(1, npt2 + 1):
                pp.xavs[j, ip] = pp.xavs[j, ip] * pressfac
                pp.xavtot[j, ip] = pp.ecoul[j, ip] + pp.xavs[j, ip]
        sums1 = pp.xavs[1, ip]
        sumtot1 = pp.xavtot[1, ip]
        x1 = xstar[1]
        x2 = xstar[2]
        for j in range(2, npt + 1):
            x2 = xstar[j]
            xdif = x2 - x1
            sums2 = pp.xavs[j, ip]
            sumtot2 = pp.xavtot[j, ip]
            pp.xavs1[j - 1, ip] = (sums2 - sums1) / xdif
            pp.xavtot1[j - 1, ip] = (sumtot2 - sumtot1) / xdif
            x1 = x2
            sums1 = sums2
            sumtot1 = sumtot2
        pp.xavs1[npt, ip] = 0.0
        #            if(ip.eq.1) write(6,*) 'generate_ep> ',ip,irefp,(xavs1(j,ip),j=1,20)
        pp.xavtot1[npt, ip] = 0.0
    # Prove that you have done it!
    #      energabs=abs_energy_ep()
    #      write(6,*) 'generate_ep> ', absolute_energy,energabs
    # Aim to adjust the size of potfac so that changes of EP absolute energy are no larger than potfacamp*energabs
    if changepotfac and pp.absolute_energy > 0.0:  # CAUTION! This doesn't work with absolute_energy - need the absolute energy of the previous potential BEFORE it was perturbed.
        if pp.potfac > r4(0.000001) and pp.potfacamp > 0.0:
            absediff = (pp.absolute_energy - energabs) / abs(pp.potfacamp * pp.absolute_energy)
            units.write(6, lw("changepotfac> ", pp.potfacamp, pp.potfac, pp.absolute_energy, energabs, absediff))
            if abs(absediff) > r4(0.1):
                absediff = sign(r4(0.1), absediff)
            pp.potfac = pp.potfac * (r4(1.0) + absediff)
            pp.potfac = max(pp.potfac, r4(0.000001))  # Smallest allowed value of potfac
            pp.potfac = min(pp.potfac, r4(1.0))  # Maximum value of potfac
    return


def _lsp_read_record(nchan, dump, rfacwrite, x, y, aread, bread, cread, xread, yread):
    """read(nchan,*,iostat=ierr) dump,dump,rfacwrite,x,dump,y,dump,aread,bread,cread

    Items are assigned in order, so on an error the items already read keep their new
    values (as in gfortran).  ``xread``/``yread`` read x and y as real*8 (lsparabola) or
    real*4 (lsparabola1).  Returns (ierr, dump, rfacwrite, x, y, aread, bread, cread).
    """
    ierr = 0
    try:
        r = units.list_read(nchan)
        dump = r.real(dump)
        dump = r.real(dump)
        rfacwrite = r.real(rfacwrite)
        x = xread(r, x)
        dump = r.real(dump)
        y = yread(r, y)
        dump = r.real(dump)
        aread = r.real(aread)
        bread = r.real(bread)
        cread = r.real(cread)
    except FortranEOF:
        ierr = -1
    except FortranIOError as e:
        ierr = e.iostat
    return ierr, dump, rfacwrite, x, y, aread, bread, cread


@routine
def lsparabola(nchan, fname, nsum, navmax, epcfrac, a, b, c, rfacwrite, aav, bav, lastqf, nextqf):
    """epsr_standalone_rev3.f:lsparabola -- running-average straight-line (quadratic
    fitting is inhibited) fit to the .erg history on unit nchan; writes fname.

    Returns (nsum, a, b, c, rfacwrite, aav, bav, lastqf, nextqf).
    """
    # character*256 fname
    # real emin,femin,rfacmin
    # real*8 x,y,sum,sumx1,sumx2,sumx3,sumx4,sumy,sumyx1,sumyx2
    # *,delta,cplus,cminus,bplus,bminus,pf1,pf2
    # real :: dump,aread,bread,cread,thisqf
    # (xm and xmax are implicit real)
    x = 0.0
    y = 0.0
    dump = f32(0.0)
    aread = f32(0.0)
    bread = f32(0.0)
    cread = f32(0.0)

    sum = 0.0
    sumx1 = 0.0
    sumx2 = 0.0
    sumx3 = 0.0
    sumx4 = 0.0
    sumy = 0.0
    sumyx1 = 0.0
    sumyx2 = 0.0
    ierr = 0
    xmax = f32(0.0)
    nav = 0
    nsum = 0
    aav = r4(0.0)
    bav = r4(0.0)
    lastqf = r4(0.0)
    nextqf = r4(0.0)
    while ierr == 0:
        # This line reads the coefficient used to fit the quadratic.
        ierr, dump, rfacwrite, x, y, aread, bread, cread = _lsp_read_record(
            nchan, dump, rfacwrite, x, y, aread, bread, cread, lambda r, v: r.dble(v), lambda r, v: r.dble(v))
        if ierr == 0:
            nsum = nsum + 1
            y = float(epcfrac) * y + float((r4(1.0) - epcfrac) * rfacwrite)
            thisqf = f32(y)
            if nsum == 1:
                nextqf = thisqf
            if x > float(xmax):
                xmax = f32(x)
            if nsum > 1:  # Ignore first point!
                lastqf = nextqf
                nextqf = thisqf
                nav = nav + 1
                if nav > navmax:
                    nav = navmax

                # form running averages so that results are not overly biassed by earlier values

                pf1 = float(r4(1.0) / f32(nav))
                pf2 = float(f32(nav - 1) / f32(nav))
                sum = pf2 * sum + pf1
                sumy = pf2 * sumy + pf1 * y
                sumx1 = pf2 * sumx1 + pf1 * x
                sumyx1 = pf2 * sumyx1 + pf1 * y * x
                xm = f32(x * x)
                sumx2 = pf2 * sumx2 + pf1 * float(xm)
                sumyx2 = pf2 * sumyx2 + pf1 * float(xm) * y
                xm = f32(float(xm) * x)
                sumx3 = pf2 * sumx3 + pf1 * float(xm)
                xm = f32(float(xm) * x)
                sumx4 = pf2 * sumx4 + pf1 * float(xm)
                aav = f32(pf2 * float(aav) + pf1 * float(aread))
                bav = f32(pf2 * float(bav) + pf1 * float(abs(bread)))

    # Backspace the file so that it is before the end of file marker

    units.backspace(nchan)

    # Now work out the coefficients the hard way! The method used here assumes all the sums are
    # positive definite

    # a) determinant of coefficients

    delta = (sum * (sumx2 * sumx4 - sumx3 * sumx3)
             + sumx1 * (sumx3 * sumx2 - sumx1 * sumx4)
             + sumx2 * (sumx1 * sumx3 - sumx2 * sumx2))
    delta = 0.0  # Inhibit quadratic fitting

    if abs(delta) > float(r4(1.0e-10)) and nsum >= navmax:

        # b) quadratic coefficient

        cplus = sumy * sumx1 * sumx3 + sumyx1 * sumx2 * sumx1 + sumyx2 * sum * sumx2
        cminus = sumy * sumx2 * sumx2 + sumyx1 * sumx3 * sum + sumyx2 * sumx1 * sumx1

        if cplus != cminus:

            c = f32((cplus - cminus) / delta)

        else:

            c = r4(0.0)

        # calculate the other coefficients - a straight line fit if c is zero

        sumy = sumy - float(c) * sumx2
        sumyx1 = sumyx1 - float(c) * sumx3
        delta = sumx1 * sumx1 - sum * sumx2
        if abs(delta) > float(r4(1.0e-10)):

            bplus = sumy * sumx1
            bminus = sumyx1 * sum
            if bplus != bminus:

                b = f32((bplus - bminus) / delta)

            else:

                b = r4(0.0)

        else:

            b = r4(0.0)

        sumy = sumy - float(b) * sumx1

        if abs(sum) > 0.0:

            a = f32(sumy / sum)

        else:

            a = r4(0.0)

    elif nsum >= 2:

        c = r4(0.0)

        # Try fitting a straight line instead

        delta = sumx1 * sumx1 - sum * sumx2
        if abs(delta) > float(r4(1.0e-10)):

            bplus = sumy * sumx1
            bminus = sumyx1 * sum
            if bplus != bminus:

                b = f32((bplus - bminus) / delta)

            else:

                b = r4(0.0)

        else:

            b = r4(0.0)

        sumy = sumy - float(b) * sumx1

        if abs(sum) > 0.0:

            a = f32(sumy / sum)

        else:

            a = r4(0.0)

    else:

        c = r4(0.0)
        b = r4(0.0)
        delta = sum

        if abs(delta) > 0.0:

            a = f32(sumy / delta)

        else:

            a = r4(0.0)

    units.write(6, lw("lsparabola> ", nsum, delta, a, b, c))

    # write a little file to show the fit

    units.open(8, fname, status="unknown")

    npt = 100
    xstep = xmax / f32(npt)

    # Calculate position and depth/height of minimum/maximum
    if c != 0.0:
        emin = -b / (r4(2.0) * c)
        femin = a + emin * (b + emin * c)
    else:
        emin = r4(0.0)
        femin = r4(0.0)
    for ix in range(1, npt + 1):

        x = float(f32(ix) * xstep)

        # format here has to match that of the .erg file so that it can be plotted with
        # the plot command.

        units.write(8, fmt(_F100_LSPARABOLA, a, b, c, x, float(a) + x * (float(b) + x * float(c)), emin, femin))

        # 100   format(7(1x,e13.6))

    units.close(8)

    return nsum, a, b, c, rfacwrite, aav, bav, lastqf, nextqf


@routine
def lsparabola1(nchan, fname, nav, navmax, epcfrac, a, b, c, rfacwrite, aav, bav, lastqf, nextqf):
    """epsr_standalone_rev3.f:lsparabola1 -- as lsparabola, but instead of using running
    averages we average only over the last navmax points.

    Returns (nav, a, b, c, rfacwrite, aav, bav, lastqf, nextqf).
    """
    # real, dimension(:), allocatable :: xval,yval,aval,bval
    # character*256 fname
    # real :: x,y,emin,femin,rfacmin,dump,aread,bread,cread,fractemp,thisqf
    # real*8 :: xm,sum,sumx1,sumx2,sumx3,sumx4,sumy,sumyx1,sumyx2
    # *,delta,cplus,cminus,bplus,bminus
    # integer nsum,ierr,i
    x = f32(0.0)
    y = f32(0.0)
    dump = f32(0.0)
    aread = f32(0.0)
    bread = f32(0.0)
    cread = f32(0.0)
    # NB: navref is not initialised in the Fortran and is written below even if no
    # point after the first was read; 0 is used (Rust potential.rs lsparabola1).
    navref = 0

    xval = real_array(navmax)
    yval = real_array(navmax)
    aval = real_array(navmax)
    bval = real_array(navmax)
    ierr = 0
    xmax = f32(0.0)
    nav = 0
    nsum = 0
    lastqf = r4(0.0)
    nextqf = r4(0.0)
    while ierr == 0:
        # This line reads the coefficient used to fit the quadratic.
        ierr, dump, rfacwrite, x, y, aread, bread, cread = _lsp_read_record(
            nchan, dump, rfacwrite, x, y, aread, bread, cread, lambda r, v: r.real(v), lambda r, v: r.real(v))
        #         fractemp=rfacwrite/(rfacwrite+y)
        #         y=epcfrac*fractemp*y+(1.0-epcfrac)*(1.0-fractemp)*rfacwrite
        if ierr == 0:
            nsum = nsum + 1
            y = epcfrac * y + (r4(1.0) - epcfrac) * rfacwrite
            thisqf = y
            if nsum == 1:
                nextqf = thisqf
            if x > xmax:
                xmax = x
            if nsum > 1:  # Ignore first point!
                lastqf = nextqf
                nextqf = thisqf
                nav = nav + 1
                navref = nav - navmax * idiv(nav - 1, navmax)
                xval[navref] = x
                yval[navref] = y
                aval[navref] = aread
                bval[navref] = bread
    units.write(6, fmt("(a,5(1x,i5))", "lsparabola1> ", nsum, nav, navref))
    # Backspace the file so that it is before the end of file marker
    units.backspace(nchan)
    nsum = min(nav, navmax)
    sum = 0.0
    sumx1 = 0.0
    sumx2 = 0.0
    sumx3 = 0.0
    sumx4 = 0.0
    sumy = 0.0
    sumyx1 = 0.0
    sumyx2 = 0.0
    aav = r4(0.0)
    bav = r4(0.0)
    # Form averages
    if nsum > 0:
        for i in range(1, nsum + 1):
            x = xval[i]
            y = yval[i]
            sum = sum + 1.0
            sumy = sumy + float(y)
            sumx1 = sumx1 + float(x)
            sumyx1 = sumyx1 + float(y * x)
            xm = float(x * x)
            sumx2 = sumx2 + xm
            sumyx2 = sumyx2 + xm * float(y)
            xm = xm * float(x)
            sumx3 = sumx3 + xm
            xm = xm * float(x)
            sumx4 = sumx4 + xm
            aav = aav + aval[i]
            bav = bav + abs(bval[i])
        aav = aav / f32(nsum)
        bav = bav / f32(nsum)
    if xval is not None:
        xval = yval = aval = bval = None
    # Now work out the coefficients the hard way! The method used here assumes all the sums are
    # positive definite

    # a) determinant of coefficients

    delta = (sum * (sumx2 * sumx4 - sumx3 * sumx3)
             + sumx1 * (sumx3 * sumx2 - sumx1 * sumx4)
             + sumx2 * (sumx1 * sumx3 - sumx2 * sumx2))
    delta = 0.0  # Inhibit quadratic fitting

    if abs(delta) > float(r4(1.0e-10)) and nsum >= navmax:

        # b) quadratic coefficient

        cplus = sumy * sumx1 * sumx3 + sumyx1 * sumx2 * sumx1 + sumyx2 * sum * sumx2
        cminus = sumy * sumx2 * sumx2 + sumyx1 * sumx3 * sum + sumyx2 * sumx1 * sumx1

        if cplus != cminus:

            c = f32((cplus - cminus) / delta)

        else:

            c = r4(0.0)

        # calculate the other coefficients - a straight line fit if c is zero

        sumy = sumy - float(c) * sumx2
        sumyx1 = sumyx1 - float(c) * sumx3
        delta = sumx1 * sumx1 - sum * sumx2
        if abs(delta) > float(r4(1.0e-10)):

            bplus = sumy * sumx1
            bminus = sumyx1 * sum
            if bplus != bminus:

                b = f32((bplus - bminus) / delta)

            else:

                b = r4(0.0)

        else:

            b = r4(0.0)

        sumy = sumy - float(b) * sumx1

        if abs(sum) > 0.0:

            a = f32(sumy / sum)

        else:

            a = r4(0.0)

    elif nav >= 2:

        c = r4(0.0)

        # Try fitting a straight line instead

        delta = sumx1 * sumx1 - sum * sumx2
        if abs(delta) > float(r4(1.0e-10)):

            bplus = sumy * sumx1
            bminus = sumyx1 * sum
            if bplus != bminus:

                b = f32((bplus - bminus) / delta)

            else:

                b = r4(0.0)

        else:

            b = r4(0.0)

        sumy = sumy - float(b) * sumx1

        if abs(sum) > 0.0:

            a = f32(sumy / sum)

        else:

            a = r4(0.0)

    else:

        c = r4(0.0)
        b = r4(0.0)
        delta = sum

        if abs(delta) > 0.0:

            a = f32(sumy / delta)

        else:

            a = r4(0.0)

    units.write(6, lw("lsparabola1> ", nsum, delta, a, b, c))

    # write a little file to show the fit

    units.open(8, fname, status="unknown")

    npt = 100
    xstep = xmax / f32(npt)

    # Calculate position and depth/height of minimum/maximum
    if c != 0.0:
        emin = -b / (r4(2.0) * c)
        femin = a + emin * (b + emin * c)
    else:
        emin = r4(0.0)
        femin = r4(0.0)
    for ix in range(1, npt + 1):

        x = f32(ix) * xstep

        # format here has to match that of the .erg file so that it can be plotted with
        # the plot command.

        units.write(8, fmt(_F100_LSPARABOLA, a, b, c, x, a + x * (b + x * c), emin, femin))

        # 100   format(7(1x,e13.6))

    units.close(8)

    return nav, a, b, c, rfacwrite, aav, bav, lastqf, nextqf


@routine
def multiply_pot_coeffs(epfact):
    """epsr_standalone_rev3.f:multiply_pot_coeffs -- multiply the EP coefficients
    pcoefp(2:ncoeffp, .) of every pair by epfact.
    """
    # real epfact
    nterms = potential_par.ncoeffp
    pcoefp = potential_par.pcoefp
    #
    # Step through the defined partials
    #
    for ip in range(1, molecom.npairs + 1):
        irefp = potential_par.irefpcof[ip]
        # The amplitude of the present EP is multiplied by the specified factor, but only
        # if the number of Gaussians defined is zero
        for ic in range(2, nterms + 1):
            pcoefp[ic, irefp] = pcoefp[ic, irefp] * epfact

    return


@routine
def multiply_pcoefex(epfact):
    """epsr_standalone_rev3.f:multiply_pcoefex -- scale the repulsive-core radius rminex
    of every pair with a repulsive core by epfact, update its symbol value and
    regenerate the reference potential.
    """
    pp = potential_par
    # real epfact

    #
    # Step through the defined partials
    #
    for ip in range(1, molecom.npairs + 1):
        irefp = pp.irefpcof[ip]
        # The amplitude of the present repulsive potential is multiplied by the specified factor
        # However the minimum amplitude is limited to 0.0001. This is to prevent the values
        # becoming absurdly small. If this still gives too much pressure then reduce the size
        # of rminex.
        if pp.pcoefex[irefp] != 0.0:
            if (abs(pp.pcoefex[irefp]) > r4(0.0001) or epfact > 1.0) and pp.rminex[irefp] > 0.0:
                #                  write(6,*) irefp,rminex(irefp),epfact
                pp.rminex[irefp] = pp.rminex[irefp] * epfact
                units.write(6, lw("multiply_pcoefex> ", irefp, pp.rminex[irefp], epfact))
                # Write the new value of this symbol, keeping the other values intact
                icount = symbols.nsymboluptodata + symbols.nfiles * symbols.nsymbolplus + irefp
                commands.sometext = fstr(symbols.symbol_value[icount], commands.nstsize)
                commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                              commands.nstsize, commands.nworddim)
                # write(symbol_value(icount),*,iostat=ierr) ...  (internal list-directed
                # write into a character*256 element)
                if commands.nwords > 1:
                    text = lw(pp.rminex[irefp],
                              substr(commands.sometext, commands.ncf[2], commands.ncl[commands.nwords]))
                else:
                    text = lw(pp.rminex[irefp])
                symbols.symbol_value[icount] = fstr(text, 256)
    generate_rp()

    return
