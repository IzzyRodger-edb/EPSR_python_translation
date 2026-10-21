"""Translation of fortran_src/ato_manipulation_rev4.f, lines 1-2343 (package E).

Module data (``module timer``, ``bragginten_mod``, ``reciprocallattice_mod``,
``diffuse_mod``) lives in ``common/<module>.py``; their module routines
(``initialise_timer``, ``initialise_diffuse_values``) are here.  Then the free routines
``ftype`` ... ``truncatecweeks`` and the molecule builder ``makemole``.

Lines 2344-end (``fmolec1``, ``get_rotgrp``, ``get_whole_rotgrp``, ``atogen``, ``press``,
``bragginten``, ``generate_mp`` ...) are translated in ``ato_manipulation_rev4_part2.py``
and are merged into this file; they are called here by bare name, as in the Fortran.
(``get_rotgrp``/``get_whole_rotgrp`` appear inside makemole only as INTERFACE blocks.)

Translation of fortran_src/ato_manipulation_rev4.f, lines 2344 to the end (package F).

fmolec1, get_direct_bonds, get_rotgrp, get_whole_rotgrp, rotatecoords, dihedral,
generate_mp, crystfrac, crystarg, bragginten, radd, sharm_diffuse, jhalf, truncconv,
change_boxsize, atogen, press, stog, lorch, abs_energy_ep, abs_energy_ep_old,
adjust_intra_wts, get_maxmolsize.

This file is merged by the integrator into ``library/ato_manipulation_rev4.py``
(lines 1-2343, package E).  Routines of the first part (``en_tethered``,
``generate_fractional_coordinates``, ``generate_cartesian_coordinates``,
``initialise_diffuse_values``) are called by bare name, as within one module; ruff
reports them as F821 here until the merge (see docs/translation_notes/F.md).

Library modules whose name equals a Fortran-module data module in ``common`` are
imported with a ``_lib`` suffix (``ergpairreference_lib`` holds the routines of
``ergpairreference.f``, ``ergpairreference`` the module data).
"""

from ..common import clustercom
from ..common import commands
from ..common import diffuse_mod
from ..common import directory_name
from ..common import epsr_dimension
from ..common import lattice
from ..common import molecom
from ..common import molecom_mod
from ..common import reciprocallattice_mod
from ..common import sears91_gudrun
from ..common import symbols
from ..common import timer
from ..common import bragginten_mod
from ..common import cellneighbourlist
from ..common import containerforces_mod
from ..common import coordconstraints_mod
from ..common import epsr16
from ..common import epsrinput_par
from ..common import ergpairreference
from ..common import potential_par
from ..common import threebody_mod
from ..common import weights
from ..common import x_rayf0
from ..fortran import Ch
from ..fortran import FortranEOF
from ..fortran import FortranIOError
from ..fortran import FortranStop
from ..fortran import PI4ATAN
from ..fortran import cosf
from ..fortran import erfcf
from ..fortran import expf
from ..fortran import f32
from ..fortran import findex
from ..fortran import findex_padded
from ..fortran import fmt
from ..fortran import fstr
from ..fortran import idiv
from ..fortran import int_array
from ..fortran import lw
from ..fortran import nint
from ..fortran import pow4
from ..fortran import powf
from ..fortran import r4
from ..fortran import read_internal
from ..fortran import real_array
from ..fortran import char_array
from ..fortran import setsub
from ..fortran import sqrtf
from ..fortran import substr
from ..fortran import units
from ..fortran import atan2f
from ..fortran import cosd_intrinsic
from ..fortran import ifix
from ..fortran import logf
from ..fortran import pow3
from ..fortran import sign
from ..fortran import sind_intrinsic
from ..fortran import sinf
from ..fortran import sinhf
from ..trace import routine
from . import ato_read_and_save
from . import auxiliary_routines
from . import file_and_directory_commands
from . import minimagec
from . import parse
from . import lngamma_mod as lngamma_mod_lib
from . import molecom_mod as molecom_mod_lib
from . import reallocation_routines
from . import molfit
from . import cellneighbourlist as cellneighbourlist_lib
from . import containerforces_mod as containerforces_mod_lib
from . import coordconstraint_mod as coordconstraint_mod_lib
from . import ergpairreference as ergpairreference_lib
from . import threebody_mod as threebody_mod_lib


# =============================================================================
# from ato_manipulation_rev4.py
# =============================================================================

# --------------------------------------------------------------------------- private helpers
def _do(start, end, step=1):
    """Iteration values of ``do i = start, end, step`` (any sign of step).

    A zero step raises ValueError (the Fortran trip count divides by the step)."""
    start = int(start)
    end = int(end)
    step = int(step)
    if step > 0:
        return range(start, end + 1, step)
    return range(start, end - 1, step)


def _iread(text, kind, old):
    """``read(text,*,iostat=ierr) v`` of one item from an internal record.

    ``kind`` is 'int' or 'real'.  Returns (v, ierr); v keeps its old value on error."""
    try:
        r = read_internal(text)
        v = r.int(old) if kind == "int" else r.real(old)
        return v, 0
    except FortranEOF:
        return old, -1
    except FortranIOError as e:
        return old, e.iostat


def _read_a(unit, old, n):
    """``read(unit,'(a)',iostat=ierr) s`` into a character*n variable: (s, ierr)."""
    try:
        return fstr(units.read_line(unit), n), 0
    except FortranEOF:
        return old, -1
    except FortranIOError as e:
        return old, e.iostat


def _read5_real(old):
    """``read(5,*,iostat=ierr) x`` for one real: (x, ierr)."""
    try:
        r = units.list_read(5)
        return r.real(old), 0
    except FortranEOF:
        return old, -1
    except FortranIOError as e:
        return old, e.iostat


# =========================================================================== module timer
@routine
def initialise_timer():
    """ato_manipulation_rev4.f:initialise_timer (module timer) -- zero the summed times."""
    timer.sumtime1 = 0
    timer.sumtime2 = 0
    return


# =========================================================================== module diffuse_mod
@routine
def initialise_diffuse_values(nfile):
    """ato_manipulation_rev4.f:initialise_diffuse_values (module diffuse_mod).

    Counts the (l,m) terms of the diffuse scattering, tabulates ln(gamma) and allocates
    the radial coefficient arrays.
    """
    # use lngamma_mod
    # use bragginten_mod
    # use reciprocallattice_mod
    # integer mdiffmax,maxlngamma

    # First count the number of terms and find the maximum ln(gamma) function required
    diffuse_mod.ndiffuse = 0
    maxlngamma = 0
    for l in _do(0, diffuse_mod.ldiffmax, diffuse_mod.ldiffstep):
        mdiffmax = l
        if diffuse_mod.mdiffstep == 0:
            mdiffmax = 0
        for m in range(0, mdiffmax + 1):
            diffuse_mod.ndiffuse = diffuse_mod.ndiffuse + 1
            # Get the maximum value of the ln(gamma) function to be called
            maxlngamma = max(maxlngamma, auxiliary_routines.maxdlm(l, m, 0), auxiliary_routines.maxdlm(l, -m, 0))
    lngamma_mod_lib.fctrln(maxlngamma)
    # Now allocate arrays and save the values
    if diffuse_mod.l1ex is not None:
        diffuse_mod.l1ex = None
        diffuse_mod.n1ex = None
    diffuse_mod.l1ex = int_array(diffuse_mod.ndiffuse)
    diffuse_mod.n1ex = int_array(diffuse_mod.ndiffuse)
    diffuse_mod.ndiffuse = 0
    for l in _do(0, diffuse_mod.ldiffmax, diffuse_mod.ldiffstep):
        mdiffmax = l
        if diffuse_mod.mdiffstep == 0:
            mdiffmax = 0
        # NB: with mdiffstep = 0 this is a DO loop with a zero step (undefined in the
        # Fortran, a division by zero in gfortran's trip count); _do raises ValueError.
        for m in _do(0, mdiffmax, diffuse_mod.mdiffstep):
            diffuse_mod.ndiffuse = diffuse_mod.ndiffuse + 1
            diffuse_mod.l1ex[diffuse_mod.ndiffuse] = l
            diffuse_mod.n1ex[diffuse_mod.ndiffuse] = m
    # Get the maximum Q of this scattering
    qmaxdiff = f32(0.0)
    for ib in range(1, bragginten_mod.nhkl + 1):
        if bragginten_mod.qhkl[ib] > qmaxdiff:
            qmaxdiff = bragginten_mod.qhkl[ib]
    diffuse_mod.qdiffstep = max(reciprocallattice_mod.ampra, reciprocallattice_mod.amprb, reciprocallattice_mod.amprc)
    diffuse_mod.nqdiff = nint(qmaxdiff / diffuse_mod.qdiffstep) + 1
    # Set up the arrays to store the radial coefficients
    if diffuse_mod.diffcoeffc is not None:
        diffuse_mod.qdiff = None
        diffuse_mod.diffcoeffc = None
        diffuse_mod.diffcoeffs = None
        diffuse_mod.ndiffcoeff = None
    diffuse_mod.qdiff = real_array(diffuse_mod.nqdiff)
    diffuse_mod.diffcoeffc = real_array(diffuse_mod.ndiffuse, diffuse_mod.nqdiff, nfile)
    diffuse_mod.diffcoeffs = real_array(diffuse_mod.ndiffuse, diffuse_mod.nqdiff, nfile)
    diffuse_mod.ndiffcoeff = int_array(diffuse_mod.nqdiff)
    for iff in range(1, nfile + 1):
        for iq in range(1, diffuse_mod.nqdiff + 1):
            diffuse_mod.qdiff[iq] = f32(iq - 1) * diffuse_mod.qdiffstep
            for id in range(1, diffuse_mod.ndiffuse + 1):
                diffuse_mod.ndiffcoeff[iq] = 0
                diffuse_mod.diffcoeffc[id, iq, iff] = 0.0
                diffuse_mod.diffcoeffs[id, iq, iff] = 0.0
    return


# =========================================================================== ftype
_F100_ftype = "(1x,'Tried to generate more than ',i3,' types of atom')"


@routine
def ftype(ia):
    """ato_manipulation_rev4.f:ftype -- assign (or create) the atom type of atom ia.

    Also sets qtype(ia) (1 for a "Q"-type atom, label starting with Q or q).
    """
    # use reallocation_routines
    # use molecom_mod
    # character*1 qtypet
    #
    # is this atom of "Q"-type?
    #
    qtypet = fstr(substr(molecom.alabel[ia], 1, 1), 1)
    if qtypet == "Q" or qtypet == "q":
        molecom.qtype[ia] = 1
    else:
        molecom.qtype[ia] = 0
    ifound = 0
    it = 0
    while it < molecom.ntype and ifound == 0:
        it = it + 1
        if molecom.alabel[ia] == molecom.type[it]:
            ifound = it
    if ifound > 0:
        molecom.itype[ia] = ifound
    elif molecom.ntype < epsr_dimension.near:
        molecom.ntype = molecom.ntype + 1
        molecom.type[molecom.ntype] = molecom.alabel[ia]
        molecom.qtypea[molecom.ntype] = molecom.qtype[ia]
        molecom.itype[ia] = molecom.ntype
    else:
        units.write(6, fmt(_F100_ftype, epsr_dimension.near))
        raise FortranStop(0)
    #      write(6,*) ntype,size(rhosphere)
    return


# =========================================================================== ctype
_F100_ctype = "(1x,'Mismatched types between atomic coordinates and ','potential values: ',a)"


@routine
def ctype(a1, ats, iexch, eread, wread, amasst, qt, rqt,
          natomsph=None, rhosph=None, radminsph=None, radmaxsph=None, pseudofacsph=None):
    """ato_manipulation_rev4.f:ctype -- store the potential parameters of atom type a1.

    natomsph, rhosph, radminsph, radmaxsph, pseudofacsph are Fortran ``optional``
    (keyword arguments, ``None`` = not present).
    """
    # use reallocation_routines
    # use molecom_mod
    # character*3 a1,ats
    # integer get_atomicnumber,n1,it
    # real, intent(in),optional :: natomsph,rhosph,radminsph,radmaxsph,pseudofacsph
    # Allocate space for the coarse-graining parameters that may be defined
    molecom_mod.natomsphere = reallocation_routines.reallocate1d_r(molecom_mod.natomsphere, molecom.ntype)
    molecom_mod.rhosphere = reallocation_routines.reallocate1d_r(molecom_mod.rhosphere, molecom.ntype)
    molecom_mod.radminsphere = reallocation_routines.reallocate1d_r(molecom_mod.radminsphere, molecom.ntype)
    molecom_mod.radmaxsphere = reallocation_routines.reallocate1d_r(molecom_mod.radmaxsphere, molecom.ntype)
    molecom_mod.pseudofacmol = reallocation_routines.reallocate1d_r(molecom_mod.pseudofacmol, molecom.ntype)
    n1 = 0
    it = 0
    while it < molecom.ntype and n1 == 0:
        it = it + 1
        if a1 == molecom.type[it]:
            n1 = it
    if n1 > 0:
        molecom.iatomicnotype[n1] = auxiliary_routines.get_atomicnumber(ats)
        molecom.atsymbol[n1] = fstr(ats, 3)
        molecom.iexchange[n1] = iexch
        molecom.epscor[n1] = eread
        molecom.sigcor[n1] = wread
        molecom.amassa[n1] = amasst
        molecom.charge[n1] = qt
        molecom.qradius[n1] = rqt
        if natomsph is not None:
            molecom_mod.natomsphere[n1] = natomsph
        else:
            molecom_mod.natomsphere[n1] = 1.0
        if rhosph is not None:
            # This value is not written to .ato file and so is not saved - calculated on reading .ato file.
            molecom_mod.rhosphere[n1] = rhosph
        else:
            molecom_mod.rhosphere[n1] = 0.0
        if radminsph is not None:
            molecom_mod.radminsphere[n1] = radminsph
        else:
            molecom_mod.radminsphere[n1] = 0.0
        if radmaxsph is not None:
            molecom_mod.radmaxsphere[n1] = radmaxsph
        else:
            molecom_mod.radmaxsphere[n1] = 0.5 * molecom.sigcor[n1]
        if pseudofacsph is not None:
            molecom_mod.pseudofacmol[n1] = pseudofacsph
        else:
            # Represents how many instances of this molecule make up a coarse-grained sphere.
            # Should be the same for all atom types in a given molecule.
            molecom_mod.pseudofacmol[n1] = 1.0
    else:
        units.write(6, fmt(_F100_ctype, Ch(a1, 3)))
    return


# =========================================================================== scatt_len
@routine
def scatt_len(chsym, massno, ifind):
    """ato_manipulation_rev4.f:scatt_len -- neutron scattering length of an isotope.

    Returns (scatt_len, ifind); ifind = 1 if the (symbol, mass number) pair was found.
    """
    # character*3 chsym,chtest
    # integer*4 massno,masstest
    # include 'sears91_gudrun.inc'
    nscatt = sears91_gudrun.mscatt
    ifind = 0
    for it in range(1, nscatt + 1):
        chtest = fstr(sears91_gudrun.chsymbol[it], 3)
        masstest = int(sears91_gudrun.chmassno[it])
        if chsym == chtest and massno == masstest:
            ifind = 1
            return sears91_gudrun.chbscat[it], ifind
    return f32(0.0), ifind


# =========================================================================== fractcoord / wholecoord
@routine
def fractcoord(a, b):
    """ato_manipulation_rev4.f:fractcoord -- a = a/b.  Returns a."""
    a = a / b
    return a


@routine
def wholecoord(a, b):
    """ato_manipulation_rev4.f:wholecoord -- a = a*b.  Returns a."""
    a = a * b
    return a


# =========================================================================== fractional / cartesian
@routine
def generate_fractional_coordinates(x, y, z):
    """ato_manipulation_rev4.f:generate_fractional_coordinates.

    Fractional coordinates relative to the lattice vectors a,b,c of the Cartesian
    (x,y,z), folded by the minimum image convention.  Returns (x, y, z).
    """
    # include 'lattice.inc'
    # real :: fa,fb,fc
    #
    # subroutine to generate the fractional coordinates to an arbitrary set of lattice vectors
    # a,b,c, for input Cartesian coordinates.
    #
    # fractional coordinates relative to (a,b,c)
    #
    fc = z / lattice.clatz
    fb = (y - fc * lattice.claty) / lattice.blaty
    fa = (x - fb * lattice.blatx - fc * lattice.clatx) / lattice.alatx
    #
    # minimum image convention
    #
    fa, fb, fc = minimagec.minimage(fa, fb, fc)
    #
    # calculate corrected coordinates
    #
    z = fc
    y = fb
    x = fa

    return x, y, z


@routine
def generate_cartesian_coordinates(x, y, z):
    """ato_manipulation_rev4.f:generate_cartesian_coordinates.

    Cartesian coordinates of the fractional (x,y,z).  Returns (x, y, z).
    """
    # include 'lattice.inc'
    # real :: fa,fb,fc
    #
    # subroutine to generate the fractional coordinates to an arbitrary set of lattice vectors
    # a,b,c, for input Cartesian coordinates.
    #
    # fractional coordinates relative to (a,b,c)
    #
    fc = z
    fb = y
    fa = x
    #
    # calculate Cartesian coordinates
    #
    z = fc * lattice.clatz
    y = fb * lattice.blaty + fc * lattice.claty
    x = fa * lattice.alatx + fb * lattice.blatx + fc * lattice.clatx

    return x, y, z


# =========================================================================== centre of mass
@routine
def centreofmassi(im):
    """ato_manipulation_rev4.f:centreofmassi -- put the centre of mass of molecule im at zero.

    For internal atom moves we only the adjust the atom coordinates for any
    change in the centre of mass.
    """
    # use molecom_mod
    #
    # determine coordinates of molecule centre of mass
    #
    x1 = f32(0.0)
    y1 = f32(0.0)
    z1 = f32(0.0)
    summass = f32(0.0)
    for ia in range(clustercom.naf[im], clustercom.nal[im] + 1):
        if molecom.qtype[ia] == 0:
            amass = molecom.amassa[molecom.itype[ia]]
            x1 = x1 + molecom.xa[ia] * amass
            y1 = y1 + molecom.ya[ia] * amass
            z1 = z1 + molecom.za[ia] * amass
            summass = summass + amass
    if summass > 0.0:
        x1 = x1 / summass
        y1 = y1 / summass
        z1 = z1 / summass
    #
    # adjust coordinates of atoms to put centre of mass at zero
    #
    for ia in range(clustercom.naf[im], clustercom.nal[im] + 1):
        molecom.xa[ia] = molecom.xa[ia] - x1
        molecom.ya[ia] = molecom.ya[ia] - y1
        molecom.za[ia] = molecom.za[ia] - z1

    return


@routine
def centreofmass(im):
    """ato_manipulation_rev4.f:centreofmass -- move the centre of mass of molecule im
    into xadd/yadd/zadd.

    For whole molecule moves or moves where part of a molecule is rotated we need
    to adjust the centre of mass coordinates as well as the atom coordinates.
    """
    # use molecom_mod
    # real coordfactor,x1,y1,z1
    #
    # determine coordinates of molecule centre of mass
    #
    x1 = f32(0.0)
    y1 = f32(0.0)
    z1 = f32(0.0)
    summass = f32(0.0)
    for ia in range(clustercom.naf[im], clustercom.nal[im] + 1):
        if molecom.qtype[ia] == 0:
            amass = molecom.amassa[molecom.itype[ia]]
            x1 = x1 + molecom.xa[ia] * amass
            y1 = y1 + molecom.ya[ia] * amass
            z1 = z1 + molecom.za[ia] * amass
            summass = summass + amass
    if summass > 0.0:
        x1 = x1 / summass
        y1 = y1 / summass
        z1 = z1 / summass
    #
    # adjust coordinates of atoms to put centre of mass at zero
    # Remember to adjust for molecular size factor if relevant!
    #
    #      coordfactor=dmaxmol/dmax
    #      x1=x1*coordfactor
    #      y1=y1*coordfactor
    #      z1=z1*coordfactor
    for ia in range(clustercom.naf[im], clustercom.nal[im] + 1):
        molecom.xa[ia] = molecom.xa[ia] - x1
        molecom.ya[ia] = molecom.ya[ia] - y1
        molecom.za[ia] = molecom.za[ia] - z1

    # Adjust the centre of mass coordinates

    clustercom.xadd[im] = clustercom.xadd[im] + x1
    clustercom.yadd[im] = clustercom.yadd[im] + y1
    clustercom.zadd[im] = clustercom.zadd[im] + z1

    return


# =========================================================================== en_tethered
@routine
def en_tethered(im, ia):
    """ato_manipulation_rev4.f:en_tethered -- squared distance of molecule im (or of its
    atom ia, if 0 < ia <= nato(im)) from its tethering position.

    With tethered molecules the energy simply derives from a harmonic force law wrt the
    tethering position.  Returns en_tethered.
    """
    # use molecom_mod
    # real en_tethered,xdif,ydif,zdif
    # integer ia,iaref

    # Get the position of this molecule relative to the tether position

    xdif = clustercom.xadd[im] - clustercom.phix[im]
    ydif = clustercom.yadd[im] - clustercom.phiy[im]
    zdif = clustercom.zadd[im] - clustercom.phiz[im]
    if ia > 0 and ia <= clustercom.nato[im]:
        iaref = int(clustercom.naf[im]) + ia - 1
        xdif = xdif + molecom.xa[iaref]
        ydif = ydif + molecom.ya[iaref]
        zdif = zdif + molecom.za[iaref]

    # Calculate the change in energy

    xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
    en_tethered = (xdif * xdif + ydif * ydif + zdif * zdif)

    return en_tethered


# =========================================================================== potential functions
@routine
def gauss6pot(r, sig, sigg):
    """ato_manipulation_rev4.f:gauss6pot -- Gaussian core minus an r**-6 term built from
    the overlapping-sphere Coulomb function.

    Returns (gauss6pot, sigg): sigg is clipped to 0.29 in the caller's variable.
    """
    # real gauss6pot,sig,sigg,r,coul,london,arg,ampgln,x,term

    x = r / sig
    coul = chsph(r4(0.5), r4(0.5), x)
    # coul**6: gfortran expands the constant power as (c*c*c)*(c*c*c) (Rust potential.rs)
    c3 = coul * coul * coul
    london = -(c3 * c3)
    if sigg > 0.29:
        sigg = r4(0.29)
    ampgln = 0.5 / (sigg * sigg)
    arg = x / sigg
    arg = 0.5 * arg * arg - ampgln
    term = f32(0.0)
    if arg < 30.0:
        term = expf(-arg)
    gauss6pot = term + london

    return gauss6pot, sigg


@routine
def chsph(sig1, sig2, r):
    """ato_manipulation_rev4.f:chsph -- Coulomb interaction of two uniformly charged
    spheres of radii sig1, sig2 at separation r.  Returns chsph."""
    ra = sig1
    rb = sig2
    if sig2 > sig1:
        ra = sig2
        rb = sig1
    if r >= (ra + rb):
        chsph = 1.0 / r
    elif r <= (ra - rb):
        if ra > 0.0:
            chsph = (3.0 * (ra * ra - 0.2 * rb * rb) - r * r) / (2.0 * ra * ra * ra)
        else:
            chsph = r4(1.0e20)
    else:
        aup1 = ra - r
        alo1 = -rb
        aup2 = +rb
        alo2 = ra - r
        chsph = sum1(r, ra, rb, aup1, alo1) + sum2(r, ra, rb, aup2, alo2)
    return chsph


@routine
def sum1(r, a, b, aup, alo):
    """ato_manipulation_rev4.f:sum1 -- first part of the overlapping-sphere integral."""
    cof = real_array(10)  # real*4 cof(10)
    a2 = a * a
    b2 = b * b
    r2 = r * r
    cof[6] = r4(1.0) / r4(6.0)
    cof[5] = 3.0 * r / (5.0)
    cof[4] = -(f32(3) * (a2 - r2) + b2) / 4.0
    cof[3] = -r * (3.0 * (a2 + b2) - r2) / 3.0
    cof[2] = 3.0 * b2 * (a2 - r2) / 2.0
    cof[1] = r * b2 * (3.0 * a2 - r2)
    fac = f32(3) / (8.0 * a2 * a * b2 * b * r)
    sumu = f32(0)
    suml = f32(0)
    for i in range(1, 7):
        ir = 7 - i
        sumu = aup * (sumu + cof[ir])
        suml = alo * (suml + cof[ir])
    sum1 = fac * (sumu - suml)
    return sum1


@routine
def sum2(r, a, b, aup, alo):
    """ato_manipulation_rev4.f:sum2 -- second part of the overlapping-sphere integral."""
    cof = real_array(10)  # real*4 cof(10)
    a2 = a * a  # noqa: F841
    b2 = b * b
    r2 = r * r  # noqa: F841
    cof[3] = -r4(1.0) / (r4(3.0))
    cof[2] = 0.0
    cof[1] = b2
    fac = f32(3) / (f32(4) * r * b2 * b)
    sumu = f32(0)
    suml = f32(0)
    for i in range(1, 4):
        ir = 4 - i
        sumu = aup * (sumu + cof[ir])
        suml = alo * (suml + cof[ir])
    sum2 = fac * (sumu - suml)
    return sum2


# =========================================================================== resetq
@routine
def resetq(im):
    """ato_manipulation_rev4.f:resetq -- resets position of Q-type atoms within molecule im
    relative to nuclei (bond-length weighted mean of the bonded atoms)."""
    # use molecom_mod
    #
    # step through atoms of molecule
    #
    for ia in range(clustercom.naf[im], clustercom.nal[im] + 1):
        nb = int(molecom.nbond[ia])
        #
        # if atom is type q, calculate revised coordinates
        #
        if molecom.qtype[ia] == 1 and nb > 0:
            wtsum = f32(0.)
            xsum = f32(0.0)
            ysum = f32(0.0)
            zsum = f32(0.0)
            for ib in range(1, nb + 1):
                ja = molecom.iabond[ia, ib]
                tb = abs(molecom.tbond[ia, ib])
                wtsum = wtsum + tb
                xsum = xsum + molecom.xa[ja] * tb
                ysum = ysum + molecom.ya[ja] * tb
                zsum = zsum + molecom.za[ja] * tb
            if wtsum > 0.0:
                molecom.xa[ia] = xsum / wtsum
                molecom.ya[ia] = ysum / wtsum
                molecom.za[ia] = zsum / wtsum
    return


# =========================================================================== truncation functions
@routine
def truncate(xx, rminpt, rmaxpt):
    """ato_manipulation_rev4.f:truncate -- potential truncation function.

    Potential is set to zero beyond rmaxpt. rminpt will represent either the decay of
    the truncation function or the point at which the truncation falls to 0.5.
    """
    # real truncate
    pi = PI4ATAN
    decay = (rmaxpt - rminpt)
    if xx >= rmaxpt:
        truncate = f32(0.0)
        return truncate
    elif xx <= rminpt:
        truncate = f32(1.0)
        return truncate
    else:
        rat = (xx - rminpt) * pi / decay
        truncate = 0.5 * (1.0 + cosf(rat))
        return truncate


@routine
def truncatec(xx, rmaxpt):
    """ato_manipulation_rev4.f:truncatec -- charged clouds truncation function according
    to Hummer et. al., JPCM, vol 6 p A141 (1994)."""
    # real truncatec
    if xx >= rmaxpt:
        truncatec = f32(0.0)
        return truncatec
    else:
        rat = xx / rmaxpt
        truncatec = pow4(1.0 - rat)
        truncatec = truncatec * (1 + rat * (1.6 + rat * 0.4))
        return truncatec


@routine
def truncatecweeks(xx, widthptc):
    """ato_manipulation_rev4.f:truncatecweeks -- erfc truncation of the charged clouds
    (Hummer et. al., JPCM, vol 6 p A141 (1994))."""
    # real truncatecweeks,xx,widthptc
    truncatecweeks = f32(0.0)
    if xx < 5.0 * widthptc:
        rat = xx / widthptc
        truncatecweeks = erfcf(rat)
    return truncatecweeks


# =========================================================================== makemole
# FORMAT statements of makemole.  A trailing '$' (no advance) is dropped from the format
# and written as advance=False.
_F100_makemole = "(a)"
_F499 = "(/'makemole> No ',a,' file specified! - routine will exit'/)"
_F445 = "(/'makemole> Unable to read files of type: ',a)"
_F444 = "(/'makemole> Could not open makemole file : ',a)"
_F447 = "(/'makemole> First line must have two numbers ','for file type: ',a)"
_F446 = "(/'makemole> File type: ',a,' Atom label add: ',i2)"
_F448 = "(/'makemole> Attempted to create more than ',i5,' atoms')"
_F111 = ("(/'makemole> Unrecognised format for mapped atom ',a"
         ",' at row ',i3,' and column ',i3//' - check the .mol file')")
_F112 = ("(/'makemole> Label referred to by mapped atom',a"
         ",' not found at row ',i3,' and column ',i3//' - check the .mol file')")
_F101 = "(/'makemole> Give the bond',' length for atom pair ',a3,1x,a3,': ')"  # ,$)
_F102 = "(/'makemole> Give the angle subtended by the triplet ',a3,1x,a3,1x,a3,': ')"  # ,$)
_F116 = "(/'makemole> Give the temperature',' (K) of the simulation: ')"
_F117 = "(/'makemole> Give the atomic number',' density of the simulation [ats/A**3]: ')"
_F115 = ("(/'makemole> Type epsilon, sigma, mass',', charge and chemical symbol '"
         ",'for component ',a3//'makemole> ? ')")
_F119 = "(i3,1x)"
_F103 = "(i1)"
_F108 = "(i1)"
_F109 = "(i2)"
_F110 = "(i3)"
_F104 = "(a,2(1x,a3),1x,f10.5)"
_F105 = "(a,3(1x,a3),1x,f10.5)"
_F121 = "(a,4(1x,i3),1x,f10.5)"
_F107 = "(a,2(1x,i3))"
_F106 = "(a,2(1x,i3),1x,f10.5)"
_F118 = "(a,2(1x,a3))"
_F113 = "(a,1x,a3,4(1x,f12.6),1x,a3,1x,e13.6,2(1x,f12.6),1x,e13.6)"
_F114 = "(a,1x,e13.6)"
_F120 = "(a,1x,f10.5,1x,f10.5)"
_F122 = "(a,6(1x,f10.5))"


@routine
def makemole():
    """ato_manipulation_rev4.f:makemole -- makes an .ato file from a .mol template file.

    7-10-2006 Modified to allow dihedral angles to be specified.
    07-01-2014 Modified to treat specified dihedral angles as rotational groups, but with
    the specified dihedral angle written into the .ato file, to be used by EPSR (and
    fmole possibly?).
    22-10-2006 Modified to allow vibrational temperature to be specified and output.
    """
    # use reallocation_routines
    # use molecom_mod
    #
    # INTERFACE get_rotgrp(ir1,ir2,labrott,nrott) / get_whole_rotgrp(im,ir1,ir2,labrott,nrott)
    #    (labrott: integer, dimension(:), allocatable, intent(out)) -- external routines of
    #    ato_manipulation_rev4_part2; they return (labrott, nrott).
    mstarlabel = 2000  # parameter (mstarlabel=2000)
    nstsize = commands.nstsize
    nworddim = commands.nworddim
    near = epsr_dimension.near
    mcomp = epsr_dimension.mcomp
    ncf = commands.ncf
    ncl = commands.ncl

    # character*256 readfile,first,fnameatm,ruler
    readfile = ""
    first = ""
    fnameatm = ""
    ruler = ""
    # character*20 ext,filetype
    ext = ""
    filetype = ""
    # character*1 iaccept
    iaccept = ""

    # integer ilenreadfile,ilenext,ilenfiletype
    # character*20 a1,offset(mstarlabel)
    a1 = ""
    offset = char_array(mstarlabel)
    # character*3 label(mstarlabel),alabelsave(mstarlabel)
    label = char_array(mstarlabel)
    alabelsave = char_array(mstarlabel)
    # character*3 la1,la2,la3,la4
    la1 = ""
    la2 = ""
    la3 = ""
    # labels for bonds defined in input stream
    lbond1 = char_array(mstarlabel)
    lbond2 = char_array(mstarlabel)
    # labels for angles defined in input stream - langle0 is the centre of a triplet
    langle1 = char_array(mstarlabel)
    langle2 = char_array(mstarlabel)
    langle3 = char_array(mstarlabel)
    # old labels to be changed
    labelold = char_array(mstarlabel)
    # labels for potential parameters and chemical symbol
    plabel = char_array(mcomp)
    cread = char_array(mcomp)
    # labels for charge radii
    qrlabel = char_array(mcomp)
    # labels to be changed on writing the .ato file
    labelnew = char_array(mstarlabel)

    # row and column coordinates of stars (bonds) and labels (atoms)
    irstar = int_array(mstarlabel)
    icstar = int_array(mstarlabel)
    irlabel = int_array(mstarlabel)
    iclabel = int_array(mstarlabel)
    # atom numbers and number of atoms for a given label
    ianumber = int_array(9, mstarlabel)
    ilabel = int_array(mstarlabel)
    # atom numbers for extra distances to be used to constrain particular atom pairs
    iaextra1 = int_array(mstarlabel)
    iaextra2 = int_array(mstarlabel)
    # atom numbers to define the axis of rotational groups
    iarot1 = int_array(mstarlabel)
    iarot2 = int_array(mstarlabel)
    iawrot1 = int_array(mstarlabel)
    iawrot2 = int_array(mstarlabel)
    iairot1 = int_array(mstarlabel)
    iairot2 = int_array(mstarlabel)
    # atom numbers to define dihedral angles
    idread = int_array(4, mstarlabel)
    # temporary rotational groups
    nrott1 = 0
    nrott2 = 0
    labrott1 = None  # integer, dimension(:), allocatable
    labrott2 = None

    # real*4 bread(mstarlabel),aread(mstarlabel),eread(mstarlabel)
    bread = real_array(mstarlabel)
    aread = real_array(mstarlabel)
    eread = real_array(mstarlabel)
    # read values for epsilon, sigma, mass, Coulomb charge, and charge radius
    pread = real_array(mcomp)
    sread = real_array(mcomp)
    mread = real_array(mcomp)
    qread = real_array(mcomp)
    qrread = real_array(mcomp)
    natomsphread = real_array(mcomp)
    radminsphread = real_array(mcomp)
    radmaxsphread = real_array(mcomp)
    cgmolfacread = real_array(mcomp)
    # dihedral angles
    dread = real_array(mstarlabel)
    # real :: simtkread,tempbread,angtempread,dihtempread
    simtkread = f32(0.0)
    tempbread = f32(0.0)
    tempangread = f32(0.0)
    tempdihread = f32(0.0)
    # real :: pseudofac !No. of real atoms to be used for each pseudo atom. Normally should
    # be 1.0, but increase for coarse graining.

    # logical unitcelldefined,found

    # external cosd,sind  -> EPSR's own auxiliary_routines.cosd/sind

    # NB: locals that the Fortran may read before setting (stack garbage); 0 is used.
    # nspace is only read from the first line of a '.mol' file; nrow, nrowmax, ncolumnmax,
    # nstar only set for '.mol'; nmap is never initialised (only incremented, never used);
    # nrowoffset/ncoloffset/irmap/icmap keep their value if an internal read fails.
    nspace = 0
    nrow = 0
    nrowmax = 0  # noqa: F841
    ncolumnmax = 0
    nstar = 0
    nmap = 0
    nrowoffset = 0
    ncoloffset = 0
    irmap = 0
    icmap = 0
    nread = 0
    vol = f32(0.0)
    #
    # initialise random number generator
    #
    #     call initran1(10)
    #
    # generate the name of the .mol file
    #
    ext = ".mol"
    ilenext = 4
    if commands.nwords < 2:
        commands.nwords = 0
        while commands.nwords <= 0:
            iaccept = fstr(file_and_directory_commands.search_for_files(ext, ilenext, iaccept), 1)
            commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 5, nstsize, nworddim)
            if commands.nwords <= 0:
                units.write(6, fmt(_F499, substr(ext, 1, ilenext)))
                return
        readfile = fstr(substr(commands.sometext, ncf[1], ncl[commands.nwords]), 256)
        ilenreadfile = int(ncl[1]) - int(ncf[1]) + 1
    else:
        readfile = fstr(substr(commands.sometext, ncf[2], ncl[commands.nwords]), 256)
        ilenreadfile = int(ncl[2]) - int(ncf[2]) + 1

    # Get the extension that was supplied - .mol files are treated differently from .gmol files

    ilenuptodot = auxiliary_routines.indexlastdot(readfile) - 1
    if ilenuptodot > 0:

        ext = fstr(substr(readfile, ilenuptodot + 1, ilenreadfile), 20)
        ilenext = ilenreadfile - ilenuptodot
        if substr(ext, 1, ilenext) != ".mol":

            units.write(6, fmt(_F445, substr(ext, 1, ilenext)))
            return

        readfile = fstr(substr(readfile, 1, ilenuptodot), 256)
        ilenreadfile = ilenuptodot

    symbols.fnameato = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                            + substr(readfile, 1, ilenreadfile) + substr(ext, 1, ilenext), 256)
    symbols.ilenfnameato = directory_name.nchardir + ilenreadfile + ilenext

    #
    # Attempt to open .mol file
    #
    try:
        units.open(10, symbols.fnameato, status="old")
        ierr = 0
    except FortranIOError as e:
        ierr = e.iostat

    # Any kind of problems and the program returns to calling routine
    # and does nothing

    if ierr != 0:

        units.write(6, fmt(_F444, substr(symbols.fnameato, 1, symbols.ilenfnameato)))
        return

    # Need to determine what file type this is

    commands.sometext, ierr = _read_a(10, commands.sometext, nstsize)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 200, nstsize, nworddim)

    if findex(substr(commands.sometext, ncf[1], ncl[1]), ".") == 1:

        filetype = fstr(substr(commands.sometext, ncf[1], ncl[1]), 20)
        ilenfiletype = int(ncl[1]) - int(ncf[1]) + 1
        if commands.nwords > 1:

            iatomlabeladd = read_internal(substr(commands.sometext, ncf[2], ncl[2])).int(0)

        else:

            units.write(6, fmt(_F447, substr(filetype, 1, ilenfiletype)))
            return

    else:

        # Assume original .mol file type

        filetype = ".mol"
        ilenfiletype = 4
        iatomlabeladd = 0
        #
        # first number represents the number of characters in a column (including the delimiting space
        # between columns)
        #
        nspace = read_internal(substr(commands.sometext, ncf[1], ncl[1])).int(nspace)

    #
    # first stage is to read the molecular template file and generate
    # the atoms and bonds reference integers
    #
    if substr(filetype, 1, ilenfiletype) == ".mol":
        #
        # read a ruler line - this can be used to show where the columns occur
        #
        ruler, ierr = _read_a(10, ruler, 256)
        nrow = 0
        nrowmax = 0  # noqa: F841
        ncolumnmax = 0
        nstar = 0

    units.write(6, fmt(_F446, substr(filetype, 1, ilenfiletype), iatomlabeladd))

    nlabel = 0
    ierr = 0
    # set the number of atom types to zero at the outset
    molecom.ntype = 0
    # no. of bond values read on input
    nbread = 0
    # no. of bond angle values read on input
    naread = 0
    # no. of extra distances read on input
    neread = 0
    # no. of rotational groups read on input
    nrread = 0
    nwrread = 0
    nirread = 0
    # no. of potential values read on input
    npread = 0
    # no. of charge radius values read on input
    nqrread = 0
    # initial temperature
    molecom.simtk = f32(300)
    simtkread = f32(0.0)
    # initial vibrational temperature
    molecom.tempb = r4(65000.0)
    tempbread = f32(0.0)
    nvread = 0  # noqa: F841
    # initial bond angle temperature (relative)
    molecom.tempangreq = r4(1.0)
    tempangread = f32(0.0)
    nangread = 0  # noqa: F841
    # initial dihedral temperature
    molecom.tempdihreq = r4(10.0)
    tempdihread = f32(0.0)
    ndihread = 0  # noqa: F841
    # initial atom number density
    rho = f32(0.0)
    # no. of labels to change before writing the .ato file
    nchread = 0
    # no. of dihedral angles on input
    ndread = 0
    # assumed values of ecore and dcore unless changed by the input.
    # These are used in the molecule generation program fmole.
    molecom.ecore = f32(0.0)
    molecom.dcore = f32(0.0)
    unitcelldefined = False
    # Default is: assume 1 real atom for each pseudo atom
    pseudofac = r4(1.0)
    ia = 0
    while ierr == 0:
        #
        # read a line of text
        #
        commands.sometext = ""
        commands.sometext, ierr = _read_a(10, commands.sometext, nstsize)
        #
        # parse the text into words. Any gaps of more than 200 spaces will
        # signal the end of the line
        #
        commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 200, nstsize, nworddim)
        nwords = commands.nwords
        sometext = commands.sometext

        def word(k):  # sometext(ncf(k):ncl(k))
            return substr(sometext, ncf[k], ncl[k])

        if nwords > 0:
            #
            # check that the first word is not a control word
            #
            first = fstr(word(1), 256)
            ilenfirst = int(ncl[1]) - int(ncf[1]) + 1
            firstw = substr(first, 1, ilenfirst)

            if firstw == "atom":

                if nwords > 6 and naread < mstarlabel:

                    if ia >= mstarlabel:

                        units.write(6, fmt(_F448, mstarlabel))
                        return

                    ia = ia + 1
                    # Ignore the atom number from ghemical
                    isref = 3
                    molecom.alabel[ia] = fstr(word(isref), 3)

                    # Store this label in case it is changed

                    alabelsave[ia] = molecom.alabel[ia]
                    #
                    # set the x and y coordinates according to position in the map
                    #
                    isref = isref + 1
                    molecom.xa[ia] = read_internal(word(isref)).real(molecom.xa[ia])
                    isref = isref + 1
                    molecom.ya[ia] = read_internal(word(isref)).real(molecom.ya[ia])
                    isref = isref + 1
                    molecom.za[ia] = read_internal(word(isref)).real(molecom.za[ia])
                    isref = isref + 1
                    molecom.nbond[ia] = read_internal(word(isref)).int(molecom.nbond[ia])
                    ib = 0
                    while ib < molecom.nbond[ia] and ib < near and isref < nwords:

                        ib = ib + 1
                        isref = isref + 1
                        molecom.iabond[ia, ib] = read_internal(word(isref)).int(molecom.iabond[ia, ib])
                        molecom.iabond[ia, ib] = molecom.iabond[ia, ib] + iatomlabeladd

                    molecom.nbond[ia] = ib

            elif firstw == "bond":
                #
                # assume the next two words are atom type labels, followed by the corresponding
                # bond length in A.
                #
                if nwords >= 4 and nbread < mstarlabel:
                    nbread = nbread + 1
                    lbond1[nbread] = fstr(word(2), 3)
                    lbond2[nbread] = fstr(word(3), 3)
                    bread[nbread] = read_internal(word(4)).real(bread[nbread])

            elif firstw == "angle":
                #
                # assume the next three words are atom type labels,
                # with the middle one being the one about which the angle subtends
                # followed by the corresponding angle in degrees.
                #
                if nwords >= 5 and naread < mstarlabel:
                    naread = naread + 1
                    langle1[naread] = fstr(word(2), 3)
                    langle2[naread] = fstr(word(3), 3)
                    langle3[naread] = fstr(word(4), 3)
                    aread[naread] = read_internal(word(5)).real(aread[naread])

            elif firstw == "extra":
                #
                # Additional distances that may needed to keep the molecule intact
                # assume the next two words are atom type labels, followed by the corresponding
                # expected separation in A.
                #
                if nwords >= 4 and neread < mstarlabel:
                    neread = neread + 1
                    iaextra1[neread] = read_internal(word(2)).int(iaextra1[neread])
                    iaextra2[neread] = read_internal(word(3)).int(iaextra2[neread])
                    eread[neread] = read_internal(word(4)).real(eread[neread])
                    iaextra1[neread] = iaextra1[neread] + iatomlabeladd
                    iaextra2[neread] = iaextra2[neread] + iatomlabeladd

            elif firstw == "rot":
                #
                # Read two atom numbers to define the axis for rotational groups
                #
                if nwords >= 3 and nrread < mstarlabel:
                    nrread = nrread + 1
                    iarot1[nrread] = read_internal(word(2)).int(iarot1[nrread])
                    iarot2[nrread] = read_internal(word(3)).int(iarot2[nrread])
                    iarot1[nrread] = iarot1[nrread] + iatomlabeladd
                    iarot2[nrread] = iarot2[nrread] + iatomlabeladd

            elif firstw == "wholerot":
                #
                # Read two atom numbers to define the axis for whole rotational groups
                # This is where a whole molecule is rotated about an axis defined by
                # two atoms within the molecule
                #
                if nwords >= 3 and nwrread < mstarlabel:
                    nwrread = nwrread + 1
                    iawrot1[nwrread] = read_internal(word(2)).int(iawrot1[nwrread])
                    iawrot2[nwrread] = read_internal(word(3)).int(iawrot2[nwrread])
                    iawrot1[nwrread] = iawrot1[nwrread] + iatomlabeladd
                    iawrot2[nwrread] = iawrot2[nwrread] + iatomlabeladd

            elif firstw == "introt":
                #
                # Read two atom numbers to define the axis for rotational groups
                #
                if nwords >= 3 and nirread < mstarlabel:
                    nirread = nirread + 1
                    iairot1[nirread] = read_internal(word(2)).int(iairot1[nirread])
                    iairot2[nirread] = read_internal(word(3)).int(iairot2[nirread])
                    iairot1[nirread] = iairot1[nirread] + iatomlabeladd
                    iairot2[nirread] = iairot2[nirread] + iatomlabeladd

            elif firstw == "dihedral":
                #
                # Read four atom numbers to define the dihedral angle plus the dihedral angle.
                # Atoms 2 and 3 define the axis about which the dihedral angle is calculated.
                #
                if nwords >= 6 and ndread < mstarlabel:
                    ndread = ndread + 1
                    for j in range(1, 5):
                        idread[j, ndread] = read_internal(word(j + 1)).int(idread[j, ndread])
                        idread[j, ndread] = idread[j, ndread] + iatomlabeladd

                    dread[ndread] = read_internal(word(6)).real(dread[ndread])

            elif firstw == "changelabel":
                #
                # Change atom labels: read old atom label and new atom label
                # these are only applied after the atom distances have been defined
                #
                if nwords > 2:
                    nchread = nchread + 1
                    labelold[nchread] = fstr(word(2), 3)
                    labelnew[nchread] = fstr(word(3), 3)

            elif firstw == "potential":
                #
                # Read epsilon, sigma, mass, charge and chemical symbol of specified label
                #
                if nwords >= 7 and npread < mcomp:
                    npread = npread + 1
                    plabel[npread] = fstr(word(2), 3)
                    pread[npread] = read_internal(word(3)).real(pread[npread])
                    sread[npread] = read_internal(word(4)).real(sread[npread])
                    mread[npread] = read_internal(word(5)).real(mread[npread])
                    qread[npread] = read_internal(word(6)).real(qread[npread])
                    cread[npread] = fstr(word(7), 3)
                    if nwords > 9:
                        natomsphread[npread] = read_internal(word(8)).real(natomsphread[npread])
                        radminsphread[npread] = read_internal(word(9)).real(radminsphread[npread])
                        radmaxsphread[npread] = read_internal(word(10)).real(radmaxsphread[npread])
                    else:
                        natomsphread[npread] = 1.0
                        radminsphread[npread] = 0.0
                        radmaxsphread[npread] = 0.5 * sread[npread]
                    if nwords > 10:
                        cgmolfacread[npread] = read_internal(word(11)).real(cgmolfacread[npread])
                    else:
                        cgmolfacread[npread] = 1.0

            # Read charge radius for specified atom type. If not specified this will
            # be set to zero

            elif firstw == "qradius":

                if nwords >= 3 and nqrread < mcomp:

                    nqrread = nqrread + 1
                    qrlabel[nqrread] = fstr(word(2), 3)
                    qrread[nqrread] = read_internal(word(3)).real(qrread[nqrread])
            #
            # set temperature of simulation
            #
            elif firstw == "temperature":

                if nwords > 1:
                    simtkread, ierr = _iread(word(2), "real", simtkread)
                    if ierr == 0 and simtkread > 0.0:
                        molecom.simtk = simtkread
                    else:
                        ierr = 0

            elif firstw == "vibtemp":

                if nwords > 1:
                    tempbread, ierr = _iread(word(2), "real", tempbread)
                    if ierr == 0 and tempbread > 0.0:
                        molecom.tempb = tempbread
                        nvread = int(ncl[2]) - int(ncf[2]) + 1  # noqa: F841
                    else:
                        ierr = 0

            elif firstw == "angtemp":

                if nwords > 1:
                    tempangread, ierr = _iread(word(2), "real", tempangread)
                    if ierr == 0 and tempangread > 0.0:
                        molecom.tempangreq = tempangread
                        nangread = int(ncl[2]) - int(ncf[2]) + 1  # noqa: F841
                    else:
                        ierr = 0

            elif firstw == "dihtemp":

                if nwords > 1:
                    tempdihread, ierr = _iread(word(2), "real", tempdihread)
                    if ierr == 0 and tempdihread > 0.0:
                        molecom.tempdihreq = tempdihread
                        ndihread = int(ncl[2]) - int(ncf[2]) + 1  # noqa: F841
                    else:
                        ierr = 0

            elif firstw == "density":

                if nwords > 1:
                    rho = read_internal(word(2)).real(rho)
            #
            # read the potential parameters used to find the molecule
            #
            elif firstw == "ecoredcore":

                if nwords > 1:
                    molecom.ecore = read_internal(word(2)).real(molecom.ecore)
                if nwords > 2:
                    molecom.dcore = read_internal(word(3)).real(molecom.dcore)
            #
            # read the unit cell parameters used to set the size and shape of the box
            #
            elif firstw == "unitcell":

                if nwords > 6:
                    unitcelldefined = True
                    r = read_internal(substr(sometext, ncf[2], ncl[7]))
                    lattice.alat = r.real(lattice.alat)
                    lattice.blat = r.real(lattice.blat)
                    lattice.clat = r.real(lattice.clat)
                    lattice.phib = r.real(lattice.phib)
                    lattice.thetac = r.real(lattice.thetac)
                    lattice.phic = r.real(lattice.phic)
                else:
                    unitcelldefined = False
            #
            # read the unit cell parameters used to set the size and shape of the box
            #
            elif firstw == "pseudofac":

                if nwords > 0:
                    pseudofac = read_internal(word(2)).real(pseudofac)
                    if pseudofac <= 0.0:
                        pseudofac = r4(1.0)  # Ignore crazy values
            #
            # otherwise treat the words as atom labels or bonds. The first word is the row
            # number and so will be ignored on input
            #
            elif nwords > 1:
                nrow = nrow + 1
                iw = 1
                while iw < nwords:
                    iw = iw + 1
                    #
                    # for each word read in, first find out if the first character is a number
                    #
                    nread, ierrread = _iread(substr(sometext, ncf[iw], ncf[iw]), "int", nread)
                    if ierrread == 0:
                        a1 = fstr(substr(sometext, ncf[iw] + 1, ncl[iw]), 20)
                    else:
                        a1 = fstr(word(iw), 20)
                        nread = 1
                    #
                    # now determine whether this character represents a bond or an atom
                    #
                    # NB: nspace is undefined in the Fortran unless the file is a '.mol' file;
                    # with the 0 used here idiv raises ZeroDivisionError (SIGFPE in the Fortran).
                    ncolumn = idiv(int(ncf[iw]) - 4 + nspace - 1, nspace)
                    if ncolumn > ncolumnmax:
                        ncolumnmax = ncolumn
                    if a1 == "*":
                        nstar = nstar + 1
                        irstar[nstar] = nrow
                        icstar[nstar] = ncolumn
                    else:
                        #
                        # set up coordinates of label and no of atoms at this label
                        #
                        nlabel = nlabel + 1
                        irlabel[nlabel] = nrow
                        iclabel[nlabel] = ncolumn
                        ilabel[nlabel] = nread
                        #
                        # column and row number of current label
                        #
                        il1 = nlabel
                        icl1 = int(iclabel[il1])
                        irl1 = int(irlabel[il1])
                        #
                        # if the first three characters of the label are 'map', this label
                        # represents a mapping to a previous label, given by the row and column offsets
                        # from the present position
                        #
                        if substr(a1, 1, 3) == "map":
                            #
                            # read the row and column offsets
                            #
                            label[nlabel] = fstr(substr(a1, 1, 3), 3)
                            nchara1 = findex_padded(a1, " ", 20) - 1
                            offset[nlabel] = fstr(substr(a1, 4, nchara1), 20)
                            nmap = nmap + 1
                            ncharoffset = nchara1 - 3
                            #
                            # test to see if this is an absolute reference
                            #
                            ilenopenbra = findex(offset[nlabel], "(")
                            ilencomma = findex(offset[nlabel], ",")
                            ilenclosebra = findex(offset[nlabel], ")")
                            ilenread1 = ilenopenbra + ilencomma + ilenclosebra
                            #
                            # to be a valid reference the sum of the positions of these three characters must
                            # be at least 1 + 3 + 5 = 9, otherwise there is no room for numbers
                            #
                            if ilenread1 > 8:
                                #
                                # absolute reference
                                #
                                irmap, ierr = _iread(substr(offset[nlabel], ilenopenbra + 1, ilencomma - 1),
                                                     "int", irmap)
                                icmap, ierr = _iread(substr(offset[nlabel], ilencomma + 1, ilenclosebra - 1),
                                                     "int", icmap)
                                if ierr != 0 or irmap <= 0 or icmap <= 0:
                                    units.write(6, fmt(_F111, substr(a1, 1, nchara1), irl1, icl1))
                                    return
                            else:
                                #
                                # alternatively try see if this is relative reference
                                #
                                ilenplus = findex(substr(offset[nlabel], 2, ncharoffset), "+")
                                ilenminus = findex(substr(offset[nlabel], 2, ncharoffset), "-")
                                ilenread = ilenplus + ilenminus
                                if ilenread > 1:
                                    ilenread = ilenread + 1
                                    #
                                    # read the offsets
                                    #
                                    nrowoffset, ierr = _iread(substr(offset[nlabel], 1, ilenread - 1),
                                                              "int", nrowoffset)
                                    ncoloffset, ierr = _iread(substr(offset[nlabel], ilenread, ncharoffset),
                                                              "int", ncoloffset)
                                    if ierr != 0:
                                        units.write(6, fmt(_F111, substr(a1, 1, nchara1), irl1, icl1))
                                        return
                                    #     write(6,*) irl1,icl1,nrowoffset,ncoloffset
                                    #
                                    # calculate the row and column of the label this map refers to
                                    #
                                    irmap = irl1 + nrowoffset
                                    icmap = icl1 + ncoloffset
                                else:
                                    units.write(6, fmt(_F111, substr(a1, 1, nchara1), irl1, icl1))
                                    return
                            #
                            # find the previous label corresponding to this label - if not found exit with an
                            # error message.
                            #
                            il = 0
                            ifound = 0
                            while il < nlabel - 1 and ifound == 0:
                                il = il + 1
                                if irlabel[il] == irmap and iclabel[il] == icmap:
                                    ifound = il
                            if ifound == 0:
                                units.write(6, fmt(_F112, substr(a1, 1, nchara1), irmap, icmap))
                                return
                            else:
                                ilabel[nlabel] = ilabel[ifound]
                                for ia1 in range(1, int(ilabel[nlabel]) + 1):
                                    ianumber[ia1, nlabel] = ianumber[ia1, ifound]
                        else:
                            #
                            # otherwise generate the atom numbers for this new atom and zero the bond array
                            #
                            label[nlabel] = fstr(a1, 3)
                            for ia1 in range(1, int(ilabel[nlabel]) + 1):
                                ia = ia + 1
                                ianumber[ia1, nlabel] = ia
                                molecom.alabel[ia] = label[nlabel]
                                #
                                # set the x and y coordinates according to position in the map
                                #
                                molecom.xa[ia] = 0.0
                                molecom.ya[ia] = 0.0
                                molecom.za[ia] = 0.0
                                #                    call ftype(ia)
                                #                    if(qtype(ia).ne.1) natomr=natomr+1
                                molecom.nbond[ia] = 0
                        #
                        # step through the existing stars to find out which occur as neighbours
                        # in a previous column on the same row, or in the previous row
                        #
                        is_ = 0
                        while is_ < nstar:
                            is_ = is_ + 1
                            irs = int(irstar[is_])
                            ics = int(icstar[is_])
                            irdif1 = int(irlabel[il1]) - irs
                            icdif1 = int(iclabel[il1]) - ics
                            if (irdif1 == 1 and abs(icdif1) < 2) or (irdif1 == 0 and icdif1 == 1):
                                #
                                # o.k. found a star in the vicinity - look for another label adjecent to this star (not including
                                # the present) one in the previous two rows and or columns
                                #
                                il2 = 1
                                while il2 < nlabel:
                                    irdif2 = abs(irs - int(irlabel[il2]))
                                    icdif2 = abs(ics - int(iclabel[il2]))
                                    if irdif2 < 2 and icdif2 < 2:
                                        #
                                        # atom found - set up the bond.
                                        #
                                        # step through all the possible pairs of atoms and identify bonds
                                        #
                                        #     write(6,*) il1,il2
                                        for ia1 in range(1, int(ilabel[il1]) + 1):
                                            ia1ref = int(ianumber[ia1, il1])
                                            for ia2 in range(1, int(ilabel[il2]) + 1):
                                                ia2ref = int(ianumber[ia2, il2])
                                                #
                                                # check that ia2ref does not already exist in the list of neighbours of ia1ref
                                                #
                                                ib = 0
                                                ifound = 0
                                                while ib < molecom.nbond[ia1ref] and ifound == 0:
                                                    ib = ib + 1
                                                    if molecom.iabond[ia1ref, ib] == ia2ref:
                                                        ifound = ia2ref
                                                if ifound == 0:
                                                    molecom.nbond[ia1ref] = molecom.nbond[ia1ref] + 1
                                                    molecom.nbond[ia2ref] = molecom.nbond[ia2ref] + 1
                                                    molecom.iabond[ia1ref, molecom.nbond[ia1ref]] = ia2ref
                                                    molecom.iabond[ia2ref, molecom.nbond[ia2ref]] = ia1ref
                                    il2 = il2 + 1
    units.close(10)
    molecom.natom = ia
    natom = ia
    #
    # write out findings so far
    #
    #      write(6,*) nrow,ncolumnmax,nstar,nlabel,ntype
    #      do it=1,ntype
    #           write(6,*) it,' ',type(it)
    #      end do
    #      do ia=1,natom
    #         write(6,*) ia
    #         write(6,*) alabel(ia),nbond(ia),(iabond(ia,ib),ib=1,nbond(ia))
    #      end do
    #
    # set the bond lengths
    #
    for ia in range(1, natom + 1):
        la1 = molecom.alabel[ia]
        ib = 0
        while ib < molecom.nbond[ia]:
            ib = ib + 1
            la2 = molecom.alabel[molecom.iabond[ia, ib]]
            #
            # search through the bonds read in (if any) and assign the appropriate bond length
            #
            ibread = 0
            ifound = 0
            while ibread < nbread and ifound == 0:
                ibread = ibread + 1
                if ((la1 == lbond1[ibread] and la2 == lbond2[ibread])
                        or (la1 == lbond2[ibread] and la2 == lbond1[ibread])):
                    ifound = ibread
            ierr = 1
            while ifound == 0 and ierr != 0:
                units.write(6, fmt(_F101, Ch(la1, 3), Ch(la2, 3)), advance=False)
                nbread = nbread + 1
                ifound = nbread
                lbond1[nbread] = la1
                lbond2[nbread] = la2
                bread[nbread], ierr = _read5_real(bread[nbread])
            molecom.tbond[ia, ib] = bread[ifound]
    clustercom.naf[1] = 1
    clustercom.nal[1] = natom
    #
    # save the number of primary bonds for this atom
    #
    for ia in range(1, natom + 1):
        molecom.nbondp[ia] = molecom.nbond[ia]
    #
    # write out findings so far
    #
    #     write(6,*) nrow,ncolumnmax,nstar,nlabel
    #     do ia=1,natom
    #        write(6,*) ia,alabel(ia),nbond(ia),(iabond(ia,ib),tbond(ia,ib),ib=1,nbond(ia))
    #     end do
    #
    # now step through the atoms and use the angles defined to create extra bonding distances
    #
    for ia2 in range(1, natom + 1):
        la2 = molecom.alabel[ia2]
        #
        # must be at least two bonds for an angle to be definable
        #
        ib1 = 0
        while ib1 < (molecom.nbondp[ia2] - 1):
            ib1 = ib1 + 1
            ia1 = int(molecom.iabond[ia2, ib1])
            la1 = molecom.alabel[ia1]
            ib2 = ib1
            while ib2 < molecom.nbondp[ia2]:
                ib2 = ib2 + 1
                ia3 = int(molecom.iabond[ia2, ib2])
                la3 = molecom.alabel[ia3]
                #
                # before we set up this bond check that a bond between these two atoms
                # does not already exist
                #
                isearch = 0
                iexist = 0
                while isearch < molecom.nbond[ia1] and iexist == 0:
                    isearch = isearch + 1
                    if ia3 == molecom.iabond[ia1, isearch]:
                        iexist = ia3
                if iexist == 0:
                    bond21 = abs(molecom.tbond[ia2, ib1])
                    bond23 = abs(molecom.tbond[ia2, ib2])
                    #
                    # search for the angle between the three atoms
                    #
                    #
                    # search through the angles read in (if any) and assign the appropriate angle
                    # to this triplet
                    #
                    iaread = 0
                    ifound = 0
                    while iaread < naread and ifound == 0:
                        iaread = iaread + 1
                        if ((la1 == langle1[iaread] and la2 == langle2[iaread] and la3 == langle3[iaread])
                                or (la1 == langle3[iaread] and la2 == langle2[iaread]
                                    and la3 == langle1[iaread])):
                            ifound = iaread
                    ierr = 1
                    while ifound == 0 and ierr != 0:
                        units.write(6, fmt(_F102, Ch(la1, 3), Ch(la2, 3), Ch(la3, 3)), advance=False)
                        naread = naread + 1
                        ifound = naread
                        langle1[naread] = la1
                        langle2[naread] = la2
                        langle3[naread] = la3
                        aread[naread], ierr = _read5_real(aread[naread])
                    #
                    # Set up distance between ib1 and ib3
                    #
                    bangle = aread[ifound]
                    b13 = (bond21 * bond21 + bond23 * bond23
                           - f32(2) * bond21 * bond23 * auxiliary_routines.cosd(bangle))
                    b13 = sqrtf(b13)
                    if molecom.nbond[ia1] < near:

                        molecom.nbond[ia1] = molecom.nbond[ia1] + 1
                        molecom.iabond[ia1, molecom.nbond[ia1]] = ia3
                        molecom.tbond[ia1, molecom.nbond[ia1]] = -b13

                    if molecom.nbond[ia3] < near:

                        molecom.nbond[ia3] = molecom.nbond[ia3] + 1
                        molecom.iabond[ia3, molecom.nbond[ia3]] = ia1
                        molecom.tbond[ia3, molecom.nbond[ia3]] = -b13
    #
    # save the number of secondary bonds
    #
    for ia in range(1, natom + 1):
        molecom.nbonds[ia] = molecom.nbond[ia]

    # Now setup any dihedral angles requested. These will be defined by means
    # of distances between atoms

    id = 0
    ir = 0
    units.write(6, lw(ndread, naread))
    while id < ndread:

        id = id + 1

        id1 = int(idread[1, id])
        id2 = int(idread[2, id])
        id3 = int(idread[3, id])
        id4 = int(idread[4, id])

        # First need to check that id4 is not already listed as a neighbour of id1.
        # If it is then nothing is done

        ifound = 0
        ib = 0
        while ib < molecom.nbond[id1] and ifound == 0:
            ib = ib + 1
            if molecom.iabond[id1, ib] == id4 and molecom.tbond[id1, ib] > 0.0:
                ifound = id4

        if ifound == 0:

            ir1 = id2
            ir2 = id3

            # Need to check whether a dihedral group for this pair of axis atoms has already been defined.

            iref = 0
            ifounddih = 0
            while iref < ir and ifounddih == 0:
                iref = iref + 1
                if (molecom_mod.rotgrp[iref] == "DIH"
                        and ((ir1 == molecom_mod.labaxi[iref, 1] and ir2 == molecom_mod.labaxi[iref, 2])
                             or (ir2 == molecom_mod.labaxi[iref, 1] and ir1 == molecom_mod.labaxi[iref, 2]))):
                    ifounddih = iref

            # Start a new rotational group:
            # Get the two possible rotational groups into temporary arrays

            if ifounddih == 0:
                labrott1, nrott1 = get_rotgrp(ir1, ir2, labrott1, nrott1)  # noqa: F821 (part 2)
                labrott2, nrott2 = get_rotgrp(ir2, ir1, labrott2, nrott2)  # noqa: F821 (part 2)
                ir = ir + 1
                molecom.ndihedralang[ir] = 1
                molecom.dihedralang[molecom.ndihedralang[ir], ir] = dread[id]
                molecom_mod_lib.testsize_rotgrps(ir)
                molecom_mod.rotgrp[ir] = "DIH"
                # Use the smallest of the two groups
                if nrott1 < nrott2:
                    molecom_mod.labaxi[ir, 1] = ir1
                    molecom_mod.labaxi[ir, 2] = ir2
                    molecom_mod.nrot[ir] = nrott1
                    molecom_mod.rotf[ir] = molecom_mod.nrottot + 1
                    molecom_mod.nrottot = molecom_mod.nrottot + int(molecom_mod.nrot[ir])
                    molecom_mod.rotl[ir] = molecom_mod.nrottot
                    molecom_mod_lib.testsize_labrot()
                    for i in range(1, int(molecom_mod.nrot[ir]) + 1):
                        labref = int(molecom_mod.rotf[ir]) + i - 1
                        molecom_mod.labrot[labref] = labrott1[i]
                    # Remember to preserve the linking to the axis atoms!!
                    molecom.labrotdih[1, molecom.ndihedralang[ir], ir] = id1
                    molecom.labrotdih[2, molecom.ndihedralang[ir], ir] = id4
                else:
                    molecom_mod.labaxi[ir, 1] = ir2
                    molecom_mod.labaxi[ir, 2] = ir1
                    molecom_mod.nrot[ir] = nrott2
                    molecom_mod.rotf[ir] = molecom_mod.nrottot + 1
                    molecom_mod.nrottot = molecom_mod.nrottot + int(molecom_mod.nrot[ir])
                    molecom_mod.rotl[ir] = molecom_mod.nrottot
                    molecom_mod_lib.testsize_labrot()
                    for i in range(1, int(molecom_mod.nrot[ir]) + 1):
                        labref = int(molecom_mod.rotf[ir]) + i - 1
                        molecom_mod.labrot[labref] = labrott2[i]
                    # Remember to preserve the linking to the axis atoms!!
                    molecom.labrotdih[1, molecom.ndihedralang[ir], ir] = id4
                    molecom.labrotdih[2, molecom.ndihedralang[ir], ir] = id1
            else:
                # Add the new dihedral angle to the existing list
                if molecom.ndihedralang[ifounddih] < near:
                    molecom.ndihedralang[ifounddih] = molecom.ndihedralang[ifounddih] + 1
                    molecom.dihedralang[molecom.ndihedralang[ifounddih], ifounddih] = dread[id]
                    if ir1 == molecom_mod.labaxi[iref, 1] and ir2 == molecom_mod.labaxi[iref, 2]:
                        molecom.labrotdih[1, molecom.ndihedralang[ifounddih], ifounddih] = id1
                        molecom.labrotdih[2, molecom.ndihedralang[ifounddih], ifounddih] = id4
                    elif ir2 == molecom_mod.labaxi[iref, 1] and ir1 == molecom_mod.labaxi[iref, 2]:
                        molecom.labrotdih[1, molecom.ndihedralang[ifounddih], ifounddih] = id4
                        molecom.labrotdih[2, molecom.ndihedralang[ifounddih], ifounddih] = id1
    clustercom.ngrot[1] = ir
    # If no dihedral groups are defined, then set ecore and dcore. This will help to stop overlaps in long chain molecules.
    #      if(ir.eq.0.and.ecore.eq.0.0) then
    #         ecore=1.0
    #         dcore=2.0
    #      end if
    #
    # (A long block of commented-out code follows in the Fortran: the old way of turning a
    # dihedral angle into an id1-id4 distance from the bond lengths 1-2, 2-3, 3-4 and the
    # angles 1-2-3, 2-3-4.)

    # finally setup any extra distances between particular atoms that have been specified

    ia1 = 0
    while ia1 < natom - 1:
        ia1 = ia1 + 1
        ia2 = ia1
        while ia2 < natom:
            ia2 = ia2 + 1
            #
            # before we set up this bond check that a bond between these two atoms
            # does not already exist
            #
            isearch = 0
            iexist = 0
            while isearch < molecom.nbond[ia1] and iexist == 0:
                isearch = isearch + 1
                if ia2 == molecom.iabond[ia1, isearch]:
                    iexist = ia2
            if iexist == 0:
                ifound = 0
                ieread = 0
                while ieread < neread and ifound == 0:
                    ieread = ieread + 1
                    if ((ia1 == iaextra1[ieread] and ia2 == iaextra2[ieread])
                            or (ia1 == iaextra2[ieread] and ia2 == iaextra1[ieread])):
                        ifound = ieread
                if ifound != 0:
                    if molecom.nbond[ia1] < near:

                        molecom.nbond[ia1] = molecom.nbond[ia1] + 1
                        molecom.iabond[ia1, molecom.nbond[ia1]] = ia2
                        molecom.tbond[ia1, molecom.nbond[ia1]] = -abs(eread[ifound])

                    if molecom.nbond[ia2] < near:

                        molecom.nbond[ia2] = molecom.nbond[ia2] + 1
                        molecom.iabond[ia2, molecom.nbond[ia2]] = ia1
                        molecom.tbond[ia2, molecom.nbond[ia2]] = -abs(eread[ifound])
    #
    # now set up the rotational groups
    #
    iref = 0
    units.write(6, lw("Getting rotational groups"))
    while iref < nrread:
        iref = iref + 1
        ir1 = int(iarot1[iref])
        ir2 = int(iarot2[iref])

        # Get the two possible rotational groups into temporary arrays

        labrott1, nrott1 = get_rotgrp(ir1, ir2, labrott1, nrott1)  # noqa: F821 (part 2)
        labrott2, nrott2 = get_rotgrp(ir2, ir1, labrott2, nrott2)  # noqa: F821 (part 2)
        ir = ir + 1
        molecom_mod_lib.testsize_rotgrps(ir)
        molecom_mod.rotgrp[ir] = "ROT"
        # Use the smallest of the two groups
        if nrott1 < nrott2:
            molecom_mod.labaxi[ir, 1] = ir1
            molecom_mod.labaxi[ir, 2] = ir2
            molecom_mod.nrot[ir] = nrott1
            molecom_mod.rotf[ir] = molecom_mod.nrottot + 1
            molecom_mod.nrottot = molecom_mod.nrottot + int(molecom_mod.nrot[ir])
            molecom_mod.rotl[ir] = molecom_mod.nrottot
            molecom_mod_lib.testsize_labrot()
            for i in range(1, int(molecom_mod.nrot[ir]) + 1):
                labref = int(molecom_mod.rotf[ir]) + i - 1
                molecom_mod.labrot[labref] = labrott1[i]
        else:
            molecom_mod.labaxi[ir, 1] = ir2
            molecom_mod.labaxi[ir, 2] = ir1
            molecom_mod.nrot[ir] = nrott2
            molecom_mod.rotf[ir] = molecom_mod.nrottot + 1
            molecom_mod.nrottot = molecom_mod.nrottot + int(molecom_mod.nrot[ir])
            molecom_mod.rotl[ir] = molecom_mod.nrottot
            molecom_mod_lib.testsize_labrot()
            for i in range(1, int(molecom_mod.nrot[ir]) + 1):
                labref = int(molecom_mod.rotf[ir]) + i - 1
                molecom_mod.labrot[labref] = labrott2[i]
        molecom.ndihedralang[ir] = 0

    # Set up the whole molecule rotational groups

    iref = 0
    while iref < nwrread:
        iref = iref + 1
        ir1 = int(iawrot1[iref])
        ir2 = int(iawrot2[iref])
        #      write(6,*) ir,iref,ir1,ir2
        # Get all the atoms in the molecule except the two atoms that form the axis.
        labrott2, nrott2 = get_whole_rotgrp(1, ir1, ir2, labrott2, nrott2)  # noqa: F821 (part 2)
        ir = ir + 1
        # NB: the Fortran assigns rotgrp(ir)='ROT' *before* testsize_rotgrps(ir) enlarges the
        # array: the store goes past the end of the old array and is lost, and the new
        # element is left uninitialised (zero bytes with gfortran/glibc, which is what the
        # reference .ato file contains; Rust makemole.rs:517).
        _nold = -1 if molecom_mod.rotgrp is None else molecom_mod.rotgrp.shape[0] - 1
        if ir <= _nold:
            molecom_mod.rotgrp[ir] = "ROT"
        molecom_mod_lib.testsize_rotgrps(ir)
        if ir > _nold:
            molecom_mod.rotgrp[ir] = "\x00\x00\x00"
        molecom_mod.labaxi[ir, 1] = ir1
        molecom_mod.labaxi[ir, 2] = ir2
        molecom_mod.nrot[ir] = nrott2
        molecom_mod.rotf[ir] = molecom_mod.nrottot + 1
        molecom_mod.nrottot = molecom_mod.nrottot + int(molecom_mod.nrot[ir])
        molecom_mod.rotl[ir] = molecom_mod.nrottot
        molecom_mod_lib.testsize_labrot()
        for i in range(1, int(molecom_mod.nrot[ir]) + 1):
            labref = int(molecom_mod.rotf[ir]) + i - 1
            molecom_mod.labrot[labref] = labrott2[i]
        molecom.ndihedralang[ir] = 0

    # Set up the internal rotational groups: this is between two atoms which are not bonded to each other and
    # aims to rotate all the atoms between them. Note that if
    # axis atoms ir1 and ir2 have more than one bond each to the group, then the group cannot be rotated, so should
    # be discarded.

    iref = 0
    while iref < nirread:
        iref = iref + 1
        ir1 = int(iairot1[iref])
        ir2 = int(iairot2[iref])
        # Get all the atoms in the molecule except the two atoms that form the axis.
        labrott1, nrott1 = get_rotgrp(ir1, ir2, labrott1, nrott1)  # noqa: F821 (part 2)
        labrott2, nrott2 = get_rotgrp(ir2, ir1, labrott2, nrott2)  # noqa: F821 (part 2)
        # Now compare the two groups and keep atoms which are common to both groups.
        found = False
        i1 = 0
        while i1 < nrott1:
            i1 = i1 + 1
            found = False
            i2 = 0
            while i2 < nrott2 and not found:
                i2 = i2 + 1
                found = bool(labrott1[i1] == labrott2[i2])
            # If this atom was not found in the second group, we remove it from the first group
            if not found:
                nrott1 = nrott1 - 1
                for i in range(i1, nrott1 + 1):
                    labrott1[i] = labrott1[i + 1]
                i1 = i1 - 1
        units.write(6, lw(ir1, ir2, nrott1))
        units.write(6, lw(*[labrott1[j] for j in range(1, nrott1 + 1)]))
        # need to check that only 1 atom of the group is a neighbour of ir1. Similarly for ir2
        # but there is no easy way to do this since the list of bonds includes bond angles and dihedral angles.
        units.write(6, lw(ir1, ir2, nrott1))
        units.write(6, lw(*[labrott1[j] for j in range(1, nrott1 + 1)]))
        ir = ir + 1
        molecom_mod_lib.testsize_rotgrps(ir)
        molecom_mod.rotgrp[ir] = "ROT"
        molecom_mod.labaxi[ir, 1] = ir1
        molecom_mod.labaxi[ir, 2] = ir2
        molecom_mod.nrot[ir] = nrott1
        molecom_mod.rotf[ir] = molecom_mod.nrottot + 1
        molecom_mod.nrottot = molecom_mod.nrottot + int(molecom_mod.nrot[ir])
        molecom_mod.rotl[ir] = molecom_mod.nrottot
        molecom_mod_lib.testsize_labrot()
        for i in range(1, int(molecom_mod.nrot[ir]) + 1):
            labref = int(molecom_mod.rotf[ir]) + i - 1
            molecom_mod.labrot[labref] = labrott1[i]
        molecom.ndihedralang[ir] = 0
    clustercom.ngrot[1] = ir
    molecom_mod.ngrrot = ir
    #
    # check that the temperature and density have been set
    #
    if molecom.simtk <= 0.0:
        units.write(6, fmt(_F116), advance=False)
        ierr = 1
        # NB: as in the Fortran, this loops for ever at end of standard input
        while ierr != 0:
            molecom.simtk, ierr = _read5_real(molecom.simtk)
    if rho <= 0.0 and not unitcelldefined:
        units.write(6, fmt(_F117), advance=False)
        ierr = 1
        while ierr != 0:
            rho, ierr = _read5_real(rho)
    #
    # write out findings so far
    #
    #     write(6,*) nrow,ncolumnmax,nstar,nlabel
    #     do ia=1,natom
    #        write(6,*) ia,alabel(ia),nbond(ia)
    #     +,(iabond(ia,ib),tbond(ia,ib),ib=1,nbond(ia))
    #     end do
    #     write(6,*) ngrot(1)
    #     do ig=1,ngrot(1)
    #        write(6,*) labaxi(ig,1),labaxi(ig,2)
    #        write(6,*) nrot(ig),(labrot(ig,j),j=1,nrot(ig))
    #     end do
    #
    # now prepare to write the .ato file
    #
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
    clustercom.tethered[1] = False
    clustercom.nato[1] = natom
    molecom.stepmi = r4(0.1)
    if clustercom.ngrot[1] > 0:
        molecom.stepri = r4(1.0)
    else:
        molecom.stepri = r4(0.0)
    molecom.stepr = r4(0.3)
    molecom.stepm = r4(1.0)

    # Initialise the atom types

    molecom.ntype = 0
    molecom.natomr = 0

    for ia in range(clustercom.naf[1], clustercom.nal[1] + 1):
        #
        # change the atom labels if this has been requested
        #
        ichread = 0
        ifound = 0
        while ichread < nchread and ifound == 0:
            ichread = ichread + 1
            if molecom.alabel[ia] == labelold[ichread]:
                ifound = ichread
        if ifound > 0:
            molecom.alabel[ia] = labelnew[ifound]
        ftype(ia)
        if molecom.qtype[ia] != 1:
            molecom.natomr = molecom.natomr + 1

    #
    # check that potential parameters have been defined for all the atom types
    #
    molecom_mod.natomsphere = reallocation_routines.reallocate1d_r(molecom_mod.natomsphere, molecom.ntype)
    molecom_mod.radminsphere = reallocation_routines.reallocate1d_r(molecom_mod.radminsphere, molecom.ntype)
    molecom_mod.radmaxsphere = reallocation_routines.reallocate1d_r(molecom_mod.radmaxsphere, molecom.ntype)
    molecom_mod.pseudofacmol = reallocation_routines.reallocate1d_r(molecom_mod.pseudofacmol, molecom.ntype)
    for it in range(1, molecom.ntype + 1):
        ipread = 0
        ifound = 0
        while ipread < npread and ifound == 0:
            ipread = ipread + 1
            if molecom.type[it] == plabel[ipread]:
                ifound = ipread
        units.write(6, lw(Ch(molecom.type[it], 3), " ", ifound))
        if ifound == 0:
            ierr = 1
            # NB: as in the Fortran, this loops for ever at end of standard input
            while ierr != 0:
                units.write(6, fmt(_F115, Ch(molecom.type[it], 3)), advance=False)
                try:
                    r = units.list_read(5)
                    pread[npread + 1] = r.real(pread[npread + 1])
                    sread[npread + 1] = r.real(sread[npread + 1])
                    mread[npread + 1] = r.real(mread[npread + 1])
                    qread[npread + 1] = r.real(qread[npread + 1])
                    cread[npread + 1] = fstr(r.word(cread[npread + 1]), 3)
                    ierr = 0
                except FortranEOF:
                    ierr = -1
                except FortranIOError as e:
                    ierr = e.iostat
            npread = npread + 1
            plabel[npread] = molecom.type[it]
            natomsphread[npread] = 1.0
            radminsphread[npread] = 0.0
            radmaxsphread[npread] = 0.5 * sread[npread]
            cgmolfacread[npread] = 1.0
            ifound = npread
        molecom.atsymbol[it] = cread[ifound]
        molecom.iexchange[it] = 0
        molecom.epscor[it] = pread[ifound]
        molecom.sigcor[it] = sread[ifound]
        molecom.amassa[it] = mread[ifound]
        molecom.charge[it] = qread[ifound]
        molecom_mod.natomsphere[it] = natomsphread[ifound]
        molecom_mod.radminsphere[it] = radminsphread[ifound]
        molecom_mod.radmaxsphere[it] = radmaxsphread[ifound]
        molecom_mod.pseudofacmol[it] = cgmolfacread[ifound]

        # Now do the same for charge radius, except that if it is not specified
        # we set it to zero.

        ifound = 0
        iqrread = 0
        while iqrread < nqrread and ifound == 0:
            iqrread = iqrread + 1
            if molecom.type[it] == qrlabel[iqrread]:
                ifound = iqrread
        if ifound == 0:
            molecom.qradius[it] = 0.0
        else:
            molecom.qradius[it] = qrread[ifound]
    #
    # set the box size
    #
    units.write(6, lw(molecom.natomr, rho, lattice.alat, lattice.blat, lattice.clat,
                      lattice.phib, lattice.thetac, lattice.phic))
    if not unitcelldefined:
        vol = f32(molecom.natomr) / rho
        clustercom.dmax = powf(vol, r4(1.0) / r4(3.0))
        lattice.alat = clustercom.dmax
        lattice.blat = clustercom.dmax
        lattice.clat = clustercom.dmax
        lattice.phib = r4(90.0)
        lattice.thetac = r4(0.0)
        lattice.phic = r4(0.0)
        lattice.cubic = True
    else:
        lattice.cubic = False
        vol = (lattice.alat * lattice.blat * auxiliary_routines.sind(lattice.phib) * lattice.clat
               * auxiliary_routines.cosd(lattice.thetac))
        clustercom.dmax = powf(vol, r4(1.0) / r4(3.0))
    units.write(6, lw(vol))
    #
    # setup molecular size factor
    #
    molecom.dmaxmol = clustercom.dmax

    #
    # write out the revised .mol file
    #
    try:
        units.open(10, symbols.fnameato, status="unknown")
        ierr = 0
    except FortranIOError as e:
        ierr = e.iostat
    if substr(filetype, 1, ilenfiletype) == ".mol":
        #
        # create a second .atm file which will show the atom
        # numbers for future reference (e.g. to set up rotational groups, etc.)
        #
        fnameatm = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                        + substr(readfile, 1, ilenreadfile) + ".atm", 256)
        #     write(6,*) fnameatm(1:50)

        try:
            units.open(11, fnameatm, status="unknown")
            ierr = 0
        except FortranIOError as e:
            ierr = e.iostat  # noqa: F841
        #
        # first number represents the number of characters in a column (including the delimiting space
        # between columns)
        #
        units.write(10, lw(nspace))
        units.write(11, lw(nspace))
        #
        # maximum number of characters in a line
        #
        ncharmax = (ncolumnmax + 1) * nspace
        if ncharmax < 80:
            ncharmax = 80
        #
        # write a blank line - this can be used for a ruler line to show where the columns occur
        #
        units.write(10, fmt(_F100_makemole, substr(ruler, 1, ncharmax)))
        units.write(11, fmt(_F100_makemole, substr(ruler, 1, ncharmax)))
        #
        # step through the rows and find the labels and stars that occur here
        #
        for ir in range(1, nrow + 1):
            #
            # the first four characters are reserved for the row number plus a space
            #
            commands.sometext = fstr(fmt(_F119, ir), nstsize)
            ncharmax = 0
            is_ = 0
            while is_ < nstar:
                is_ = is_ + 1
                if irstar[is_] == ir:
                    icfirst = (int(icstar[is_]) - 1) * nspace + 2 + 4
                    commands.sometext = setsub(commands.sometext, icfirst, icfirst, "*", nstsize)
                    if icfirst > ncharmax:
                        ncharmax = icfirst
            il = 0
            while il < nlabel:
                il = il + 1
                if irlabel[il] == ir:
                    icfirst = (int(iclabel[il]) - 1) * nspace + 2 + 4
                    if ilabel[il] > 1:
                        commands.sometext = setsub(commands.sometext, icfirst - 1, icfirst - 1,
                                                   fmt(_F103, ilabel[il]), nstsize)
                    commands.sometext = setsub(commands.sometext, icfirst, icfirst + 2, label[il], nstsize)
                    if label[il] != "map":
                        if icfirst + 2 > ncharmax:
                            ncharmax = icfirst + 2
                    else:
                        ilenoffset = findex_padded(offset[il], " ", 20) - 1
                        commands.sometext = setsub(commands.sometext, icfirst + 3, icfirst + 2 + ilenoffset,
                                                   substr(offset[il], 1, ilenoffset), nstsize)
                        if (icfirst + 2 + ilenoffset) > ncharmax:
                            ncharmax = icfirst + 2 + ilenoffset
            units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
            #
            # write the atom numbers
            #
            il = 0
            while il < nlabel:
                il = il + 1
                if irlabel[il] == ir:
                    icfirst = (int(iclabel[il]) - 1) * nspace + 2 + 4
                    if ianumber[1, il] < 10:
                        commands.sometext = setsub(commands.sometext, icfirst, icfirst + 2,
                                                   fmt(_F108, ianumber[1, il]), nstsize)
                    elif ianumber[1, il] < 100:
                        commands.sometext = setsub(commands.sometext, icfirst, icfirst + 2,
                                                   fmt(_F109, ianumber[1, il]), nstsize)
                    else:
                        commands.sometext = setsub(commands.sometext, icfirst - 1, icfirst + 1,
                                                   fmt(_F110, ianumber[1, il]), nstsize)
                    if icfirst + 2 > ncharmax:
                        ncharmax = icfirst + 2
            units.write(11, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
        units.close(11)

    elif substr(filetype, 1, ilenfiletype) == ".gmol":

        units.write(10, lw(substr(filetype, 1, ilenfiletype), iatomlabeladd))

        #       For gmol file, write the atoms

        ia1 = 0
        while ia1 < natom:

            ia1 = ia1 + 1
            units.write(10, lw("atom ", ia1 - iatomlabeladd, " ", Ch(alabelsave[ia1], 3), " ",
                               molecom.xa[ia1], molecom.ya[ia1], molecom.za[ia1], molecom.nbondp[ia1],
                               *[int(molecom.iabond[ia1, j]) - iatomlabeladd
                                 for j in range(1, int(molecom.nbondp[ia1]) + 1)]))
    #
    # now write the bond lengths, angles, extra distances, rotational groups, and potential parameters
    #
    ibread = 0
    while ibread < nbread:
        ibread = ibread + 1
        commands.sometext = ""
        commands.sometext = fstr(fmt(_F104, "bond", Ch(lbond1[ibread], 3), Ch(lbond2[ibread], 3),
                                     bread[ibread]), nstsize)
        ncharmax = 4 + 2 * 4 + 1 + 10
        units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
    iaread = 0
    while iaread < naread:
        iaread = iaread + 1
        commands.sometext = ""
        commands.sometext = fstr(fmt(_F105, "angle", Ch(langle1[iaread], 3), Ch(langle2[iaread], 3),
                                     Ch(langle3[iaread], 3), aread[iaread]), nstsize)
        ncharmax = 5 + 3 * 4 + 1 + 10
        units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
    id = 0
    while id < ndread:
        id = id + 1
        commands.sometext = ""
        commands.sometext = fstr(fmt(_F121, "dihedral", *[int(idread[j, id]) - iatomlabeladd for j in range(1, 5)],
                                     dread[id]), nstsize)
        ncharmax = 8 + 4 * 4 + 1 + 10
        units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
    irread = 0
    while irread < nrread:
        irread = irread + 1
        commands.sometext = ""
        commands.sometext = fstr(fmt(_F107, "rot", int(iarot1[irread]) - iatomlabeladd,
                                     int(iarot2[irread]) - iatomlabeladd), nstsize)
        ncharmax = 3 + 2 * 4
        units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
    irread = 0
    while irread < nwrread:
        irread = irread + 1
        commands.sometext = ""
        commands.sometext = fstr(fmt(_F107, "wholerot", int(iawrot1[irread]) - iatomlabeladd,
                                     int(iawrot2[irread]) - iatomlabeladd), nstsize)
        ncharmax = 8 + 2 * 4
        units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
    irread = 0
    while irread < nirread:
        irread = irread + 1
        commands.sometext = ""
        commands.sometext = fstr(fmt(_F107, "introt", int(iairot1[irread]) - iatomlabeladd,
                                     int(iairot2[irread]) - iatomlabeladd), nstsize)
        ncharmax = 6 + 2 * 4
        units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
    ieread = 0
    while ieread < neread:
        ieread = ieread + 1
        commands.sometext = ""
        commands.sometext = fstr(fmt(_F106, "extra", int(iaextra1[ieread]) - iatomlabeladd,
                                     int(iaextra2[ieread]) - iatomlabeladd, eread[ieread]), nstsize)
        ncharmax = 5 + 2 * 4 + 1 + 10
        units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
    #
    # write the change labels
    #
    ichread = 0
    while ichread < nchread:
        ichread = ichread + 1
        commands.sometext = ""
        commands.sometext = fstr(fmt(_F118, "changelabel", Ch(labelold[ichread], 3), Ch(labelnew[ichread], 3)),
                                 nstsize)
        ncharmax = 11 + 8
        units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
    commands.sometext = fstr(fmt("(a,1x,f8.2)", "pseudofac", pseudofac), nstsize)
    units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, len(commands.sometext))))
    ipread = 0
    while ipread < npread:
        ipread = ipread + 1
        commands.sometext = ""
        commands.sometext = fstr(fmt(_F113, "potential", Ch(plabel[ipread], 3),
                                     pread[ipread], sread[ipread],
                                     mread[ipread], qread[ipread], Ch(cread[ipread], 3),
                                     natomsphread[ipread], radminsphread[ipread], radmaxsphread[ipread],
                                     cgmolfacread[ipread]), nstsize)
        ncharmax = len(commands.sometext)
        units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
    iqrread = 0
    while iqrread < nqrread:
        iqrread = iqrread + 1
        commands.sometext = ""
        commands.sometext = fstr(fmt(_F113, "qradius", Ch(qrlabel[iqrread], 3), qrread[iqrread]), nstsize)
        ncharmax = 8 + 4 + 2 * 11
        units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
    #
    # write temperature and density
    #
    commands.sometext = fstr(fmt(_F114, "temperature", molecom.simtk), nstsize)
    ncharmax = 11 + 14
    units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
    commands.sometext = fstr(fmt(_F114, "vibtemp", molecom.tempb), nstsize)
    ncharmax = 7 + 14
    units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
    commands.sometext = fstr(fmt(_F114, "angtemp", molecom.tempangreq), nstsize)
    ncharmax = 7 + 14
    units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
    commands.sometext = fstr(fmt(_F114, "dihtemp", molecom.tempdihreq), nstsize)
    ncharmax = 7 + 14
    units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
    #
    # write ecore and dcore
    #
    commands.sometext = fstr(fmt(_F120, "ecoredcore", molecom.ecore, molecom.dcore), nstsize)
    ncharmax = 10 + 22
    units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
    if unitcelldefined:
        commands.sometext = fstr(fmt(_F122, "unitcell", lattice.alat, lattice.blat, lattice.clat,
                                     lattice.phib, lattice.thetac, lattice.phic), nstsize)
        ncharmax = 8 + 66
        units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
    else:
        commands.sometext = fstr(fmt(_F114, "density", rho), nstsize)
        ncharmax = 7 + 14
        units.write(10, fmt(_F100_makemole, substr(commands.sometext, 1, ncharmax)))
    units.close(10)

    #
    # write the .ato file
    #
    symbols.fnameato = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                            + substr(readfile, 1, ilenreadfile) + ".ato", 256)
    symbols.ilenfnameato = directory_name.nchardir + ilenreadfile + 4

    #
    # define the number of molecule types and their template files
    #
    clustercom.nmoltype = 1
    clustercom.moltemplate[clustercom.nmoltype] = fstr(substr(readfile, 1, ilenreadfile), 256)
    units.write(6, lw(substr(clustercom.moltemplate[clustercom.nmoltype], 1, ilenreadfile)))
    units.write(6, lw(substr(symbols.fnameato, 1, symbols.ilenfnameato)))
    clustercom.interntrans[1] = 0.0
    clustercom.internrot[1] = 0.0
    clustercom.wholerot[1] = 0.0
    clustercom.wholetrans[1] = 0.0

    ato_read_and_save.moloutc(symbols.fnameato)

    return


# =============================================================================
# from ato_manipulation_rev4_part2.py
# =============================================================================

# --------------------------------------------------------------------------- private helpers
# OpenMP runtime: the reference runs with OMP_NUM_THREADS=1 (TRANSLATION_GUIDE 2.11).
def _omp_get_max_threads():
    return 1


def _omp_get_num_threads():
    return 1


def _atan2d_intrinsic(y, x):
    """gfortran's ATAN2D intrinsic for real*4: ``atan2f(y, x) * (180/pi)`` with the
    constant folded as real(180)/real(pi) = 57.2957764 (checked against gfortran 15)."""
    return atan2f(y, x) * (r4(180.0) / PI4ATAN)


# --------------------------------------------------------------------------- fmolec1
_F101_fmolec1 = "(/'fmolec1> Average no. of neighbours per molecule = ',i0)"
_F102_fmolec1 = "('fmolec1> Average triplet energy per molecule = ',e12.5)"
_F103_fmolec1 = "('fmolec1> Average tethered energy per molecule = ',e12.5)"
_F104_fmolec1 = "('fmolec1> Average coordination constraint energy per molecule = ',e12.5/)"


@routine
def fmolec1(ierr, num_threads_input, ftest, ftestb, ftethered, fcoordconstraints):
    """ato_manipulation_rev4.f:fmolec1 -- compute and store all molecular energies.

    Intramolecular energies go first into each molecule's energy-pair list, then the
    intermolecular pair energies from the cell neighbour lists, three-body, coordination
    constraint and container-force terms.
    Returns (ierr, num_threads_input, ftest, ftestb, ftethered, fcoordconstraints).
    """
    # Set the number of threads for parallel sections.
    # If num_threads is zero, we determine the maximum possible number of threads and use half of this value
    if num_threads_input <= 0:
        num_threads = idiv(_omp_get_max_threads(), 2)
        num_threads_input = num_threads  # Do this to prevent this value being reset if fmolec1 is called more than once.
    else:
        num_threads = num_threads_input
    if num_threads > _omp_get_max_threads():
        num_threads = _omp_get_max_threads()
    # call OMP_SET_NUM_THREADS(num_threads)  -- serial translation: no effect
    units.write(6, lw('fmolec1> Number of parallel threads to be used: ', num_threads,
                      _omp_get_num_threads(), _omp_get_max_threads()))

    ierr = 0

    cellneighbourlist_lib.get_date_and_time('Start initialise_ergpair#', 0)
    ergpairreference_lib.initialise_ergpair(clustercom.nmol)
    cellneighbourlist_lib.get_date_and_time('Finish initialise_ergpair#', 1)

    # Setup the triplets of atoms if three-body terms are requested
    if threebody_mod.nthreeb > 0:
        units.write(6, fmt('(a,i3,1x,i3)', 'fmolec1> No. of three-body terms = ',
                           threebody_mod.nthreeb, threebody_mod.nthreebangsmax))
        threebody_mod_lib.find_triplets(clustercom.nmol, molecom.natom, clustercom.xadd, clustercom.yadd,
                                        clustercom.zadd, molecom.xa, molecom.ya, molecom.za, molecom.mref,
                                        molecom.itype, molecom.ntype, molecom.type)
    if coordconstraints_mod.ncoordconstraints > 0:
        coordconstraint_mod_lib.initialise_coords(molecom.natom)
        for id_ in range(1, coordconstraints_mod.ncoordconstraints + 1):
            coordconstraints_mod.coordconstraints[id_] = coordconstraint_mod_lib.find_constraints(
                coordconstraints_mod.coordconstraints[id_], clustercom.nmol, molecom.natom, clustercom.naf,
                clustercom.nal, clustercom.xadd, clustercom.yadd, clustercom.zadd, molecom.xa, molecom.ya,
                molecom.za, molecom.mref, molecom.itype, molecom.ntype, molecom.type, molecom.cfrac,
                clustercom.dmax)
        coordconstraint_mod_lib.writecoordcons(molecom.natom, molecom.itype, molecom.type)
    #
    # step through atom pairs and calculate separation.
    # If within the specified distance range (rmaxpt) then add to the neighbour list
    # and calculate inter- or intra-molecular energy as appropriate
    #
    ftest = f32(0.0)
    ftestb = f32(0.0)
    # Make sure the intra-molecular energy weights are set
    molecom.tempangwt = molecom.tempangreq
    molecom.tempdihwt = molecom.tempdihreq
    molecom.tempang = molecom.tempangwt
    molecom.tempdih = molecom.tempdihwt * molecom.tempb
    cellneighbourlist_lib.get_date_and_time('fmolec1> Start energy calculation#', 0)
    # NB: sumevib, sumeang, sumedih are not initialised in the Fortran; they are only
    # accumulated by en_mole_intra and never used here, so 0 is as good as any value.
    sumevib = f32(0.0)
    sumeang = f32(0.0)
    sumedih = f32(0.0)
    # Do the intramolecular terms first so that they go in the first element of each energy pair list
    nmol = clustercom.nmol
    for im in range(1, nmol + 1):
        # Get the intramolecular energy of this molecule
        energy, sumevib, sumeang, sumedih = molfit.en_mole_intra(im, sumevib, sumeang, sumedih)
        ftestb = ftestb + energy
        jm = im
        ergpairreference_lib.set_ergpair(im, jm, energy)
    sumetriplets = f32(0.0)
    sumecoordcon = f32(0.0)
    if threebody_mod.nthreeb > 0:
        for id_ in range(1, threebody_mod.nthreeb + 1):
            threebody_mod.sume_triplet_dist[id_] = 0.0
    for ic in range(1, cellneighbourlist.ncell + 1):
        if cellneighbourlist.nmolsincell[ic] > 0:
            ergpairreference.nnbref = cellneighbourlist_lib.make_neighbour_list_from_cells(
                ic, ergpairreference.nnbref, ergpairreference.nbref)
            for iref in range(1, int(cellneighbourlist.nmolsincell[ic]) + 1):
                im = int(cellneighbourlist.molsincell[iref, ic])
                nnbref = ergpairreference.nnbref
                nbref = ergpairreference.nbref
                ergold = ergpairreference.ergold
                # !$omp parallel do private(inn,jm) shared(ergold)
                for inn in range(1, nnbref + 1):
                    jm = int(nbref[inn])
                    ergold[inn] = 0.0
                    if jm > im:  # Avoid double counting pairs of molecules
                        ergold[inn] = molfit.en_mole_pair(im, jm)
                # !$omp end parallel do
                # NB (FORTRAN_BUGS.md section 2): the store loop increments inn before
                # use, so it starts at inn=2 and the first neighbour's pair energy is never
                # stored.  Always active in the Fortran; reproduced.
                ergold = ergpairreference.ergold
                inn = 1
                while inn < nnbref:
                    inn = inn + 1
                    jm = int(nbref[inn])
                    if ergold[inn] != 0.0 and jm > im:
                        ergpairreference_lib.set_ergpair(im, jm, ergold[inn])
                        ftest = ftest + ergold[inn]
                # Get the triplet energy for this molecule
                if threebody_mod.nthreeb > 0:
                    threebody_mod.ergold_triplet[im] = threebody_mod_lib.en_mole_triplets(
                        im, clustercom.xadd, clustercom.yadd, clustercom.zadd, clustercom.phix, clustercom.phiy,
                        clustercom.phiz, molecom.xa, molecom.ya, molecom.za, molecom.mref, clustercom.tethered)
                    sumetriplets = sumetriplets + threebody_mod.ergold_triplet[im]
                    for id_ in range(1, threebody_mod.nthreeb + 1):
                        threebody_mod.sume_triplet_dist[id_] = (threebody_mod.sume_triplet_dist[id_]
                                                                + threebody_mod.en_mole_triplet_dist[id_])
                    ftest = ftest + threebody_mod.ergold_triplet[im]
                if coordconstraints_mod.ncoordconstraints > 0:
                    coordconstraints_mod.ergold_coord[im] = coordconstraint_mod_lib.en_mole_coord(
                        im, clustercom.naf, clustercom.nal, clustercom.xadd, clustercom.yadd, clustercom.zadd,
                        molecom.xa, molecom.ya, molecom.za, molecom.mref, molecom.itype, molecom.qtype,
                        molecom.npairref, molecom.amassr, clustercom.tethered)
                    ftest = ftest + coordconstraints_mod.ergold_coord[im]
                    sumecoordcon = sumecoordcon + coordconstraints_mod.ergold_coord[im]
                if containerforces_mod.numcontainerforces > 0:
                    containerforces_mod.ergold_container[im] = containerforces_mod_lib.en_containerforce(
                        im, clustercom.naf[im], clustercom.nal[im], clustercom.xadd, clustercom.yadd,
                        clustercom.zadd, molecom.xa, molecom.ya, molecom.za, molecom.itype, molecom.mref,
                        molecom.npairref, potential_par.irefpcof, potential_par.rmaxex, potential_par.pcoefexmax,
                        epsrinput_par.expecf, epsrinput_par.expecfpwr, clustercom.dmax, potential_par.harmonic)
                    ftest = ftest + containerforces_mod.ergold_container[im]
    fcoordconstraints = sumecoordcon / f32(nmol)
    units.write(6, lw(ftest / f32(nmol), ftestb / f32(nmol)))
    # For tethered molecules, need to get their tethering energies
    ftethered = f32(0.0)
    for im in range(1, nmol + 1):
        if clustercom.tethered[im]:
            ftethered = ftethered + en_tethered(im, clustercom.tetheredatom[im])  # noqa: F821 (part 1)
    ftethered = ftethered / f32(nmol)
    cellneighbourlist_lib.get_date_and_time('fmolec1> Finish energy calculation#', 1)

    # calculate the average number of neighbours per molecule

    neighbsum = 0  # NB: set by get_neighbsum
    neighbsum = ergpairreference_lib.get_neighbsum(neighbsum)
    units.write(6, fmt(_F101_fmolec1, neighbsum))
    units.write(6, fmt(_F102_fmolec1, sumetriplets / f32(nmol)))
    if threebody_mod.nthreeb > 0:
        units.write(6, lw('fmolec1> Total energy for each triplet: ',
                          *[threebody_mod.sume_triplet_dist[id_] / f32(nmol)
                            for id_ in range(1, threebody_mod.nthreeb + 1)]))
    units.write(6, fmt(_F103_fmolec1, ftethered))
    units.write(6, fmt(_F104_fmolec1, fcoordconstraints))

    # normalise energies per molecule

    ftest = ftest / f32(nmol)
    ftestb = ftestb / f32(nmol)
    return ierr, num_threads_input, ftest, ftestb, ftethered, fcoordconstraints


# --------------------------------------------------------------------------- get_direct_bonds
@routine
def get_direct_bonds():
    """ato_manipulation_rev4.f:get_direct_bonds -- count the direct bonds of each atom.

    Gets the direct bonds between atoms (as opposed to those involving an intermediary
    such as a bond angle).  The assumption is that for a triangle of three atoms, the
    longest bond of the three is not a direct bond, hence this method would not work if
    the bondangle were less than ~60 degrees (which is unlikely in practice).
    """
    nbond = molecom.nbond
    tbond = molecom.tbond
    iabond = molecom.iabond
    nbondp = molecom.nbondp
    # First step is to sort bonds for each atom in order of increasing bond length.
    # This is because to use get_rotgrp we require only the number of direct bonds for each atom - their atom numbers
    # are assumed to be at the beginning of the list, as generated with makemole
    for ia in range(1, molecom.natom + 1):
        ib1 = 0
        ib2 = ib1 + 1
        while ib2 < nbond[ia]:
            ib1 = ib1 + 1
            ib2 = ib1 + 1
            if abs(tbond[ia, ib2]) < abs(tbond[ia, ib1]):
                iabondsave = iabond[ia, ib2]
                tbondsave = tbond[ia, ib2]
                iabond[ia, ib2] = iabond[ia, ib1]
                tbond[ia, ib2] = tbond[ia, ib1]
                iabond[ia, ib1] = iabondsave
                tbond[ia, ib1] = tbondsave
                # Start again to be sure
                ib1 = 0
                ib2 = ib1 + 1

    for ia in range(1, molecom.natom + 1):
        nbondp[ia] = 0
        irefb = 0
        while irefb < nbond[ia]:
            irefb = irefb + 1
            # Only direct bonds are positive
            if tbond[ia, irefb] > 0.0:
                nbondp[ia] = nbondp[ia] + 1
    return


# --------------------------------------------------------------------------- get_rotgrp
@routine
def get_rotgrp(ir1, ir2, labrott, nrott):
    """ato_manipulation_rev4.f:get_rotgrp (line 2579) -- atoms rotated about ir1-ir2.

    labrott is ``allocatable, intent(out)`` in the Fortran (deallocated on entry and
    grown with reallocate1d_i), so the new array is returned, following the
    reallocation_routines convention (docs/SIGNATURE_EXCEPTIONS.md).
    Returns (labrott, nrott).
    """
    # intent(out) allocatable: deallocated on entry
    labrott = None

    # For each of these atoms check how many atoms would be involved in the rotational
    # group. Choose the one with the smallest number of rotatable atoms

    nbondp = molecom.nbondp
    iabond = molecom.iabond
    nrott = 0
    nrottmax = 0
    iselect = ir1
    ir2test = False  # Used to find out if ir2 is joined indirectly to ir1 via direct bonds between other atoms.
    # If so then atoms directly bonded to ir2 need to included in the rotational group.
    nrotnew = 0
    while nrotnew <= nrott:
        #
        # step through the primary bonded atoms associated with the current atom iselect.
        # If it is not one of the axis atoms and does not already exist in the rotational
        # group add it to the group.
        #
        iref = 0
        while iref < nbondp[iselect]:
            iref = iref + 1
            iaref = int(iabond[iselect, iref])
            if not ir2test:
                ir2test = iselect != ir1 and iaref == ir2  # We only test this if it has not already tested positive.

            #
            # check that this atom is not one of the axis atoms
            #
            if iaref != ir1 and iaref != ir2:
                #
                # now check that it does not already exist within the rotational group
                #
                itest = 0
                ifound = 0
                while itest < nrott and ifound == 0:
                    itest = itest + 1
                    if iaref == labrott[itest]:
                        ifound = iaref
                if ifound == 0:
                    #
                    # add this atom to the group
                    #
                    nrott = nrott + 1
                    if nrott > nrottmax:
                        nrottmax = nrottmax + 10
                        labrott = reallocation_routines.reallocate1d_i(labrott, nrottmax)
                    labrott[nrott] = iaref
        nrotnew = nrotnew + 1
        if nrotnew <= nrott:
            iselect = int(labrott[nrotnew])
    # Finally add in any direct bonds associated with ir2, which might have been missed so far
    if ir2test:
        iref = 0
        while iref < nbondp[ir2]:
            iref = iref + 1
            iaref = int(iabond[ir2, iref])
            #
            # check that this atom is not one of the axis atoms
            #
            if iaref != ir1 and iaref != ir2:
                #
                # now check that it does not already exist within the rotational group
                #
                itest = 0
                ifound = 0
                while itest < nrott and ifound == 0:
                    itest = itest + 1
                    if iaref == labrott[itest]:
                        ifound = iaref
                if ifound == 0:
                    #
                    # add this atom to the group
                    #
                    nrott = nrott + 1
                    if nrott > nrottmax:
                        nrottmax = nrottmax + 10
                        labrott = reallocation_routines.reallocate1d_i(labrott, nrottmax)
                    labrott[nrott] = iaref
    units.write(6, lw(ir1, ir2, nrott, *[labrott[i] for i in range(1, nrott + 1)]))

    return labrott, nrott


@routine
def get_whole_rotgrp(im, ir1, ir2, labrott, nrott):
    """ato_manipulation_rev4.f:get_whole_rotgrp (line 2682) -- all atoms of molecule im
    except the two that form the axis.

    labrott: as for get_rotgrp (the reallocated array is returned).
    Returns (labrott, nrott).
    """
    # intent(out) allocatable: deallocated on entry
    labrott = None

    nrott = 0
    nrottmax = 0
    for ia in range(int(clustercom.naf[im]), int(clustercom.nal[im]) + 1):
        if ia != ir1 and ia != ir2:
            #
            # add this atom to the group
            #
            nrott = nrott + 1
            if nrott > nrottmax:
                nrottmax = nrottmax + 10
                labrott = reallocation_routines.reallocate1d_i(labrott, nrottmax)
            labrott[nrott] = ia

    return labrott, nrott


# --------------------------------------------------------------------------- rotatecoords
@routine
def rotatecoords(czphi, szphi, cztheta, sztheta, czchi, szchi, x, y, z):
    """ato_manipulation_rev4.f:rotatecoords -- rotate x,y,z by the Euler angles.

    The rotation zphi about the z-axis is performed first, followed by a rotation ztheta
    about the new y axis, followed by a further rotation zchi about the z-axis.
    Returns (x, y, z).
    """
    # Perform the initial rotation about the z-axis

    xdum = x * czphi + y * szphi
    ydum = y * czphi - x * szphi
    zdum = z

    # Perform the rotation ztheta about the y-axis

    x = xdum * cztheta - zdum * sztheta
    y = ydum
    z = zdum * cztheta + xdum * sztheta

    # Perform final the rotation zchi about the z-axis

    xdum = x * czchi + y * szchi
    ydum = y * czchi - x * szchi
    zdum = z

    x = xdum
    y = ydum
    z = zdum

    return x, y, z


# --------------------------------------------------------------------------- dihedral
@routine
def dihedral(im1, lb1, lb2, ld1, ld2):
    """ato_manipulation_rev4.f:dihedral -- dihedral angle (degrees) about lb1-lb2.

    Calculates the dihedral angle about the axis defined by lb1 and lb2 (in the same
    molecule) by atoms ld1 (bonded to lb1) and ld2 (bonded to lb2).
    Returns dihedral.
    """
    cd = real_array(3)
    co1 = real_array(3)
    co2 = real_array(3)
    cv1 = real_array(3)
    cv2 = real_array(3)
    cpar = real_array(3)
    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za

    pi = PI4ATAN
    ia1 = int(clustercom.naf[im1]) + lb1 - 1
    ia2 = int(clustercom.naf[im1]) + lb2 - 1
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
    dmag = sqrtf(auxiliary_routines.scapro(cd, cd))
    cd[1] = cd[1] / dmag
    cd[2] = cd[2] / dmag
    cd[3] = cd[3] / dmag
    #
    # now determine the orientation vectors co for base atom 1 and base atom 2, using
    # the values of ltyp2 and ltyp3
    #
    io1 = int(clustercom.naf[im1]) + ld1 - 1
    io2 = int(clustercom.naf[im1]) + ld2 - 1
    co1[1] = xa[io1] - x1
    co1[2] = ya[io1] - y1
    co1[3] = za[io1] - z1
    co2[1] = xa[io2] - x2
    co2[2] = ya[io2] - y2
    co2[3] = za[io2] - z2
    #
    # define new vectors normal to both cd and co
    #
    auxiliary_routines.vecpro(co1, cd, cv1)
    auxiliary_routines.vecpro(co2, cd, cv2)
    #
    # renormalize these to unit vectors
    #
    dmag = sqrtf(auxiliary_routines.scapro(cv1, cv1))
    cv1[1] = cv1[1] / dmag
    cv1[2] = cv1[2] / dmag
    cv1[3] = cv1[3] / dmag
    dmag = sqrtf(auxiliary_routines.scapro(cv2, cv2))
    cv2[1] = cv2[1] / dmag
    cv2[2] = cv2[2] / dmag
    cv2[3] = cv2[3] / dmag
    # calculate the cosine of the angle between cv1 and cv2
    cosphi = auxiliary_routines.scapro(cv1, cv2)
    if cosphi > 1.0:
        cosphi = f32(1.0)
    if cosphi < -1.0:
        cosphi = f32(-1.0)
    #
    # calculate the vector product of cv1 and cv2 - this should be either parallel
    # (phi positive) or antiparallel to cd (phi negative)
    #
    auxiliary_routines.vecpro(cv1, cv2, cpar)
    sinphi = auxiliary_routines.scapro(cpar, cd)
    result = atan2f(sinphi, cosphi) * f32(180) / pi

    return result


# --------------------------------------------------------------------------- generate_mp
@routine
def generate_mp():
    """ato_manipulation_rev4.f:generate_mp -- WCA repulsive potential to prevent
    molecular overlaps (xav, xav1), and the container-force pair lists.
    """
    x = epsr16.x
    xstar = epsr16.xstar
    npt = epsr16.npt
    xav = potential_par.xav
    xav1 = potential_par.xav1
    wfac = potential_par.wfac
    power = molecom.power

    pi = PI4ATAN  # noqa: F841 (unused in the Fortran too)
    if power > 0.0:
        ratt = powf(r4(1.0e10), r4(1.0) / power)  # noqa: F841 (unused in the Fortran too)
    #
    # maximum value of a sum - this it to prevent overflows
    #
    summax = expf(r4(25.0))  # noqa: F841 (unused in the Fortran too)
    nterms = potential_par.ncoeffp  # noqa: F841
    decayp = potential_par.psigma2  # noqa: F841
    #
    # set up WCA between pairs.
    #
    # first define prefactor for the well depth
    #
    rat = r4(6.0) / power
    rat1 = sqrtf(rat)
    rat2 = r4(1.0) / rat1
    power1 = (power + 6) / (power - 6)
    rat3 = powf(rat1, power1)
    prefac = (rat2 - rat1) * rat3
    prefac = r4(1.0) / prefac
    # Determine minimum of LJ potential - this is maximum non-zero value of WCA
    rat4 = r4(1.0) / (power - r4(6.0))
    rmaxwca = (molecom.dcore / molecom.sizefactor) * powf(power / r4(6.0), rat4)
    # NB: it and jt are never set in generate_mp (uninitialised locals in the Fortran);
    # 0 reproduces the reference: no atoms match, so the container-force lists are empty
    # (FORTRAN_BUGS.md section 2; Rust potential.rs:561-563).
    it = 0
    jt = 0
    for ic in range(1, molecom.npairs + 1):
        ect = molecom.ecore
        dc = molecom.dcore / molecom.sizefactor
        x1 = x[1]
        x2 = x[2]
        for j in range(2, npt + 1):
            x2 = x[j]
            #
            # modified L-J potential n-6 where n>6
            #
            ddif = dc / x2
            if ddif > 0.0 and ect > 0.0 and x2 < rmaxwca:
                fac = prefac * ect
                xx = logf(ddif)
                pwr1 = r4(0.5) * (power + r4(6.0)) * xx
                if abs(pwr1) > 20.0:
                    pwr1 = sign(r4(20.0), pwr1)
                pwr2 = r4(0.5) * (power - r4(6.0)) * xx
                if abs(pwr2) > 40.0:
                    pwr2 = sign(r4(40.0), pwr2)
                sumlj = r4(2.0) * fac * sinhf(pwr2) * expf(pwr1) + ect
            else:
                sumlj = f32(0.0)
            xav[j, ic] = sumlj * wfac[j]
            x1 = x2
        xav[1, ic] = xav[2, ic]
        #
        # set up gradients of these potentials
        #
        sumlj1 = xav[1, ic]
        x1 = xstar[1]
        x2 = xstar[2]
        for j in range(2, npt + 1):
            x2 = xstar[j]
            xdif = x2 - x1
            sumlj2 = xav[j, ic]
            xav1[j - 1, ic] = (sumlj2 - sumlj1) / xdif
            x1 = x2
            sumlj1 = sumlj2
        xav1[npt, ic] = 0.0
        irefp = int(potential_par.irefpcof[ic])
        # Get a list of atom pairs for this type. If rminex(irefp).ge.0.0 then the list is empty.
        containerforces_mod_lib.get_containerforce_list(clustercom.nmol, ic, potential_par.rmaxex[irefp], it, jt,
                                                        molecom.itype, molecom.qtype, molecom.natom,
                                                        molecom.qtypea, molecom.ntype)

    return


# --------------------------------------------------------------------------- crystfrac / crystarg
@routine
def crystfrac(x, y, z, fa, fb, fc):
    """ato_manipulation_rev4.f:crystfrac -- fractional coordinates of x,y,z in terms of
    the crystallographic axes a,b,c, with the minimum image convention applied.
    Returns (fa, fb, fc).
    """
    if lattice.cubic:
        fc = z / lattice.clatz
        fb = y / lattice.blaty
        fa = x / lattice.alatx
        # Apply minimum image convention
        fa, fb, fc = minimagec.minimage(fa, fb, fc)
    else:
        fc = z / lattice.clatz
        fb = (y - fc * lattice.claty) / lattice.blaty
        fa = (x - fb * lattice.blatx - fc * lattice.clatx) / lattice.alatx
        # Apply minimum image convention
        fa, fb, fc = minimagec.minimage(fa, fb, fc)
    return fa, fb, fc


@routine
def crystarg(fa, fb, fc, h, k, l):  # noqa: E741 (Fortran name)
    """ato_manipulation_rev4.f:crystarg -- 2 pi (fa h + fb k + fc l).
    Returns crystarg.
    """
    return lattice.twopi * (fa * h + fb * k + fc * l)


# --------------------------------------------------------------------------- bragginten
_F101_bragginten = "(1x,a,4(1x,e12.5))"
_F102_bragginten = "('(1x,\"# \",i0,3(1x,i3),1x,f10.5,',i3.3,'(1x,e12.5))')"
_F103_bragginten = "('(1x,f10.5,',i3.3,'(1x,e12.5))')"


@routine
def bragginten(qstep, qmax, qwidthin, fwhmq, multa, multb, multc, usea, useb, usec, hklqmin, nbragg, qbragg,
               intenbragg):
    """ato_manipulation_rev4.f:bragginten -- Bragg intensities up to qmax.

    Calculates Bragg intensities up to the specified Q max assuming a unit cell which
    a/multa, b/multb and c/multc.  intenbragg(npairs,nbragg) is filled in place and the
    <run>.bragg summary file is written.
    """
    B = bragginten_mod
    R = reciprocallattice_mod
    npairs = molecom.npairs
    ntype = molecom.ntype
    natom = molecom.natom

    units.write(6, lw('bragginten> num_threads_to be used: ', _omp_get_max_threads()))
    pi = PI4ATAN
    lattice.twopi = r4(2.0) * pi
    twopi = lattice.twopi
    # Output file name
    outputfile = fstr(substr(symbols.fname_run, 1, symbols.ilenfname_run), 256)
    ncharoutput = symbols.ilenfname_run
    # Generate the reciprocal lattice vectors based on supplied information
    auxiliary_routines.vecpro(lattice.blatc, lattice.clatc, R.recipra)
    auxiliary_routines.vecpro(lattice.clatc, lattice.alatc, R.reciprb)
    auxiliary_routines.vecpro(lattice.alatc, lattice.blatc, R.reciprc)
    vol = auxiliary_routines.scapro(lattice.alatc, R.recipra)
    auxiliary_routines.divvec(R.recipra, vol)
    auxiliary_routines.divvec(R.reciprb, vol)
    auxiliary_routines.divvec(R.reciprc, vol)
    # Correct for number of unit cells along each direction
    auxiliary_routines.mulvec(R.recipra, f32(multa) * twopi)
    auxiliary_routines.mulvec(R.reciprb, f32(multb) * twopi)
    auxiliary_routines.mulvec(R.reciprc, f32(multc) * twopi)
    sqrmoda = auxiliary_routines.scapro(R.recipra, R.recipra)
    sqrmodb = auxiliary_routines.scapro(R.reciprb, R.reciprb)
    sqrmodc = auxiliary_routines.scapro(R.reciprc, R.reciprc)
    # Volume of unit cell
    vol = vol / f32(multa * multb * multc)
    # Get the volume of the reciprocal lattice
    auxiliary_routines.vecpro(R.reciprb, R.reciprc, R.rbragg)
    R.volrecipr = auxiliary_routines.scapro(R.recipra, R.rbragg)
    # Get the amplitude of each reciprocal lattice vector
    R.ampra = sqrtf(sqrmoda)
    R.amprb = sqrtf(sqrmodb)
    R.amprc = sqrtf(sqrmodc)
    recipra = R.recipra
    reciprb = R.reciprb
    reciprc = R.reciprc
    rbragg = R.rbragg
    # Get maximum number of steps in each reciprocal lattice direction
    qwidth = abs(qwidthin) + qmax * fwhmq
    qwidthfac = r4(10.0)
    hmax = usea * (nint((qmax + qwidthfac * qwidth) / abs(recipra[1])) + 1)
    kmax = useb * (nint((qmax + qwidthfac * qwidth) / abs(reciprb[2])) + 1)
    lmax = usec * (nint((qmax + qwidthfac * qwidth) / abs(reciprc[3])) + 1)
    # Maximum number of hkl values
    nhklmax = (hmax + 1) * (2 * kmax + 1) * (2 * lmax + 1)
    units.write(6, lw('Maximum number of hkl values ', nhklmax))
    # Allocate arrays
    B.hihkl = int_array(nhklmax)
    B.kihkl = int_array(nhklmax)
    B.lihkl = int_array(nhklmax)
    B.qhkl = real_array(nhklmax)
    hihkl = B.hihkl
    kihkl = B.kihkl
    lihkl = B.lihkl
    qhkl = B.qhkl
    # Now run over the allowed h,k,l values and calculate scattering amplitudes and intensities
    # (NB: more than nhklmax reflections would write out of bounds in the Fortran; here
    # NumPy raises IndexError.)
    ib = 0
    hi = 0  # Only calculate positive values of hi and use inversion symmetry to get negative hi
    cellneighbourlist_lib.get_date_and_time('Start setting up hkl values and fractional coordinates#', 0)
    while hi <= hmax:
        kmax = useb * (nint((qmax + qwidthfac * qwidth - f32(hi) * recipra[2]) / abs(reciprb[2])) + 1)
        if kmax > 0:
            kmin = -kmax
        else:
            kmin = 0
        ki = kmin
        while ki <= kmax:
            lmax = usec * (nint((qmax + qwidthfac * qwidth - f32(hi) * recipra[3] - f32(ki) * reciprb[3])
                                / abs(reciprc[3])) + 1)
            if lmax > 0:
                lmin = -lmax
            else:
                lmin = 0
            li = lmin
            while li <= lmax:
                for it in range(1, 4):
                    rbragg[it] = f32(hi) * recipra[it] + f32(ki) * reciprb[it] + f32(li) * reciprc[it]
                qhkltest = sqrtf(auxiliary_routines.scapro(rbragg, rbragg))
                if qhkltest <= (qmax + qwidthfac * qwidth):
                    ib = ib + 1
                    hihkl[ib] = hi
                    kihkl[ib] = ki
                    lihkl[ib] = li
                    units.write(6, lw(ib, hi, ki, li))
                    qhkl[ib] = qhkltest
                li = li + 1
            ki = ki + 1
        hi = hi + 1
    B.nhkl = ib
    nhkl = B.nhkl
    units.write(6, lw('Total number of hkl values ', nhkl))
    # (fa is a local allocatable, never allocated on entry: the deallocate is skipped)
    fa = real_array(natom)
    fb = real_array(natom)
    fc = real_array(natom)
    B.ramp = real_array(nhkl, ntype)
    B.iamp = real_array(nhkl, ntype)
    B.intenhkl = real_array(nhkl, npairs)
    B.iqvalues = int_array(nhkl)
    B.qunique = real_array(nhkl)
    B.nqunique = int_array(nhkl)
    natyperef = int_array(ntype)
    iatyperef = int_array(natom, ntype)
    qunique = B.qunique
    nqunique = B.nqunique
    iqvalues = B.iqvalues
    ramp = B.ramp
    iamp = B.iamp
    intenhkl = B.intenhkl
    # Set up the unique q values
    B.ntotalqunique = 0
    ntotalqunique = 0
    nquniquemax = 0
    for ib in range(1, nhkl + 1):
        qhkltest = qhkl[ib]
        iq = 0
        found = False
        while iq < ntotalqunique and not found:
            iq = iq + 1
            found = abs(qhkltest - qunique[iq]) < 1.0e-5
        if not found:
            ntotalqunique = ntotalqunique + 1
            qunique[ntotalqunique] = qhkltest
            nqunique[ntotalqunique] = 1
            nquniquemax = max(nquniquemax, int(nqunique[ntotalqunique]))
        else:
            nqunique[iq] = nqunique[iq] + 1
            nquniquemax = max(nquniquemax, int(nqunique[iq]))
    B.ntotalqunique = ntotalqunique
    B.ibvalues = int_array(nquniquemax, ntotalqunique)
    ibvalues = B.ibvalues
    # Now assign the ib values for each unique Q value.
    # Set counters to zero
    for iq in range(1, ntotalqunique + 1):
        nqunique[iq] = 0
    # Set up the revised Q scale and intensities
    # Store the ib reference numbers for each unique Bragg peak
    for ib in range(1, nhkl + 1):
        qhkltest = qhkl[ib]
        iq = 0
        found = False
        while iq < ntotalqunique and not found:
            iq = iq + 1
            found = abs(qhkltest - qunique[iq]) < 1.0e-5
        if found:
            nqunique[iq] = nqunique[iq] + 1
            ibvalues[nqunique[iq], iq] = ib
            iqvalues[ib] = iq
    # Set up fractional coordinates
    for ia in range(1, natom + 1):
        im = int(molecom.mref[ia])
        xuse = molecom.xa[ia] + clustercom.xadd[im]
        yuse = molecom.ya[ia] + clustercom.yadd[im]
        zuse = molecom.za[ia] + clustercom.zadd[im]
        fa[ia], fb[ia], fc[ia] = crystfrac(xuse, yuse, zuse, fa[ia], fb[ia], fc[ia])
    cellneighbourlist_lib.get_date_and_time('Finished setting up hkl values and fractional coordinates#', 1)
    cellneighbourlist_lib.get_date_and_time('Start calculating Bragg intensities#', 0)
    ramp[:] = 0.0
    iamp[:] = 0.0
    natyperef[:] = 0
    for ia in range(1, natom + 1):
        it = int(molecom.itype[ia])
        natyperef[it] = natyperef[it] + 1
        iatyperef[natyperef[it], it] = ia
    # !$omp parallel do private(it,ir,ia,ib,hi,ki,li,huse,kuse,luse,arg,ramptemp,iamptemp) shared(ramp,iamp)
    for it in range(1, ntype + 1):
        for ir in range(1, int(natyperef[it]) + 1):
            ia = int(iatyperef[ir, it])
            for ib in range(1, nhkl + 1):
                hi = int(hihkl[ib])
                ki = int(kihkl[ib])
                li = int(lihkl[ib])
                huse = f32(hi * multa)
                kuse = f32(ki * multb)
                luse = f32(li * multc)
                arg = twopi * (fa[ia] * huse + fb[ia] * kuse + fc[ia] * luse)
                ramptemp = cosf(arg)
                iamptemp = -sinf(arg)
                ramp[ib, it] = ramp[ib, it] + ramptemp
                iamp[ib, it] = iamp[ib, it] + iamptemp
    # !$omp end parallel do
    # Bragg intensity for each pair of atom types, averaged for the number of unit cells,
    # and per atom
    arg = r4(1.0) / f32(molecom.natomr)
    # Sum the real and imaginary parts
    # !$omp parallel do private(ip,it,jt,ib) shared(intenhkl)
    for ip in range(1, npairs + 1):
        it = int(molecom.n1extr[ip])
        jt = int(molecom.n2extr[ip])
        for ib in range(1, nhkl + 1):
            intenhkl[ib, ip] = ramp[ib, it] * ramp[ib, jt] + iamp[ib, it] * iamp[ib, jt]
            intenhkl[ib, ip] = intenhkl[ib, ip] * arg
            if it != jt:
                intenhkl[ib, ip] = intenhkl[ib, ip] * r4(2.0)
    # !$omp end parallel do
    # Calculate the total Bragg intensity at each unique Q value
    B.intenqunique = real_array(npairs, ntotalqunique)
    intenqunique = B.intenqunique
    # !$omp parallel do private(ip,iq,iref,ib) shared(intenqunique)
    for ip in range(1, npairs + 1):
        for iq in range(1, ntotalqunique + 1):
            for iref in range(1, int(nqunique[iq]) + 1):
                ib = int(ibvalues[iref, iq])
                if hihkl[ib] > 0:
                    # Double intensity to capture inversion symmetry
                    intenqunique[ip, iq] = intenqunique[ip, iq] + r4(2.0) * intenhkl[ib, ip]
                else:
                    intenqunique[ip, iq] = intenqunique[ip, iq] + intenhkl[ib, ip]
    # !$omp end parallel do
    cellneighbourlist_lib.get_date_and_time('Finished calculating Bragg intensities#', 1)
    # Zero intensities and get scattering levels for each file
    for ip in range(1, npairs + 1):
        for ib in range(1, nbragg + 1):
            intenbragg[ip, ib] = 0.0
    cellneighbourlist_lib.get_date_and_time('Start broadening Bragg intensities#', 0)
    # Set up the Bragg peaks, combining hkl which have the same Q value
    # Factor outside Gaussian
    gfactor = R.volrecipr / (r4(2.0) * powf(twopi, r4(1.5)))
    gfactor0 = R.volrecipr / (powf(twopi, r4(1.5)))
    lfactor = R.volrecipr / (pi * pi)
    # Step through hkl values and add broadened intensities to structure factor
    for iq in range(1, ntotalqunique + 1):
        qhkltemp = qunique[iq]
        # Q dependent peak width
        qwidth = abs(qwidthin) + qhkltemp * fwhmq
        if qwidthin < 0.0:
            # For Lorentzian peak intensity rises fast near Q=0 so need to include all the Q values up to the peak position.
            nfirst = 1
        else:
            nfirst = max(1, nint((qhkltemp - r4(15.0) * qwidth) / qstep) - 1)
        nlast = min(nbragg, nint((qhkltemp + r4(15.0) * qwidth) / qstep) + 1)
        if qhkltemp < hklqmin:
            nfirst = 1
            nlast = 0
        it = nfirst - 1
        while it < nlast:
            it = it + 1
            sum_ = f32(0.0)
            if qwidthin > 0.0:
                # Gaussian peak shape
                if qhkltemp > 0.0 and qbragg[it] > 0.0:
                    diff = (qhkltemp - qbragg[it]) / qwidth
                    add = (qhkltemp + qbragg[it]) / qwidth
                    if abs(diff) < 8.0:
                        sum_ = expf(-(r4(0.5) * diff * diff))
                    if abs(add) < 8.0:
                        sum_ = sum_ - expf(-(r4(0.5) * add * add))
                    sum_ = sum_ * gfactor / (qhkltemp * qbragg[it] * qwidth)
                else:
                    add = qhkltemp / qwidth
                    diff = qbragg[it] / qwidth
                    sum_ = r4(0.5) * (diff * diff + add * add)
                    if sum_ < 32.0:
                        sum_ = gfactor0 * expf(-sum_) / pow3(qwidth)
                    else:
                        sum_ = f32(0.0)
            else:
                # Lorentzian peak shape
                diff = (qhkltemp - qbragg[it]) / qwidth
                add = (qhkltemp + qbragg[it]) / qwidth
                sum_ = lfactor / ((r4(1.0) + diff * diff) * (r4(1.0) + add * add) * pow3(qwidth))
            for ip in range(1, npairs + 1):
                intenbragg[ip, it] = intenbragg[ip, it] + intenqunique[ip, iq] * sum_
    cellneighbourlist_lib.get_date_and_time('Finished broadening Bragg intensities#', 1)
    # Subtract the single atom scattering and renormalise the scattering to per atom
    cfraccg = molecom.cfraccg
    ip = 0
    for it in range(1, ntype + 1):
        for jt in range(it, ntype + 1):
            ip = ip + 1
            if it == jt:
                cprod = cfraccg[it] * cfraccg[jt]
                if cprod > 0.0:
                    for ib in range(1, nbragg + 1):
                        intenbragg[ip, ib] = (intenbragg[ip, ib] - cfraccg[it]) / cprod
                else:
                    for ib in range(1, nbragg + 1):
                        intenbragg[ip, ib] = 0.0
            else:
                cprod = r4(2.0) * cfraccg[it] * cfraccg[jt]
                if cprod > 0.0:
                    for ib in range(1, nbragg + 1):
                        intenbragg[ip, ib] = intenbragg[ip, ib] / cprod
                else:
                    for ib in range(1, nbragg + 1):
                        intenbragg[ip, ib] = 0.0
    # Create a filename for the output
    fname = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                 + substr(outputfile, 1, ncharoutput) + '.bragg', 256)
    units.open(10, fname, status='unknown')
    alatc = lattice.alatc
    blatc = lattice.blatc
    clatc = lattice.clatc
    units.write(10, lw('# Multiplicities along a, b and c: ', multa, multb, multc))
    units.write(10, fmt(_F101_bragginten, '# Volume of unit cell: ', vol))
    units.write(10, fmt(_F101_bragginten, '# Volume of reciprocal lattice: ', R.volrecipr))
    units.write(10, fmt(_F101_bragginten, '# Lattice vectors a: ', alatc[1], alatc[2], alatc[3],
                        sqrtf(auxiliary_routines.scapro(alatc, alatc))))
    units.write(10, fmt(_F101_bragginten, '# Lattice vectors b: ', blatc[1], blatc[2], blatc[3],
                        sqrtf(auxiliary_routines.scapro(blatc, blatc))))
    units.write(10, fmt(_F101_bragginten, '# Lattice vectors c: ', clatc[1], clatc[2], clatc[3],
                        sqrtf(auxiliary_routines.scapro(clatc, clatc))))
    units.write(10, fmt(_F101_bragginten, '# Reciprocal lattice vectors a: ', recipra[1], recipra[2], recipra[3],
                        R.ampra))
    units.write(10, fmt(_F101_bragginten, '# Reciprocal lattice vectors b: ', reciprb[1], reciprb[2], reciprb[3],
                        R.amprb))
    units.write(10, fmt(_F101_bragginten, '# Reciprocal lattice vectors c: ', reciprc[1], reciprc[2], reciprc[3],
                        R.amprc))
    units.write(10, fmt('(a,i0)', '# Number of Bragg refections calculated: ', nhkl))
    # Prepare the output format statement
    formatout = fstr(fmt(_F102_bragginten, npairs), 256)
    units.write(6, fmt('(a)', formatout))
    units.write(10, fmt('(1x,"# ",a)', 'Actual Bragg Q values'))
    units.write(10, fmt('(1x,"# ",a)', 'Unique Bragg Q values and multiplicities'))
    for iq in range(1, ntotalqunique + 1):
        units.write(10, fmt('(1x,"#",i6,1x,f10.5,1x,i0)', iq, qunique[iq], nqunique[iq]))
    formatout = fstr(fmt(_F103_bragginten, npairs), 256)
    units.write(6, lw(substr(formatout, 1, findex_padded(formatout, ' ', 256))))
    units.close(10)

    B.ramp = None
    B.iamp = None
    if not diffuse_mod.diffusedefined:
        B.hihkl = None
        B.kihkl = None
        B.lihkl = None
        B.qhkl = None
        B.intenhkl = None
    # deallocate(fa,fb,fc,qunique,nqunique,intenqunique,ibvalues,iqvalues,natyperef,iatyperef)
    B.qunique = None
    B.nqunique = None
    B.intenqunique = None
    B.ibvalues = None
    B.iqvalues = None

    return


@routine
def radd(rsum, rterm):
    """ato_manipulation_rev4.f:radd -- rsum + rterm.
    Returns radd.
    """
    return rsum + rterm


# --------------------------------------------------------------------------- sharm_diffuse
_F333_sharm_diffuse = "(a12,i3.3,a9)"


def _intenhkl_ref(ip, ib):
    """Element ``intenhkl(ip,ib)`` of the array allocated as ``intenhkl(nhkl,npairs)``.

    sharm_diffuse indexes intenhkl with its dimensions swapped; without bounds checking
    gfortran addresses element (ip-1) + (ib-1)*nhkl of the column-major storage.  This
    returns the (row, column) of that element in the (nhkl, npairs) array.
    """
    nhkl = bragginten_mod.intenhkl.shape[0] - 1
    k = (ip - 1) + (ib - 1) * nhkl
    return k % nhkl + 1, k // nhkl + 1


@routine
def sharm_diffuse(qstep, nq, nsums, nfile):
    """ato_manipulation_rev4.f:sharm_diffuse -- spherical-harmonic coefficients of the
    diffuse (non-Bragg) scattering, written to <run>dNN.DIFFUSE.h01, adding to existing
    diffuse scattering if accumulating.

    Uses the reference build's patched run-time FORMAT (tools/build_fortran_reference.sh):
    ``'(1h(,'//il1fmt...`` instead of ``'(1h('//il1fmt...``.
    Returns nsums.
    """
    B = bragginten_mod
    R = reciprocallattice_mod
    D = diffuse_mod
    il1fmt = ''
    in1fmt = ''

    if not D.diffusedefined or B.nhkl <= 0:  # Get out of here if no diffuse scattering has been defined or can be calculated
        return nsums
    pi = PI4ATAN  # noqa: F841 (unused in the Fortran too)
    npairs = molecom.npairs
    nhkl = B.nhkl
    # NB: intenhkl(nhkl,npairs) is indexed intenhkl(ip,ib) below.  When nhkl > npairs the
    # Fortran reads and writes past the end of the array (undefined behaviour); that
    # cannot be reproduced, so stop here as the Rust port does (bragg.rs sharm_diffuse).
    if npairs > 0 and _intenhkl_ref(npairs, nhkl)[1] > npairs:
        raise RuntimeError('sharm_diffuse: intenhkl(ip,ib) runs past the end of intenhkl(nhkl,npairs) '
                           '(nhkl > npairs) - undefined behaviour in the Fortran')
    intenhkl = B.intenhkl
    # Normalise intensities to atomic fractions
    ip = 0
    for it in range(1, molecom.ntype + 1):
        for jt in range(it, molecom.ntype + 1):
            ip = ip + 1
            cprod = molecom.cfrac[it] * molecom.cfrac[jt]
            if it != jt:
                cprod = r4(2.0) * cprod
            for ib in range(1, nhkl + 1):
                r, c = _intenhkl_ref(ip, ib)
                intenhkl[r, c] = intenhkl[r, c] / cprod
    # Get the number of terms
    initialise_diffuse_values(nfile)  # noqa: F821 (part 1, module diffuse_mod)
    intenfile = real_array(nfile)
    hihkl = B.hihkl
    kihkl = B.kihkl
    lihkl = B.lihkl
    qhkl = B.qhkl
    recipra = R.recipra
    reciprb = R.reciprb
    reciprc = R.reciprc
    rbragg = R.rbragg
    diffcoeffc = D.diffcoeffc
    diffcoeffs = D.diffcoeffs
    ndiffcoeff = D.ndiffcoeff
    # Step through the Bragg peaks, picking out those that match the requested diff patterns
    for ib in range(1, nhkl + 1):
        hi = int(hihkl[ib])
        ki = int(kihkl[ib])
        li = int(lihkl[ib])
        # Calculate contribution to sharm coefficients, avoiding Bragg peaks
        if (D.multdiffa * idiv(hi, D.multdiffa) != hi or D.multdiffb * idiv(ki, D.multdiffb) != ki
                or D.multdiffc * idiv(li, D.multdiffc) != li):
            iqdiff = nint(qhkl[ib] / D.qdiffstep) + 1
            for it in range(1, 4):
                rbragg[it] = f32(hi) * recipra[it] + f32(ki) * reciprb[it] + f32(li) * reciprc[it]
            xbragg = rbragg[1]
            ybragg = rbragg[2]
            zbragg = rbragg[3]
            if qhkl[ib] > 0.0:
                ctl = zbragg / qhkl[ib]
            else:
                ctl = f32(1.0)
            # atan2d: gfortran intrinsic (not declared external in this unit)
            phil = _atan2d_intrinsic(ybragg, xbragg)
            # Weight the intensities for each data file
            for iff in range(1, nfile + 1):
                ix = int(weights.xfile[iff])
                intenfile[iff] = 0.0
                if ix > 0:
                    # If X-ray data, need to apply the Q-dependent weights
                    iq = nint(qhkl[ib] / qstep) + 1
                    if iq > 0 and iq <= nq:
                        for ip in range(1, npairs + 1):
                            r, c = _intenhkl_ref(ip, ib)
                            # weightx is real*8: the sum is done in double, then stored as real
                            intenfile[iff] = float(intenfile[iff]) + x_rayf0.weightx[iq, ip, ix] * float(intenhkl[r, c])
                else:
                    for ip in range(1, npairs + 1):
                        r, c = _intenhkl_ref(ip, ib)
                        # cwt is real*8
                        intenfile[iff] = float(intenfile[iff]) + weights.cwt[iff, ip] * float(intenhkl[r, c])
            for id_ in range(1, D.ndiffuse + 1):
                l = int(D.l1ex[id_])  # noqa: E741 (Fortran name)
                m = int(D.n1ex[id_])
                factor = f32(2 * l + 1)
                dlmcoeff = factor * auxiliary_routines.dlm(l, m, 0, ctl)
                # cosd/sind: gfortran intrinsics (not declared external in this unit)
                dlmcoeffc = cosd_intrinsic(f32(m) * phil) * dlmcoeff
                dlmcoeffs = sind_intrinsic(f32(m) * phil) * dlmcoeff
                # If h > 0, add on the term for which (h,k,l) --> (-h,-k,-l) (Inversion image, theta --> pi-theta, phi --> phi + pi)
                if hi > 0:
                    dlmcoeff = factor * auxiliary_routines.dlm(l, m, 0, -ctl)
                    dlmcoeffc = dlmcoeffc + cosd_intrinsic(f32(m) * (phil + r4(180.0))) * dlmcoeff
                    dlmcoeffs = dlmcoeffs + sind_intrinsic(f32(m) * (phil + r4(180.0))) * dlmcoeff
                for iff in range(1, nfile + 1):
                    diffcoeffc[id_, iqdiff, iff] = diffcoeffc[id_, iqdiff, iff] + dlmcoeffc * intenfile[iff]
                    diffcoeffs[id_, iqdiff, iff] = diffcoeffs[id_, iqdiff, iff] + dlmcoeffs * intenfile[iff]
            if abs(hi) > 0:
                ndiffcoeff[iqdiff] = ndiffcoeff[iqdiff] + 2  # Account for negative hi
            else:
                ndiffcoeff[iqdiff] = ndiffcoeff[iqdiff] + 1
    nblock = 2 * D.ndiffuse + 1
    # Set up a title for the first line
    formattitle = fstr(fmt(_F333_sharm_diffuse, '(1h#,1x,a14,', nblock, '(1x,a29))'), 256)  # noqa: F841
    ifmt1 = 'i1'
    ifmt2 = 'i2'
    ifmt3 = 'i3'
    ifmt4 = 'i4'
    titledat = char_array(D.ndiffuse)
    for id_ in range(1, D.ndiffuse + 1):
        il1 = int(D.l1ex[id_])
        in1 = int(D.n1ex[id_])
        if il1 < 10:
            il1fmt = ifmt1
        elif il1 < 100:
            il1fmt = ifmt2
        elif il1 < 1000:
            il1fmt = ifmt3
        elif il1 < 10000:
            il1fmt = ifmt4
        if in1 < 10:
            in1fmt = ifmt1
        elif in1 < 100:
            in1fmt = ifmt2
        elif in1 < 1000:
            in1fmt = ifmt3
        elif in1 < 10000:
            in1fmt = ifmt4
        # Patched as in the reference build: '(1h(' -> '(1h(,' (tools/build_fortran_reference.sh)
        formatout = fstr('(1h(,' + il1fmt + ',1h,,' + ifmt1 + ',1h,,' + il1fmt + ',1h;,' + in1fmt + ',1h,,' + ifmt1
                         + ',1h))', 256)
        titledat[id_] = fstr(fmt(formatout, il1, 0, il1, in1, 0), 29)
    titler = '   Q'  # axis title
    # Set up the files from and to which the diffuse scattering will be read/written
    nqdiff = D.nqdiff
    ndiffuse = D.ndiffuse
    D.diffcoefftc = real_array(nqdiff, ndiffuse)
    D.diffcoeffts = real_array(nqdiff, ndiffuse)
    D.diffcoeffoldc = real_array(nqdiff, ndiffuse)
    D.diffcoeffolds = real_array(nqdiff, ndiffuse)
    tempdiff = real_array(nqdiff)
    outputfile = fstr(substr(symbols.fname_run, 1, symbols.ilenfname_run), 256)
    ncharoutput = findex(outputfile, '.EPSR') - 1
    for iff in range(1, nfile + 1):
        D.diffusename = fstr(fmt('("d",i2.2)', iff), 3)
        fname = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                     + substr(outputfile, 1, ncharoutput) + substr(D.diffusename, 1, 3) + '.DIFFUSE.h01', 256)
        for iq in range(1, nqdiff + 1):
            if ndiffcoeff[iq] > 0:  # If possible average coefficients to number of contributions
                for id_ in range(1, ndiffuse + 1):
                    diffcoeffc[id_, iq, iff] = diffcoeffc[id_, iq, iff] / f32(ndiffcoeff[iq])
                    diffcoeffs[id_, iq, iff] = diffcoeffs[id_, iq, iff] / f32(ndiffcoeff[iq])
        # Smooth the coefficients prior to generating maps. Don't do this very much - it distorts image drastically
        for id_ in range(1, ndiffuse + 1):
            for iq in range(1, nqdiff + 1):
                tempdiff[iq] = diffcoeffc[id_, iq, iff]
            auxiliary_routines.smood(1, nqdiff, 1, tempdiff)
            for iq in range(1, nqdiff + 1):
                D.diffcoefftc[iq, id_] = tempdiff[iq]
            for iq in range(1, nqdiff + 1):
                tempdiff[iq] = diffcoeffs[id_, iq, iff]
            auxiliary_routines.smood(1, nqdiff, 1, tempdiff)
            for iq in range(1, nqdiff + 1):
                D.diffcoeffts[iq, id_] = tempdiff[iq]
        ext = 'h'
        from . import read_and_save_files  # imported here: avoids an import cycle at module load
        nsums = read_and_save_files.outgofzs_sharm(fname, ext, nqdiff, ndiffuse, nsums, D.qdiff, D.diffcoefftc,
                                                   D.diffcoeffts, D.diffcoeffoldc, D.diffcoeffolds, titler,
                                                   titledat, nqdiff, ndiffuse)
    # deallocate(intenfile,tempdiff,diffcoefftc,diffcoeffts,diffcoeffoldc,diffcoeffolds)
    D.diffcoefftc = None
    D.diffcoeffts = None
    D.diffcoeffoldc = None
    D.diffcoeffolds = None
    if B.intenhkl is not None:
        B.qhkl = None
        B.intenhkl = None
    if D.diffcoeffc is not None:
        D.qdiff = None
        D.diffcoeffc = None
        D.diffcoeffs = None
        D.ndiffcoeff = None
    return nsums


# --------------------------------------------------------------------------- jhalf / truncconv
@routine
def jhalf(x):
    """ato_manipulation_rev4.f:jhalf -- sin(x)/x (1 at x = 0).
    Returns jhalf.
    """
    if abs(x) > 0.0:
        result = sinf(x) / x
    else:
        result = f32(1.0)
    return result


@routine
def truncconv(q, q0, rm, qstep, pi):
    """ato_manipulation_rev4.f:truncconv -- truncation convolution function for the
    peak at q0 when the maximum radius is rm.
    Returns truncconv.
    """
    result = f32(0.0)
    if abs(q) > 0.0:
        if abs(q0) > 0.0:
            factor = q0 * rm * qstep / (pi * q)
            arg1 = (q0 - q) * rm
            arg2 = (q0 + q) * rm
            result = factor * (jhalf(arg1) - jhalf(arg2))
        else:
            factor = 2 * q0 * q0 * rm * qstep / (pi * q * q)
            arg = q * rm
            result = factor * (jhalf(arg) - cosf(arg))
    elif abs(q0) > 0.0:
        factor = 2 * rm * qstep / (pi)
        arg = q0 * rm
        result = factor * (jhalf(arg) - cosf(arg))
    else:
        # truncconv=2*(rm*rm)*rm/(pi*3.0)
        result = f32(0.0)
    return result


# --------------------------------------------------------------------------- change_boxsize
@routine
def change_boxsize(afac, bfac, cfac, newphib, newthetac, newphic):
    """ato_manipulation_rev4.f:change_boxsize -- multiply the lattice vectors by the
    specified factors and relocate molecules relative to the new axes.

    cosd/sind are not declared external here: gfortran's intrinsics are used.
    """
    L = lattice
    C = clustercom

    # Need to place each molecule in same relative position in the box as at present
    # First generate the current fractional coordinates

    for im in range(1, C.nmol + 1):
        C.xadd[im], C.yadd[im], C.zadd[im] = generate_fractional_coordinates(  # noqa: F821 (part 1)
            C.xadd[im], C.yadd[im], C.zadd[im])
        C.phix[im], C.phiy[im], C.phiz[im] = generate_fractional_coordinates(  # noqa: F821 (part 1)
            C.phix[im], C.phiy[im], C.phiz[im])

    # set up lattice parameters for applying the minimum image convention

    L.alat = L.alat * afac
    L.blat = L.blat * bfac
    L.clat = L.clat * cfac
    L.phib = newphib
    L.thetac = newthetac
    L.phic = newphic
    L.alatx = L.alat
    L.blatx = L.blat * cosd_intrinsic(L.phib)
    L.blaty = L.blat * sind_intrinsic(L.phib)
    L.clatx = L.clat * sind_intrinsic(L.thetac) * cosd_intrinsic(L.phic)
    L.claty = L.clat * sind_intrinsic(L.thetac) * sind_intrinsic(L.phic)
    L.clatz = L.clat * cosd_intrinsic(L.thetac)

    # Now convert fractional coordinates back to Cartesian coordinates

    for im in range(1, C.nmol + 1):
        C.xadd[im], C.yadd[im], C.zadd[im] = generate_cartesian_coordinates(  # noqa: F821 (part 1)
            C.xadd[im], C.yadd[im], C.zadd[im])
        C.phix[im], C.phiy[im], C.phiz[im] = generate_cartesian_coordinates(  # noqa: F821 (part 1)
            C.phix[im], C.phiy[im], C.phiz[im])

    #
    # define the Cartesian coordinates of the lattice vectors
    #
    L.alatc[1] = L.alatx
    L.alatc[2] = 0.0
    L.alatc[3] = 0.0
    L.blatc[1] = L.blatx
    L.blatc[2] = L.blaty
    L.blatc[3] = 0.0
    L.clatc[1] = L.clatx
    L.clatc[2] = L.claty
    L.clatc[3] = L.clatz
    L.bxax = L.blatx / L.alatx
    L.cxax = L.clatx / L.alatx
    L.cyby = L.claty / L.blaty
    #
    # calculate the perpendicular distance between the faces of the lattice
    #
    perpab = auxiliary_routines.perpd(L.alatc, L.blatc, L.clatc)
    perpbc = auxiliary_routines.perpd(L.blatc, L.clatc, L.alatc)
    perpca = auxiliary_routines.perpd(L.clatc, L.alatc, L.blatc)
    units.write(6, lw(perpab, perpbc, perpca))

    # set up maximum radius to be used for Fourier transforms and to truncate
    # the potential. This is defined as the minimum perpendicular distance between
    # planes.
    #
    # G(r) and potentials will be calculated to rmax, or maximum allowed distance
    # based on number of points and step size, whichever is smaller.
    #
    C.dmax = min(perpab, perpbc, perpca)
    epsr16.cell = C.dmax
    if L.alat != L.blat or L.blat != L.clat or L.clat != L.alat:
        L.cubic = False
    elif L.phib != 90.0 or L.thetac != 0.0 or L.phic != 0.0:
        L.cubic = False
    else:
        L.cubic = True

    return


# --------------------------------------------------------------------------- atogen
@routine
def atogen(wav, ngrsamples):
    """ato_manipulation_rev4.f:atogen -- histogram the configuration into the
    intramolecular (wav(1,..)) and intermolecular (wav(2,..)) pair distributions,
    normalised by the number of pairs, and count close contacts (potex).

    wav is real wav(2,mpz,mpairs), filled in place.
    Returns ngrsamples.
    """
    # Local aliases for speed (none of these arrays is reallocated inside atogen,
    # except molecom_mod.wavinter which is fetched after its allocation).
    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za
    itype = molecom.itype
    npairref = molecom.npairref
    xadd = clustercom.xadd
    yadd = clustercom.yadd
    zadd = clustercom.zadd
    naf = clustercom.naf
    nal = clustercom.nal
    mtype = clustercom.mtype
    potex = potential_par.potex
    potexmax = potential_par.potexmax
    irefpcof = potential_par.irefpcof
    rminex = potential_par.rminex
    npt1 = epsr16.npt1
    npairs = molecom.npairs
    nmol = clustercom.nmol
    sizefactor = molecom.sizefactor
    dmax = clustercom.dmax
    coarsegraining = molecom_mod.coarsegraining
    mic = minimagec.minimagec

    sumwavself = real_array(epsr_dimension.mpairs)

    pi = PI4ATAN
    pi2 = r4(2.0) * pi
    spi2 = sqrtf(pi2)  # noqa: F841 (unused in the Fortran too)
    #
    # set weights to zero
    #
    for iz in range(1, npt1 + 1):
        for in_ in range(1, npairs + 1):
            wav[1, iz, in_] = 0.0
            wav[2, iz, in_] = 0.0
    if npairs > 0:
        for ic in range(1, npairs + 1):
            potex[ic] = 0.0
            potexmax[ic] = 0.0
    if coarsegraining:
        molecom_mod.wavinter = real_array(npt1, npairs, molecom_mod.maxnscalefactype,
                                          molecom_mod.maxnscalefactype)
        molecom_mod.wavinter[:] = 0.0
    wavinter = molecom_mod.wavinter
    stepcl = epsr16.cellst / dmax
    #
    # generate intramolecular spacings
    #
    for im in range(1, nmol + 1):
        if clustercom.nato[im] > 1:
            for ia in range(int(naf[im]), int(nal[im])):
                x1 = xa[ia]
                y1 = ya[ia]
                z1 = za[ia]
                for ja in range(ia + 1, int(nal[im]) + 1):
                    it = itype[ia]
                    jt = itype[ja]
                    icref = npairref[it, jt]
                    x2 = xa[ja]
                    y2 = ya[ja]
                    z2 = za[ja]
                    xx = (x1 - x2) * sizefactor
                    yy = (y1 - y2) * sizefactor
                    zz = (z1 - z2) * sizefactor
                    #
                    # apply minimum image convention
                    #
                    xx, yy, zz = mic(xx, yy, zz)
                    dtest = sqrtf(xx * xx + yy * yy + zz * zz)
                    # find the x value this corresponds to
                    iz1 = nint(dtest / stepcl) + 2
                    #
                    # accumulate the partial g(r) for this pair
                    #
                    if iz1 > 1 and iz1 <= npt1:
                        wav[1, iz1, icref] = wav[1, iz1, icref] + r4(2.0)
    #
    # generate intermolecular spacings
    #
    csampled = real_array(molecom.ntype)
    iz1values = int_array(molecom.natom)
    icrefvalues = int_array(molecom.natom)
    molsampled = int_array(clustercom.nmoltype)
    molfraction = real_array(clustercom.nmoltype)
    # Define probability for using a molecule
    if ngrsamples <= 0:
        ngrsamples = 1000
    if ngrsamples > nmol:
        ngrsamples = nmol
    probaccept = f32(ngrsamples) / f32(nmol)
    units.write(6, lw('atogen> Acceptance probability = ', probaccept))
    for it in range(1, molecom.ntype + 1):
        csampled[it] = 0.0
    # Get the molecular fraction of each molecule type
    for it in range(1, clustercom.nmoltype + 1):
        molsampled[it] = 0
        molfraction[it] = 0.0
    for im in range(1, nmol + 1):
        molfraction[mtype[im]] = molfraction[mtype[im]] + r4(1.0)
    for it in range(1, clustercom.nmoltype + 1):
        molfraction[it] = molfraction[it] / f32(nmol)
    im1 = 1
    im2 = nmol
    addwav = f32(1.0)
    if probaccept >= 1.0:
        im1 = 1
        im2 = nmol - 1
        addwav = f32(2.0)
        # If probaccept.ge.1.0 this means we only access nmol-1 first molecules. Hence we need to count the last molecule when calculating the number of
        # first molecules sampled.
        for ia in range(int(naf[nmol]), int(nal[nmol]) + 1):
            it = itype[ia]
            csampled[it] = csampled[it] + r4(1.0)
    for im in range(im1, im2 + 1):
        # Always select molecules that are weakly represented
        # (ran1() is the first operand, so it is called for every im)
        if (auxiliary_routines.ran1() < probaccept or molsampled[mtype[im]] < 1
                or molfraction[mtype[im]] < probaccept):
            molsampled[mtype[im]] = molsampled[mtype[im]] + 1
            xc1 = xadd[im]
            yc1 = yadd[im]
            zc1 = zadd[im]
            jm1 = 1
            jm2 = nmol
            if probaccept >= 1.0:
                jm1 = im + 1
                jm2 = nmol
            for ia in range(int(naf[im]), int(nal[im]) + 1):
                it = itype[ia]
                csampled[it] = csampled[it] + 1
            for ia in range(int(naf[im]), int(nal[im]) + 1):
                x1 = xa[ia] + xc1
                y1 = ya[ia] + yc1
                z1 = za[ia] + zc1
                it = itype[ia]
                # !$omp parallel do private(jm,xc2,yc2,zc2,ja,jaref,jt,icref,x2,y2,z2,xx,yy,zz,dtest,dtesttest,iz1,irefp,scalefac) shared(iz1values,icrefvalues,potex)
                for jm in range(jm1, jm2 + 1):
                    if jm != im:
                        xc2 = xadd[jm]
                        yc2 = yadd[jm]
                        zc2 = zadd[jm]
                        for ja in range(int(naf[jm]), int(nal[jm]) + 1):
                            jt = itype[ja]
                            icref = npairref[it, jt]
                            x2 = xa[ja] + xc2
                            y2 = ya[ja] + yc2
                            z2 = za[ja] + zc2
                            xx = x1 - x2
                            yy = y1 - y2
                            zz = z1 - z2
                            #
                            # apply minimum image convention
                            #
                            xx, yy, zz = mic(xx, yy, zz)
                            dtest = sqrtf(xx * xx + yy * yy + zz * zz)
                            # find the x value this corresponds to
                            iz1 = nint(dtest / stepcl) + 2
                            #
                            # accumulate the partial g(r) for this pair
                            #
                            if iz1 > 1 and iz1 <= npt1:
                                iz1values[ja] = iz1
                                icrefvalues[ja] = icref
                            else:
                                iz1values[ja] = 0
                            #
                            # now check whether this pair of atoms is below a specified minimum distance or above a specified maximum distance
                            #
                            dtesttest = dtest * dmax
                            scalefac = f32(1.0)
                            if coarsegraining:
                                scalefac = r4(0.5) * (molecom_mod.scalefacmol[im] + molecom_mod.scalefacmol[jm])
                            irefp = irefpcof[icref]
                            # if(.not.check_coord_bonded(ia,ja,mref,itype)) then
                            # if(pcoefex(irefp).ge.0.0.and..not.check_coord_bonded(ia,ja,mref,itype)) then
                            # if(pcoefex(irefp).ge.0.0) then
                            if dtesttest < scalefac * abs(rminex[irefp]):
                                potex[icref] = potex[icref] + r4(1.0)
                            # end if
                    else:
                        for ja in range(int(naf[jm]), int(nal[jm]) + 1):
                            iz1values[ja] = 0
                # !$omp end parallel do
                for ja in range(int(naf[jm1]), int(nal[jm2]) + 1):
                    iz1 = iz1values[ja]
                    if iz1 > 0:
                        icref = icrefvalues[ja]
                        wav[2, iz1, icref] = wav[2, iz1, icref] + addwav
                # Repeat for coarse graining
                if coarsegraining:
                    ibinia = molecom_mod.assignatomscalefac[ia]
                    for ja in range(int(naf[jm1]), int(nal[jm2]) + 1):
                        iz1 = iz1values[ja]
                        if iz1 > 0:
                            jt = itype[ja]
                            ibinja = molecom_mod.assignatomscalefac[ja]
                            icref = icrefvalues[ja]
                            if it <= jt:  # Pairs are referred in the order it <= jt
                                wavinter[iz1, icref, ibinia, ibinja] = wavinter[iz1, icref, ibinia, ibinja] + addwav
                            else:
                                wavinter[iz1, icref, ibinja, ibinia] = wavinter[iz1, icref, ibinja, ibinia] + addwav
        # Check overlaps if container forces are in use
        if containerforces_mod.numcontainerforces > 0.0:  # Deal with container forces
            containerforces_mod_lib.dist_containerforce(im, naf[im], nal[im], xadd, yadd, zadd, xa, ya, za,
                                                        itype, molecom.mref, npairref, irefpcof,
                                                        potential_par.rmaxex, potexmax, dmax, r4(1.0))
    #
    # normalize sums
    #
    sumwavself[:] = 0
    cfrac = molecom.cfrac
    ycoh = epsr16.ycoh
    for ic in range(1, npairs + 1):
        it = int(molecom.n1extr[ic])
        jt = int(molecom.n2extr[ic])
        cprod = cfrac[it] * cfrac[jt] * f32(molecom.natomr)
        cprod2 = csampled[it] * cfrac[jt]
        if it != jt:
            cprod = cprod * 2
            cprod2 = cprod2 + csampled[jt] * cfrac[it]
        if cprod > 0.0:
            for iz in range(1, npt1 + 1):
                wav[1, iz, ic] = wav[1, iz, ic] / cprod
                sumwavself[ic] = sumwavself[ic] + wav[1, iz, ic]
        else:
            for iz in range(1, npt1 + 1):
                wav[1, iz, ic] = 0.0
        if cprod2 > 0.0:
            for iz in range(1, npt1 + 1):
                wav[2, iz, ic] = wav[2, iz, ic] / cprod2
            # Normalise potex as well, to give the coordination number within rminex from atom type
            potex[ic] = potex[ic] * cfrac[jt] / cprod2
            if coarsegraining:
                for iz in range(1, npt1 + 1):
                    for ibinia in range(1, int(molecom_mod.nsffbintype[it]) + 1):
                        for ibinja in range(1, int(molecom_mod.nsffbintype[jt]) + 1):
                            wavinter[iz, ic, ibinia, ibinja] = wavinter[iz, ic, ibinia, ibinja] / cprod2
        else:
            for iz in range(1, npt1 + 1):
                wav[2, iz, ic] = ycoh[iz]
            if coarsegraining:
                ycohfactor = r4(1.0) / f32(int(molecom_mod.nsffbintype[it]) * int(molecom_mod.nsffbintype[jt]))
                for iz in range(1, npt1 + 1):
                    ycohmod = ycoh[iz] * ycohfactor
                    for ibinia in range(1, int(molecom_mod.nsffbintype[it]) + 1):
                        for ibinja in range(1, int(molecom_mod.nsffbintype[jt]) + 1):
                            wavinter[iz, ic, ibinia, ibinja] = ycohmod

    return ngrsamples


# --------------------------------------------------------------------------- press
@routine
def press(energy, pressc, energf, pressf, energabs, energthreeb, prthreebout):
    """ato_manipulation_rev4.f:press -- energy and pressure from intermolecular
    interactions.
    Returns (energy, pressc, energabs, energthreeb, prthreebout).
    """
    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za
    itype = molecom.itype
    mref = molecom.mref
    npairref = molecom.npairref
    xadd = clustercom.xadd
    yadd = clustercom.yadd
    zadd = clustercom.zadd
    naf = clustercom.naf
    nal = clustercom.nal
    xstar = epsr16.xstar
    npt1 = epsr16.npt1
    epot = potential_par.epot
    epot1 = potential_par.epot1
    ecoul = potential_par.ecoul
    ecoul1 = potential_par.ecoul1
    xavs = potential_par.xavs
    xavs1 = potential_par.xavs1
    nmol = clustercom.nmol
    mic = minimagec.minimagec
    check_coord_bonded = coordconstraint_mod_lib.check_coord_bonded
    # check_coord_bonded returns .false. without looking at anything when there are
    # no coordination constraints; skipping the call then gives the same result.
    use_ccb = coordconstraints_mod.ncoordconstraints > 0

    stepcl = epsr16.cellst / clustercom.dmax
    sumc = f32(0.)
    presst = 0.
    #
    # absolute energy of EP coefficients
    #
    energabs = abs_energy_ep()
    for im in range(1, nmol):
        # Get neighbour list
        ergpairreference_lib.make_neighbour_list_from_mols(im)
        nnbref = ergpairreference.nnbref
        nbref = ergpairreference.nbref
        ergold = ergpairreference.ergold
        ergnew = ergpairreference.ergnew
        xc1 = xadd[im]
        yc1 = yadd[im]
        zc1 = zadd[im]
        # Zero the energy and pressure accumulators
        for jref in range(1, nnbref + 1):
            ergold[jref] = 0
            ergnew[jref] = 0
        # !$omp parallel do private(...) shared(ergold,ergnew)
        for jref in range(1, nnbref + 1):
            jm = int(nbref[jref])
            # Avoid double counting pairs
            if jm > im:
                for ia in range(int(naf[im]), int(nal[im]) + 1):
                    it = itype[ia]
                    x1 = xa[ia] + xc1
                    y1 = ya[ia] + yc1
                    z1 = za[ia] + zc1
                    for ja in range(int(naf[jm]), int(nal[jm]) + 1):
                        xc2 = xadd[jm]
                        yc2 = yadd[jm]
                        zc2 = zadd[jm]
                        jt = itype[ja]
                        icref = npairref[it, jt]
                        #
                        # set up separation of atom pair
                        #
                        x2 = xa[ja] + xc2
                        y2 = ya[ja] + yc2
                        z2 = za[ja] + zc2
                        xdif = x2 - x1
                        ydif = y2 - y1
                        zdif = z2 - z1
                        #
                        # apply minimum image convention
                        #
                        xdif, ydif, zdif = mic(xdif, ydif, zdif)
                        #
                        # separation of centres (for pressure calculation)
                        #
                        xcdif = xc2 - xc1
                        ycdif = yc2 - yc1
                        zcdif = zc2 - zc1
                        xcdif, ycdif, zcdif = mic(xcdif, ycdif, zcdif)
                        dsq = xdif * xdif + ydif * ydif + zdif * zdif
                        dtest = sqrtf(dsq)
                        # find the r value this corresponds to
                        rat = dtest / stepcl
                        iz = ifix(rat) + 2
                        if iz <= npt1:
                            dtdif = dtest - xstar[iz]
                            ipz = iz
                            if dtest > 0.0:
                                prfac = (xcdif * xdif + ycdif * ydif + zcdif * zdif) / dtest
                            else:
                                prfac = f32(0.0)
                            # sum reference potential
                            # (prpot is real*8: the real*4 product is converted on assignment)
                            if use_ccb and check_coord_bonded(ia, ja, mref, itype):
                                enpot = epot[iz, icref] + epot1[iz, icref] * dtdif
                                prpot = float(r4(0.5) * (epot1[ipz, icref] + epot1[ipz - 1, icref]) * prfac)
                            else:
                                enpot = ecoul[iz, icref] + ecoul1[iz, icref] * dtdif
                                prpot = float(r4(0.5) * (ecoul1[ipz, icref] + ecoul1[ipz - 1, icref]) * prfac)
                            #
                            # sum intermolecular energy and pressure from empirical potential
                            #
                            enwt = xavs[iz, icref] + xavs1[iz, icref] * dtdif
                            # For absolute empirical energy do not include the repulsive term (but this is no longer included
                            # in xavs, so we can use xavs).
                            enabswt = enwt  # noqa: F841
                            prwt = r4(0.5) * (xavs1[ipz, icref] + xavs1[ipz - 1, icref]) * prfac
                            ergold[jref] = ergold[jref] + enwt + enpot
                            # (ergnew + prwt) in real*4, + prpot in real*8, stored as real*4
                            ergnew[jref] = float(ergnew[jref] + prwt) + prpot
        # !$omp end parallel do
        # Accumulate values
        for jref in range(1, nnbref + 1):
            sumc = sumc + ergold[jref]
            presst = presst + float(ergnew[jref])
    # Add on the triple body terms if present
    energthreeb = f32(0.0)
    pressthreeb = 0.0
    if threebody_mod.nthreeb > 0:
        threebody_mod_lib.find_triplets(nmol, molecom.natom, xadd, yadd, zadd, xa, ya, za, mref, itype,
                                        molecom.ntype, molecom.type)
        for im in range(1, nmol + 1):
            energthreeb = energthreeb + threebody_mod_lib.en_mole_triplets(
                im, xadd, yadd, zadd, clustercom.phix, clustercom.phiy, clustercom.phiz, xa, ya, za, mref,
                clustercom.tethered)
            pressthreeb = pressthreeb + float(threebody_mod_lib.pr_mole_triplets(im, xadd, yadd, zadd, xa, ya, za,
                                                                                 mref))
        sumc = sumc + energthreeb
        presst = presst + pressthreeb
        energthreeb = energthreeb / f32(nmol)
        prthreebout = f32(-(100.0 * (pressthreeb / 3.0)) / float(r4(6.0222) * clustercom.volume))
    if coordconstraints_mod.ncoordconstraints > 0:
        for id_ in range(1, coordconstraints_mod.ncoordconstraints + 1):
            coordconstraints_mod.coordconstraints[id_] = coordconstraint_mod_lib.find_constraints(
                coordconstraints_mod.coordconstraints[id_], nmol, molecom.natom, naf, nal, xadd, yadd, zadd, xa,
                ya, za, mref, itype, molecom.ntype, molecom.type, molecom.cfrac, clustercom.dmax)
        for im in range(1, nmol + 1):
            sumc = sumc + coordconstraint_mod_lib.en_mole_coord(im, naf, nal, xadd, yadd, zadd, xa, ya, za, mref,
                                                                itype, molecom.qtype, npairref, molecom.amassr,
                                                                clustercom.tethered)
    # energy in kJ/mole
    energy = sumc / f32(nmol)
    # pressure in kbar
    pressc = f32(((float(f32(nmol) * molecom.simt) - presst / 3.0) * 100.0)
                 / float(r4(6.0222) * clustercom.volume))
    wsum = potential_par.wsum
    for if_ in range(1, epsr16.nfile + 1):
        energf[if_] = energf[if_] / (f32(nmol) * wsum[if_])
        pressf[if_] = ((-pressf[if_] / r4(3.0)) * r4(100.0)) / (r4(6.0222) * clustercom.volume * wsum[if_])

    return energy, pressc, energabs, energthreeb, prthreebout


# --------------------------------------------------------------------------- stog / lorch
@routine
def stog(rho, q, sofq, r, gofr, nqmin, nqmax, nr, npairs, nqdim, nrdim, npdim, xwt):
    """ato_manipulation_rev4.f:stog -- direct transform of sofq to r space, using a
    Lorch function to remove excess truncation oscillations.

    gofr(nrdim,npdim) is filled in place.
    """
    pi = PI4ATAN
    qstep = q[4] - q[3]
    sfac = r4(0.5) * qstep / pi / pi / rho
    sum_ = real_array(nqdim, npairs)
    sum_[:] = 0
    for j in range(1, nr + 1):
        rr = r[j]
        if rr > 0.0:
            # !$omp parallel do private(in,qmax,i,qq,qr,prod) shared(q,sum,nqmin,nqmax,pi)
            for in_ in range(1, npairs + 1):
                qmax = pi / q[nqmax[in_]]
                for i in range(int(nqmin[in_]), int(nqmax[in_]) + 1):
                    qq = q[i]
                    qr = qq * rr
                    prod = qq * sinf(qr) * lorch(qmax * qq, pi)
                    sum_[i, in_] = sofq[i, in_] * prod
            # !$omp end parallel do
            for in_ in range(1, npairs + 1):
                sumi = f32(0.0)
                for i in range(int(nqmin[in_]), int(nqmax[in_]) + 1):
                    sumi = sumi + sum_[i, in_]
                gofr[j, in_] = xwt[in_] + sfac * sumi / rr
        else:
            # !$omp parallel do private(in,qmax,i,qq,prod) shared(q,sum,nqmin,nqmax,pi)
            for in_ in range(1, npairs + 1):
                qmax = pi / q[nqmax[in_]]
                for i in range(int(nqmin[in_]), int(nqmax[in_]) + 1):
                    qq = q[i]
                    prod = qq * qq * lorch(qmax * qq, pi)
                    sum_[i, in_] = sofq[i, in_] * prod
            # !$omp end parallel do
            for in_ in range(1, npairs + 1):
                sumi = f32(0.0)
                for i in range(int(nqmin[in_]), int(nqmax[in_]) + 1):
                    sumi = sumi + sum_[i, in_]
                gofr[j, in_] = xwt[in_] + sfac * sumi
    return


@routine
def lorch(x, pi):
    """ato_manipulation_rev4.f:lorch -- Lorch window sin(x)/x for x < pi, else 0.
    Returns lorch.
    """
    if x == 0:
        result = f32(1.0)
    elif x < pi:
        result = sinf(x) / x
    else:
        result = f32(0.0)
    return result


# --------------------------------------------------------------------------- abs_energy_ep
@routine
def abs_energy_ep():
    """ato_manipulation_rev4.f:abs_energy_ep -- largest (max - min) range of the EP
    coefficients over all pairs.
    Returns abs_energy_ep.
    """
    #
    # Start at 0.5 Angstrom
    #
    ifirst = ifix(r4(0.5) / potential_par.pdstep) + 1  # noqa: F841 (unused in the Fortran too)
    #
    # absolute energy of EP
    #
    energabs = f32(0.0)
    # Sum the absolute values of the potentials directly
    #
    pcoefp = potential_par.pcoefp
    for ip in range(1, potential_par.npairspcof + 1):
        energabspairmin = f32(0.0)
        energabspairmax = f32(0.0)
        for i in range(2, potential_par.ncoeffp + 1):
            term = pcoefp[i, ip]
            energabspairmin = min(energabspairmin, term)
            energabspairmax = max(energabspairmax, term)
        units.write(6, fmt('(a,1x,i3,1x,a3,1x,a3,1x,f12.4)', 'abs_energy_ep>', ip,
                           Ch(molecom.type[molecom.n1extr[ip]], 3), Ch(molecom.type[molecom.n2extr[ip]], 3),
                           energabspairmax - energabspairmin))
        energabs = max(energabs, (energabspairmax - energabspairmin))
    result = energabs

    return result


@routine
def abs_energy_ep_old():
    """ato_manipulation_rev4.f:abs_energy_ep_old -- old version of abs_energy_ep (xavs).

    NB: the Fortran assigns ``abs_energy_ep=energabs`` (an implicit local variable here),
    so the function result abs_energy_ep_old is never set (garbage).  0.0 is returned.
    Returns abs_energy_ep_old.
    """
    #
    # absolute energy of EP
    #
    energabs = f32(0.0)
    # Sum the absolute values of the potentials directly
    #
    xavs = potential_par.xavs
    for ip in range(1, molecom.npairs + 1):
        irefp = potential_par.irefpcof[ip]  # noqa: F841
        # cfracprod=cfrac(n1extr(ip))*cfrac(n2extr(ip))
        energabspairmin = f32(0.0)
        energabspairmax = f32(0.0)
        for i in range(1, epsr16.npt1 + 1):
            # term=ycoh(i)*cfracprod*abs(xavs(i,irefp))
            # term=cfracprod*abs(xavs(i,ip))
            # term=ycoh(i)*xavs(i,ip)
            term = xavs[i, ip]
            energabspairmin = min(energabspairmin, term)
            energabspairmax = max(energabspairmax, term)
        units.write(6, fmt('(a,1x,i3,1x,a3,1x,a3,1x,f12.4)', 'abs_energy_ep>', ip,
                           Ch(molecom.type[molecom.n1extr[ip]], 3), Ch(molecom.type[molecom.n2extr[ip]], 3),
                           energabspairmax - energabspairmin))
        energabs = max(energabs, (energabspairmax - energabspairmin))
    abs_energy_ep_local = energabs  # noqa: F841 (abs_energy_ep=energabs: a local in the Fortran)

    return f32(0.0)


# --------------------------------------------------------------------------- adjust_intra_wts
@routine
def adjust_intra_wts():
    """ato_manipulation_rev4.f:adjust_intra_wts -- set the bond-angle and dihedral
    energy weights (in this version simply reset to the requested values).
    """
    # Adjust the bond angle energy weight
    # First we try to set the intra-molecular weights to give the required energy
    # ratios between vibrational, angle and dihedral
    sumevib = f32(0.0)
    sumeang = f32(0.0)
    sumedih = f32(0.0)
    molecom.tempang = f32(1.0)
    molecom.tempdih = f32(10.0)
    for im in range(1, clustercom.nmol + 1):
        # Get the intramolecular energy of this molecule
        # NB: SIGNATURES.md lists only sumedih as returned by en_mole_intra, although
        # sumevib and sumeang are accumulated too (docs/SIGNATURE_EXCEPTIONS.md).
        energy, sumevib, sumeang, sumedih = molfit.en_mole_intra(im, sumevib, sumeang, sumedih)
    units.write(6, lw('adjust_intra_wts> Vibrational, bond angle and dihedral energies: ', sumevib, sumeang,
                      sumedih))
    if sumevib > 0.0:
        ratioang = sumeang / sumevib
        ratiodih = sumedih / sumevib
    else:
        ratioang = f32(1.0)
        ratiodih = f32(1.0)
    units.write(6, lw('adjust_intra_wts> Relative bond angle and dihedral energies: ', ratioang, ratiodih))
    # (the commented-out relative-weight adjustment of the Fortran is omitted)
    molecom.tempangwt = molecom.tempangreq
    molecom.tempang = molecom.tempangwt
    molecom.tempdihwt = molecom.tempdihreq
    molecom.tempdih = molecom.tempdihwt * molecom.tempb
    units.write(6, lw('adjust_intra_wts> Bond angle and dihedral energy weights: ', molecom.tempang,
                      molecom.tempdih))
    for imoltype in range(1, clustercom.nmoltype + 1):
        units.write(6, lw('adjust_intra_wts> Internal step sizes: ', imoltype, clustercom.interntrans[imoltype],
                          clustercom.internrot[imoltype]))
    return


# --------------------------------------------------------------------------- get_maxmolsize
@routine
def get_maxmolsize():
    """ato_manipulation_rev4.f:get_maxmolsize -- twice the largest molecular diameter."""
    molecom.maxmolsize = f32(0.0)
    for im in range(1, clustercom.nmol + 1):
        if clustercom.moldiam[im] > molecom.maxmolsize:
            molecom.maxmolsize = clustercom.moldiam[im]
    molecom.maxmolsize = r4(2.0) * molecom.maxmolsize
