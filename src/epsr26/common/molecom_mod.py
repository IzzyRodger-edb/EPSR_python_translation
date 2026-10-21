"""Fortran module data of ``module molecom_mod`` (fortran_src/molecom_mod.f90).

Hand-written (TRANSLATION_GUIDE 5.2).  The routines of the module
(``testsize_rotgrps``, ``testsize_labrot``) are in ``library/molecom_mod.py``.

Allocatable arrays start as ``None``; uninitialised module scalars are 0 / ``.false.``
(gfortran puts module variables in zero-initialised static storage).
"""

import sys

from ..fortran.fmath import f32, r4
from . import register


def reset():
    """(Re)create the module variables of molecom_mod."""
    global ngrrot, ngrrotmax, rotgrp, nrot, labaxi, rotf, rotl, mrot
    global nrottot, nrottotmax, labrot
    global coarsegraining, qscaleinputfactor, selfbroadpwr, sphereformfactor
    global intraformfactor, intraformfactorself, scalefacmol, scalefacmoluse, rfractmol
    global natomsphere, rhosphere, radminsphere, radmaxsphere
    global scalefactypemin, scalefactypemax, scalefacratio, nsffbintype, assignatomscalefac
    global maxnscalefactype, ntermssff, fracntermssff, wavinter, prodssff
    global ntermsintra, ntermsintraself, pseudofacmol, natomcg, natomrcg, rhocg
    ngrrot = 0  # Number of rotational groups
    ngrrotmax = 0  # and current maximum number
    rotgrp = None  # character(len=3)(:) Rotational group label (currently only ROT or DIH)
    nrot = None  # integer(:) Number of rotatable atoms in each rotational group
    labaxi = None  # integer(:,:) Pair of atoms about which rotations are made
    rotf = None  # integer(:) First in list of rotatable atoms for each group
    rotl = None  # integer(:) Last in list of rotatable atoms for each group
    mrot = None  # integer(:) Molecule number for each rotational group

    nrottot = 0  # Total number of rotatable atoms in all groups
    nrottotmax = 0
    labrot = None  # integer(:) List of rotatable atoms for all groups
    # Coarse-graining variables
    coarsegraining = False  # True if coarse graining is in place
    qscaleinputfactor = f32(0.0)  # Scale factor on input Q scale - needed for coarse graining
    selfbroadpwr = f32(0.0)  # Broadening power for self terms - attempts to mimic molecular structure
    sphereformfactor = None  # real(:,:,:) Used if coarse graining is on.
    intraformfactor = None  # real(:,:) Used if coarse graining is on.
    intraformfactorself = None  # real(:,:)
    scalefacmol = None  # real(:) Scale factor assigned to size of each molecule
    scalefacmoluse = None  # real(:)
    rfractmol = None  # real(:) Relative broadening for each molecule type.
    natomsphere = None  # real(:) Atomic number density, and minimum and maximum radii for coarse graining
    rhosphere = None  # real(:)
    radminsphere = None  # real(:)
    radmaxsphere = None  # real(:)
    scalefactypemin = None  # real(:) Maximum and minimum scale factor for each molecule type
    scalefactypemax = None  # real(:)
    scalefacratio = r4(0.2)  # Used to determine the number of distinct spherical form factors for a given molecule type
    nsffbintype = None  # integer(:) Number of spherical form factor bins for a given atom type
    assignatomscalefac = None  # integer(:) Assigned spherical form factor bin for each atom
    maxnscalefactype = 0  # Maximum number of distinct spherical form factors
    ntermssff = None  # integer(:,:) No. of terms for each sphere form factor bin and each atom type
    fracntermssff = None  # real(:,:) Fractional no. of terms for each sphere form factor bin and each atom type
    wavinter = None  # real(:,:,:,:) Used for distinguishing spherical form factors of different scale factors when coarse graining
    prodssff = None  # real(:,:) Temporary form factor product used by gensq
    ntermsintra = None  # integer(:) No. of terms for intra correlations
    ntermsintraself = None  # integer(:) No. of terms for intra self correlations
    pseudofacmol = None  # real(:) No. of real atoms to be used for each pseudo atom. Normally should be 1.0, but increase for coarse graining.
    natomcg = f32(0.0)  # Counts total number of coarse grain atoms
    natomrcg = f32(0.0)
    rhocg = f32(0.0)


register(sys.modules[__name__])
reset()
