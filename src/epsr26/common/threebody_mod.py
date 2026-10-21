"""Module data of ``module threebody_mod`` (fortran_src/threebody_mod.f90).

Hand-written (TRANSLATION_GUIDE 5.2).  The routines of the module are in
``library/threebody_mod.py``.  Allocatable arrays start as ``None``; module
integers are static in gfortran, i.e. zero before first assignment.
"""


def reset():
    """(Re)create the module variables of threebody_mod."""
    global nthreeb, nthreebangsmax, nthreebneighbmolmax
    global threebgrp, threebamp, threebrmax, threebrwidth, threebawidth
    global nthreebangs, threebangs, threebamps
    global nthreebneighbmol, threebatom, threebneighbmol, threebatomneighb, threebsamemol
    global ltyp1, ltyp2, ltyp3, atom1list, atom2list, atom3list, natom1list, natom2list, natom3list
    global ntripletsinmol, tripletlist, firsttripletinmol, ntripletsinmolmax
    global ergold_triplet, ergnew_triplet, sume_triplet_dist, en_mole_triplet_dist, consthreeb
    # Number of three-body groups, maximum number of three-body atoms and angles
    nthreeb = 0
    nthreebangsmax = 0
    nthreebneighbmolmax = 0
    # character(len=3), dimension(:,:): Three atom types to form the group (3 per group).
    threebgrp = None
    # Three-body potential amplitude
    threebamp = None
    # Maximum distances for atom pairs (2 per group). Initially the second will be set equal to the first.
    threebrmax = None
    # Nominal width of three body potential function in radius and angel (deg.)
    threebrwidth = None
    threebawidth = None
    # Number of angles in group
    nthreebangs = None
    # List of requested angles and amplitudes in group.
    threebangs = None
    threebamps = None

    # Number of neighbouring molecules, for which three-body energies can be defined
    nthreebneighbmol = None
    # List of atom numbers that can form a three-body potential
    threebatom = None
    # List of neighbours for each molecule in a triplet
    threebneighbmol = None
    # List of acceptable neighbours of each atom in threebatom
    threebatomneighb = None
    # Array of 3 1s or 0s to say whether atoms 1 and 2, 2 and 3, and 1 and 3 are on same molecule.
    # 1 means they are on the same molecule.
    threebsamemol = None
    # The following arrays are used to make lists of atoms of the specified types
    ltyp1 = 0
    ltyp2 = 0
    ltyp3 = 0
    atom1list = None
    atom2list = None
    atom3list = None
    natom1list = 0
    natom2list = 0
    natom3list = 0
    # The following arrays list the number of triplets for each molecule, the list of triplets,
    # and the first entry in this list for each molecule
    ntripletsinmol = None
    tripletlist = None
    firsttripletinmol = None
    ntripletsinmolmax = 0
    # The following arrays store the triplet energies of each molecule
    ergold_triplet = None
    ergnew_triplet = None
    sume_triplet_dist = None
    en_mole_triplet_dist = None
    # Stores the constant for the potential for each triplet distribution
    consthreeb = None


from . import register  # noqa: E402
import sys  # noqa: E402

register(sys.modules[__name__])
reset()
