"""Fortran module data of ``module molecule_coordinates_mod``
(fortran_src/molecule_coordinates_mod.f90).

Hand-written (TRANSLATION_GUIDE 5.2).  The routine of the module
(``get_coords_from_xyz``) is in ``library/molecule_coordinates_mod.py``.

Allocatable arrays start as ``None``; uninitialised module scalars are 0 / ''
(gfortran puts module variables in zero-initialised static storage).
"""

import sys

from . import register


def reset():
    """(Re)create the module variables of molecule_coordinates_mod."""
    global chsymbolt, centralmolecule, rt, phit, thetat
    global nofatoms, nread, lenatomsequence, xyzfilename, atomsequence
    chsymbolt = None  # character(len=3)(:)
    centralmolecule = None  # integer(:)
    rt = None  # real(:)
    phit = None  # real(:)
    thetat = None  # real(:)
    nofatoms = 0
    nread = 0
    lenatomsequence = 0
    xyzfilename = ''  # character(len=256)
    atomsequence = ''  # character(len=256)


register(sys.modules[__name__])
reset()
