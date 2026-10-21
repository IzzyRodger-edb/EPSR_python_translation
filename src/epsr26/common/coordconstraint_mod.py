"""Fortran module ``coordconstraint_mod`` (fortran_src/coordconstraint_mod.f90): the
derived type ``type(coordconstraint)``.

The module declares no variables, only the type; its routines (``init_constraint``,
``find_constraints``, ``find_type``, ``en_coordcon``, ``en_coordcon_total``) are in
``library/coordconstraint_mod.py``.

``coordconstraint()`` creates one object with every component in the state gfortran
gives a fresh ``type(coordconstraint)`` variable: allocatable components unallocated
(``None``), the others uninitialised in the Fortran (here 0, '' or .false.).  Arrays
follow the usual convention (index 0 unused; ``iabondcc[j, ia]`` is ``iabondcc(j,ia)``).
"""

from ..fortran.fmath import f32


class coordconstraint:
    """``type coordconstraint`` (coordconstraint_mod.f90:8-32)."""

    def __init__(self):
        self.atc = ''  # character(len=3): atom type at the centre of the constraint
        self.ats = None  # character(len=3), allocatable(:): possible secondary atom types
        self.nc = 0  # expected overall coordination number
        self.alternatenc = 0  # alternate coordination number for the coordinating species
        self.cctemp = f32(0.0)  # temperature factor associated with this constraint
        self.ccgaussweight = f32(0.0)
        self.ccgausswidth = f32(0.0)
        self.ccdistort = f32(0.0)  # distortion on coordination potential
        self.ns = 0  # number of secondary atoms associated with this constraint
        self.atclistsize = 0
        self.atclist = None  # integer, allocatable(:): atom numbers of type atc
        self.nbondcc = None  # integer, allocatable(:)
        self.nearnbondcc = None  # integer, allocatable(:)
        self.expecnbondcc = None  # integer, allocatable(:)
        self.iabondcc = None  # integer, allocatable(:,:): secondary atoms for each atom
        self.wtcc = None  # integer, allocatable(:,:)
        self.ds = None  # real, allocatable(:): expected bond distances for types ats
        self.dsmax = None  # real, allocatable(:)
        self.tbondcc = None  # real, allocatable(:,:): expected bonds for each neighbour
        self.adistcc = None  # real, allocatable(:,:): actual distances for each neighbour
        self.hastb = False  # calculation of tb terms is to be performed
        self.hastbintra = False
        self.tbweight = f32(0.0)
        self.tbgaussweight = f32(0.0)
        self.tbintraweight = f32(0.0)
        self.tbintragaussweight = f32(0.0)
        self.tbdistort = f32(0.0)
        self.tbgausswidth = f32(0.0)
        self.tbintradistort = f32(0.0)
        self.tbintragausswidth = f32(0.0)
        self.ntbangles = 0
        self.ntbintraangles = 0
        self.tbangles = None  # real, allocatable(:)
        self.tbcosines = None  # real, allocatable(:)
        self.tbintraangles = None  # real, allocatable(:)
        self.tbintracosines = None  # real, allocatable(:)
        self.tbintraats = None  # character(len=3), allocatable(:): secondary types for intra angles

    def __repr__(self):
        return f"coordconstraint(atc={self.atc!r}, nc={self.nc}, ns={self.ns})"


def reset():
    """The module has no variables."""
    pass


from . import register  # noqa: E402
import sys  # noqa: E402

register(sys.modules[__name__])
reset()
