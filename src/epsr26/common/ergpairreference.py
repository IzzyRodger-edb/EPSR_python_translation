"""Data of Fortran ``module ergpairreference`` (fortran_src/ergpairreference.f).

Hand-written (TRANSLATION_GUIDE 5.2).  The routines of the module are in
``library/ergpairreference.py``.  Module integers are static in gfortran and start at
0; allocatable arrays start unallocated (``None``).
"""


def reset():
    """(Re)create the module variables of ergpairreference."""
    global ergpair, nergpair, nergpairmax, molergpair, ergpairref, nergpairm
    global nnbrefmax, ergold, ergnew, nbref, nnbref, nmolhere
    # Stores and manipulates the arrays to store and reference the set of pair energies
    # real, dimension(:), allocatable :: ergpair
    ergpair = None
    # integer :: nergpair,nergpairmax
    nergpair = 0
    nergpairmax = 0
    # integer, dimension(:,:), allocatable :: molergpair,ergpairref
    molergpair = None
    ergpairref = None
    # integer, dimension(:), allocatable :: nergpairm
    nergpairm = None
    # integer :: nnbrefmax
    nnbrefmax = 0
    # real, dimension(:), allocatable :: ergold,ergnew
    ergold = None
    ergnew = None
    # integer, dimension(:), allocatable :: nbref
    nbref = None
    # integer :: nnbref,nmolhere
    nnbref = 0
    nmolhere = 0


from . import register  # noqa: E402
import sys  # noqa: E402

register(sys.modules[__name__])
reset()
