"""Data of Fortran ``module ergpairstore`` (fortran_src/ergpairreference.f).

Hand-written (TRANSLATION_GUIDE 5.2).  The routines of the module are in
``library/ergpairreference.py``.  (The module is not used by any linked program.)
"""


def reset():
    """(Re)create the module variables of ergpairstore."""
    global ergold, ergnew, nbref, nnbref
    # real, dimension(:), allocatable :: ergold,ergnew
    ergold = None
    ergnew = None
    # integer, dimension(:), allocatable :: nbref
    nbref = None
    # integer :: nnbref
    nnbref = 0


from . import register  # noqa: E402
import sys  # noqa: E402

register(sys.modules[__name__])
reset()
