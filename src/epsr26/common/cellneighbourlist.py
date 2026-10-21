"""Data of Fortran ``module cellneighbourlist`` (fortran_src/cellneighbourlist.f).

Hand-written (TRANSLATION_GUIDE 5.2).  The routines of the module are in
``library/cellneighbourlist.py``.  Module integers are static in gfortran and start at
0; allocatable arrays start unallocated (``None``).
"""


def reset():
    """(Re)create the module variables of cellneighbourlist."""
    global ncelldiva, ncelldivb, ncelldivc, ncelldivsq, ncell, ncellneighb, maxcellno, nmolcell
    global nmolsincell, molcellno, iacellneighb, ibcellneighb, iccellneighb, molsincell
    # integer :: ncelldiva,ncelldivb,ncelldivc,ncelldivsq,ncell,ncellneighb,maxcellno,nmolcell
    ncelldiva = 0
    ncelldivb = 0
    ncelldivc = 0
    ncelldivsq = 0
    ncell = 0
    ncellneighb = 0
    maxcellno = 0
    nmolcell = 0
    # integer, dimension(:), allocatable :: nmolsincell,molcellno,iacellneighb,ibcellneighb,iccellneighb
    nmolsincell = None
    molcellno = None
    iacellneighb = None
    ibcellneighb = None
    iccellneighb = None
    # integer, dimension(:,:), allocatable :: molsincell
    molsincell = None


from . import register  # noqa: E402
import sys  # noqa: E402

register(sys.modules[__name__])
reset()
