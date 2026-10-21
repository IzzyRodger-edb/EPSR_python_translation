"""Module data of ``module containerforces_mod`` (fortran_src/containerforces_mod.f90).

Hand-written (TRANSLATION_GUIDE 5.2).  The routines of the module are in
``library/containerforces_mod.py``.  Allocatable arrays start as ``None``; module
integers are static in gfortran, i.e. zero before first assignment.
"""


def reset():
    """(Re)create the module variables of containerforces_mod."""
    global numcontainerforces, numpairs, numcentraltype, numsecondarytype
    global pairrefstore, centraltype, secondarytype, ergold_container, ergnew_container
    numcontainerforces = 0
    # integer, dimension(:), allocatable :: numpairs,numcentraltype,numsecondarytype
    numpairs = None
    numcentraltype = None
    numsecondarytype = None
    # integer, dimension(:,:), allocatable :: pairrefstore
    pairrefstore = None
    # List of atoms (q-type) at centres of "containers"
    centraltype = None
    # List of atom types to be contained (not q-type)
    secondarytype = None
    # real, dimension(:), allocatable :: ergold_container,ergnew_container
    ergold_container = None
    ergnew_container = None


from . import register  # noqa: E402
import sys  # noqa: E402

register(sys.modules[__name__])
reset()
