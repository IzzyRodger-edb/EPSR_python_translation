"""Fortran module ``coordconstraints_mod`` (fortran_src/coordconstraint_mod.f90:1243-1453):
module data.

* ``coordconstraints`` -- ``type(coordconstraint), dimension(:), allocatable``: ``None``
  while unallocated, otherwise a Python **list** of ``coordconstraint`` objects with
  index 0 unused (``coordconstraints[ic]`` is ``coordconstraints(ic)``; ``size`` is
  ``len(coordconstraints) - 1``).  A list (not a NumPy array) so that
  ``reallocate1d_cc`` can resize it in place.
* ``ncoordconstraints`` -- integer (static storage, 0 at program start).
* ``ergold_coord``, ``ergnew_coord`` -- ``real, dimension(:), allocatable``.

The routines of the module are in ``library/coordconstraint_mod.py``.
"""


def reset():
    """(Re)create the variables of module coordconstraints_mod."""
    global coordconstraints, ncoordconstraints, ergold_coord, ergnew_coord
    coordconstraints = None
    ncoordconstraints = 0
    ergold_coord = None
    ergnew_coord = None


from . import register  # noqa: E402
import sys  # noqa: E402

register(sys.modules[__name__])
reset()
