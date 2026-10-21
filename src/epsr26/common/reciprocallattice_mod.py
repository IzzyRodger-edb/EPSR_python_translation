"""Module data of ``module reciprocallattice_mod`` (fortran_src/ato_manipulation_rev4.f,
lines 31-36)."""

import sys

from ..fortran.arrays import real_array
from ..fortran.fmath import f32
from . import register


def reset():
    """(Re)create the module variables of ``module reciprocallattice_mod``."""
    global recipra, reciprb, reciprc, rbragg, volrecipr, ampra, amprb, amprc
    # real, dimension(3) :: recipra,reciprb,reciprc,rbragg !Components of the reciprocal lattice vector
    recipra = real_array(3)
    reciprb = real_array(3)
    reciprc = real_array(3)
    rbragg = real_array(3)
    # real :: volrecipr,ampra,amprb,amprc
    volrecipr = f32(0.0)
    ampra = f32(0.0)
    amprb = f32(0.0)
    amprc = f32(0.0)


register(sys.modules[__name__])
reset()
