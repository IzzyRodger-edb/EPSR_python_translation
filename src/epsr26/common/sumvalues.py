"""COMMON /sumvalues/ (hand-written; TRANSLATION_GUIDE 5.1).

Declared directly inside three routines of fortran_src/epsr_standalone_rev3.f
(``runepsr17``, ``sumgr`` and ``sumsq``) as::

      real*8 fc1,fc2,add1(2),add2(2),sumtot1(2),sumtot2(2),errort(2)
      common/sumvalues/fc1,fc2,add1,add2,sumtot1,sumtot2,errort

``runepsr17`` sets the accumulation factors fc1, fc2 before each call of sumsq/sumgr,
which use them (and the other members as work space).  Access as ``sumvalues.fc1``.
"""

import sys

from ..fortran.arrays import dble_array


def reset():
    """(Re)create the COMMON storage of /sumvalues/."""
    global fc1, fc2, add1, add2, sumtot1, sumtot2, errort
    fc1 = 0.0
    fc2 = 0.0
    add1 = dble_array(2)
    add2 = dble_array(2)
    sumtot1 = dble_array(2)
    sumtot2 = dble_array(2)
    errort = dble_array(2)


from . import register  # noqa: E402

register(sys.modules[__name__])
reset()
