"""Module data of ``module timer`` (fortran_src/ato_manipulation_rev4.f, lines 1-18).

The module routine ``initialise_timer`` is in ``library/ato_manipulation_rev4.py``.
"""

import sys

from . import register


def reset():
    """(Re)create the module variables of ``module timer``."""
    global sumtime1, sumtime2, time1, time2
    # integer(kind=8) :: sumtime1,sumtime2,time1,time2
    sumtime1 = 0
    sumtime2 = 0
    time1 = 0
    time2 = 0


register(sys.modules[__name__])
reset()
