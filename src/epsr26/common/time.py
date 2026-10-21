"""``common /time/`` declared inside ``get_date_and_time`` (fortran_src/cellneighbourlist.f).

Hand-written (TRANSLATION_GUIDE 5.1, a COMMON written directly inside a routine).
Imported as ``from ..common import time as time_common`` (the name is the COMMON block
name; it does not shadow the standard library because imports are absolute).
"""


def reset():
    """(Re)create the COMMON storage of /time/."""
    global startyear, startabsolutetime, startmusecs
    # common /time/startyear,startabsolutetime,startmusecs   (all integer)
    startyear = 0
    startabsolutetime = 0
    startmusecs = 0


from . import register  # noqa: E402
import sys  # noqa: E402

register(sys.modules[__name__])
reset()
