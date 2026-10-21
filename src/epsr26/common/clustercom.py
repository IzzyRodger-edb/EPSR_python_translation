"""Generated from fortran_src/clustercom.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``clustercom.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401
from .epsr_dimension import mcomp, mmol  # noqa: F401

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of clustercom.inc."""
    global moltype, moltemplate, nmol, nato, ngrot, naf, nal, nrf
    global nrl, mtype, nmoltype, nmovei, ninter, nwhole, nmovem, nmoveit
    global nintert, nwholet, nmovemt, nstepreject, tetheredatom, moltetheredatom, tethered, moltethered
    global xadd, yadd, zadd, phix, phiy, phiz, dmax, volume
    global amasst, moldiam, amasstsqrt, interntrans, internrot, wholerot, wholetrans
    # common /clusterc/
    moltype = char_array(mcomp)
    moltemplate = char_array(mcomp)
    # common /clusteri/
    nmol = 0
    nato = int_array(mmol)
    ngrot = int_array(mmol)
    naf = int_array(mmol)
    nal = int_array(mmol)
    nrf = int_array(mmol)
    nrl = int_array(mmol)
    mtype = int_array(mmol)
    nmoltype = 0
    nmovei = int_array(2, mcomp)
    ninter = int_array(2, mcomp)
    nwhole = int_array(2, mcomp)
    nmovem = int_array(2, mcomp)
    nmoveit = int_array(2)
    nintert = int_array(2)
    nwholet = int_array(2)
    nmovemt = int_array(2)
    nstepreject = int_array(mmol)
    tetheredatom = int_array(mmol)
    moltetheredatom = int_array(mcomp)
    # common /clusterl/
    tethered = logical_array(mmol)
    moltethered = logical_array(mcomp)
    # common /clusterr/
    xadd = real_array(mmol)
    yadd = real_array(mmol)
    zadd = real_array(mmol)
    phix = real_array(mmol)
    phiy = real_array(mmol)
    phiz = real_array(mmol)
    dmax = f32(0.0)
    volume = f32(0.0)
    amasst = real_array(mmol)
    moldiam = real_array(mmol)
    amasstsqrt = real_array(mmol)
    interntrans = real_array(mcomp)
    internrot = real_array(mcomp)
    wholerot = real_array(mcomp)
    wholetrans = real_array(mcomp)


reset()
