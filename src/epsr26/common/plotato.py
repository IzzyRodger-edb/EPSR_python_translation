"""Generated from fortran_src/plotato.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``plotato.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401
from .epsr_dimension import matom, mcomp, mpairs  # noqa: F401

mclass = 120

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of plotato.inc."""
    global fnamesplot, classlab, classlab1, classlab2, cbstr1, cbstr2, classrad, classrgb
    global classbond1, classbond2, bondrad, xat, yat, zat, unit, posfrac
    global ptscaler, scaleatoms, xangle, yangle, dstr1, dstr2, radstr, bondrgb
    global boxrgb, boxlwidth, nclassread, classtype, classcolour, classcoloura, nclassbond, atomref
    global nclplt, nplt, iplot, ilenfnamesplot, ibstr1, ibstr2, nclass, iplotclass
    global ntypplot, itypplot, nmolcentre, appendxyz
    # common /plotato_c/
    fnamesplot = char_array((0, mcomp))
    classlab = char_array(mclass)
    classlab1 = char_array(mpairs)
    classlab2 = char_array(mpairs)
    cbstr1 = char_array(mpairs)
    cbstr2 = char_array(mpairs)
    # common /plotato_r/
    classrad = real_array(mclass)
    classrgb = real_array(3, mclass)
    classbond1 = real_array(mpairs)
    classbond2 = real_array(mpairs)
    bondrad = real_array(mpairs)
    xat = real_array(matom)
    yat = real_array(matom)
    zat = real_array(matom)
    unit = f32(0.0)
    posfrac = f32(0.0)
    ptscaler = f32(0.0)
    scaleatoms = f32(0.0)
    xangle = f32(0.0)
    yangle = f32(0.0)
    dstr1 = real_array(mpairs)
    dstr2 = real_array(mpairs)
    radstr = real_array(mpairs)
    bondrgb = real_array(3)
    boxrgb = real_array(3)
    boxlwidth = f32(0.0)
    # common /plotato_i/
    nclassread = 0
    classtype = int_array(mclass)
    classcolour = int_array(mclass)
    classcoloura = int_array(mclass)
    nclassbond = 0
    atomref = int_array(matom)
    nclplt = 0
    nplt = int_array(mcomp)
    iplot = int_array(matom)
    ilenfnamesplot = int_array((0, mcomp))
    ibstr1 = int_array(mpairs)
    ibstr2 = int_array(mpairs)
    nclass = 0
    iplotclass = int_array(mcomp)
    ntypplot = 0
    itypplot = int_array(mcomp)
    nmolcentre = 0
    # common /plotato_l/
    appendxyz = False


reset()
