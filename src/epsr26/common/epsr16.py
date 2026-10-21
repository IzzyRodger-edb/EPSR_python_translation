"""Generated from fortran_src/epsr16.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``epsr16.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401
from .epsr_dimension import mfile, mpairs, mpts, mpz  # noqa: F401

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of epsr16.inc."""
    global x, xstar, y, ycoh, entwt, foqols, xsav, ysav
    global xavbar, rlog, soqinter, sq2inter, soqintra, sq2intra, add, cellst
    global sigma1, sigma2, addcon, rhofac, xwt, cell, datfc2, sm1sav
    global rmaxcell, cellmax, upol, rho, amubar, amubar2, rejrate, qminm
    global qminp, gorinter, gr2inter, gorintra, gr2intra, energyref, energyemp, energytot
    global npt, npt1, npt2, nmove, nfile, ndist, nqmin, nqmax
    global nqpts, nptmin, nptmax, entsig
    # common /dista/
    x = real_array(mpz)
    xstar = real_array(mpz)
    y = real_array(mpz)
    ycoh = real_array(mpz)
    entwt = real_array(mpz)
    foqols = real_array(mpz, mpts)
    xsav = real_array(mpz)
    ysav = real_array(mpz)
    xavbar = f32(0.0)
    rlog = f32(0.0)
    soqinter = real_array(mpts, mpairs)
    sq2inter = real_array(mpts, mpairs)
    soqintra = real_array(mpts, mpairs)
    sq2intra = real_array(mpts, mpairs)
    add = f32(0.0)
    cellst = f32(0.0)
    sigma1 = real_array(mfile)
    sigma2 = real_array(mfile)
    addcon = real_array(mfile)
    rhofac = real_array(mfile)
    xwt = real_array(mfile)
    cell = f32(0.0)
    datfc2 = f32(0.0)
    sm1sav = real_array(mpz)
    rmaxcell = f32(0.0)
    cellmax = f32(0.0)
    upol = f32(0.0)
    rho = f32(0.0)
    amubar = f32(0.0)
    amubar2 = f32(0.0)
    rejrate = f32(0.0)
    qminm = f32(0.0)
    qminp = f32(0.0)
    # common /xaver/
    gorinter = real_array(mpz, mpairs)
    gr2inter = real_array(mpz, mpairs)
    gorintra = real_array(mpz, mpairs)
    gr2intra = real_array(mpz, mpairs)
    energyref = real_array(mpairs)
    energyemp = real_array(mpairs)
    energytot = real_array(mpairs)
    # common /idista/
    npt = 0
    npt1 = 0
    npt2 = 0
    nmove = 0
    nfile = 0
    ndist = 0
    nqmin = int_array(mfile)
    nqmax = int_array(mfile)
    nqpts = int_array(mpts)
    nptmin = 0
    nptmax = 0
    # common /xaver8/
    entsig = real_array(mpz, mpairs)


reset()
