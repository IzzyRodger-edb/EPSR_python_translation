"""Generated from fortran_src/potential_par.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``potential_par.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401
from .epsr_dimension import mfile, mpairs, mpz  # noqa: F401

mcoeff = 1000
mpgauss = 5

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of potential_par.inc."""
    global tlab1, tlab2, addpotential, ncoeffp, npcoeff, npits, npitss, npfreq
    global nzero, nnonzero, npairspcof, irefpcof, nergref, nneighb, npgauss, pdecay
    global pstepd, pstepp, pwt, pwtold, expgr, expsq, pdmin, pdstep
    global pdmax, paccept, psigma2, psigma2read, psigmafac, pcoefd, pcoefp, rminex
    global coordex, rmaxex, pcoefex, pcoefexmax, potex, potexmax, weight, uwt
    global uwait, xav, xav1, wsum, xavs, xavs1, wfac, wfacc
    global epot, epot1, ecoul, ecoul1, xavtot, xavtot1, feedback, feedbackold
    global efacm, efacold, efacp, rpgauss, wpgauss, hpgauss, absolute_energy, potfacamp
    global potfac, refpotfac, ereqstep, erequnit, coulfac, gaussian, harmonic
    # common /potential_c/
    tlab1 = char_array(mpairs)
    tlab2 = char_array(mpairs)
    addpotential = ''
    # common /potential_i/
    ncoeffp = 0
    npcoeff = int_array(mfile)
    npits = 0
    npitss = 0
    npfreq = 0
    nzero = 0
    nnonzero = 0
    npairspcof = 0
    irefpcof = int_array(mpairs)
    nergref = 0
    nneighb = 0
    npgauss = int_array(mpairs)
    # common /potential_r/
    pdecay = real_array(mfile)
    pstepd = real_array(mfile)
    pstepp = real_array(mfile)
    pwt = real_array(mcoeff)
    pwtold = real_array(mcoeff)
    expgr = real_array(3000000)
    expsq = real_array(3000000)
    pdmin = f32(0.0)
    pdstep = f32(0.0)
    pdmax = f32(0.0)
    paccept = f32(0.0)
    psigma2 = f32(0.0)
    psigma2read = f32(0.0)
    psigmafac = f32(0.0)
    pcoefd = real_array(mcoeff, mfile)
    pcoefp = real_array(mcoeff, mpairs)
    rminex = real_array(mpairs)
    coordex = real_array(mpairs)
    rmaxex = real_array(mpairs)
    pcoefex = real_array(mpairs)
    pcoefexmax = real_array(mpairs)
    potex = real_array(mpairs)
    potexmax = real_array(mpairs)
    weight = real_array(2, mfile, mpairs)
    uwt = real_array(mpairs)
    uwait = real_array(2, 850, mpairs)
    xav = real_array(mpz, mpairs)
    xav1 = real_array(mpz, mpairs)
    wsum = real_array(mfile)
    xavs = real_array(mpz, mpairs)
    xavs1 = real_array(mpz, mpairs)
    wfac = real_array(mpz)
    wfacc = real_array(mpz)
    epot = real_array(mpz, mpairs)
    epot1 = real_array(mpz, mpairs)
    ecoul = real_array(mpz, mpairs)
    ecoul1 = real_array(mpz, mpairs)
    xavtot = real_array(mpz, mpairs)
    xavtot1 = real_array(mpz, mpairs)
    feedback = f32(0.0)
    feedbackold = f32(0.0)
    efacm = real_array(mfile)
    efacold = real_array(mfile)
    efacp = real_array(mpairs)
    rpgauss = real_array(mpgauss, mpairs)
    wpgauss = real_array(mpgauss, mpairs)
    hpgauss = real_array(mpgauss, mpairs)
    absolute_energy = f32(0.0)
    potfacamp = f32(0.0)
    potfac = f32(0.0)
    refpotfac = f32(0.0)
    ereqstep = f32(0.0)
    erequnit = f32(0.0)
    coulfac = real_array(mpairs)
    # common /potential_l/
    gaussian = False
    harmonic = False


reset()
