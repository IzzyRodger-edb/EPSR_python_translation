"""Fortran module data of ``module sharm_mod`` (fortran_src/sharm_routines_rev1.f, lines 1-30).

Hand-written (TRANSLATION_GUIDE 5.2).  The routines of sharm_routines_rev1.f
(``get_l_and_n_values``, ``form_cgcoeffs``, ``form_dlmn``, ...) are in
``library/sharm_routines_rev1.py``.  Used by the programs sharm and sdf.
"""

import sys

from . import register


def reset():
    """(Re)create the module variables of sharm_mod."""
    global nll1v, nll2v, nllv, n1step, n2step, ll1v, ll2v, llv
    global nlmax, ndist, ldistinct, ll1, ll2, ll, n1, n2, alsec
    global mdlmn, mz, mdist, nm, il1rf, il2rf, ilrf, illrf
    global phim, cthetm, chim, xcent, ycent, zcent, philab, dlmnc, dlmnmc, chiphi
    global mcg, termcg
    global phimt1, cthetmt1, chimt1, phimt2, cthetmt2, chimt2, dlmnmct, chiphit
    global z, gofz, gofzs, gzt, gzt2, titledat
    global performtest, outputaveragemol
    # Specification of l1, l2 and l values
    nll1v = 0
    nll2v = 0
    nllv = 0
    n1step = 0
    n2step = 0
    ll1v = None  # integer, dimension(:), allocatable
    ll2v = None
    llv = None
    # Arrays of l1,l2,l,n1,n2 values
    nlmax = 0
    ndist = 0
    ldistinct = 0
    ll1 = None  # integer, dimension(:), allocatable
    ll2 = None
    ll = None
    n1 = None
    n2 = None
    alsec = None  # real, dimension(:), allocatable
    # Arrays of euler angles and coefficients
    mdlmn = 0
    mz = 0
    mdist = 0
    nm = None  # integer, dimension(:), allocatable: first or second molecule
    il1rf = None  # reference integers
    il2rf = None
    ilrf = None
    illrf = None
    phim = None  # real, dimension(:), allocatable
    cthetm = None
    chim = None
    xcent = None
    ycent = None
    zcent = None
    philab = None
    dlmnc = None  # real, dimension(:,:), allocatable
    dlmnmc = None  # real, dimension(:,:,:), allocatable
    chiphi = None
    # Array of Clebsch-Gordan coefficients
    mcg = 0
    termcg = None
    # Arrays of angles and coefficients for testing
    phimt1 = None
    cthetmt1 = None
    chimt1 = None
    phimt2 = None
    cthetmt2 = None
    chimt2 = None
    dlmnmct = None
    chiphit = None
    # Input and output arrays
    z = None
    gofz = None  # real, dimension(:,:), allocatable
    gofzs = None
    gzt = None
    gzt2 = None
    titledat = None  # character(len=29), dimension(:), allocatable
    # Logical variables to perform a test calculation and to output average or single
    # molecules for plotting purposes
    performtest = False
    outputaveragemol = False


register(sys.modules[__name__])
reset()
