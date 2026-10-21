"""Generated from fortran_src/molecom.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``molecom.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401
from .epsr_dimension import matom, mcomp, mpairs, near  # noqa: F401

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of molecom.inc."""
    global tol, tolout, step, tempf, tempb, tempc, stepf, stepd
    global steps, tempdih, tempang, tempdihreq, tempangreq, tempdihwt, tempangwt, stepr
    global stepri, stepm, stepmi, xa, ya, za, sxa, sya
    global sza, amassa, amassr, tbond, bound, eatomt, eatomv, xs
    global ys, zs, fold, trial, trialt, trialv, simt, simtk
    global charge, qradius, qcor, dcor, sigcor, epscor, res, sig1m
    global sig2m, scaleb, distob, dscor, ecor, xbar, ybar, zbar
    global dscorm, x2bar, y2bar, z2bar, sig1b, sig2b, sig3b, xid
    global yid, zid, scalef, distor, sig1b2, sig2b2, sig3b2, scal2b
    global dist2b, cfrac, cfraccg, dcore, ecore, dcorm, dcorsq, dcorsm
    global power, sizefactor, decline_rate, dmaxmol, refine_intra_fac, coul_intra_fac, sumcharge, dihedralang
    global maxmolsize, alabel, type, atsymbol, typefile, labrotdih, ndihedralang, nearn
    global nbondp, nbonds, nbond, iabond, ntype, itype, npair, iadd
    global ntry, ntryo, nbar, iok, natom, natomr, mref, qtype
    global qtypea, npairs, npairref, n1extr, n2extr, iexchange, ntypefile, iatomicnotype
    global natomsmoved, atomsmoved, imageatom, refine_intra, coul_intra
    # common /molecom_r/
    tol = f32(0.0)
    tolout = f32(0.0)
    step = f32(0.0)
    tempf = f32(0.0)
    tempb = f32(0.0)
    tempc = f32(0.0)
    stepf = f32(0.0)
    stepd = f32(0.0)
    steps = f32(0.0)
    tempdih = f32(0.0)
    tempang = f32(0.0)
    tempdihreq = f32(0.0)
    tempangreq = f32(0.0)
    tempdihwt = f32(0.0)
    tempangwt = f32(0.0)
    stepr = f32(0.0)
    stepri = f32(0.0)
    stepm = f32(0.0)
    stepmi = f32(0.0)
    xa = real_array(matom)
    ya = real_array(matom)
    za = real_array(matom)
    sxa = real_array(matom)
    sya = real_array(matom)
    sza = real_array(matom)
    amassa = real_array(mcomp)
    amassr = real_array(mpairs)
    tbond = real_array(matom, near)
    bound = real_array(matom)
    eatomt = real_array(matom)
    eatomv = real_array(matom)
    xs = f32(0.0)
    ys = f32(0.0)
    zs = f32(0.0)
    fold = f32(0.0)
    trial = f32(0.0)
    trialt = real_array(matom)
    trialv = real_array(matom)
    simt = f32(0.0)
    simtk = f32(0.0)
    charge = real_array(mcomp)
    qradius = real_array(mcomp)
    qcor = real_array(mpairs)
    dcor = real_array(mpairs)
    sigcor = real_array(mcomp)
    epscor = real_array(mcomp)
    res = f32(0.0)
    sig1m = f32(0.0)
    sig2m = f32(0.0)
    scaleb = f32(0.0)
    distob = f32(0.0)
    dscor = real_array(mpairs)
    ecor = real_array(mpairs)
    xbar = real_array(matom)
    ybar = real_array(matom)
    zbar = real_array(matom)
    dscorm = real_array(mpairs)
    x2bar = real_array(matom)
    y2bar = real_array(matom)
    z2bar = real_array(matom)
    sig1b = f32(0.0)
    sig2b = f32(0.0)
    sig3b = f32(0.0)
    xid = real_array(matom)
    yid = real_array(matom)
    zid = real_array(matom)
    scalef = f32(0.0)
    distor = f32(0.0)
    sig1b2 = f32(0.0)
    sig2b2 = f32(0.0)
    sig3b2 = f32(0.0)
    scal2b = f32(0.0)
    dist2b = f32(0.0)
    cfrac = real_array(mcomp)
    cfraccg = real_array(mcomp)
    dcore = f32(0.0)
    ecore = f32(0.0)
    dcorm = f32(0.0)
    dcorsq = f32(0.0)
    dcorsm = f32(0.0)
    power = f32(0.0)
    sizefactor = f32(0.0)
    decline_rate = f32(0.0)
    dmaxmol = f32(0.0)
    refine_intra_fac = f32(0.0)
    coul_intra_fac = f32(0.0)
    sumcharge = f32(0.0)
    dihedralang = real_array(near, matom)
    maxmolsize = f32(0.0)
    # common /molecom_c/
    alabel = char_array(matom)
    type = char_array(mcomp)
    atsymbol = char_array(mcomp)
    typefile = char_array(mcomp)
    # common /molecom_i/
    labrotdih = int_array(2, near, matom)
    ndihedralang = int_array(matom)
    nearn = int_array(matom)
    nbondp = int_array(matom)
    nbonds = int_array(matom)
    nbond = int_array(matom)
    iabond = int_array(matom, near)
    ntype = 0
    itype = int_array(matom)
    npair = 0
    iadd = int_array(matom)
    ntry = 0
    ntryo = 0
    nbar = 0
    iok = int_array(matom)
    natom = 0
    natomr = 0
    mref = int_array(matom)
    qtype = int_array(matom)
    qtypea = int_array(mcomp)
    npairs = 0
    npairref = int_array(mcomp, mcomp)
    n1extr = int_array(mpairs)
    n2extr = int_array(mpairs)
    iexchange = int_array(mcomp)
    ntypefile = 0
    iatomicnotype = int_array(mcomp)
    natomsmoved = 0
    atomsmoved = int_array(near)
    imageatom = int_array(matom)
    # common /molecom_l/
    refine_intra = False
    coul_intra = False


reset()
