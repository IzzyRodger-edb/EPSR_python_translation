"""Translation of fortran_src/sharm_routines_rev1.f: the routines shared by the programs
sharm and sdf (spherical harmonic coefficients of the orientational pair correlation
function), plus form_eulerangles, which sdfcube uses too.

The data of ``module sharm_mod`` (lines 1-30 of the Fortran file) is in
``common/sharm_mod.py`` (TRANSLATION_GUIDE 5.2).  dlm, maxdlm, cgcoef, maxcgcoef,
scapro, vecpro and ran1 are EPSR's functions in auxiliary_routines.f.

Routines: get_l_and_n_values, form_cgcoeffs, getrefint, getrefintnon, form_eulerangles,
get_max_size_dlmn, form_dlmn, form_sum_prod_dlmn_cgcoeff, populate_disk.
See docs/translation_notes/Q.md.
"""

from ..common import lattice, sharm, sharm_mod, sharm_sdf  # noqa: F401
from ..fortran import (
    acosf,
    atan2f,
    cosf,
    f32,
    fmt,
    idiv,
    ifix,
    int_array,
    logf,
    lw,
    r4,
    real_array,
    sinf,
    sqrtf,
    units,
)
from ..trace import routine
from . import auxiliary_routines, minimagec, reallocation_routines
from . import lngamma_mod as lngamma_mod_lib

_F46 = "(6(1x,i5))"
_F500 = "(1x,i6,5(1x,i3),3(1x,i5),1x,e12.5)"
_F193 = "(1x,'rho = ',e13.6)"


def _do(a, b, s):
    """The iteration values of the Fortran ``do i=a,b,s`` (s > 0)."""
    return range(a, b + 1, s)


@routine
def get_l_and_n_values():
    """sharm_routines_rev1.f:get_l_and_n_values -- set up the l1,l2,l,n1,n2 values of
    the coefficients (ll1, ll2, ll, n1, n2 of sharm_mod), then form the Clebsch-Gordan
    coefficients (form_cgcoeffs)."""
    S = sharm_mod
    # set up l and n values
    # First count the number of terms
    nterm = 0
    # l1 values
    for ill1 in range(1, S.nll1v + 1):
        il1 = int(S.ll1v[ill1])
        # l2 values
        for ill2 in range(1, S.nll2v + 1):
            il2 = int(S.ll2v[ill2])
            if (S.ldistinct == 0 and il2 <= il1) or S.ldistinct == 1:
                # l values - use closure relation on CG coefficients to determine range
                jlmax = abs(il1 + il2)
                jlmin = abs(il1 - il2)
                for ill in range(1, S.nllv + 1):
                    il = int(S.llv[ill])
                    if il >= jlmin and il <= jlmax:
                        # n1 values
                        if S.n1step > 0:
                            n1max = il1 + 1
                            n1s = S.n1step
                        else:
                            n1max = 1
                            n1s = 1
                        for jn1 in _do(1, n1max, n1s):
                            in1 = jn1 - 1
                            # n2 values
                            if S.n2step > 0:
                                nts = idiv(il2, S.n2step)
                                if in1 > 0:  # Avoid redundant coefficients
                                    n2min = -nts * S.n2step
                                else:
                                    n2min = 0
                                n2max = nts * S.n2step
                                n2s = S.n2step
                            else:
                                n2min = 0
                                n2max = 0
                                n2s = 1
                            for in2 in _do(n2min, n2max, n2s):
                                nterm = nterm + 1
    S.ndist = nterm
    S.nlmax = nterm
    if S.ll1 is not None:
        S.ll1 = S.ll2 = S.ll = S.n1 = S.n2 = S.il1rf = S.il2rf = S.ilrf = S.illrf = S.alsec = None
    S.ll1 = int_array(S.nlmax)
    S.ll2 = int_array(S.nlmax)
    S.ll = int_array(S.nlmax)
    S.n1 = int_array(S.nlmax)
    S.n2 = int_array(S.nlmax)
    S.il1rf = int_array(S.nlmax)
    S.il2rf = int_array(S.nlmax)
    S.ilrf = int_array(S.nlmax)
    S.illrf = int_array(S.nlmax)
    S.alsec = real_array(S.nlmax)
    nterm = 0
    # save them!
    # l1 values
    for ill1 in range(1, S.nll1v + 1):
        il1 = int(S.ll1v[ill1])
        # l2 values
        for ill2 in range(1, S.nll2v + 1):
            il2 = int(S.ll2v[ill2])
            if (S.ldistinct == 0 and il2 <= il1) or S.ldistinct == 1:
                # l values - use closure relation on CG coefficients to determine range
                jlmax = abs(il1 + il2)
                jlmin = abs(il1 - il2)
                for ill in range(1, S.nllv + 1):
                    il = int(S.llv[ill])
                    if il >= jlmin and il <= jlmax:
                        # n1 values
                        if S.n1step > 0:
                            n1max = il1 + 1
                            n1s = S.n1step
                        else:
                            n1max = 1
                            n1s = 1
                        for jn1 in _do(1, n1max, n1s):
                            in1 = jn1 - 1
                            # n2 values
                            if S.n2step > 0:
                                nts = idiv(il2, S.n2step)
                                if in1 > 0:  # Avoid redundant coefficients
                                    n2min = -nts * S.n2step
                                else:
                                    n2min = 0
                                n2max = nts * S.n2step
                                n2s = S.n2step
                            else:
                                n2min = 0
                                n2max = 0
                                n2s = 1
                            for in2 in _do(n2min, n2max, n2s):
                                if nterm < S.nlmax:
                                    nterm = nterm + 1
                                    S.ll1[nterm] = il1
                                    S.ll2[nterm] = il2
                                    S.ll[nterm] = il
                                    S.n1[nterm] = in1
                                    S.n2[nterm] = in2
    for i in range(1, S.ndist + 1):
        units.write(6, fmt(_F46, i, S.ll1[i], S.ll2[i], S.ll[i], S.n1[i], S.n2[i]))
    # Form Clebsch-Gordan coefficients ln+ gamma values
    form_cgcoeffs()
    return


@routine
def form_cgcoeffs():
    """sharm_routines_rev1.f:form_cgcoeffs -- the reference integers of the Dlmn arrays
    and the Clebsch-Gordan coefficients (termcg) of every coefficient; alsec(id) is the
    sum of their squares."""
    S = sharm_mod
    for id_ in range(1, S.ndist + 1):
        S.alsec[id_] = f32(0.0)
    #
    # set up the required CG coefficients and determine maximum number of ln gamma values required
    #
    # First count how many needed
    ic = 0
    mgammat = 0
    for id_ in range(1, S.ndist + 1):
        il1 = int(S.ll1[id_])
        in1 = int(S.n1[id_])
        il2 = int(S.ll2[id_])
        in2 = int(S.n2[id_])
        il = int(S.ll[id_])
        #
        # initialize reference integer for molecule 1
        #
        S.il1rf[id_] = 1 + getrefint(il1, in1, S.ll1v, S.nll1v, S.n1step) + il1
        #
        # initialize reference integer for molecule 2
        #
        S.il2rf[id_] = 1 + getrefint(il2, in2, S.ll2v, S.nll2v, S.n2step) + il2
        #
        # initialize reference integer for laboratory function
        #
        S.ilrf[id_] = 1 + getrefintnon(il, S.llv, S.nllv) + il
        m1max = il1
        m1min = -il1
        for m1 in range(m1min, m1max + 1):
            mgammat = max(mgammat, auxiliary_routines.maxdlm(il1, m1, in1))  # Include the dlm's in the calculation
            mgammat = max(mgammat, auxiliary_routines.maxdlm(il1, m1, -in1))  # Include the dlm's in the calculation
            m2max = min(il2, il - m1)
            m2min = max(-il2, -il - m1)
            for m2 in range(m2min, m2max + 1):
                ic = ic + 1
                m = m1 + m2
                mgammat = max(mgammat, auxiliary_routines.maxcgcoef(il1, il2, il, m1, m2))
                mgammat = max(mgammat, auxiliary_routines.maxdlm(il2, m2, in2))  # Include the dlm's in the calculation
                mgammat = max(mgammat, auxiliary_routines.maxdlm(il2, m2, -in2))  # Include the dlm's in the calculation
                mgammat = max(mgammat, auxiliary_routines.maxdlm(il, m, 0))  # Include the dlm's in the calculation
        #         mmin=-il
        #         mmax=il
        #         do m=mmin,mmax
        #            mgammat=max(mgammat,maxdlm(il,m,0))
        #         end do
    units.write(6, lw('form_cgcoeffs> Creating ', mgammat, ' log gamma values'))
    lngamma_mod_lib.fctrln(mgammat)
    S.mcg = ic
    if S.termcg is not None:
        S.termcg = None
    S.termcg = real_array(S.mcg)
    ic = 0
    for id_ in range(1, S.ndist + 1):
        il1 = int(S.ll1[id_])
        in1 = int(S.n1[id_])
        il2 = int(S.ll2[id_])
        in2 = int(S.n2[id_])
        il = int(S.ll[id_])
        m1max = il1
        m1min = -il1
        for m1 in range(m1min, m1max + 1):
            m2max = min(il2, il - m1)
            m2min = max(-il2, -il - m1)
            for m2 in range(m2min, m2max + 1):
                ic = ic + 1
                S.termcg[ic] = auxiliary_routines.cgcoef(il1, il2, il, m1, m2)
                S.alsec[id_] = S.alsec[id_] + S.termcg[ic] * S.termcg[ic]
        #         write(6,*) id,il1,in1,il1rf(id),il2,in2,il2rf(id),il,ilrf(id)
    units.write(6, lw('form_cgcoeffs> ', S.mcg, ' Clebsch-Gordan coefficients created'))
    return


@routine
def getrefint(lreq, nreq, llv, nllv, nstep):
    """sharm_routines_rev1.f:getrefint -- gets the reference integers for the Dlmn for
    specified l and n values.  Returns getrefint."""
    sum_ = 0
    # Get the number of terms up to, but not including the current value of id
    id_ = 1
    while id_ < nllv and llv[id_] < lreq:
        l = int(llv[id_])  # noqa: E741
        if nstep > 0:
            nts = idiv(l, nstep) + 1
        else:
            nts = 1
        # For each l and n value there will be 2l+1 m values. However since -n is treated at
        # same time as +n, we treat -n same as +n
        sum_ = sum_ + nts * (2 * l + 1)
        id_ = id_ + 1
    # Add on the n values for the current id
    l = int(llv[id_])  # noqa: E741
    ntmax = abs(nreq)
    if nstep > 0:
        nts = idiv(ntmax, nstep)
        sum_ = sum_ + nts * (2 * l + 1)
    getrefint = sum_
    return getrefint


@routine
def getrefintnon(lreq, llv, nllv):
    """sharm_routines_rev1.f:getrefintnon -- gets the reference integers for the Dlmn
    for specified l with n values set to zero.  Returns getrefintnon."""
    sum_ = 0
    # Get the number of terms up to, but not including the current value of id
    id_ = 1
    while id_ < nllv and llv[id_] < lreq:
        l = int(llv[id_])  # noqa: E741
        # For each l and n value there will be 2l+1 m values. However since -n is treated at
        # same time as +n, we treat -n same as +n
        sum_ = sum_ + 2 * l + 1
        id_ = id_ + 1
    getrefintnon = sum_
    return getrefintnon


@routine
def form_eulerangles(imol1, xt, yt, zt, naxis1, naxis2, ncz, ncy, lztyp, lytyp, alab, blab, clab,
                     phim, cthetm, chim):
    """sharm_routines_rev1.f:form_eulerangles -- the Euler angles of a molecule from
    the atoms defining its two axes.

    Vectors alab, blab and clab must be defined before calling this routine.
    Returns (phim, cthetm, chim).
    """
    scapro = auxiliary_routines.scapro
    vecpro = auxiliary_routines.vecpro
    # real, dimension(3) :: pt,at,bt,ct,a1,b1,c1,perp1 (locals; alab,blab,clab are dummies)
    # NB: at and bt are not set in the Fortran when the secondary axis is undefined (opt=0);
    # they are then used only with oa1/ob1/oc1 = 0, which makes the result independent of
    # them.  They start as zeros here.
    pt = real_array(3)
    at = real_array(3)
    bt = real_array(3)
    ct = real_array(3)
    a1 = real_array(3)
    b1 = real_array(3)
    c1 = real_array(3)
    perp1 = real_array(3)

    # Get the vectors defining the two specified axes

    ierr = 0  # noqa: F841
    ifindz = 0
    ifindy = 0
    czx = f32(0.0)
    czy = f32(0.0)
    czz = f32(0.0)
    cyx = f32(0.0)
    cyy = f32(0.0)
    cyz = f32(0.0)
    for iz in range(1, ncz + 1):
        ia1 = int(lztyp[iz])
        ifindz = ifindz + 1
        czx = czx + xt[ia1]
        czy = czy + yt[ia1]
        czz = czz + zt[ia1]
    for iy in range(1, ncy + 1):
        ia1 = int(lytyp[iy])
        ifindy = ifindy + 1
        cyx = cyx + xt[ia1]
        cyy = cyy + yt[ia1]
        cyz = cyz + zt[ia1]
    if ifindz > 0:
        ct[1] = czx / f32(ifindz)
        ct[2] = czy / f32(ifindz)
        ct[3] = czz / f32(ifindz)
    else:
        ct[1] = 0.0
        ct[2] = 0.0
        ct[3] = 1.0
    if ifindy > 0:
        pt[1] = cyx / f32(ifindy)
        pt[2] = cyy / f32(ifindy)
        pt[3] = cyz / f32(ifindy)
    else:
        pt[1] = 0.0
        pt[2] = 0.0
        pt[3] = 0.0
    # c     write(6,991) p1(1),p1(2),p1(3)
    # c     write(6,991) c1(1),c1(2),c1(3)
    # renormalize these coordinates to unit vectors
    #
    oct_ = sqrtf(scapro(ct, ct))
    #
    # if oct is zero it means molecule is spherically symmetric, so no thetam can
    # be calculated for this molecule
    #
    if oct_ > 0.0:
        for i in range(1, 4):
            ct[i] = ct[i] / oct_
        oct_ = f32(1.0)
    else:
        oct_ = f32(0.0)
    opt = sqrtf(scapro(pt, pt))
    # if opt is zero it means molecule is linear so no chi dependence can be calculated.
    if opt > 0.0:
        for i in range(1, 4):
            pt[i] = pt[i] / opt
        opt = f32(1.0)
        # vector product of pt with ct (z1) gives the at axis
        vecpro(pt, ct, at)
        oat = sqrtf(scapro(at, at))
        # renormalize to give unit vector
        for i in range(1, 4):
            at[i] = at[i] / oat
        oat = f32(1.0)
        # vector product of ct to at gives  bt
        vecpro(ct, at, bt)
        obt = sqrtf(scapro(bt, bt))
        # renormalize to give unit vector
        for i in range(1, 4):
            bt[i] = bt[i] / obt
        obt = f32(1.0)
    else:
        oat = f32(0.0)
        obt = f32(0.0)
    oa1 = f32(0.0)
    ob1 = f32(0.0)
    oc1 = f32(0.0)
    # put axis coordinates into correct arrays, depending on values of naxis1,naxis2
    if naxis1 == 'z' and naxis2 == 'y':
        for i in range(1, 4):
            a1[i] = at[i]
            b1[i] = bt[i]
            c1[i] = ct[i]
        oa1 = oat
        ob1 = obt
        oc1 = oct_
    elif naxis1 == '-z' and naxis2 == 'y':
        for i in range(1, 4):
            a1[i] = at[i]
            b1[i] = bt[i]
            c1[i] = -ct[i]
        oa1 = oat
        ob1 = obt
        oc1 = oct_
    elif naxis1 == 'y' and naxis2 == 'z':
        for i in range(1, 4):
            a1[i] = -at[i]
            b1[i] = ct[i]
            c1[i] = bt[i]
        oa1 = oat
        ob1 = oct_
        oc1 = obt
    elif naxis1 == 'y' and naxis2 == '-z':
        for i in range(1, 4):
            a1[i] = -at[i]
            b1[i] = ct[i]
            c1[i] = -bt[i]
        oa1 = oat
        ob1 = oct_
        oc1 = obt
    elif naxis1 == 'y' and naxis2 == 'x':
        for i in range(1, 4):
            a1[i] = bt[i]
            b1[i] = ct[i]
            c1[i] = at[i]
        oa1 = obt
        ob1 = oct_
        oc1 = oat
    elif naxis1 == 'x' and naxis2 == 'y':
        for i in range(1, 4):
            a1[i] = ct[i]
            b1[i] = bt[i]
            c1[i] = -at[i]
        oa1 = oct_
        ob1 = obt
        oc1 = oat
    elif naxis1 == 'z' and naxis2 == 'x':
        for i in range(1, 4):
            a1[i] = bt[i]
            b1[i] = -at[i]
            c1[i] = ct[i]
        oa1 = obt
        ob1 = oat
        oc1 = oct_
    elif naxis1 == '-z' and naxis2 == 'x':
        for i in range(1, 4):
            a1[i] = bt[i]
            b1[i] = -at[i]
            c1[i] = -ct[i]
        oa1 = obt
        ob1 = oat
        oc1 = oct_
    elif naxis1 == 'x' and naxis2 == 'z':
        for i in range(1, 4):
            a1[i] = ct[i]
            b1[i] = at[i]
            c1[i] = bt[i]
        oa1 = oct_
        ob1 = oat
        oc1 = obt
    elif naxis1 == 'x' and naxis2 == '-z':
        for i in range(1, 4):
            a1[i] = ct[i]
            b1[i] = at[i]
            c1[i] = -bt[i]
        oa1 = oct_
        ob1 = oat
        oc1 = obt
    #
    # now evaluate theta value for this molecule
    #
    if oc1 > 0.0:
        cthetm = scapro(clab, c1)
        #
        # vector product of clab to c1 gives the y' axis about which rotation thetam was
        # made
        #
        vecpro(clab, c1, perp1)
        #
        # renormalize to give a unit vector
        #
        operp1 = sqrtf(scapro(perp1, perp1))
        if operp1 > 0.0:
            for i in range(1, 4):
                perp1[i] = perp1[i] / operp1
            operp1 = f32(1.0)
            sphim = -scapro(perp1, alab)
            cphim = scapro(perp1, blab)
            phim = atan2f(sphim, cphim)
        else:
            phim = f32(0.0)
        # now calculate chim
        #
        # scalar product of perp1 with a1 and b1 give chim rotation of molecule
        # about molecules c-axis
        #
        if oa1 > 0.0 and ob1 > 0.0:
            schim = scapro(perp1, a1)
            cchim = scapro(perp1, b1)
            chim = atan2f(schim, cchim)
        else:
            chim = f32(0.0)
    else:
        cthetm = f32(0.0)
        phim = f32(0.0)
        chim = f32(0.0)
    return phim, cthetm, chim


@routine
def get_max_size_dlmn():
    """sharm_routines_rev1.f:get_max_size_dlmn -- mdlmn, the size of the per-molecule
    and laboratory Dlmn arrays."""
    S = sharm_mod
    ilm = 0
    for i in range(1, S.nll1v + 1):
        il = int(S.ll1v[i])
        ntstep = S.n1step
        # n values
        if ntstep > 0:
            ntmax = il
            nts = ntstep
        else:
            ntmax = 0
            nts = 1
        for in_ in _do(0, ntmax, nts):
            ilm = ilm + 2 * il + 1
    S.mdlmn = ilm
    units.write(6, lw(S.mdlmn, S.nll1v, ilm))
    ilm = 0
    for i in range(1, S.nll2v + 1):
        il = int(S.ll2v[i])
        ntstep = S.n2step
        # n values
        if ntstep > 0:
            ntmax = il
            nts = ntstep
        else:
            ntmax = 0
            nts = 1
        for in_ in _do(0, ntmax, nts):
            ilm = ilm + 2 * il + 1
    S.mdlmn = max(S.mdlmn, ilm)
    units.write(6, lw(S.mdlmn, S.nll2v, ilm))
    ilm = 0
    for ilt in range(1, S.nllv + 1):
        il = int(S.llv[ilt])
        ilm = ilm + 2 * il + 1
    S.mdlmn = max(S.mdlmn, ilm)
    units.write(6, lw(S.mdlmn, S.nllv, ilm))
    return


@routine
def form_dlmn(imol1, ierr):
    """sharm_routines_rev1.f:form_dlmn -- calculate the D coefficients (dlmnmc, chiphi)
    of molecule imol1.  Returns ierr."""
    S = sharm_mod
    dlm = auxiliary_routines.dlm
    lltt = None  # integer, dimension(:), allocatable
    #
    # calculate the D coefficients for this molecule
    #
    ierr = 0
    lltt = reallocation_routines.reallocate1d_i(lltt, max(S.nll1v, S.nll2v))
    if S.nm[imol1] == 1:
        nllt = S.nll1v
        for i in range(1, S.nll1v + 1):
            lltt[i] = S.ll1v[i]
        ntstep = S.n1step
    else:
        nllt = S.nll2v
        for i in range(1, S.nll2v + 1):
            lltt[i] = S.ll2v[i]
        ntstep = S.n2step
    #
    # l values
    #
    nm1 = int(S.nm[imol1])
    ilmnlm = 0
    for ilt in range(1, nllt + 1):
        il = int(lltt[ilt])
        # range of m values
        mmin = -il
        mmax = il
        # n values
        if ntstep > 0:
            ntmax = il
            nts = ntstep
        else:
            ntmax = 0
            nts = 1
        for in_ in _do(0, ntmax, nts):
            for m in range(mmin, mmax + 1):
                ilmnlm = ilmnlm + 1
                # c     write(6,501) il,m,in,ilmnlm
                # error exit if too many coeficients demanded
                if ilmnlm > S.mdlmn:
                    ierr = ilmnlm - S.mdlmn
                    units.write(6, lw('form_dlmn> Array reference ', ilmnlm, ' larger than allowed ', S.mdlmn))
                    return ierr
                if S.performtest:
                    S.chiphi[1, imol1, ilmnlm] = f32(m) * S.phimt1[imol1] + f32(in_) * S.chimt1[imol1]
                    S.dlmnmc[1, imol1, ilmnlm] = dlm(il, m, in_, S.cthetmt1[imol1])
                    S.chiphit[1, imol1, ilmnlm] = f32(m) * S.phimt2[imol1] + f32(in_) * S.chimt2[imol1]
                    S.dlmnmct[1, imol1, ilmnlm] = dlm(il, m, in_, S.cthetmt2[imol1])
                    if in_ != 0 and (nm1 == 2 or (nm1 == 1 and S.ldistinct == 0)):
                        S.chiphi[2, imol1, ilmnlm] = f32(m) * S.phimt1[imol1] - f32(in_) * S.chimt1[imol1]
                        S.dlmnmc[2, imol1, ilmnlm] = dlm(il, m, -in_, S.cthetmt1[imol1])
                        S.chiphit[2, imol1, ilmnlm] = f32(m) * S.phimt2[imol1] - f32(in_) * S.chimt2[imol1]
                        S.dlmnmct[2, imol1, ilmnlm] = dlm(il, m, -in_, S.cthetmt2[imol1])
                else:
                    S.chiphi[1, imol1, ilmnlm] = f32(m) * S.phim[imol1] + f32(in_) * S.chim[imol1]
                    S.dlmnmc[1, imol1, ilmnlm] = dlm(il, m, in_, S.cthetm[imol1])
                    if in_ != 0 and (nm1 == 2 or (nm1 == 1 and S.ldistinct == 0)):
                        S.chiphi[2, imol1, ilmnlm] = f32(m) * S.phim[imol1] - f32(in_) * S.chim[imol1]
                        S.dlmnmc[2, imol1, ilmnlm] = dlm(il, m, -in_, S.cthetm[imol1])
    if S.performtest:
        units.write(6, fmt('(1x,i5,6(1x,e12.5))', imol1, S.phimt1[imol1], S.cthetmt1[imol1], S.chimt1[imol1],
                           S.phimt2[imol1], S.cthetmt2[imol1], S.chimt2[imol1]))
    else:
        units.write(6, fmt('(1x,i5,3(1x,e12.5))', imol1, S.phim[imol1], S.cthetm[imol1], S.chim[imol1]))
    return ierr


@routine
def form_sum_prod_dlmn_cgcoeff(nmolsharm, ierr):
    """sharm_routines_rev1.f:form_sum_prod_dlmn_cgcoeff -- sum the products of the
    molecular and laboratory Dlmn and the Clebsch-Gordan coefficients over all pairs of
    molecules (gofz, gofzs), then normalise them.  Returns ierr."""
    S = sharm_mod
    scapro = auxiliary_routines.scapro
    vecpro = auxiliary_routines.vecpro
    dlm = auxiliary_routines.dlm
    pdif = sharm_sdf.pdif
    perp1 = sharm_sdf.perp1
    alab = sharm_sdf.alab
    blab = sharm_sdf.blab
    clab = sharm_sdf.clab
    ierr = 0
    nsec = 0
    itimestotal = 0
    gfac = r4(4.0) * sharm_sdf.pi / r4(3.0)
    nm = S.nm
    #
    # initialize counters for number of molecules of each kind in the system
    #
    if S.ldistinct == 0:
        if nm[nmolsharm] == 1 or nm[nmolsharm] == 2:
            lsum1 = 1
            lsum2 = 1
        else:
            lsum1 = 0
            lsum2 = 0
    else:
        lsum1 = 0
        lsum2 = 0
        if nm[nmolsharm] == 1:
            lsum1 = 1
        if nm[nmolsharm] == 2:
            lsum2 = 1
    # Count how many type 1 molecules there are
    notype1 = 0
    notype2 = 0
    for imol1 in range(1, nmolsharm + 1):
        if S.ldistinct == 0:
            if nm[imol1] == 1 or nm[imol1] == 2:
                notype1 = notype1 + 1
                notype2 = notype2 + 1
        else:
            if nm[imol1] == 1:
                notype1 = notype1 + 1
            if nm[imol1] == 2:
                notype2 = notype2 + 1
    if S.ldistinct == 0:
        ntimestotal = notype1 * (notype1 - 1)
    else:
        ntimestotal = notype1 * notype2
    units.write(6, lw(nmolsharm, nm[nmolsharm], lsum1, lsum2))
    # Local names for the arrays of the innermost loops (the same objects)
    nz = sharm.nz
    zstep = sharm.zstep
    ndist = S.ndist
    nllv = S.nllv
    llv = S.llv
    ll1 = S.ll1
    ll2 = S.ll2
    ll = S.ll
    n2 = S.n2
    il1rf = S.il1rf
    il2rf = S.il2rf
    ilrf = S.ilrf
    xcent = S.xcent
    ycent = S.ycent
    zcent = S.zcent
    philab = S.philab
    dlmnc = S.dlmnc
    termcg = list(S.termcg)  # termcg is not changed below; list of its numpy.float32 elements
    gofz = S.gofz
    gofzs = S.gofzs
    performtest = S.performtest
    ll_ndist = ll[ndist]
    for imol1 in range(1, nmolsharm - 1 + 1):
        if nm[imol1] != 0:
            # c     call system_clock(time3,rate,max)
            #
            # increment counters for number of molecules of types 1 and 2 in the system
            #
            if S.ldistinct == 0:
                if nm[imol1] == 1 or nm[imol1] == 2:
                    lsum1 = lsum1 + 1
                    lsum2 = lsum2 + 1
            else:
                if nm[imol1] == 1:
                    lsum1 = lsum1 + 1
                if nm[imol1] == 2:
                    lsum2 = lsum2 + 1
            # NB: ntimestotal = 0 (no pairs) is an integer division by zero in the Fortran too
            units.write(6, fmt('(1x,i10,1x,i10,1x,i8,"% done",a)', itimestotal, ntimestotal,
                               (itimestotal * 100) // ntimestotal, chr(13)), advance=False)
            xc1 = xcent[imol1]
            yc1 = ycent[imol1]
            zc1 = zcent[imol1]
            #
            # calc now look at second molecules
            #
            nm2f = imol1 + 1
            for imol2 in range(nm2f, nmolsharm + 1):
                if imol1 != imol2 and nm[imol2] != 0:
                    iflag1 = 0
                    iflag2 = 0
                    nm1 = nm[imol1]
                    nm2 = nm[imol2]
                    if S.ldistinct == 0 and nm1 == 1 and nm2 == 1:
                        iflag1 = 1
                        iflag2 = 1
                    if S.ldistinct == 1 and nm1 == 1 and nm2 == 2:
                        iflag1 = 1
                        iflag2 = 0
                    if S.ldistinct == 1 and nm1 == 2 and nm2 == 1:
                        iflag1 = 0
                        iflag2 = 1
                    if (iflag1 + iflag2) != 0:
                        itimestotal = itimestotal + iflag1 + iflag2
                        xc2 = xcent[imol2]
                        yc2 = ycent[imol2]
                        zc2 = zcent[imol2]
                        xdif = xc2 - xc1
                        ydif = yc2 - yc1
                        zdif = zc2 - zc1
                        # apply periodic boundary conditions
                        xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                        # c                xc2=xc2+anint(xc1-xc2)
                        # c calculate separation of molecular centres
                        pdif[1] = xdif
                        pdif[2] = ydif
                        pdif[3] = zdif
                        # renormalize to unit vector
                        opdif = sqrtf(scapro(pdif, pdif))
                        d12 = opdif
                        for i in range(1, 4):
                            pdif[i] = pdif[i] / opdif
                        opdif = f32(1.0)
                        iref = ifix(d12 / zstep) + 1
                        # check that distance is within bounds prescribed for this distribution
                        if iref <= nz:
                            # apply cosine rule to determine cos(theta) for intermolecular axis
                            cz = scapro(clab, pdif)
                            # vector perpendicular to c1 and pdif - gives the rotated y-axis
                            vecpro(clab, pdif, perp1)
                            # normalize to give a unit vector
                            operp1 = sqrtf(scapro(perp1, perp1))
                            if operp1 > 0.0:
                                for i in range(1, 4):
                                    perp1[i] = perp1[i] / operp1
                                operp1 = f32(1.0)
                                # scalar product with a1 to give -sin(phi)
                                sphil = -scapro(alab, perp1)
                                # scalar product with b1 to give cos(phi)
                                cphil = scapro(blab, perp1)
                                # angle phi in lab coordinate system:
                                phil = atan2f(sphil, cphil)
                            else:
                                phil = f32(0.0)
                            if performtest:
                                #  call populate_disk(pi,phillab,thetallab,0.0,phinew,thetanew,chinew,thetaspread,phispread)
                                phil = sharm_sdf.phillab
                                cz = cosf(sharm_sdf.thetallab)
                            #
                            # calculate coefficients for laboratory coordinates
                            #
                            # l values
                            #
                            ilm = 0
                            for ilt in range(1, nllv + 1):
                                il = int(llv[ilt])
                                # Don't bother with values disallowed by the CG coefficients
                                if il <= ll_ndist:
                                    # range of m values
                                    mmin = -il
                                    mmax = il
                                    for m in range(mmin, mmax + 1):
                                        ilm = ilm + 1
                                        # error exit if too many coeficients demanded
                                        philab[ilm] = f32(m) * phil
                                        dlmnc[1, ilm] = dlm(il, m, 0, cz)
                                        # (-1)**m is an integer power: +-1 also for m < 0
                                        dlmnc[2, ilm] = f32(1 if m % 2 == 0 else -1) * dlm(il, m, 0, -cz)
                            #
                            # FINALLY we sum the products of coefficients
                            #
                            if iflag1 == 1:  # Molecule 1 is at the origin
                                _sum_products(imol1, imol2, 1, iref, ndist, ll1, ll2, ll, n2, il1rf, il2rf, ilrf,
                                              philab, dlmnc, termcg, gofz, gofzs, performtest)
                                nsec = nsec + 1
                            if iflag2 == 1 and not performtest:  # Molecule 2 is at the origin
                                _sum_products(imol2, imol1, 2, iref, ndist, ll1, ll2, ll, n2, il1rf, il2rf, ilrf,
                                              philab, dlmnc, termcg, gofz, gofzs, performtest)
                                nsec = nsec + 1
            # c     call system_clock(time4,rate,max)
            # c     write(6,*) rate,max,time4-time3,nsec
    units.write(6, fmt('(1x,i10,1x,i10,1x,i8,"% done")', itimestotal, ntimestotal, (itimestotal * 100) // ntimestotal))
    # 502   format(1x,i10,2(1x,i5))
    for id_ in range(1, ndist + 1):
        units.write(6, fmt(_F500, id_, ll1[id_], ll2[id_], ll[id_], S.n1[id_], n2[id_], il1rf[id_], il2rf[id_],
                           ilrf[id_], S.alsec[id_]))
    # renormalize g(r)'s to total number of atoms
    # density of second atoms in sphere
    rho = f32(lsum1) * f32(lsum2) / sharm.vol
    units.write(6, fmt(_F193, rho))
    gfac1 = gfac * rho
    z = S.z
    for id_ in range(1, ndist + 1):
        z0 = f32(0.0)
        z03 = f32(0.0)
        for iz in range(1, nz + 1):
            zm = z[iz]
            zdiff = zm - z0
            z1 = zm + zdiff
            z13 = z1 * z1 * z1
            # average number of particles in this volume element
            rat = gfac1 * (z13 - z03)
            gofz[iz, id_] = f32((2 * int(ll1[id_]) + 1) * (2 * int(ll2[id_]) + 1)) * gofz[iz, id_] / rat
            gofzs[iz, id_] = f32((2 * int(ll1[id_]) + 1) * (2 * int(ll2[id_]) + 1)) * gofzs[iz, id_] / rat
            if id_ == 1:
                gofz[iz, id_] = gofz[iz, id_] - r4(1.0)
            z0 = z1
            z03 = z13
            # c           write(6,996) iz,zm,gofz(iz,id)
            # 996   format(1x,i5,2(1x,e12.5))
    return ierr


def _sum_products(ma, mb, irow, iref, ndist, ll1, ll2, ll, n2, il1rf, il2rf, ilrf, philab, dlmnc, termcg,
                  gofz, gofzs, performtest):
    """The two "FINALLY we sum the products of coefficients" blocks of
    form_sum_prod_dlmn_cgcoeff (sharm_routines_rev1.f lines 845-926 and 927-997).

    The blocks are identical except for which molecule is at the origin: the first block
    (molecule 1 at the origin) is ma=imol1, mb=imol2, irow=1 (dlmnc(1,ilmn)); the second
    (molecule 2 at the origin, only without performtest) is ma=imol2, mb=imol1, irow=2.
    ``chiphi(1,ma,ilmn1)``, ``dlmnmc(1,ma,ilmn1)`` belong to the molecule at the origin and
    ``chiphi(k,mb,ilmn2)``, ``dlmnmc(k,mb,ilmn2)`` (k=2 if n2 < 0, else 1) to the other one;
    with performtest the latter come from chiphit, dlmnmct.  The arithmetic is the
    Fortran's, term by term in the same order.
    """
    S = sharm_mod
    # The rows used here as Python lists of the same numpy.float32 elements (list indexing is
    # faster than ndarray indexing; the values and the arithmetic are unchanged)
    chiphi_a = list(S.chiphi[1, ma, :])
    dlmnmc_a = list(S.dlmnmc[1, ma, :])
    if performtest:
        chiphi_b1 = list(S.chiphit[1, mb, :])
        chiphi_b2 = list(S.chiphit[2, mb, :])
        dlmnmc_b1 = list(S.dlmnmct[1, mb, :])
        dlmnmc_b2 = list(S.dlmnmct[2, mb, :])
    else:
        chiphi_b1 = list(S.chiphi[1, mb, :])
        chiphi_b2 = list(S.chiphi[2, mb, :])
        dlmnmc_b1 = list(S.dlmnmc[1, mb, :])
        dlmnmc_b2 = list(S.dlmnmc[2, mb, :])
    dlmnc_r = list(dlmnc[irow, :])
    philab = list(philab)
    _cosf = cosf
    _sinf = sinf
    ic = 0
    for id_ in range(1, ndist + 1):
        il1 = int(ll1[id_])
        il2 = int(ll2[id_])
        in2 = int(n2[id_])
        il = int(ll[id_])
        #
        # initialize reference integer for molecule 1
        #
        il1ref = int(il1rf[id_])
        #
        # initialize reference integer for molecule 2
        #
        il2ref = int(il2rf[id_])
        #
        # initialize reference integer for laboratory function
        #
        ilref = int(ilrf[id_])
        if in2 < 0:
            chiphi_b = chiphi_b2
            dlmnmc_b = dlmnmc_b2
        else:
            chiphi_b = chiphi_b1
            dlmnmc_b = dlmnmc_b1
        gz = gofz[iref, id_]
        gzs = gofzs[iref, id_]
        m1max = il1
        m1min = -il1
        for m1 in range(m1min, m1max + 1):
            m2max = min(il2, il - m1)
            m2min = max(-il2, -il - m1)
            for m2 in range(m2min, m2max + 1):
                m = m1 + m2
                ic = ic + 1
                ilmn1 = il1ref + m1
                ilmn2 = il2ref + m2
                ilmn = ilref + m
                phisum = philab[ilmn] - chiphi_a[ilmn1] - chiphi_b[ilmn2]
                cosphi = _cosf(phisum)
                sinphi = _sinf(phisum)
                prod = dlmnmc_a[ilmn1] * dlmnmc_b[ilmn2]
                # 997   format(1x,6(1x,i4),3(1x,e12.5))
                factor = prod * dlmnc_r[ilmn] * termcg[ic]
                add = factor * cosphi
                adds = factor * sinphi
                # gofz(iref,id)=gofz(iref,id)+add (accumulated in gz, stored after the loops)
                gz = gz + add
                gzs = gzs + adds
        gofz[iref, id_] = gz
        gofzs[iref, id_] = gzs


@routine
def populate_disk(pi, phi0, theta0, chi0, phinew, thetanew, chinew, thetaspread, phispread):
    """sharm_routines_rev1.f:populate_disk -- randomly populates a tangential disk at
    (theta0,phi0) with a Gaussian-like distribution at constant radius, giving a new
    (theta,phi).  Generates random chi about chi0.

    Returns (phinew, thetanew, chinew).
    """
    ran1 = auxiliary_routines.ran1
    ctheta0 = cosf(theta0)
    stheta0 = sinf(theta0)
    x0 = stheta0 * cosf(phi0)
    y0 = stheta0 * sinf(phi0)
    z0 = ctheta0
    z02 = z0 * z0
    cons = f32(-1.0)
    deltax = f32(0.0)
    deltay = f32(0.0)
    u2 = f32(0.0)
    while cons < 0.0:
        u1 = ran1()
        u2 = ran1()
        amp = sqrtf(r4(-2.0) * logf(u1))
        deltax = amp * cosf(r4(2.0) * pi * u2) * phispread
        deltay = amp * sinf(r4(2.0) * pi * u2) * thetaspread
        cons = z02 - (deltax * (deltax + r4(2.0) * x0) + deltay * (deltay + r4(2.0) * y0))
    # Choose the root closest to z0
    cthetanew = sqrtf(cons)
    root1 = acosf(cthetanew)
    root2 = acosf(-cthetanew)
    dif1 = abs(theta0 - root1) / thetaspread
    dif2 = abs(theta0 - root2) / thetaspread
    # If either or these differences is greater than 5 deviations, choose the other, otherwise
    # toss a coin.
    if (dif2 - dif1) > 5.0:
        thetanew = root1
    elif (dif1 - dif2) > 5.0:
        thetanew = root2
    else:  # Toss a coin to choose between the two roots
        if ran1() < 0.5:
            thetanew = root1
        else:
            thetanew = root2
    sthetanew = sqrtf(r4(1.0) - cthetanew * cthetanew)
    cphinew = (x0 + deltax) / sthetanew
    sphinew = (y0 + deltay) / sthetanew
    phinew = atan2f(sphinew, cphinew)
    u1 = u2
    u2 = ran1()
    chinew = chi0 + sqrtf(f32(-2) * logf(u1)) * cosf(f32(2) * pi * u2) * phispread
    #      write(6,*) x0,y0,z0,deltax,deltay,root1,root2,dif1,dif2
    #      write(6,*) phi0,theta0,chi0,phinew,thetanew,chinew
    return phinew, thetanew, chinew
