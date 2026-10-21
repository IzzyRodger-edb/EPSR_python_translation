"""Translation of fortran_src/plotsharm_routines.f: smoothing and interpolation helpers
of the spherical-harmonic plotting programs (mapgr, mapgrentropy): smoosharm,
quadinterp2d, quadinterp, sumoverm2, suminterp."""

from ..fortran import cosd_intrinsic, f32, isign, r4, sind_intrinsic
from ..trace import routine
from . import auxiliary_routines


@routine
def smoosharm(nqmin, nqmax, nsmoo, diff):
    """plotsharm_routines.f:smoosharm -- 3-point smoothing keeping end points fixed.

    diff (real, dimension(*)) is smoothed in place.
    """
    n1 = nqmin + 1
    n2 = nqmax - 1
    diff[nqmax] = 0.0  # Last value is set to zero exactly
    if (n2 - n1) < 2 or nsmoo == 0:
        return
    for it in range(1, nsmoo + 1):
        sum1 = r4(0.25) * (diff[nqmin] + diff[n1])
        sm = diff[n1]
        for i in range(n1, n2 + 1):
            j = i + 1
            sn = diff[j]
            sum2 = r4(0.25) * (sm + sn)
            diff[i] = sum1 + sum2
            sm = sn
            sum1 = sum2
    return


@routine
def quadinterp2d(grw, rrat, thrat, nro, nrp, nrm, ntho, nthp, nthm, ntdim, nrdim):
    """plotsharm_routines.f:quadinterp2d -- quadratic interpolation of grw(ntdim,nrdim)
    along r, then in theta.  Returns quadinterp2d."""
    #
    # quadratic interpolation along r
    #
    delta = rrat
    dp = r4(1.0) + delta
    dm = r4(1.0) - delta
    dpm = dm * dp
    delta2 = r4(0.5) * delta
    am = delta2 * (dp * grw[nthm, nrp] - dm * grw[nthm, nrm]) + dpm * grw[nthm, nro]
    ao = delta2 * (dp * grw[ntho, nrp] - dm * grw[ntho, nrm]) + dpm * grw[ntho, nro]
    ap = delta2 * (dp * grw[nthp, nrp] - dm * grw[nthp, nrm]) + dpm * grw[nthp, nro]
    #
    # now interpolate in theta
    #
    delta = thrat
    dp = r4(1.0) + delta
    dm = r4(1.0) - delta
    dpm = dm * dp
    delta2 = r4(0.5) * delta
    quadinterp2d = delta2 * (dp * ap - dm * am) + dpm * ao
    return quadinterp2d


@routine
def quadinterp(grw, rrat, thrat, phrat, nro, nrp, nrm, ntho, nthp, nthm, npho, nphp, nphm, ntdim, npdim, nrdim):
    """plotsharm_routines.f:quadinterp -- quadratic interpolation of grw(ntdim,npdim,nrdim)
    along r, then theta, then phi.  Returns quadinterp."""
    #
    # quadratic interpolation along r
    #
    delta = rrat
    dp = r4(1.0) + delta
    dm = r4(1.0) - delta
    dpm = dm * dp
    delta2 = r4(0.5) * delta
    amm = delta2 * (dp * grw[nthm, nphm, nrp] - dm * grw[nthm, nphm, nrm]) + dpm * grw[nthm, nphm, nro]
    amo = delta2 * (dp * grw[nthm, npho, nrp] - dm * grw[nthm, npho, nrm]) + dpm * grw[nthm, npho, nro]
    amp = delta2 * (dp * grw[nthm, nphp, nrp] - dm * grw[nthm, nphp, nrm]) + dpm * grw[nthm, nphp, nro]
    aom = delta2 * (dp * grw[ntho, nphm, nrp] - dm * grw[ntho, nphm, nrm]) + dpm * grw[ntho, nphm, nro]
    aoo = delta2 * (dp * grw[ntho, npho, nrp] - dm * grw[ntho, npho, nrm]) + dpm * grw[ntho, npho, nro]
    aop = delta2 * (dp * grw[ntho, nphp, nrp] - dm * grw[ntho, nphp, nrm]) + dpm * grw[ntho, nphp, nro]
    apm = delta2 * (dp * grw[nthp, nphm, nrp] - dm * grw[nthp, nphm, nrm]) + dpm * grw[nthp, nphm, nro]
    apo = delta2 * (dp * grw[nthp, npho, nrp] - dm * grw[nthp, npho, nrm]) + dpm * grw[nthp, npho, nro]
    app = delta2 * (dp * grw[nthp, nphp, nrp] - dm * grw[nthp, nphp, nrm]) + dpm * grw[nthp, nphp, nro]
    #
    # now interpolate in theta
    #
    delta = thrat
    dp = r4(1.0) + delta
    dm = r4(1.0) - delta
    dpm = dm * dp
    delta2 = r4(0.5) * delta
    am = delta2 * (dp * apm - dm * amm) + dpm * aom
    ao = delta2 * (dp * apo - dm * amo) + dpm * aoo
    ap = delta2 * (dp * app - dm * amp) + dpm * aop
    #
    # now interpolate in phi
    #
    delta = phrat
    dp = r4(1.0) + delta
    dm = r4(1.0) - delta
    dpm = dm * dp
    delta2 = r4(0.5) * delta
    quadinterp = delta2 * (dp * ap - dm * am) + dpm * ao
    return quadinterp


@routine
def sumoverm2(l1, l2, l, n1, m2f, m2l, n2, phil, ctl, phi, cz, chi, phs, sumc, sums):
    """plotsharm_routines.f:sumoverm2 -- sum over m2 of the rotational invariant terms.

    cosd/sind are not declared external here, so gfortran (>= 14) uses its intrinsics
    (fmath.cosd_intrinsic/sind_intrinsic), not EPSR's own (Rust sharm_routines.rs).
    Returns (sumc, sums).
    """
    # integer l1,l2,l,n1,m2f,m2l,n2,m2
    # real cgcoef,dlm,phil,ctl,phi,cz,chi,phs,sumc,sums,phisum,cg,dl2m2n2,dlm0
    # For negative n1, we need the complex conjugate (See Gray and Gubbins, eq 2.25). Sharm gives us only +ve n1,
    # but all n2
    phssine = -f32(isign(1, n1))
    m1 = n1  # Imposed by requirement that orientation of molecule at origin is (000)
    for m2 in range(m2f, m2l + 1):
        m = m1 + m2
        cg = auxiliary_routines.cgcoef(l1, l2, l, m1, m2)
        dl2m2n2 = auxiliary_routines.dlm(l2, m2, n2, cz)
        dlm0 = auxiliary_routines.dlm(l, m, 0, ctl)
        phisum = f32(m2) * phi + f32(n2) * chi - f32(m) * phil  # Reverse signs because of complex conjugate
        prod = phs * cg * dl2m2n2 * dlm0
        sumc = sumc + prod * cosd_intrinsic(phisum)
        sums = sums + phssine * prod * sind_intrinsic(phisum)
    return sumc, sums


@routine
def suminterp(grw, rrat, thrat, nro, nrp, nrm, ntho, nthp, nthm, np0, nphp, nphm, ntdim, npdim, nrdim):
    """plotsharm_routines.f:suminterp -- average over phi of the quadratic interpolation
    of grw(ntdim,npdim,nrdim) in r and theta.  Returns (suminterp, nphp, nphm)."""
    sum_ = f32(0.0)
    nph = nphp - nphm + 1
    if nph > np0:
        nphm = 1
        nphp = np0
        nph = np0
    for npht in range(nphm, nphp + 1):
        npho = npht
        # Apply cyclic boundary conditions on phi
        while npho > np0:
            npho = npho - np0
        while npho < 1:
            npho = npho + np0
        #
        # quadratic interpolation along r
        #
        delta = rrat
        dp = r4(1.0) + delta
        dm = r4(1.0) - delta
        dpm = dm * dp
        delta2 = r4(0.5) * delta
        amo = delta2 * (dp * grw[nthm, npho, nrp] - dm * grw[nthm, npho, nrm]) + dpm * grw[nthm, npho, nro]
        aoo = delta2 * (dp * grw[ntho, npho, nrp] - dm * grw[ntho, npho, nrm]) + dpm * grw[ntho, npho, nro]
        apo = delta2 * (dp * grw[nthp, npho, nrp] - dm * grw[nthp, npho, nrm]) + dpm * grw[nthp, npho, nro]
        #
        # now interpolate in theta
        #
        delta = thrat
        dp = r4(1.0) + delta
        dm = r4(1.0) - delta
        dpm = dm * dp
        delta2 = r4(0.5) * delta
        ao = delta2 * (dp * apo - dm * amo) + dpm * aoo
        sum_ = sum_ + ao
    suminterp = sum_ / f32(nph)
    return suminterp, nphp, nphm
