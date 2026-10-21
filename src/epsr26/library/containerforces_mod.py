"""Translation of fortran_src/containerforces_mod.f90: the routines of module containerforces_mod.

A pair potential line in the ``.pcof`` file may give a maximum distance ``rmaxex`` and
amplitude ``pcoefexmax`` for a (Q-type, real) atom pair; every atom of the real type is
then pulled towards the nearest Q atom of that type once it is further than ``rmaxex``.

The module data is in ``common/containerforces_mod.py``.
"""

from ..common import containerforces_mod as cf  # Fortran module data
from ..fortran import expf, f32, logf, sqrtf
from ..trace import routine
from . import minimagec, reallocation_routines


@routine
def get_containerforce_list(nmol, ic, rmaxex, it, jt, itype, qtype, natom, qtypea, ntype):
    """containerforces_mod.f90:get_containerforce_list -- list the central (Q-type) atoms
    of container force ic (pair types it, jt) if rmaxex > 0.
    """
    # NB: typesecondary/typecentral are not initialised in the Fortran when no type
    # matches; 0 is used (Rust containerforces.rs:43).
    typesecondary = 0
    typecentral = 0

    # initialise numcontainerforces if ic .eq. 1
    if ic == 1:
        cf.numcontainerforces = 0
        cf.ergold_container = reallocation_routines.reallocate1d_r(cf.ergold_container, nmol)
        cf.ergnew_container = reallocation_routines.reallocate1d_r(cf.ergnew_container, nmol)
        for im in range(1, nmol + 1):
            cf.ergold_container[im] = 0.0
            cf.ergnew_container[im] = 0.0
    cf.numpairs = reallocation_routines.reallocate1d_i(cf.numpairs, ic)
    cf.numpairs[ic] = 0
    if rmaxex <= 0.0:
        return  # Will only count types if rmaxex .gt. 0.0)
    cf.numcontainerforces = cf.numcontainerforces + 1
    # First determine the type of the central and secondary atoms
    for nt in range(1, ntype + 1):
        if qtypea[nt] == 0:
            if nt == it:
                typesecondary = it
            elif nt == jt:
                typesecondary = jt
        else:
            if nt == it:
                typecentral = it
            elif nt == jt:
                typecentral = jt
    # Count the number of atoms of types it and jt
    numcentral = 0
    numsecondary = 0
    for ia in range(1, natom + 1):
        if qtype[ia] == 1 and (itype[ia] == it or itype[ia] == jt):
            numcentral = numcentral + 1
        elif qtype[ia] == 0 and (itype[ia] == it or itype[ia] == jt):
            numsecondary = numsecondary + 1
    cf.numpairs[ic] = numcentral * numsecondary
    if cf.numpairs[ic] > 0:
        n = cf.numcontainerforces
        cf.pairrefstore = reallocation_routines.reallocate2d_i(cf.pairrefstore, 2, n)
        cf.numcentraltype = reallocation_routines.reallocate1d_i(cf.numcentraltype, n)
        cf.numsecondarytype = reallocation_routines.reallocate1d_i(cf.numsecondarytype, n)
        cf.centraltype = reallocation_routines.reallocate2d_i(cf.centraltype, numcentral, n)
        # Store the secondary atom type number for this atom pair - this will be used to identify secondary atoms
        cf.pairrefstore[1, n] = typecentral
        cf.pairrefstore[2, n] = typesecondary
        cf.numcentraltype[n] = numcentral
        cf.numsecondarytype[n] = numsecondary
        # Now assign the values
        numcentral = 0
        for ia in range(1, natom + 1):
            if qtype[ia] == 1 and (itype[ia] == it or itype[ia] == jt):
                numcentral = numcentral + 1
                cf.centraltype[numcentral, n] = ia
    else:
        # Just in case no pairs were found in fact - unlikely.
        cf.numpairs[ic] = 0
        if cf.numcontainerforces > 0:
            cf.numcontainerforces = cf.numcontainerforces - 1


@routine
def en_containerforce(im, naf, nal, xadd, yadd, zadd, xa, ya, za, itype, mref, npairref, irefpcof, rmaxex,
                      pcoefexmax, expecf, expecfpwr, dmax, harmonic):
    """containerforces_mod.f90:en_containerforce -- container energy of molecule im.

    Calculates a quadratic energy if a specified atom pair has a minimum separation greater
    than specified rmaxex.  Returns en_containerforce.
    """
    en_containerforce = f32(0.0)
    if cf.numcontainerforces == 0:
        return en_containerforce
    pairrefstore = cf.pairrefstore
    numcentraltype = cf.numcentraltype
    centraltype = cf.centraltype
    sume = f32(0.0)
    # Limits the gradient on the container force to prevent the energy become too large
    gradlimit = expecfpwr
    # If gradlimit is zero or less, then make the gradient limit very large - this is to make
    # program compatible with previous versions
    if gradlimit > 0.0:
        loggradlimit = logf(gradlimit)
    else:
        gradlimit = f32(1.0)
        loggradlimit = f32(0.0)
    for ic in range(1, cf.numcontainerforces + 1):
        for ia in range(naf, nal + 1):
            if itype[ia] == pairrefstore[2, ic]:
                icref = npairref[itype[ia], pairrefstore[1, ic]]
                irefp = irefpcof[icref]
                expecfuse = expecf * sqrtf(abs(rmaxex[irefp]))
                mindist = f32(0.0)
                x1 = xadd[im] + xa[ia]
                y1 = yadd[im] + ya[ia]
                z1 = zadd[im] + za[ia]
                # Need the minimum distance between ia and all the "central" atoms
                for jc in range(1, numcentraltype[ic] + 1):
                    ja = centraltype[jc, ic]
                    jm = mref[ja]
                    x2 = xadd[jm] + xa[ja]
                    y2 = yadd[jm] + ya[ja]
                    z2 = zadd[jm] + za[ja]
                    xdif = x2 - x1
                    ydif = y2 - y1
                    zdif = z2 - z1
                    # apply minimum image convention
                    xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                    dsq = xdif * xdif + ydif * ydif + zdif * zdif
                    if mindist == 0.0 or dsq < mindist:
                        mindist = dsq
                dtest = sqrtf(mindist) * dmax
                if not harmonic:
                    rtest = rmaxex[irefp]
                    ddiff = (dtest - rtest) / expecfuse
                    if gradlimit > 0.0 and ddiff >= loggradlimit:
                        # Put linear term to prevent divergent exponential values
                        sume = sume + abs(pcoefexmax[irefp]) * gradlimit * (1.0 + ddiff - loggradlimit)
                    else:
                        sume = sume + abs(pcoefexmax[irefp]) * expf(ddiff)
                else:
                    rtest = rmaxex[irefp] - expecfuse
                    ddiff = (dtest - rtest) / expecfuse
                    if ddiff > 0.0:
                        if gradlimit > 0.0 and ddiff > gradlimit:
                            sume = sume + 0.5 * abs(pcoefexmax[irefp]) * gradlimit * (2 * ddiff - gradlimit)
                        else:
                            sume = sume + 0.5 * abs(pcoefexmax[irefp]) * ddiff * ddiff
                # if(isnan(sume)) then
                #    write(6,*) 'en_containerforce> ',im,itype(ia),irefp,rtest,ddiff,gradlimit,pcoefexmax(irefp)
                #    stop
                # end if
    en_containerforce = sume
    return en_containerforce


@routine
def dist_containerforce(im, naf, nal, xadd, yadd, zadd, xa, ya, za, itype, mref, npairref, irefpcof, rmaxex,
                        potexmax, dmax, addwav):
    """containerforces_mod.f90:dist_containerforce -- for each requested container force,
    calculates whether atom pair exceeds maximum distance.  If so, increments weighting on
    this term (potexmax, in place).
    """
    if cf.numcontainerforces == 0:
        return
    pairrefstore = cf.pairrefstore
    numcentraltype = cf.numcentraltype
    centraltype = cf.centraltype
    for ic in range(1, cf.numcontainerforces + 1):
        # addpcoef=addwav/real(numsecondarytype(ic))
        addpcoef = addwav
        for ia in range(naf, nal + 1):
            if itype[ia] == pairrefstore[2, ic]:
                icref = npairref[itype[ia], pairrefstore[1, ic]]
                irefp = irefpcof[icref]
                mindist = f32(0.0)
                x1 = xadd[im] + xa[ia]
                y1 = yadd[im] + ya[ia]
                z1 = zadd[im] + za[ia]
                # Need the minimum distance between ia and all the "central" atoms
                for jc in range(1, numcentraltype[ic] + 1):
                    ja = centraltype[jc, ic]
                    jm = mref[ja]
                    x2 = xadd[jm] + xa[ja]
                    y2 = yadd[jm] + ya[ja]
                    z2 = zadd[jm] + za[ja]
                    xdif = x2 - x1
                    ydif = y2 - y1
                    zdif = z2 - z1
                    # apply minimum image convention
                    xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                    dsq = xdif * xdif + ydif * ydif + zdif * zdif
                    if mindist == 0.0 or dsq < mindist:
                        mindist = dsq
                dtest = sqrtf(mindist) * dmax
                if dtest > rmaxex[irefp]:
                    potexmax[icref] = potexmax[icref] + addpcoef
