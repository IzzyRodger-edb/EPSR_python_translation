"""Translation of fortran_src/molfit.f: one Monte Carlo cycle and the energy routines.

``molfit`` visits every molecule once (cell by cell) and gives it a whole-molecule rotation
and translation, accepted with the Metropolis criterion on the change of intermolecular
energy (``en_change``).  Every ``intra`` cycles the cycle instead "shakes" the molecules
(``movei`` + ``en_change_intra``) and rotates their rotational groups (``rotati``).  Step
sizes are adapted per molecule type towards the acceptance ratio ``rejrate``.

These are the innermost loops of EPSR.  Module/COMMON arrays are aliased into locals for
speed where no reallocation can happen inside the routine; the arithmetic is the
Fortran's, statement by statement, in real*4.
"""

import math
import time

from ..common import cellneighbourlist as cnl  # Fortran module data (cellneighbourlist)
from ..common import clustercom, epsr16, epsrinput_par, molecom, potential_par
from ..common import containerforces_mod as cf_d  # Fortran module data (containerforces_mod)
from ..common import coordconstraints_mod as cc_d  # Fortran module data (coordconstraints_mod)
from ..common import ergpairreference as erg_d  # Fortran module data (ergpairreference)
from ..common import molecom_mod  # Fortran module data
from ..common import threebody_mod as tb_d  # Fortran module data (threebody_mod)
from ..common import timer  # Fortran module data (module timer, ato_manipulation_rev4.f)
from ..fortran import PI4ATAN, cosf, expf, f32, idiv, ifix, powf, r4, real_array, sign, sinf, sqrtf
from ..trace import routine
from . import (
    ato_manipulation_rev4,
    auxiliary_routines,
    containerforces_mod,
    coordconstraint_mod,
    ergpairreference,
    minimagec,
    threebody_mod,
)

# real, parameter :: lim=15.0,limit=exp(lim) (en_morse_mod1, en_morse_mod2) and exp(30.0)
# (en_gauss, en_harmonic): gfortran folds exp() of a constant at compile time (correctly
# rounded; equal to glibc expf for these arguments).
_LIM = r4(15.0)
_LIMIT = f32(math.exp(15.0))
_EXP30 = f32(math.exp(30.0))


def _system_clock8():
    """``call system_clock(count)`` with an integer(kind=8) count: gfortran gives the
    CLOCK_MONOTONIC time in nanoseconds (count_rate 1000000000).  Timing only."""
    return time.monotonic_ns()


@routine
def molfit(fadds, faddbs, nisum, intra, rotfreq, ireset):
    """molfit.f:molfit -- one Monte Carlo cycle.

    If nisum >= intra, intramolecular moves and internal rotations of every molecule;
    otherwise whole-molecule rotations and translations, cell by cell (with internal
    rotations every rotfreq cycles).  Then adapts the step sizes.
    Returns (fadds, faddbs, nisum).
    """
    # real, dimension(:), allocatable :: xr,yr,zr
    xr = None
    yr = None
    zr = None
    # NB: sumevib, sumeang, sumedih and actualstep are uninitialised locals in the Fortran;
    # en_change_intra / rotati set them before they are used.
    sumevib = f32(0.0)
    sumeang = f32(0.0)
    sumedih = f32(0.0)
    actualstep = f32(0.0)

    # COMMON arrays (never reallocated)
    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za
    xadd = clustercom.xadd
    yadd = clustercom.yadd
    zadd = clustercom.zadd
    naf = clustercom.naf
    nato = clustercom.nato
    ngrot = clustercom.ngrot
    nrf = clustercom.nrf
    nrl = clustercom.nrl
    mtype = clustercom.mtype
    nmovei = clustercom.nmovei
    ninter = clustercom.ninter
    nwhole = clustercom.nwhole
    nmovem = clustercom.nmovem
    nstepreject = clustercom.nstepreject
    tetheredatom = clustercom.tetheredatom
    interntrans = clustercom.interntrans
    internrot = clustercom.internrot
    wholerot = clustercom.wholerot
    wholetrans = clustercom.wholetrans
    ran1 = auxiliary_routines.ran1

    temp = r4(1.0) / molecom.simt
    # tol=tolout/temp
    molecom.tol = molecom.tolout / clustercom.dmax
    for j in range(1, 3):
        for imoltype in range(1, clustercom.nmoltype + 1):
            nmovei[j, imoltype] = 0
            ninter[j, imoltype] = 0
            nwhole[j, imoltype] = 0
            nmovem[j, imoltype] = 0

    # Set the maximum translational step size

    trstepmax = r4(1.0) / f32(molecom.natomr)
    trstepmax = powf(trstepmax, r4(1.0) / r4(3.0))
    # Limit the smallest step size to that which will make a difference within the available precision.
    trstepmin = max(0.0001 * trstepmax, r4(1.0e-6))
    rotstepmax = r4(1.0)
    rotstepmin = r4(1.0e-5)
    fadd = f32(0.0)
    faddb = f32(0.0)
    # Ensure correct weights on angle and dihedral terms are set
    molecom.tempang = molecom.tempangwt
    molecom.tempdih = molecom.tempdihwt * molecom.tempb
    #
    # if step size is non-zero, then do internal molecule structure refinement. This is to
    # prevent the build up of round-off errors from allowing the molecule to
    # become too distorted.
    #
    if nisum >= intra:

        for im in range(1, clustercom.nmol + 1):
            ergpairreference.make_neighbour_list_from_mols(im)
            n = int(nato[im])
            if xr is None:
                xr = real_array(n)
                yr = real_array(n)
                zr = real_array(n)
            elif xr.shape[0] - 1 < n:
                xr = None
                yr = None
                zr = None
                xr = real_array(n)
                yr = real_array(n)
                zr = real_array(n)
            mt = int(mtype[im])
            if nato[im] > 1 and molecom.stepmi > 0.0:
                # write(6,*) im,nergpairm(im),molergpair(1,im),ergpair(ergpairref(1,im))
                # first save the coordinates of this molecule
                stepmithis = interntrans[mt]
                if molecom.stepmi > 0.0 and stepmithis == 0.0:
                    stepmithis = trstepmin
                interntrans[mt] = stepmithis
                molecom.xs = xadd[im]
                molecom.ys = yadd[im]
                molecom.zs = zadd[im]
                nat = naf[im] - 1
                for ia in range(1, n + 1):
                    nat = nat + 1
                    xr[ia] = xa[nat]
                    yr[ia] = ya[nat]
                    zr[ia] = za[nat]
                fdiff = f32(0.0)  # noqa: F841
                fdiffb = f32(0.0)

                # move particles and calculate new chisq
                # do not move the first atom of each molecule if this is the only
                # atom in the molecule!
                nmovei[1, mt] = nmovei[1, mt] + 1
                movei(im, stepmithis)
                # now adjust the coordinates so that the centre of mass has not changed
                # call centreofmassi(im) !This step appears to produce bad overlaps - use after
                # accepting move, otherwise it could cause errors!
                fdiffb, sumevib, sumeang, sumedih = en_change_intra(im, sumevib, sumeang, sumedih)
                ediff = fdiffb * temp
                if ediff <= 0.0:
                    etest = r4(10.0)
                elif ediff < 20.:
                    etest = expf(-ediff)
                else:
                    etest = r4(0.)
                rtest = ran1()
                if rtest < etest:
                    faddb = faddb + fdiffb
                    ergpairreference.save_energies()
                    # now adjust the coordinates so that the centre of mass has not changed
                    ato_manipulation_rev4.centreofmassi(im)
                else:
                    nmovei[2, mt] = nmovei[2, mt] + 1
                    nat = naf[im] - 1
                    for ia in range(1, n + 1):
                        nat = nat + 1
                        xa[nat] = xr[ia]
                        ya[nat] = yr[ia]
                        za[nat] = zr[ia]
                # now try rotating internal groups if no of rotational groups is greater than 0
                # This doesn't work because you need to include the interactions with neighbouring
                # molecules
                # Save all the relevant pair energies
                ergpairreference.save_erg_pairs(im)
            # Now do the internal rotations
            if molecom.stepri > 0.0 and ngrot[im] > 0:
                # if(.false.) then !Eliminate this step for the time being.
                if internrot[mt] == 0.0:
                    internrot[mt] = rotstepmin
                #
                # perform internal rotations
                #
                ninter[1, mt] = ninter[1, mt] + 1
                # Save the coordinates of this group
                nat = naf[im] - 1
                for ia in range(1, n + 1):
                    nat = nat + 1
                    xr[ia] = xa[nat]
                    yr[ia] = ya[nat]
                    zr[ia] = za[nat]
                molecom.xs = xadd[im]
                molecom.ys = yadd[im]
                molecom.zs = zadd[im]
                fdiff = f32(0.0)
                fdiffb = f32(0.0)
                for ig in range(nrf[im], nrl[im] + 1):
                    actualstep = rotati(im, ig, internrot[mt], actualstep)
                # Check that change in tethering energy is acceptable
                if check_tethered(im, tetheredatom[im]):
                    fdiffb, sumevib, sumeang, sumedih = en_change_intra(im, sumevib, sumeang, sumedih)
                    fdiff = en_change(im)
                    ediff = (fdiff + fdiffb) * temp
                    if ediff <= 0.0:
                        etest = r4(10.0)
                    elif ediff < 20.:
                        etest = expf(-ediff)
                    else:
                        etest = r4(0.)
                    if ran1() < etest:
                        nstepreject[im] = 0
                        fadd = fadd + fdiff
                        faddb = faddb + fdiffb
                        ergpairreference.save_energies()
                        threebody_mod.save_threeb_energies(im)
                        coordconstraint_mod.save_coord_energies(im)
                        # Adjust the coordinates of centre of mass to correct for any change
                    else:
                        # Move rejected - restore coordinates
                        nstepreject[im] = nstepreject[im] + 1
                        # Restore the coordinates of this group
                        nat = naf[im] - 1
                        for ia in range(1, n + 1):
                            nat = nat + 1
                            xa[nat] = xr[ia]
                            ya[nat] = yr[ia]
                            za[nat] = zr[ia]
                        ninter[2, mt] = ninter[2, mt] + 1
                    # Save all the relevant pair energies
                    ergpairreference.save_erg_pairs(im)
                else:
                    # Move rejected - restore coordinates
                    nstepreject[im] = nstepreject[im] + 1
                    # Restore the coordinates of this group
                    nat = naf[im] - 1
                    for ia in range(1, n + 1):
                        nat = nat + 1
                        xa[nat] = xr[ia]
                        ya[nat] = yr[ia]
                        za[nat] = zr[ia]
                    xadd[im] = molecom.xs
                    yadd[im] = molecom.ys
                    zadd[im] = molecom.zs
                    ninter[2, mt] = ninter[2, mt] + 1
            else:
                if ngrot[im] == 0:
                    internrot[mt] = 0.0
            # write(6,*) 'molfit> Vibrational, bond angle and dihedral energies: ',im,sumevib,sumeang,sumedih
    else:
        # Step through cells and molecules
        nmolsincell = cnl.nmolsincell
        molsincell = cnl.molsincell
        for ic in range(1, cnl.ncell + 1):
            # write(6,*) 'molfit-rot> Cell ',ic,fadd,faddb
            if nmolsincell[ic] > 0:
                for imref in range(1, nmolsincell[ic] + 1):
                    im = int(molsincell[imref, ic])
                    n = int(nato[im])
                    if xr is None:
                        xr = real_array(n)
                        yr = real_array(n)
                        zr = real_array(n)
                    elif xr.shape[0] - 1 < n:
                        xr = None
                        yr = None
                        zr = None
                        xr = real_array(n)
                        yr = real_array(n)
                        zr = real_array(n)
                    timer.time1 = _system_clock8()
                    ergpairreference.make_neighbour_list_from_mols(im)
                    timer.time2 = _system_clock8()
                    timer.sumtime1 = timer.sumtime1 + idiv(timer.time2 - timer.time1, 1000)
                    mt = int(mtype[im])
                    # Only do this every rotfreq iterations
                    if molecom.stepri > 0.0 and ngrot[im] > 0 and nisum == rotfreq * idiv(nisum, rotfreq):
                        if internrot[mt] == 0.0:
                            internrot[mt] = rotstepmin
                        #
                        # perform internal rotations
                        #
                        ninter[1, mt] = ninter[1, mt] + 1
                        # Save the coordinates of this group
                        nat = naf[im] - 1
                        for ia in range(1, n + 1):
                            nat = nat + 1
                            xr[ia] = xa[nat]
                            yr[ia] = ya[nat]
                            zr[ia] = za[nat]
                        molecom.xs = xadd[im]
                        molecom.ys = yadd[im]
                        molecom.zs = zadd[im]
                        fdiff = f32(0.0)
                        fdiffb = f32(0.0)
                        for ig in range(nrf[im], nrl[im] + 1):
                            actualstep = rotati(im, ig, internrot[mt], actualstep)
                        # Check that change in tethering energy is acceptable
                        if check_tethered(im, tetheredatom[im]):
                            fdiffb, sumevib, sumeang, sumedih = en_change_intra(im, sumevib, sumeang, sumedih)
                            fdiff = en_change(im)
                            ediff = (fdiff + fdiffb) * temp
                            if ediff <= 0.0:
                                etest = r4(10.0)
                            elif ediff < 20.:
                                etest = expf(-ediff)
                            else:
                                etest = r4(0.)
                            if ran1() < etest:
                                nstepreject[im] = 0
                                fadd = fadd + fdiff
                                faddb = faddb + fdiffb
                                ergpairreference.save_energies()
                                threebody_mod.save_threeb_energies(im)
                                coordconstraint_mod.save_coord_energies(im)
                            else:
                                # Move rejected - restore coordinates
                                nstepreject[im] = nstepreject[im] + 1
                                # Restore the coordinates of this group
                                nat = naf[im] - 1
                                for ia in range(1, n + 1):
                                    nat = nat + 1
                                    xa[nat] = xr[ia]
                                    ya[nat] = yr[ia]
                                    za[nat] = zr[ia]
                                ninter[2, mt] = ninter[2, mt] + 1
                            # Save all the relevant pair energies
                            ergpairreference.save_erg_pairs(im)
                        else:
                            # Move rejected - restore coordinates
                            nstepreject[im] = nstepreject[im] + 1
                            # Restore the coordinates of this group
                            nat = naf[im] - 1
                            for ia in range(1, n + 1):
                                nat = nat + 1
                                xa[nat] = xr[ia]
                                ya[nat] = yr[ia]
                                za[nat] = zr[ia]
                            xadd[im] = molecom.xs
                            yadd[im] = molecom.ys
                            zadd[im] = molecom.zs
                            ninter[2, mt] = ninter[2, mt] + 1
                    else:
                        if ngrot[im] == 0:
                            internrot[mt] = 0.0
                    # Rotate and translate whole molecules and calculate energy change. Only rotate
                    # molecules that have more than one atom
                    # First save the coordinates of this molecule
                    nat = naf[im] - 1
                    for ia in range(1, n + 1):
                        nat = nat + 1
                        xr[ia] = xa[nat]
                        yr[ia] = ya[nat]
                        zr[ia] = za[nat]
                    molecom.xs = xadd[im]
                    molecom.ys = yadd[im]
                    molecom.zs = zadd[im]
                    fdiff = f32(0.0)
                    fdiffb = f32(0.0)
                    if molecom.stepr > 0.0 and nato[im] > 1:
                        if wholerot[mt] <= 0.0:
                            wholerot[mt] = rotstepmin
                        nwhole[1, mt] = nwhole[1, mt] + 1
                        rotatm(im, wholerot[mt])
                    else:
                        wholerot[mt] = 0.0
                    if molecom.stepm > 0.0:
                        # now translate each molecule whole
                        if wholetrans[mt] <= 0.0:
                            wholetrans[mt] = trstepmin
                        nmovem[1, mt] = nmovem[1, mt] + 1
                        # NB (FORTRAN_BUGS.md section 2, reproduced): movem doubles its stepsz
                        # dummy until the move changes the position; the actual argument is
                        # wholetrans(mtype(im)), so the doubled step is stored permanently.
                        wholetrans[mt] = movem(im, wholetrans[mt])
                    else:
                        wholetrans[mt] = 0.0
                    # Check that change in tethering energy is acceptable
                    if check_tethered(im, tetheredatom[im]):
                        # fdiffb=en_change_intra(im,sumevib,sumeang,sumedih)
                        fdiffb = f32(0.0)
                        fdiff = en_change(im)
                        ediff = (fdiff + fdiffb) * temp
                        if ediff <= 0.0:
                            etest = r4(10.0)
                        elif ediff < 20.:
                            etest = expf(-ediff)
                        else:
                            etest = r4(0.)
                        if ran1() < etest:
                            # write(6,*) im,fdiffb,fdiff,ediff,etest
                            nstepreject[im] = 0
                            fadd = fadd + fdiff
                            faddb = faddb + fdiffb
                            ergpairreference.save_energies()
                            threebody_mod.save_threeb_energies(im)
                        else:
                            # Move rejected - restore coordinates
                            nstepreject[im] = nstepreject[im] + 1
                            nat = naf[im] - 1
                            for ia in range(1, n + 1):
                                nat = nat + 1
                                xa[nat] = xr[ia]
                                ya[nat] = yr[ia]
                                za[nat] = zr[ia]
                            xadd[im] = molecom.xs
                            yadd[im] = molecom.ys
                            zadd[im] = molecom.zs
                            if molecom.stepr > 0:
                                nwhole[2, mt] = nwhole[2, mt] + 1
                            if molecom.stepm > 0:
                                nmovem[2, mt] = nmovem[2, mt] + 1
                        # Save all the relevant pair energies
                        ergpairreference.save_erg_pairs(im)
                    else:
                        # Move rejected - restore coordinates
                        nstepreject[im] = nstepreject[im] + 1
                        nat = naf[im] - 1
                        for ia in range(1, n + 1):
                            nat = nat + 1
                            xa[nat] = xr[ia]
                            ya[nat] = yr[ia]
                            za[nat] = zr[ia]
                        xadd[im] = molecom.xs
                        yadd[im] = molecom.ys
                        zadd[im] = molecom.zs
                        if molecom.stepr > 0:
                            nwhole[2, mt] = nwhole[2, mt] + 1
                        if molecom.stepm > 0:
                            nmovem[2, mt] = nmovem[2, mt] + 1
        # Set the weighting on the dihedral energy
        # (commented-out code in the Fortran: recompute sumevib, sumeang, sumedih over all
        # molecules and print the relative bond angle and dihedral energies)

    # do im=1,nmol
    #     call minimagec(xadd(im),yadd(im),zadd(im))
    #     if(tethered(im)) call minimagec(phix(im),phiy(im),phiz(im))
    # end do

    if nisum >= intra:
        nisum = 1
    else:
        nisum = nisum + 1

    # Adjust step sizes to give required rejection rate

    if ireset != 1:

        rejrate = epsr16.rejrate
        for imoltype in range(1, clustercom.nmoltype + 1):

            if nmovei[2, imoltype] > 0:
                rat = rejrate * f32(nmovei[1, imoltype]) / f32(nmovei[2, imoltype])
                if rat > 1.1:
                    rat = r4(1.1)
                if rat < 0.9:
                    rat = r4(0.9)
            else:
                if nmovei[1, imoltype] > 0:
                    rat = r4(1.1)
                else:
                    rat = r4(1.0)
            interntrans[imoltype] = interntrans[imoltype] * rat
            if interntrans[imoltype] > trstepmax:
                interntrans[imoltype] = trstepmax
            if interntrans[imoltype] < trstepmin and interntrans[imoltype] != 0.0:
                interntrans[imoltype] = trstepmin

            if ninter[2, imoltype] > 0:
                rat = rejrate * f32(ninter[1, imoltype]) / f32(ninter[2, imoltype])
                if rat > 1.1:
                    rat = r4(1.1)
                if rat < 0.9:
                    rat = r4(0.9)
            else:
                if ninter[1, imoltype] > 0:
                    rat = r4(1.1)
                else:
                    rat = r4(1.0)
            internrot[imoltype] = internrot[imoltype] * rat
            if internrot[imoltype] > rotstepmax:
                internrot[imoltype] = rotstepmax
            if internrot[imoltype] < rotstepmin and internrot[imoltype] != 0.0:
                internrot[imoltype] = rotstepmin

            if nwhole[2, imoltype] > 0:
                rat = rejrate * f32(nwhole[1, imoltype]) / f32(nwhole[2, imoltype])
                if rat > 1.1:
                    rat = r4(1.1)
                if rat < 0.9:
                    rat = r4(0.9)
            else:
                if nwhole[1, imoltype] > 0:
                    rat = r4(1.1)
                else:
                    rat = r4(1.0)
            wholerot[imoltype] = wholerot[imoltype] * rat
            if wholerot[imoltype] > rotstepmax:
                wholerot[imoltype] = rotstepmax
            if wholerot[imoltype] < rotstepmin and wholerot[imoltype] != 0.0:
                wholerot[imoltype] = rotstepmin

            if nmovem[2, imoltype] > 0:
                rat = rejrate * f32(nmovem[1, imoltype]) / f32(nmovem[2, imoltype])
                if rat > 1.1:
                    rat = r4(1.1)
                if rat < 0.9:
                    rat = r4(0.9)
            else:
                if nmovem[1, imoltype] > 0:
                    rat = r4(1.1)
                else:
                    rat = r4(1.0)
            wholetrans[imoltype] = wholetrans[imoltype] * rat
            if wholetrans[imoltype] > trstepmax:
                wholetrans[imoltype] = trstepmax
            if wholetrans[imoltype] < trstepmin and wholetrans[imoltype] != 0.0:
                wholetrans[imoltype] = trstepmin

    fadds = fadds + fadd / f32(clustercom.nmol)
    faddbs = faddbs + faddb / f32(clustercom.nmol)

    return fadds, faddbs, nisum


@routine
def movei(im, stepmithis):
    """molfit.f:movei -- move every atom of molecule im by a random amount."""
    #
    # move the atoms
    #
    for ia in range(clustercom.naf[im], clustercom.nal[im] + 1):
        movea(ia, stepmithis)
    return


@routine
def movea(ia, stepmithis):
    """molfit.f:movea -- move atom ia by random amounts (z, then y, then x)."""
    ran1 = auxiliary_routines.ran1
    # amass=sqrt(amassa(itype(ia)))
    # stepsz=stepmithis/amass
    stepsz = stepmithis
    # move atoms by random amounts
    molecom.za[ia] = molecom.za[ia] + stepsz * (r4(2.) * ran1() - 1)
    molecom.ya[ia] = molecom.ya[ia] + stepsz * (r4(2.) * ran1() - 1)
    molecom.xa[ia] = molecom.xa[ia] + stepsz * (r4(2.) * ran1() - 1)
    return


@routine
def interm(im, fdiff, fdiffb):
    """molfit.f:interm -- empty (its body is commented out in the Fortran)."""
    #
    # perform internal rotations
    #
    # do ig=nrf(im),nrl(im)
    #    call rotati(ig,internrot(mtype(im)))
    # end do
    #
    # calculate change in energy
    #
    # fdiff=0.0
    # fdiffb=0.0
    # call en_change_intra(im,fdiffb)
    # call en_change(im,fdiff)
    return


@routine
def rotati(im, ig, steprithis, actualstep):
    """molfit.f:rotati -- rotate the atoms of rotational group ig about its axis by a
    random angle.  Returns actualstep (the angle used).
    """
    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za
    labaxi = molecom_mod.labaxi
    labrot = molecom_mod.labrot
    # pvector is the vector about which the rotation occurs
    pvector = real_array(3)
    # xvector is the vector from centre of rotation to an individual atom
    # xpvector is the component of this vector along the axis of rotation
    # xrvector is the component of this vector perpendicular to the axis of rotation
    # xrpvector is the vector perpendicular to both axis of rotation and xrvector
    xvector = real_array(3)
    xpvector = real_array(3)
    xrvector = real_array(3)
    xrpvector = real_array(3)
    pi = PI4ATAN

    # define random new orientation
    phi = pi * steprithis * (2 * auxiliary_routines.ran1() - 1)
    actualstep = phi
    cphi = cosf(phi)
    # sphi=sin(phi)
    sphi = sign(sqrtf(1.0 - cphi * cphi), phi)
    # define length of axis between the two points - the direction 1 -> 2 is
    # treated as the positive going transformed z axis
    pvector[1] = xa[labaxi[ig, 2]] - xa[labaxi[ig, 1]]
    pvector[2] = ya[labaxi[ig, 2]] - ya[labaxi[ig, 1]]
    pvector[3] = za[labaxi[ig, 2]] - za[labaxi[ig, 1]]
    rab2 = auxiliary_routines.scapro(pvector, pvector)
    if rab2 > 0:
        rab = sqrtf(rab2)
        # Unit vector along axis of rotation
        for i in range(1, 4):
            pvector[i] = pvector[i] / rab
        # step through atoms in group
        for ia in range(molecom_mod.rotf[ig], molecom_mod.rotl[ig] + 1):
            nat = labrot[ia]
            # first calculate vector from centre of rotation to each atom in group
            nax = labaxi[ig, 1]
            xvector[1] = xa[nat] - xa[nax]
            xvector[2] = ya[nat] - ya[nax]
            xvector[3] = za[nat] - za[nax]
            # Vector parallel to axis of rotation
            frac = auxiliary_routines.scapro(pvector, xvector)
            for i in range(1, 4):
                xpvector[i] = frac * pvector[i]
            # Vector perpendicular to axis of rotation
            for i in range(1, 4):
                xrvector[i] = xvector[i] - xpvector[i]
            # Vector perpendicular to plane defined by axis of rotation and xrvector
            auxiliary_routines.vecpro(pvector, xrvector, xrpvector)
            # Generate rotated coordinates. The vector component parallel to the axis of rotation
            # is not rotated
            xa[nat] = xa[nax] + xpvector[1]
            ya[nat] = ya[nax] + xpvector[2]
            za[nat] = za[nax] + xpvector[3]
            xa[nat] = xa[nat] + xrvector[1] * cphi + xrpvector[1] * sphi
            ya[nat] = ya[nat] + xrvector[2] * cphi + xrpvector[2] * sphi
            za[nat] = za[nat] + xrvector[3] * cphi + xrpvector[3] * sphi
    return actualstep


@routine
def wholem(im):
    """molfit.f:wholem -- rotate entire molecule im."""
    #
    # rotate entire molecule
    #
    rotatm(im, clustercom.wholerot[clustercom.mtype[im]])
    #
    # Calculate the change in energy
    #
    # fdiff=0.0
    # fdiffb=0.0
    # call en_change(im,fdiff)
    return


@routine
def rotatm(im, stepsz):
    """molfit.f:rotatm -- routine to rotate entire molecule (about z, then y, then x)."""
    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za
    ran1 = auxiliary_routines.ran1
    pi = PI4ATAN

    iref1 = clustercom.naf[im]
    iref2 = clustercom.nal[im]
    if iref2 > iref1:
        # define rotations about x,y or z-axes
        # rotate about z axis
        phi = pi * stepsz * (2 * ran1() - 1)
        cphi = cosf(phi)
        # sphi=sin(phi)
        sphi = sign(sqrtf(1.0 - cphi * cphi), phi)
        for ia in range(iref1, iref2 + 1):
            xasave = xa[ia] * cphi - ya[ia] * sphi
            ya[ia] = xa[ia] * sphi + ya[ia] * cphi
            xa[ia] = xasave
        # rotate about y-axis
        phi = pi * stepsz * (2 * ran1() - 1)
        cphi = cosf(phi)
        # sphi=sin(phi)
        sphi = sign(sqrtf(1.0 - cphi * cphi), phi)
        for ia in range(iref1, iref2 + 1):
            zasave = za[ia] * cphi - xa[ia] * sphi
            xa[ia] = za[ia] * sphi + xa[ia] * cphi
            za[ia] = zasave
        # rotate about x-axis
        phi = pi * stepsz * (2 * ran1() - 1)
        cphi = cosf(phi)
        # sphi=sin(phi)
        sphi = sign(sqrtf(1.0 - cphi * cphi), phi)
        for ia in range(iref1, iref2 + 1):
            yasave = ya[ia] * cphi - za[ia] * sphi
            za[ia] = ya[ia] * sphi + za[ia] * cphi
            ya[ia] = yasave
    return


@routine
def transm(im):
    """molfit.f:transm -- move molecule im."""
    #
    # move molecule
    #
    mt = clustercom.mtype[im]
    clustercom.wholetrans[mt] = movem(im, clustercom.wholetrans[mt])
    #
    # Calculate the change in energy
    #
    # fdiff=0.0
    # fdiffb=0.0
    # call en_change(im,fdiff)

    return


@routine
def movem(im, stepsz):
    """molfit.f:movem -- translate molecule im by random amounts (z, then y, then x).

    As in the Fortran, stepsz (a dummy passed by reference) is doubled until the move
    actually changes the position (FORTRAN_BUGS.md section 2, reproduced).
    Returns stepsz.
    """
    xadd = clustercom.xadd
    yadd = clustercom.yadd
    zadd = clustercom.zadd
    ran1 = auxiliary_routines.ran1
    # store present coordinates
    molecom.xs = xadd[im]
    molecom.ys = yadd[im]
    molecom.zs = zadd[im]
    # move atoms by random amounts
    zadd[im] = zadd[im] + stepsz * (r4(2.) * ran1() - 1.0)
    yadd[im] = yadd[im] + stepsz * (r4(2.) * ran1() - 1.0)
    xadd[im] = xadd[im] + stepsz * (r4(2.) * ran1() - 1.0)
    nochange = molecom.xs == xadd[im] and molecom.ys == yadd[im] and molecom.zs == zadd[im]
    while nochange:
        stepsz = stepsz * 2.0
        zadd[im] = zadd[im] + stepsz * (r4(2.) * ran1() - 1.0)
        yadd[im] = yadd[im] + stepsz * (r4(2.) * ran1() - 1.0)
        xadd[im] = xadd[im] + stepsz * (r4(2.) * ran1() - 1.0)
        nochange = molecom.xs == xadd[im] and molecom.ys == yadd[im] and molecom.zs == zadd[im]
    return stepsz


@routine
def en_change_intra(im, sumevib, sumeang, sumedih):
    """molfit.f:en_change_intra -- change in the intramolecular energy of molecule im
    relative to ergold(1).  Returns (en_change_intra, sumevib, sumeang, sumedih).
    """
    timer.time1 = _system_clock8()
    sumevib = f32(0.0)
    sumeang = f32(0.0)
    sumedih = f32(0.0)
    # Calculate the change in energy as a result of moving all the atoms in the molecule
    # NB: per SIGNATURES.md en_mole_intra returns only (value, sumedih); see
    # docs/SIGNATURE_EXCEPTIONS.md (sumevib, sumeang are not passed back).
    e, sumevib, sumeang, sumedih = en_mole_intra(im, sumevib, sumeang, sumedih)
    erg_d.ergnew[1] = e
    en_change_intra = erg_d.ergnew[1] - erg_d.ergold[1]
    # write(6,*) im,molergpair(1,im),ergpair(ergpairref(1,im)),ergold(1),ergnew(1),fdiffb
    timer.time2 = _system_clock8()
    timer.sumtime2 = timer.sumtime2 + idiv(timer.time2 - timer.time1, 1000)

    return en_change_intra, sumevib, sumeang, sumedih


@routine
def en_change_intra_rot(im, ig, sumevib, sumeang, sumedih):
    """molfit.f:en_change_intra_rot -- calculates the intra-molecular energy change of a
    molecule after rotating group ig, using the current coordinates.
    Returns (en_change_intra_rot, sumevib, sumeang, sumedih).
    """
    timer.time1 = _system_clock8()
    sumevib = f32(0.0)
    sumeang = f32(0.0)
    sumedih = f32(0.0)
    # Calculate the change in energy as a result of moving the atoms in rotational group ig in the molecule
    # The safest way - may not be the quickest - is simply calculate the full intramolecular energy of the
    # revised molecule
    e, sumevib, sumeang, sumedih = en_mole_intra(im, sumevib, sumeang, sumedih)
    erg_d.ergnew[1] = e
    en_change_intra_rot = erg_d.ergnew[1] - erg_d.ergold[1]
    timer.time2 = _system_clock8()
    timer.sumtime2 = timer.sumtime2 + idiv(timer.time2 - timer.time1, 1000)
    return en_change_intra_rot, sumevib, sumeang, sumedih


@routine
def en_change(im):
    """molfit.f:en_change -- calculates the change in energy as a result of moving all the
    atoms in molecule im (its neighbours nbref(2..nnbref), three-body, container and
    coordination-constraint terms).  Fills ergnew.  Returns en_change.
    """
    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za
    itype = molecom.itype
    qtype = molecom.qtype
    npairref = molecom.npairref
    mref = molecom.mref
    xadd = clustercom.xadd
    yadd = clustercom.yadd
    zadd = clustercom.zadd
    naf = clustercom.naf
    nal = clustercom.nal
    tethered = clustercom.tethered
    xstar = epsr16.xstar
    npt2 = epsr16.npt2
    # ergpairreference arrays: reallocated only by make_neighbour_list_from_mols
    ergnew = erg_d.ergnew
    ergold = erg_d.ergold
    nbref = erg_d.nbref
    nnbref = erg_d.nnbref
    coarsegraining = molecom_mod.coarsegraining
    scalefacmol = molecom_mod.scalefacmol
    ncoordconstraints = cc_d.ncoordconstraints
    check_coord_bonded = coordconstraint_mod.check_coord_bonded

    stepcl = epsr16.cellst / clustercom.dmax
    timer.time1 = _system_clock8()
    # !$omp parallel do private(inn)
    for inn in range(2, nnbref + 1):
        ergnew[inn] = 0.0
    # !$omp end parallel do
    # !$omp parallel do private(...) shared(...)  (serial with OMP_NUM_THREADS=1)
    # Step through list of neighbouring molecules
    naf_im = naf[im]
    nal_im = nal[im]
    for inn in range(2, nnbref + 1):
        jm = int(nbref[inn])
        # Scale factor on intermolecular spacings
        scalefac = f32(1.0)
        if coarsegraining:
            scalefac = 0.5 * (scalefacmol[im] + scalefacmol[jm])
        # Don't bother if both molecules are tethered and both atoms are Q-type.
        bothtethered = tethered[im] and tethered[jm]
        sume = f32(0.0)
        naf_jm = naf[jm]
        nal_jm = nal[jm]
        # Step through moved atoms on molecule im
        for ia in range(naf_im, nal_im + 1):
            it = itype[ia]
            x1 = xa[ia] + xadd[im]
            y1 = ya[ia] + yadd[im]
            z1 = za[ia] + zadd[im]
            qtypeia = qtype[ia] == 1 and bothtethered
            for ja in range(naf_jm, nal_jm + 1):
                bothqatoms = qtype[ja] == 1 and qtypeia
                # Don't calculate energy if both molecules are tethered and both atoms are Q types,
                # or the atoms form a coordination constraint.
                if not bothqatoms:
                    jt = itype[ja]
                    icref = npairref[it, jt]
                    x2 = xa[ja] + xadd[jm]
                    y2 = ya[ja] + yadd[jm]
                    z2 = za[ja] + zadd[jm]
                    # (check_coord_bonded is .false. without side effects when
                    # ncoordconstraints is 0; the test is skipped then, for speed)
                    if ncoordconstraints > 0 and check_coord_bonded(ia, ja, mref, itype):
                        sume = sume + en_atom_inter_scale_coordextra(x1, y1, z1, x2, y2, z2,
                                                                     stepcl, xstar, icref, npt2, scalefac)
                    else:
                        sume = sume + en_atom_inter_scale(x1, y1, z1, x2, y2, z2,
                                                          stepcl, xstar, icref, npt2, scalefac)
        ergnew[inn] = ergnew[inn] + sume
    # !$omp end parallel do
    fdiff = f32(0.0)
    for inn in range(2, nnbref + 1):
        fdiff = fdiff + ergnew[inn] - ergold[inn]
    if tb_d.nthreeb > 0:
        # Get the new triplet energy
        tb_d.ergnew_triplet[im] = threebody_mod.en_mole_triplets(im, xadd, yadd, zadd, clustercom.phix,
                                                                 clustercom.phiy, clustercom.phiz, xa, ya, za,
                                                                 mref, tethered)
        fdiff = fdiff + tb_d.ergnew_triplet[im] - tb_d.ergold_triplet[im]
        # IMPORTANT! Include the neighbours of im which may have been affected by the movement of im!
        if tb_d.nthreebneighbmol[im] > 0:
            for it in range(1, tb_d.nthreebneighbmol[im] + 1):
                imref = tb_d.threebneighbmol[it, im]
                tb_d.ergnew_triplet[imref] = threebody_mod.en_mole_triplets(imref, xadd, yadd, zadd,
                                                                            clustercom.phix, clustercom.phiy,
                                                                            clustercom.phiz, xa, ya, za,
                                                                            mref, tethered)
                fdiff = fdiff + tb_d.ergnew_triplet[imref] - tb_d.ergold_triplet[imref]
    # Get the new container force energy if requested
    if cf_d.numcontainerforces > 0:
        cf_d.ergnew_container[im] = containerforces_mod.en_containerforce(
            im, int(naf[im]), int(nal[im]), xadd, yadd, zadd, xa, ya, za,
            itype, mref, npairref, potential_par.irefpcof, potential_par.rmaxex, potential_par.pcoefexmax,
            epsrinput_par.expecf, epsrinput_par.expecfpwr, clustercom.dmax, potential_par.harmonic)
        # write(6,*) im,ergnew_container(im)
        fdiff = fdiff + cf_d.ergnew_container[im] - cf_d.ergold_container[im]
    if cc_d.ncoordconstraints > 0:
        cc_d.ergnew_coord[im] = coordconstraint_mod.en_mole_coord(im, naf, nal, xadd, yadd, zadd, xa, ya, za,
                                                                  mref, itype, qtype, npairref, molecom.amassr,
                                                                  tethered)
        fdiff = fdiff + cc_d.ergnew_coord[im] - cc_d.ergold_coord[im]
    en_change = fdiff
    timer.time2 = _system_clock8()
    timer.sumtime2 = timer.sumtime2 + idiv(timer.time2 - timer.time1, 1000)

    return en_change


@routine
def en_atom_intra(ia, ja, x1, y1, z1, x2, y2, z2, sumevib, sumeang):
    """molfit.f:en_atom_intra -- calculates the internal energy of atom pair ia,ja in a
    molecule (it is assumed this pair have been preselected not to be the same atom).
    Returns (en_atom_intra, sumevib, sumeang).
    """
    en_atom_intra = f32(0.0)
    # If ja is an image atom of ia, do nothing: image atoms are ones which mimic the atom they relate
    # to and should lie on top of them
    imageatom = molecom.imageatom
    if ja == imageatom[ia]:
        return en_atom_intra, sumevib, sumeang
    itype = molecom.itype
    qtype = molecom.qtype
    nbond = molecom.nbond
    iabond = molecom.iabond
    tbond = molecom.tbond
    xstar = epsr16.xstar
    npt2 = epsr16.npt2
    sume = f32(0)
    # Calculate the internal energy for this atom pair
    stepcl = epsr16.cellst / molecom.dmaxmol
    # Set the atom type numbers
    it = itype[ia]
    jt = itype[ja]
    icref = molecom.npairref[it, jt]
    notqtype = qtype[ia] == 0 and qtype[ja] == 0
    # calculate separation of atoms within this molecule
    xdifs = x2 - x1
    ydifs = y2 - y1
    zdifs = z2 - z1
    xdif = xdifs
    ydif = ydifs
    zdif = zdifs
    # Apply minimum image convention
    xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
    # first sum intramolecular harmonic potentials if needed
    jref = 0
    iflag = 0
    image = ia
    while jref < nbond[image] and iflag == 0:
        jref = jref + 1
        if ja == iabond[image, jref]:
            iflag = jref
    # If a bond was not found, then check the image atom for bonds
    if iflag == 0:
        image = 0
        if imageatom[ia] > 0:
            image = int(imageatom[ia])
            jref = 0
            while jref < nbond[image] and iflag == 0:
                jref = jref + 1
                if ja == iabond[image, jref]:
                    iflag = jref
    if iflag > 0 and image > 0:
        dsq = xdif * xdif + ydif * ydif + zdif * zdif
        dtest = sqrtf(dsq)
        # find the r value this corresponds to
        rat = dtest / stepcl
        iz = ifix(rat) + 2
        wtms = molecom.amassr[icref]
        tb = max(stepcl, abs(tbond[image, jref]))
        wtm = wtms / tb
        if tbond[image, jref] < 0.0:
            wtm = wtm * molecom.tempang  # Emphasize bond angles if requested.
        ddif = tb - dtest
        # sume=en_morse(ddif,sqrt(wtm),tempang)*tempb*simt
        sume = wtm * ddif * ddif * molecom.tempb * molecom.simt
        if tbond[image, jref] > 0.0:
            sumevib = sumevib + sume
        if tbond[image, jref] < 0.0:
            sumeang = sumeang + sume  # Emphasize bond angles if requested.
        if molecom.refine_intra and iz <= npt2 and notqtype:
            dtdif = dtest - xstar[iz]
            sume = sume + molecom.refine_intra_fac * (potential_par.xavs[iz, icref]
                                                      + potential_par.xavs1[iz, icref] * dtdif)
    else:
        # Check whether ja is listed as forming a dihedral angle with ia
        im = molecom.mref[ia]
        found = False
        # Step through the rotational groups for this molecule and pick out the dihedrals -
        # Any atom pairs that are part of a dihedral are not used to calculate interatomic energies from
        # the reference potential
        ndihedralang = molecom.ndihedralang
        labrotdih = molecom.labrotdih
        ir = clustercom.nrf[im] - 1
        nrl_im = clustercom.nrl[im]
        while ir < nrl_im and not found:
            ir = ir + 1
            # Step through all the dihedral angles and check whether these two atoms form a dihedral. If so
            # their energy is calculated separately - not here.
            kd = 0
            while kd < ndihedralang[ir] and not found:
                kd = kd + 1
                io1 = labrotdih[1, kd, ir]
                io2 = labrotdih[2, kd, ir]
                found = (io1 == ia and io2 == ja) or (io1 == ja and io2 == ia)
        # if no bonds were found for this atom pair, and also no dihedrals are defined for this pair apply
        # molecular and empirical potential in the usual way.
        # Q-type atoms are not included in this calculation since their positions are
        # defined by reference to other atoms in the molecule.
        if not found:
            # apply minimum image convention if atoms are not connected by bonds
            sizefactor = molecom.sizefactor
            xdif = xdifs * sizefactor
            ydif = ydifs * sizefactor
            zdif = zdifs * sizefactor
            xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
            dsq = xdif * xdif + ydif * ydif + zdif * zdif
            dtest = sqrtf(dsq)
            # find the r value this corresponds to
            rat = dtest / stepcl
            iz = ifix(rat) + 2
            if iz <= npt2:
                dtdif = dtest - xstar[iz]
                sume = potential_par.xav[iz, icref] + potential_par.xav1[iz, icref] * dtdif
                if molecom.refine_intra and notqtype:
                    sume = sume + molecom.refine_intra_fac * (potential_par.xavs[iz, icref]
                                                              + potential_par.xavs1[iz, icref] * dtdif)
        else:
            sume = f32(0.0)
    en_atom_intra = sume
    return en_atom_intra, sumevib, sumeang


@routine
def en_morse(ddif, wtm, wtms, distort):
    """molfit.f:en_morse -- distort controls the amount of assymmetry in the bond. Larger
    means less symmetric, smaller more symmetric."""
    en = f32(0.0)
    arg = wtm * ddif * distort
    # write(6,*) ddif,wtm,arg
    if distort > 0.0:
        en = expf(arg) - 1.0
        en = en * en
    en_morse = wtms * en
    return en_morse


@routine
def en_morse_mod(ddif, wtm, wtms, distort):
    """molfit.f:en_morse_mod -- distort controls the amount of assymmetry in the bond."""
    en = f32(0.0)
    arg = wtm * ddif * distort
    if distort > 0.0:
        en = (expf(arg) - 1.0) / distort
    else:
        en = wtm * ddif
    en = en * en
    arg = 0.5 * arg * arg
    if arg < 30.0:
        en = en - expf(-arg)
    en_morse_mod = wtms * en
    return en_morse_mod


@routine
def en_morse_mod1(ddif, wtm, distort):
    """molfit.f:en_morse_mod1 -- distort controls the amount of assymmetry in the bond."""
    lim = _LIM
    limit = _LIMIT
    en = f32(0.0)
    if distort > 0.0:
        arg = wtm * ddif * distort
        if arg > 0.0:
            if arg < lim:
                en = expf(arg) - 1.0
            else:
                en = (arg / lim) * limit - 1.0
        else:
            en = abs(arg)
        # if(arg.gt.-1.0) then !Linear region for arg < -1.0
        #    en=en*en
        # end if
        if en < limit or arg > 0.0:
            en = en * en
        else:
            en = limit * en
    en_morse_mod1 = en
    return en_morse_mod1


@routine
def en_morse_mod2(ddif, wtm, wtms, distort):
    """molfit.f:en_morse_mod2 -- distort controls the amount of assymmetry in the bond."""
    lim = _LIM
    limit = _LIMIT
    en = f32(0.0)
    args = wtm * ddif
    if distort > 0.0:
        arg = args * distort
        if arg > 0.0:
            if arg < lim:
                en = expf(arg) - 1.0
            else:
                en = (arg / lim) * limit - 1.0
        else:
            en = abs(arg)
        # if(arg.gt.-1.0) then !Linear region for arg < -1.0
        #    en=en*en
        # end if
        if en < limit or arg > 0.0:
            en = en * en
        else:
            en = limit * en
    # args=0.5*args*args
    # if(args.lt.30.0) en=en-wtms*exp(-args)
    en_morse_mod2 = wtms * en
    return en_morse_mod2


@routine
def en_gauss(ddif, wtm, wtms, distort, gausswt, gausswidth):
    """molfit.f:en_gauss -- distort controls the amount of assymmetry in the bond."""
    en = f32(0.0)
    arg = wtm * ddif
    args = arg / distort
    # if(arg.lt.0.0) then
    #    args=arg*gausswidth
    # end if
    if args > 30.0:
        en = (args / 30.0) * _EXP30 - 1.0
    elif args < -30.0:
        en = r4(-1.0)
    else:
        en = expf(args) - 1.0
    en = wtms * (en * en - 1.0)
    if gausswt != 0.0 and gausswidth > 0.0:
        args = arg / gausswidth
        # if(args.lt.0.0) args=arg*gausswidth
        args = 0.5 * args * args
        if args < 30.0:
            en = en - gausswt * expf(-args)
        # en=en+gausswt*args
    en_gauss = en
    return en_gauss


@routine
def en_harmonic(ddif, wtm, wtms, distort, gausswt, gausswidth):
    """molfit.f:en_harmonic -- distort controls the amount of assymmetry in the bond."""
    en = f32(0.0)
    arg = wtm * ddif
    args = arg / distort
    # if(arg.lt.0.0) then
    #    args=arg*gausswidth
    # end if
    if args > 30.0:
        en = (args / 30.0) * _EXP30 - 1.0
    elif args < -30.0:
        en = r4(-1.0)
    else:
        en = expf(args) - 1.0
    en = wtms * (en * en - 1.0)
    if gausswt != 0.0 and gausswidth > 0.0:
        args = arg / gausswidth
        # if(args.lt.0.0) args=arg*gausswidth
        # args=0.5*args*args
        # if(args.lt.30.0) then
        #    en=en-gausswt*exp(-args)/gausswidth
        # end if
        if abs(args) < 1.0:
            en = en + gausswt * args * args
        else:
            en = en + gausswt * abs(args)
    en_harmonic = en
    return en_harmonic


@routine
def en_atom_inter(ia, ja, x1, y1, z1, x2, y2, z2):
    """molfit.f:en_atom_inter -- calculates the intermolecular energy between atoms ia and
    ja on distinct molecules im and jm (preselected to occur on distinct molecules)."""
    sume = f32(0.0)
    stepcl = epsr16.cellst / clustercom.dmax
    # write(6,*) ia,ja,x1,y1,z1,x2,y2,z2,stepcl
    # Set the atom type numbers
    it = molecom.itype[ia]
    jt = molecom.itype[ja]
    icref = molecom.npairref[it, jt]
    # write(6,*) ia,ja,it,jt,icref
    xdif = x2 - x1
    ydif = y2 - y1
    zdif = z2 - z1
    # apply minimum image convention
    xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
    dsq = xdif * xdif + ydif * ydif + zdif * zdif
    dtest = sqrtf(dsq)
    # write(6,*) ia,ja,dtest
    # find the r value this corresponds to
    rat = dtest / stepcl
    iz = ifix(rat) + 2
    # write(6,*) ia,ja,it,jt,iz,npt2
    if iz <= epsr16.npt2:
        # sum the reference and empirical potentials
        dtdif = dtest - epsr16.xstar[iz]
        sume = potential_par.xavtot[iz, icref] + potential_par.xavtot1[iz, icref] * dtdif
    else:
        sume = f32(0.0)
    en_atom_inter = sume
    return en_atom_inter


@routine
def en_atom_inter_test(x1, y1, z1, x2, y2, z2, stepcl, xstar, icref, npt2):
    """molfit.f:en_atom_inter_test -- intermolecular energy between two atoms on distinct
    molecules from the total (reference + empirical) potential."""
    sume = f32(0.0)
    xdif = x2 - x1
    ydif = y2 - y1
    zdif = z2 - z1
    # apply minimum image convention
    xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
    dsq = xdif * xdif + ydif * ydif + zdif * zdif
    dtest = sqrtf(dsq)
    # find the r value this corresponds to
    rat = dtest / stepcl
    iz = ifix(rat) + 2
    if iz <= npt2:
        dtdif = dtest - xstar[iz]
        # sum the reference and empirical potentials
        sume = potential_par.xavtot[iz, icref] + potential_par.xavtot1[iz, icref] * dtdif
    else:
        sume = f32(0.0)
    en_atom_inter_test = sume
    return en_atom_inter_test


@routine
def en_atom_inter_scale(x1, y1, z1, x2, y2, z2, stepcl, xstar, icref, npt2, scalefac):
    """molfit.f:en_atom_inter_scale -- intermolecular energy between two atoms on distinct
    molecules: empirical potential at the distance, reference potential at the distance
    scaled by scalefac."""
    sume = f32(0.0)
    xdif = x2 - x1
    ydif = y2 - y1
    zdif = z2 - z1
    # apply minimum image convention
    xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
    dsq = xdif * xdif + ydif * ydif + zdif * zdif
    dtest = sqrtf(dsq)
    # find the r value this corresponds to
    rat = dtest / stepcl
    iz = ifix(rat) + 2
    if iz <= npt2:
        dtdif = dtest - xstar[iz]
        # sum the empirical potentials
        sume = potential_par.xavs[iz, icref] + potential_par.xavs1[iz, icref] * dtdif
    # The reference potential uses the scaled distance
    dtest = dtest * scalefac
    rat = dtest / stepcl
    iz = ifix(rat) + 2
    if iz <= npt2:
        dtdif = dtest - xstar[iz]
        # sum the reference potentials
        sume = sume + potential_par.ecoul[iz, icref] + potential_par.ecoul1[iz, icref] * dtdif
    en_atom_inter_scale = sume
    return en_atom_inter_scale


@routine
def en_atom_inter_scale_coord(dtest, tb, wtm, wtms, distort, gausswt, gausswidth):
    """molfit.f:en_atom_inter_scale_coord -- coordination-constraint pair energy (Gaussian
    well for tb < 0, harmonic otherwise)."""
    ddif = abs(tb) - dtest
    if tb < 0.0:
        en_atom_inter_scale_coord = en_gauss(ddif, sqrtf(wtm), wtms, distort, gausswt, gausswidth)
    else:
        en_atom_inter_scale_coord = en_harmonic(ddif, sqrtf(wtm), wtms, distort, gausswt, gausswidth)
    # en_atom_inter_scale_coord=wtm*ddif*ddif
    return en_atom_inter_scale_coord


@routine
def en_atom_inter_scale_coordextra(x1, y1, z1, x2, y2, z2, stepcl, xstar, icref, npt2, scalefac):
    """molfit.f:en_atom_inter_scale_coordextra -- as en_atom_inter_scale for an atom pair
    forming a coordination constraint (reference potential epot instead of ecoul)."""
    sume = f32(0.0)
    xdif = x2 - x1
    ydif = y2 - y1
    zdif = z2 - z1
    # apply minimum image convention
    xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
    dsq = xdif * xdif + ydif * ydif + zdif * zdif
    dtest = sqrtf(dsq)
    # find the r value this corresponds to
    rat = dtest / stepcl
    iz = ifix(rat) + 2
    if iz <= npt2:
        dtdif = dtest - xstar[iz]
        # sum the empirical potentials
        sume = sume + potential_par.xavs[iz, icref] + potential_par.xavs1[iz, icref] * dtdif
    # The reference potential uses the scaled distance
    dtest = dtest * scalefac
    rat = dtest / stepcl
    iz = ifix(rat) + 2
    if iz <= npt2:
        dtdif = dtest - xstar[iz]
        # sum the reference potentials
        sume = sume + potential_par.epot[iz, icref] + potential_par.epot1[iz, icref] * dtdif
    en_atom_inter_scale_coordextra = sume
    return en_atom_inter_scale_coordextra


@routine
def en_mole_intra(im, sumevib, sumeang, sumedih):
    """molfit.f:en_mole_intra -- calculates the intra-molecular energy of a molecule using
    the current coordinates.

    Returns (en_mole_intra, sumevib, sumeang, sumedih): all three sums are updated, the
    first two through en_atom_intra (docs/SIGNATURE_EXCEPTIONS.md).
    """
    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za
    naf_im = clustercom.naf[im]
    nal_im = clustercom.nal[im]
    en_mole_intra = f32(0.0)
    for ia in range(naf_im, nal_im):
        for ja in range(ia + 1, nal_im + 1):
            sume, sumevib, sumeang = en_atom_intra(ia, ja, xa[ia], ya[ia], za[ia], xa[ja], ya[ja], za[ja],
                                                   sumevib, sumeang)
            en_mole_intra = en_mole_intra + sume
    if clustercom.ngrot[im] > 0 and molecom.tempdih > 0.0:
        sume = molecom.tempdih * molecom.simt * en_mole_dihedral(im)
        en_mole_intra = en_mole_intra + sume
        sumedih = sumedih + sume
    return en_mole_intra, sumevib, sumeang, sumedih


@routine
def en_mole_dihedral(im):
    """molfit.f:en_mole_dihedral -- dihedral energy of the DIH rotational groups of
    molecule im, using the current coordinates."""
    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za
    labaxi = molecom_mod.labaxi
    rotgrp = molecom_mod.rotgrp
    labrotdih = molecom.labrotdih
    ndihedralang = molecom.ndihedralang
    dihedralang = molecom.dihedralang
    cd = real_array(3)
    co1 = real_array(3)
    co2 = real_array(3)

    en_mole_dihedral = f32(0.0)
    # Step through the rotational groups for this molecule and pick out the dihedrals
    for ir in range(clustercom.nrf[im], clustercom.nrl[im] + 1):
        if rotgrp[ir] == 'DIH':
            #
            # define the separation vector cd between the two base atoms and normalize it to
            # unit value
            #
            ia1 = labaxi[ir, 1]
            ia2 = labaxi[ir, 2]
            x1 = xa[ia1]
            y1 = ya[ia1]
            z1 = za[ia1]
            x2 = xa[ia2]
            y2 = ya[ia2]
            z2 = za[ia2]
            cd[1] = (x2 - x1)
            cd[2] = (y2 - y1)
            cd[3] = (z2 - z1)
            dmag = sqrtf(auxiliary_routines.scapro(cd, cd))
            cd[1] = cd[1] / dmag
            cd[2] = cd[2] / dmag
            cd[3] = cd[3] / dmag
            # write(6,*) ir,ia1,ia2,ndihedralang(ir)
            # Step through all the dihedral angles and calculate their contribution to the energy
            kd = 0
            while kd < ndihedralang[ir]:
                kd = kd + 1
                #
                # now determine the orientation vectors co for base atom 1 and base atom 2, using
                # the values of ltyp2 and ltyp3
                #
                io1 = labrotdih[1, kd, ir]
                io2 = labrotdih[2, kd, ir]
                # write(6,*) kd,io1,io2
                co1[1] = xa[io1] - x1
                co1[2] = ya[io1] - y1
                co1[3] = za[io1] - z1
                co2[1] = xa[io2] - x2
                co2[2] = ya[io2] - y2
                co2[3] = za[io2] - z2
                en_mole_dihedral = en_mole_dihedral + en_dihedral(co1, co2, cd, dihedralang[kd, ir])
    return en_mole_dihedral


@routine
def en_change_dihedral(im, ig, xr, yr, zr):
    """molfit.f:en_change_dihedral -- change of the dihedral energy of group ig of molecule
    im between the saved coordinates (xr, yr, zr) and the current ones."""
    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za
    labaxi = molecom_mod.labaxi
    labrotdih = molecom.labrotdih
    ndihedralang = molecom.ndihedralang
    dihedralang = molecom.dihedralang
    cd = real_array(3)
    co1 = real_array(3)
    co2 = real_array(3)
    # cv1(3),cv2(3),cpar(3) are declared but not used

    en_change_dihedral = f32(0.0)
    # Step through the rotational groups for this molecule and pick out the dihedrals
    ir = ig
    iaadd = 1 - clustercom.naf[im]
    if molecom_mod.rotgrp[ir] == 'DIH':
        #
        # define the separation vector cd between the two base atoms and normalize it to
        # unit value
        #
        ia1 = labaxi[ir, 1]
        ia2 = labaxi[ir, 2]
        x1 = xa[ia1]
        y1 = ya[ia1]
        z1 = za[ia1]
        x2 = xa[ia2]
        y2 = ya[ia2]
        z2 = za[ia2]
        cd[1] = (x2 - x1)
        cd[2] = (y2 - y1)
        cd[3] = (z2 - z1)
        dmag = sqrtf(auxiliary_routines.scapro(cd, cd))
        cd[1] = cd[1] / dmag
        cd[2] = cd[2] / dmag
        cd[3] = cd[3] / dmag
        # write(6,*) ir,ia1,ia2,ndihedralang(ir)
        # Step through all the dihedral angles and calculate their contribution to the energy
        kd = 0
        while kd < ndihedralang[ir]:
            kd = kd + 1
            #
            # now determine the orientation vectors co for base atom 1 and base atom 2, using
            # the values of ltyp2 and ltyp3
            #
            # Energy of current dihedral
            io1 = labrotdih[1, kd, ir]
            io2 = labrotdih[2, kd, ir]
            co1[1] = xa[io1] - x1
            co1[2] = ya[io1] - y1
            co1[3] = za[io1] - z1
            co2[1] = xa[io2] - x2
            co2[2] = ya[io2] - y2
            co2[3] = za[io2] - z2
            en_change_dihedral = en_change_dihedral + en_dihedral(co1, co2, cd, dihedralang[kd, ir])
            # Energy of old dihedral
            ir1 = io1 + iaadd
            ir2 = io2 + iaadd
            co1[1] = xr[ir1] - x1
            co1[2] = yr[ir1] - y1
            co1[3] = zr[ir1] - z1
            co2[1] = xr[ir2] - x2
            co2[2] = yr[ir2] - y2
            co2[3] = zr[ir2] - z2
            en_change_dihedral = en_change_dihedral - en_dihedral(co1, co2, cd, dihedralang[kd, ir])

    return en_change_dihedral


@routine
def en_dihedral(co1, co2, cd, dihedralang):
    """molfit.f:en_dihedral -- dihedral energy for two vectors, co1 and co2, with respect to
    the axis vector cd and specified angle dihedralang:
    (cos phi - cos phi0)**2 + (sin phi - sin phi0)**2."""
    cv1 = real_array(3)
    cv2 = real_array(3)
    cpar = real_array(3)

    # define new vectors normal to both cd and co
    auxiliary_routines.vecpro(co1, cd, cv1)
    auxiliary_routines.vecpro(co2, cd, cv2)
    # renormalize these to unit vectors
    dmag = sqrtf(auxiliary_routines.scapro(cv1, cv1))
    cv1[1] = cv1[1] / dmag
    cv1[2] = cv1[2] / dmag
    cv1[3] = cv1[3] / dmag
    dmag = sqrtf(auxiliary_routines.scapro(cv2, cv2))
    cv2[1] = cv2[1] / dmag
    cv2[2] = cv2[2] / dmag
    cv2[3] = cv2[3] / dmag
    # calculate the cosine of the angle between cv1 and cv2
    cosphi = auxiliary_routines.scapro(cv1, cv2)
    if cosphi > 1.0:
        cosphi = r4(1.0)
    if cosphi < -1.0:
        cosphi = r4(-1.0)
    # calculate the vector product of cv1 and cv2 - this should be either parallel
    # (phi positive) or antiparallel to cd (phi negative)
    auxiliary_routines.vecpro(cv1, cv2, cpar)
    sinphi = auxiliary_routines.scapro(cpar, cd)
    cosdih = cosf(dihedralang)
    sindih = sinf(dihedralang)
    difc = cosdih - cosphi
    difs = sindih - sinphi
    diff = difc * difc + difs * difs
    # write(6,*) ir,kd,cosdih,cosphi,sindih,sinphi
    en_dihedral = diff
    return en_dihedral


@routine
def en_mole_pair(im, jm):
    """molfit.f:en_mole_pair -- calculates the inter-molecular energy of a pair of molecules
    using the current coordinates."""
    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za
    itype = molecom.itype
    qtype = molecom.qtype
    npairref = molecom.npairref
    mref = molecom.mref
    xadd = clustercom.xadd
    yadd = clustercom.yadd
    zadd = clustercom.zadd
    xstar = epsr16.xstar
    npt2 = epsr16.npt2
    ncoordconstraints = cc_d.ncoordconstraints
    check_coord_bonded = coordconstraint_mod.check_coord_bonded

    stepcl = epsr16.cellst / clustercom.dmax
    # Don't calculate pair energies if both molecules are tethered and both atoms are Q-type
    bothtethered = clustercom.tethered[im] and clustercom.tethered[jm]
    scalefac = f32(1.0)
    if molecom_mod.coarsegraining:
        scalefac = 0.5 * (molecom_mod.scalefacmol[im] + molecom_mod.scalefacmol[jm])
    sume = f32(0.0)
    naf_jm = clustercom.naf[jm]
    nal_jm = clustercom.nal[jm]
    for ia in range(clustercom.naf[im], clustercom.nal[im] + 1):
        it = itype[ia]
        x1 = xa[ia] + xadd[im]
        y1 = ya[ia] + yadd[im]
        z1 = za[ia] + zadd[im]
        qtypeia = qtype[ia] == 1 and bothtethered
        for ja in range(naf_jm, nal_jm + 1):
            bothqatoms = qtype[ja] == 1 and qtypeia
            if not bothqatoms:
                jt = itype[ja]
                icref = npairref[it, jt]
                x2 = xa[ja] + xadd[jm]
                y2 = ya[ja] + yadd[jm]
                z2 = za[ja] + zadd[jm]
                # (check_coord_bonded is .false. without side effects when
                # ncoordconstraints is 0; the test is skipped then, for speed)
                if ncoordconstraints > 0 and check_coord_bonded(ia, ja, mref, itype):
                    sume = sume + en_atom_inter_scale_coordextra(x1, y1, z1, x2, y2, z2,
                                                                 stepcl, xstar, icref, npt2, scalefac)
                else:
                    sume = sume + en_atom_inter_scale(x1, y1, z1, x2, y2, z2,
                                                      stepcl, xstar, icref, npt2, scalefac)
    en_mole_pair = sume
    return en_mole_pair


@routine
def test_en_change_tethered(im, ia, temp):
    """molfit.f:test_en_change_tethered -- Metropolis test on the change in tethering
    energy."""
    ediff = en_change_tethered(im, ia) * temp
    etest = r4(10.0)
    if ediff <= 0.0:
        etest = r4(10.0)
    elif ediff < 30.0:
        etest = expf(-ediff)
    else:
        etest = r4(0.0)
    test_en_change_tethered = bool(auxiliary_routines.ran1() < etest)
    return test_en_change_tethered


@routine
def en_change_tethered(im, ia):
    """molfit.f:en_change_tethered -- with tethered molecules, for whole molecule moves or
    moves where part of a molecule is rotated we need to calculate the change in energy that
    derives from the change in position of centre of mass (relative to xs, ys, zs)."""
    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za
    nato_im = clustercom.nato[im]

    # Get the old energy of this molecule relative to the tether position

    xdif = molecom.xs - clustercom.phix[im]
    ydif = molecom.ys - clustercom.phiy[im]
    zdif = molecom.zs - clustercom.phiz[im]
    if ia > 0 and ia <= nato_im:
        iaref = clustercom.naf[im] + ia - 1
        xdif = xdif + xa[iaref]
        ydif = ydif + ya[iaref]
        zdif = zdif + za[iaref]
    xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
    term = xdif * xdif + ydif * ydif + zdif * zdif

    # Get the new position of this molecule relative to the tether position

    xdif = clustercom.xadd[im] - clustercom.phix[im]
    ydif = clustercom.yadd[im] - clustercom.phiy[im]
    zdif = clustercom.zadd[im] - clustercom.phiz[im]
    if ia > 0 and ia <= nato_im:
        iaref = clustercom.naf[im] + ia - 1
        xdif = xdif + xa[iaref]
        ydif = ydif + ya[iaref]
        zdif = zdif + za[iaref]

    # Calculate the change in energy

    xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
    en_change_tethered = (xdif * xdif + ydif * ydif + zdif * zdif - term)
    return en_change_tethered


@routine
def check_tethered(im, ia):
    """molfit.f:check_tethered -- with tethered molecules, check that the molecule has not
    moved more than the allowed distance (tol) from the tethering point; if it has, accept
    only if the tethering energy decreases.  Always .true. for untethered molecules."""
    # Get the position of this molecule relative to the tether position

    check_tethered = not clustercom.tethered[im]
    if check_tethered:
        return True
    xa = molecom.xa
    ya = molecom.ya
    za = molecom.za
    iaref = clustercom.naf[im] + ia - 1
    xdif = clustercom.xadd[im] - clustercom.phix[im]
    ydif = clustercom.yadd[im] - clustercom.phiy[im]
    zdif = clustercom.zadd[im] - clustercom.phiz[im]
    if ia > 0 and ia <= clustercom.nato[im]:
        iaref = clustercom.naf[im] + ia - 1
        xdif = xdif + xa[iaref]
        ydif = ydif + ya[iaref]
        zdif = zdif + za[iaref]
    xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
    dist = sqrtf(xdif * xdif + ydif * ydif + zdif * zdif)
    check_tethered = bool(dist <= molecom.tol)
    # If outsice the allowed range, only accept the move if the change in tethering energy is negative -
    # particle gets closer to tethering position
    if not check_tethered:
        check_tethered = bool(en_change_tethered(im, ia) <= 0.0)
        # (commented-out code in the Fortran: pull the molecule back to distance tol)
    return check_tethered
