"""Translation of fortran_src/threebody_mod.f90: the routines of module threebody_mod.

Three-body (triplet) potentials.  A ``tb`` line of the ``.pcof`` file defines a group of
triplets A1-A2-A3 (A2 = centre); ``find_triplets`` lists them per centre molecule and
``en_mole_triplets`` gives the energy of the triplets centred on a molecule.

The module data is in ``common/threebody_mod.py``.

FORTRAN_BUGS.md 1.1 (angles corrupted to radians after every save) is fixed entirely in
``read_and_save_files`` (package B, the Rust ``write_pcof_threebody``); no part of that
fix lives in the module routines translated here.
"""

import math

from ..common import threebody_mod as td  # Fortran module data
from ..fortran import PI4ATAN, Ch, acosf, expf, f32, fmt, lw, r4, real_array, sinf, sqrtf, units
from ..trace import routine
from . import minimagec, reallocation_routines

# gfortran folds exp() of a constant at compile time (correctly rounded); equal to expf here.
_EXP30 = f32(math.exp(30.0))

_F100 = "(/'find_triplets> atoms not matched to those in .ato file for triplet',i3)"


def _size(array):
    """size(array) of a 1-D allocatable (0 if not allocated: undefined in the Fortran)."""
    return 0 if array is None else array.shape[0] - 1


@routine
def find_triplets(nmol, natom, xadd, yadd, zadd, xa, ya, za, mref, itype, ntype, type):
    """threebody_mod.f90:find_triplets -- list the triplets of every three-body group,
    per centre molecule, and the triplet neighbours of each molecule.
    """
    # Initialise certain arrays
    # NB: ergold_triplet/ergnew_triplet are freshly allocated (uninitialised) in the
    # Fortran; zero here.
    if td.ergold_triplet is not None:
        td.ergold_triplet = None
        td.ergnew_triplet = None
    td.ergold_triplet = real_array(nmol)
    td.ergnew_triplet = real_array(nmol)

    if td.nthreeb <= 0:
        return

    td.ntripletsinmol = reallocation_routines.reallocate2d_i(td.ntripletsinmol, nmol, td.nthreeb)
    td.firsttripletinmol = reallocation_routines.reallocate2d_i(td.firsttripletinmol, nmol, td.nthreeb)
    td.sume_triplet_dist = reallocation_routines.reallocate1d_r(td.sume_triplet_dist, td.nthreeb)
    td.en_mole_triplet_dist = reallocation_routines.reallocate1d_r(td.en_mole_triplet_dist, td.nthreeb)
    ntripletsmax = 10000
    td.tripletlist = reallocation_routines.reallocate3d_i(td.tripletlist, 3, ntripletsmax, td.nthreeb)
    # Arrays which will list the triplet neighbours of each molecule
    td.nthreebneighbmolmax = 4
    td.nthreebneighbmol = reallocation_routines.reallocate1d_i(td.nthreebneighbmol, nmol)
    td.threebneighbmol = reallocation_routines.reallocate2d_i(td.threebneighbmol, td.nthreebneighbmolmax, nmol)
    # Set the initial values
    nthreebneighbmol = td.nthreebneighbmol
    for im in range(1, nmol + 1):
        nthreebneighbmol[im] = 0

    ntripletsinmol = td.ntripletsinmol
    firsttripletinmol = td.firsttripletinmol
    threebgrp = td.threebgrp
    threebsamemol = td.threebsamemol
    threebrmax = td.threebrmax

    # First identify the type numbers of the three specified atoms

    for id in range(1, td.nthreeb + 1):
        for im1 in range(1, nmol + 1):
            ntripletsinmol[im1, id] = 0
            firsttripletinmol[im1, id] = 0
        ntriplets = 0
        ifound1 = 0
        ifound2 = 0
        ifound3 = 0
        it = 0
        while it < ntype and (ifound1 == 0 or ifound2 == 0 or ifound3 == 0):
            it = it + 1
            units.write(6, fmt('(1x,i3,1x,a,3(1x,a),9(1x,e12.5))', it, Ch(type[it], 3),
                               Ch(threebgrp[1, id], 3), Ch(threebgrp[2, id], 3), Ch(threebgrp[3, id], 3),
                               *[td.threebangs[j, id] for j in range(1, td.nthreebangs[id] + 1)]))
            if type[it] == threebgrp[1, id]:
                ifound1 = it
            if type[it] == threebgrp[2, id]:
                ifound2 = it
            if type[it] == threebgrp[3, id]:
                ifound3 = it
        if ifound1 == 0 or ifound2 == 0 or ifound3 == 0:
            units.write(6, fmt(_F100, id))
            return
        td.ltyp1 = ifound1
        td.ltyp2 = ifound2
        td.ltyp3 = ifound3
        ltyp1 = td.ltyp1
        ltyp2 = td.ltyp2
        ltyp3 = td.ltyp3

        # Make a list of all the atoms of type ltyp1, ltype2, and ltyp3

        td.natom1list = 0
        td.natom2list = 0
        td.natom3list = 0
        for ia2 in range(1, natom + 1):
            if itype[ia2] == ltyp1:
                td.natom1list = td.natom1list + 1
            if itype[ia2] == ltyp2:
                td.natom2list = td.natom2list + 1
            if itype[ia2] == ltyp3:
                td.natom3list = td.natom3list + 1
        units.write(6, fmt('(4(i2,1x),3(i5,1x))', id, ltyp1, ltyp2, ltyp3, td.natom1list, td.natom2list,
                           td.natom3list))
        td.atom1list = reallocation_routines.reallocate1d_i(td.atom1list, td.natom1list)
        units.write(6, fmt('(2(i2,1x),2(i5,1x))', id, ltyp1, td.natom1list, _size(td.atom1list)))
        td.atom2list = reallocation_routines.reallocate1d_i(td.atom2list, td.natom2list)
        units.write(6, fmt('(2(i2,1x),2(i5,1x))', id, ltyp1, td.natom2list, _size(td.atom2list)))
        td.atom3list = reallocation_routines.reallocate1d_i(td.atom3list, td.natom3list)
        units.write(6, fmt('(2(i2,1x),2(i5,1x))', id, ltyp1, td.natom3list, _size(td.atom3list)))
        atom1list = td.atom1list
        atom2list = td.atom2list
        atom3list = td.atom3list
        natom1list = td.natom1list
        natom2list = td.natom2list
        natom3list = td.natom3list
        # Repeat this saving the atom numbers
        iref1 = 0
        iref2 = 0
        iref3 = 0
        for ia2 in range(1, natom + 1):
            if itype[ia2] == ltyp1:
                iref1 = iref1 + 1
                atom1list[iref1] = ia2
            if itype[ia2] == ltyp2:
                iref2 = iref2 + 1
                atom2list[iref2] = ia2
            if itype[ia2] == ltyp3:
                iref3 = iref3 + 1
                atom3list[iref3] = ia2
        units.write(6, lw(id, natom1list, natom2list, natom3list))

        s1 = threebsamemol[1, id]
        s2 = threebsamemol[2, id]
        s3 = threebsamemol[3, id]
        # For each type 2 atom we need to determine which type 1 and type 3 atoms are within the specified distances,
        # on the understanding that the type 1 and type 3 atoms cannot be the same atom as type 2, nor can
        # the atom of type 3 be the same as the atom of type 1
        for iref2 in range(1, natom2list + 1):
            ia2 = atom2list[iref2]
            im2 = mref[ia2]
            # write(6,*) 'Molecule: ',im2
            # Get the coordinates of this atom
            x2 = xa[ia2] + xadd[im2]
            y2 = ya[ia2] + yadd[im2]
            z2 = za[ia2] + zadd[im2]
            iref1 = 0
            # Search for all possible pairs of type 1 and type 3 atoms that satisfy the distance constraints
            while iref1 < natom1list:
                # Search type 1 atoms
                foundia1 = False
                while iref1 < natom1list and not foundia1:
                    iref1 = iref1 + 1
                    ia1 = atom1list[iref1]
                    im1 = mref[ia1]
                    if ia1 != ia2 and ((s1 == 1 and im1 == im2) or (s1 == 0 and im1 != im2)):
                        # Get the coordinates of this atom
                        x1 = xa[ia1] + xadd[im1]
                        y1 = ya[ia1] + yadd[im1]
                        z1 = za[ia1] + zadd[im1]
                        xdif = x1 - x2
                        ydif = y1 - y2
                        zdif = z1 - z2
                        # apply periodic boundary conditions
                        xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                        dsq = xdif * xdif + ydif * ydif + zdif * zdif
                        d12 = sqrtf(dsq)
                        # add12=threebrmax(1,id)*threebrwidth(id)
                        foundia1 = bool(d12 < threebrmax[1, id])
                # No point in searching for type 3 atoms if no type 1 was found
                if foundia1:
                    # write(6,*) 'Molecule: ',im2,im1,ia1
                    iref3 = 0
                    foundia3 = False
                    for iref3 in range(1, natom3list + 1):
                        ia3 = atom3list[iref3]
                        im3 = mref[ia3]
                        # Need to ensure unique combinations of ia3 and ia1 are selected
                        if (((s2 == 1 and im2 == im3) or (s2 == 0 and im2 != im3))
                                and ((s3 == 1 and im1 == im3) or (s3 == 0 and im1 != im3))
                                and ((ltyp3 == ltyp1 and ia3 > ia1) or (ltyp3 != ltyp1 and ia3 != ia1))):
                            # Get the coordinates of this atom
                            x3 = xa[ia3] + xadd[im3]
                            y3 = ya[ia3] + yadd[im3]
                            z3 = za[ia3] + zadd[im3]
                            xdif = x3 - x2
                            ydif = y3 - y2
                            zdif = z3 - z2
                            # apply periodic boundary conditions
                            xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                            dsq = xdif * xdif + ydif * ydif + zdif * zdif
                            d32 = sqrtf(dsq)
                            # add32=threebrmax(2,id)*threebrwidth(id)
                            foundia3 = bool(d32 < threebrmax[2, id])
                            # .and.d32.gt.(threebrmax(2,id)-add32)
                            if foundia3:
                                ntriplets = ntriplets + 1
                                if ntriplets > ntripletsmax:
                                    ntripletsmax = ntriplets + 10000
                                    td.tripletlist = reallocation_routines.reallocate3d_i(td.tripletlist, 3, ntripletsmax,
                                                                                          td.nthreeb)
                                tripletlist = td.tripletlist
                                tripletlist[1, ntriplets, id] = ia1
                                tripletlist[2, ntriplets, id] = ia2
                                tripletlist[3, ntriplets, id] = ia3
                                ntripletsinmol[im2, id] = ntripletsinmol[im2, id] + 1
                                if ntripletsinmol[im2, id] == 1:
                                    firsttripletinmol[im2, id] = ntriplets
                                # write(6,*) im2,ia1,ia2,ia3
                                foundia3 = False
                                # Need to check that im1 and im3 are listed as triplet neighbours of im2,
                                # assuming they are not the same molecule
                                if im1 != im2:
                                    it = 0
                                    found = False
                                    while it < nthreebneighbmol[im2] and not found:
                                        it = it + 1
                                        found = im1 == td.threebneighbmol[it, im2]
                                    if not found:
                                        nthreebneighbmol[im2] = nthreebneighbmol[im2] + 1
                                        if nthreebneighbmol[im2] > td.nthreebneighbmolmax:
                                            td.nthreebneighbmolmax = td.nthreebneighbmolmax + 4
                                            td.threebneighbmol = reallocation_routines.reallocate2d_i(
                                                td.threebneighbmol, td.nthreebneighbmolmax, nmol)
                                        td.threebneighbmol[nthreebneighbmol[im2], im2] = im1
                                if im3 != im2:
                                    it = 0
                                    found = False
                                    while it < nthreebneighbmol[im2] and not found:
                                        it = it + 1
                                        found = im3 == td.threebneighbmol[it, im2]
                                    if not found:
                                        nthreebneighbmol[im2] = nthreebneighbmol[im2] + 1
                                        if nthreebneighbmol[im2] > td.nthreebneighbmolmax:
                                            td.nthreebneighbmolmax = td.nthreebneighbmolmax + 4
                                            td.threebneighbmol = reallocation_routines.reallocate2d_i(
                                                td.threebneighbmol, td.nthreebneighbmolmax, nmol)
                                        td.threebneighbmol[nthreebneighbmol[im2], im2] = im3
        # do im2=1,min(20,nmol)
        #    write(6,'(a,i3,1x,i7,1x,i9,1x,i9)') 'find_triplets> No. of triplets in mol: ' &
        #       ,id,im2,ntripletsinmol(im2,id),firsttripletinmol(im2,id)
        # end do
    get_triplet_pot_cons()
    # do im=1,nmol
    #    write(6,*) im,nthreebneighbmol(im),(threebneighbmol(it,im),it=1,nthreebneighbmol(im))
    # end do


@routine
def get_triplet_pot_cons():
    """threebody_mod.f90:get_triplet_pot_cons -- potential constant of each triplet
    distribution (printed only).
    """
    pi = PI4ATAN
    if td.nthreeb <= 0:
        return
    td.consthreeb = reallocation_routines.reallocate1d_r(td.consthreeb, td.nthreeb)
    for id in range(1, td.nthreeb + 1):
        sumangle = f32(0.0)
        sumlorentz = f32(0.0)
        nangle = 1001
        for ia in range(1, nangle + 1):
            angle = f32(ia - 1) * pi / f32(nangle - 1)
            sineangle = sinf(angle)
            sumangle = sumangle + sineangle
            for ib in range(1, td.nthreebangs[id] + 1):
                # NB: subtracts the stored cosines from an angle in radians, as the Fortran
                # (FORTRAN_BUGS.md section 3; print-only value)
                arg = (angle - td.threebangs[ib, id]) / td.threebawidth[id]
                arg = 0.5 * arg * arg
                if arg < 30.:
                    sumlorentz = sumlorentz + td.threebamps[ib, id] * expf(-arg)
                # sumlorentz=sumlorentz+threebamps(ib,id)*sineangle/(1.0+arg*arg)
        td.consthreeb[id] = sumlorentz / sumangle
    units.write(6, lw('get_triplet_pot_cons> Potential constants: ',
                      *[td.consthreeb[id] for id in range(1, td.nthreeb + 1)]))
    return


@routine
def en_mole_triplets(imol, xadd, yadd, zadd, phix, phiy, phiz, xa, ya, za, mref, tethered):
    """threebody_mod.f90:en_mole_triplets -- three-body energy associated with molecule imol
    (also sets en_mole_triplet_dist).  Returns en_mole_triplets.
    """
    # Range of triplets to be calculated

    en_mole_triplets = f32(0.0)
    if td.nthreeb <= 0:
        return en_mole_triplets
    threebamp = td.threebamp
    threebangs = td.threebangs
    threebamps = td.threebamps
    threebawidth = td.threebawidth
    threebrmax = td.threebrmax
    threebrwidth = td.threebrwidth
    nthreebangs = td.nthreebangs
    tripletlist = td.tripletlist
    ntripletsinmol = td.ntripletsinmol
    firsttripletinmol = td.firsttripletinmol
    en_mole_triplet_dist = td.en_mole_triplet_dist
    ensum = f32(0.0)
    for id in range(1, td.nthreeb + 1):
        en_mole_triplet_dist[id] = 0.0
        if abs(threebamp[id]) > 0.0:
            # OK now step through each set of pairs and find the angle value which is closest to the specified angle
            entrisum = f32(0.0)
            itf = firsttripletinmol[imol, id]
            itl = itf + ntripletsinmol[imol, id] - 1
            for it in range(itf, itl + 1):
                # Get coordinates of atoms in triplet
                ia1 = tripletlist[1, it, id]
                ia2 = tripletlist[2, it, id]
                ia3 = tripletlist[3, it, id]
                im1 = mref[ia1]
                im2 = mref[ia2]
                im3 = mref[ia3]
                x1 = xa[ia1] + xadd[im1]
                x2 = xa[ia2] + xadd[im2]
                x3 = xa[ia3] + xadd[im3]
                y1 = ya[ia1] + yadd[im1]
                y2 = ya[ia2] + yadd[im2]
                y3 = ya[ia3] + yadd[im3]
                z1 = za[ia1] + zadd[im1]
                z2 = za[ia2] + zadd[im2]
                z3 = za[ia3] + zadd[im3]
                distmax12 = threebrmax[1, id] * (1.0 + threebrwidth[id])  # noqa: F841
                distmax32 = threebrmax[2, id] * (1.0 + threebrwidth[id])  # noqa: F841
                # Get the distances 1 - 2, 2 - 3
                xdif = x1 - x2
                ydif = y1 - y2
                zdif = z1 - z2
                xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                d12sq = xdif * xdif + ydif * ydif + zdif * zdif
                d12 = sqrtf(d12sq)
                xdif = x3 - x2
                ydif = y3 - y2
                zdif = z3 - z2
                xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                d32sq = xdif * xdif + ydif * ydif + zdif * zdif
                d32 = sqrtf(d32sq)
                # Get distance 1 - 3
                xdif = x1 - x3
                ydif = y1 - y3
                zdif = z1 - z3
                xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                d13sq = xdif * xdif + ydif * ydif + zdif * zdif
                cost = (d12sq + d32sq - d13sq) / (2.0 * d12 * d32)
                # gfortran max/min: the first argument unless the second compares greater/less
                cost = max(r4(-1.0), min(r4(1.0), cost))
                # cosa=acos(cost)
                entrimin = f32(0.0)
                ib = 1
                ibmin = ib
                arg = cost - threebangs[ib, id]
                while ib < nthreebangs[id]:
                    ib = ib + 1
                    argtest = cost - threebangs[ib, id]
                    if abs(argtest) < abs(arg):
                        arg = argtest
                        ibmin = ib
                arg = arg / threebawidth[id]
                term = (0.5 * arg * arg)
                entrimin = entrimin + threebamps[ibmin, id] * term
                truncatefactor = f32(1.0)
                entrisum = entrisum + truncatefactor * entrimin
                en_mole_triplet_dist[id] = en_mole_triplet_dist[id] + truncatefactor * entrimin
                # if(im2.lt.21) then
                #    write(6,*) id,im1,im2,im3,ia1,ia2,ia3,d12,d32,cosa,threebangs(imin,id),entrisum
                # end if
            if ntripletsinmol[imol, id] > 0:
                en_mole_triplet_dist[id] = threebamp[id] * en_mole_triplet_dist[id]
                ensum = ensum + threebamp[id] * entrisum
    # if(imol.lt.21) write(6,*) imol,ensum,(id,en_mole_triplet_dist(id),id=1,nthreeb)
    en_mole_triplets = ensum
    return en_mole_triplets


@routine
def en_mole_triplets_old(imol, xadd, yadd, zadd, phix, phiy, phiz, xa, ya, za, mref, tethered):
    """threebody_mod.f90:en_mole_triplets_old -- older version of en_mole_triplets (not called).

    The Fortran declares a *local* ``en_mole_triplets`` and assigns that; the function
    result ``en_mole_triplets_old`` is never set.  Returns en_mole_triplets_old.
    """
    # NB: the result en_mole_triplets_old is never assigned in the Fortran (undefined);
    # 0.0 is returned.  The routine is not called anywhere.
    en_mole_triplets_old = f32(0.0)
    ia2ia3pairs = None  # integer, dimension(:,:), allocatable

    # Range of triplets to be calculated

    en_mole_triplets = f32(0.0)
    if td.nthreeb <= 0:
        return en_mole_triplets_old
    ensum = f32(0.0)
    for id in range(1, td.nthreeb + 1):
        td.en_mole_triplet_dist[id] = 0.0
        if td.ntripletsinmol[imol, id] > 0 and abs(td.threebamp[id]) > 0.0:
            ia2ia3pairs = reallocation_routines.reallocate2d_i(ia2ia3pairs, 2, td.ntripletsinmol[imol, id])
            itf = td.firsttripletinmol[imol, id]
            itl = itf + td.ntripletsinmol[imol, id] - 1
            # Step through listed triplets.
            # It will noted that if nato(imol)>1, there may be more than one triplet associated with the atom pair ia2,ia3.
            # We will need to identify the one of these pairs that is closest to the specified bond angle.
            nia2ia3pairs = 0
            for it in range(itf, itl + 1):
                # Get coordinates of atoms in triplet
                ia1 = td.tripletlist[1, it, id]
                ia2 = td.tripletlist[2, it, id]
                ia3 = td.tripletlist[3, it, id]
                found = False
                iia2ia3pairs = 0
                while iia2ia3pairs < nia2ia3pairs and not found:
                    iia2ia3pairs = iia2ia3pairs + 1
                    found = ia2 == ia2ia3pairs[1, iia2ia3pairs] and ia3 == ia2ia3pairs[2, iia2ia3pairs]
                if not found:
                    nia2ia3pairs = nia2ia3pairs + 1
                    ia2ia3pairs[1, nia2ia3pairs] = ia2
                    ia2ia3pairs[2, nia2ia3pairs] = ia3
            # OK now step through each set of pairs and find the angle value which is closest to the specified angle
            for ipairs in range(1, nia2ia3pairs + 1):
                icount = 0
                imin = 0  # noqa: F841
                enmin = f32(0)  # noqa: F841
                itmin = 0  # noqa: F841
                d12min = f32(0)  # noqa: F841
                d32min = f32(0)  # noqa: F841
                entrisum = f32(0)
                for it in range(itf, itl + 1):
                    # Get coordinates of atoms in triplet
                    ia1 = td.tripletlist[1, it, id]
                    ia2 = td.tripletlist[2, it, id]
                    ia3 = td.tripletlist[3, it, id]
                    if ia2 == ia2ia3pairs[1, ipairs] and ia3 == ia2ia3pairs[2, ipairs]:
                        im1 = mref[ia1]
                        im2 = mref[ia2]
                        im3 = mref[ia3]
                        x1 = xa[ia1] + xadd[im1]
                        x2 = xa[ia2] + xadd[im2]
                        x3 = xa[ia3] + xadd[im3]
                        y1 = ya[ia1] + yadd[im1]
                        y2 = ya[ia2] + yadd[im2]
                        y3 = ya[ia3] + yadd[im3]
                        z1 = za[ia1] + zadd[im1]
                        z2 = za[ia2] + zadd[im2]
                        z3 = za[ia3] + zadd[im3]
                        distmax12 = td.threebrmax[1, id] * (1.0 + td.threebrwidth[id])  # noqa: F841
                        distmax32 = td.threebrmax[2, id] * (1.0 + td.threebrwidth[id])  # noqa: F841
                        # Get the distances 1 - 2, 2 - 3
                        xdif = x1 - x2
                        ydif = y1 - y2
                        zdif = z1 - z2
                        xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                        d12sq = xdif * xdif + ydif * ydif + zdif * zdif
                        d12 = sqrtf(d12sq)
                        xdif = x3 - x2
                        ydif = y3 - y2
                        zdif = z3 - z2
                        xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                        d32sq = xdif * xdif + ydif * ydif + zdif * zdif
                        d32 = sqrtf(d32sq)
                        # Get distance 1 - 3
                        xdif = x1 - x3
                        ydif = y1 - y3
                        zdif = z1 - z3
                        xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                        d13sq = xdif * xdif + ydif * ydif + zdif * zdif
                        cost = (d12sq + d32sq - d13sq) / (2.0 * d12 * d32)
                        cosa = acosf(cost)
                        entrimin = f32(0)
                        # Get the closest value from the list of requested angles. Also save the value of
                        # ia1 that corresponds to this minimum.
                        for ib in range(1, td.nthreebangs[id] + 1):
                            icount = icount + 1
                            test = abs(cosa - td.threebangs[ib, id])
                            if ib == 1:
                                entrimin = test
                            elif test < entrimin:
                                entrimin = test
                        truncatefactor = f32(1.0)
                        # if(distmax12.gt.0.0.and.im1.ne.im2) then
                        #    d12star=d12/distmax12
                        #    truncatefactor=truncatefactor*truncate_threeb(d12star)
                        # end if
                        # if(distmax32.gt.0.0.and.im2.ne.im3) then
                        #    d32star=d32/distmax32
                        #    truncatefactor=truncatefactor*truncate_threeb(d32star)
                        # end if
                        # Use a harmonic three-body potential to emphasize the angular part.

                        entrisum = entrisum + entrimin * entrimin * truncatefactor
                # (commented-out code in the Fortran: tethering of ia3 between ia1 and ia2 when the
                # requested angle is zero)
                ensum = ensum + td.threebamp[id] * entrisum
                td.en_mole_triplet_dist[id] = td.en_mole_triplet_dist[id] + td.threebamp[id] * entrisum
    en_mole_triplets = ensum  # noqa: F841
    return en_mole_triplets_old


@routine
def pr_mole_triplets(imol, xadd, yadd, zadd, xa, ya, za, mref):
    """threebody_mod.f90:pr_mole_triplets -- three-body pressure associated with molecule
    imol (real*8).  Returns pr_mole_triplets.
    """
    # Range of triplets to be calculated

    pr_mole_triplets = 0.0
    if td.nthreeb <= 0:
        return pr_mole_triplets
    threebamp = td.threebamp
    threebangs = td.threebangs
    threebamps = td.threebamps
    threebawidth = td.threebawidth
    threebrmax = td.threebrmax
    nthreebangs = td.nthreebangs
    tripletlist = td.tripletlist
    ntripletsinmol = td.ntripletsinmol
    firsttripletinmol = td.firsttripletinmol
    prsum = 0.0
    for id in range(1, td.nthreeb + 1):
        if ntripletsinmol[imol, id] > 0 and abs(threebamp[id]) > 0.0:
            itf = firsttripletinmol[imol, id]
            itl = itf + ntripletsinmol[imol, id] - 1
            # Step through listed triplets
            for it in range(itf, itl + 1):
                # Get coordinates of atoms in triplet
                ia1 = tripletlist[1, it, id]
                ia2 = tripletlist[2, it, id]
                ia3 = tripletlist[3, it, id]
                im1 = mref[ia1]
                im2 = mref[ia2]
                im3 = mref[ia3]
                x1 = xa[ia1] + xadd[im1]
                x2 = xa[ia2] + xadd[im2]
                x3 = xa[ia3] + xadd[im3]
                y1 = ya[ia1] + yadd[im1]
                y2 = ya[ia2] + yadd[im2]
                y3 = ya[ia3] + yadd[im3]
                z1 = za[ia1] + zadd[im1]
                z2 = za[ia2] + zadd[im2]
                z3 = za[ia3] + zadd[im3]
                # Get the distances 1 - 2, 2 - 3
                xdif = x1 - x2
                ydif = y1 - y2
                zdif = z1 - z2
                xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                d12sq = xdif * xdif + ydif * ydif + zdif * zdif
                d12 = sqrtf(d12sq)
                xdif = x3 - x2
                ydif = y3 - y2
                zdif = z3 - z2
                xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                d32sq = xdif * xdif + ydif * ydif + zdif * zdif
                d32 = sqrtf(d32sq)
                if d12 < threebrmax[1, id] and d32 < threebrmax[2, id]:
                    # Get relative distances
                    d12star = d12 / threebrmax[1, id]  # noqa: F841
                    d32star = d32 / threebrmax[2, id]  # noqa: F841
                    # Get distance 1 - 3
                    xdif = x1 - x3
                    ydif = y1 - y3
                    zdif = z1 - z3
                    xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                    d13sq = xdif * xdif + ydif * ydif + zdif * zdif
                    # d13=sqrt(d13sq)
                    cost = (d12sq + d32sq - d13sq) / (2.0 * d12 * d32)
                    xcdif = xadd[im1] - xadd[im3]
                    ycdif = yadd[im1] - yadd[im3]
                    zcdif = zadd[im1] - zadd[im3]
                    xcdif, ycdif, zcdif = minimagec.minimagec(xcdif, ycdif, zcdif)
                    if d12 > 0.0 and d32 > 0.0:
                        # Don't forget d(cost)/d13
                        prfac = -(xcdif * xdif + ycdif * ydif + zcdif * zdif) / (d12 * d32)
                    else:
                        prfac = f32(0.0)
                    entrimin = f32(0.0)  # noqa: F841
                    ib = 1
                    ibmin = ib
                    arg = cost - threebangs[ib, id]
                    while ib < nthreebangs[id]:
                        ib = ib + 1
                        argtest = cost - threebangs[ib, id]
                        if abs(argtest) < abs(arg):
                            arg = argtest
                            ibmin = ib
                    pot = prfac * threebamps[ibmin, id] * arg / (threebawidth[id] * threebawidth[id])
                    prsum = prsum + float(threebamp[id] * pot)
    pr_mole_triplets = prsum
    return pr_mole_triplets


@routine
def truncate_threeb(x):
    """threebody_mod.f90:truncate_threeb -- smooth cut-off exp(1/(x-1)+1) for x < 1."""
    truncate_threeb = f32(0.0)
    if x >= 1.0:
        return truncate_threeb
    arg = r4(1.0) / (x - 1.0) + 1.0
    if arg > -30.0:
        truncate_threeb = expf(arg)
    else:
        truncate_threeb = f32(0.0)
    return truncate_threeb


@routine
def mod_morse(x, sigma, width):
    """threebody_mod.f90:mod_morse."""
    mod_morse = f32(0.0)
    if width > 0.0:
        arg = (sigma - x) / width
        if arg > 30.0:
            mod_morse = (arg - 1.0) * _EXP30
        elif arg < -30.0:
            mod_morse = f32(0.0)
        else:
            mod_morse = (arg - 1.0) * expf(arg)
    return mod_morse


@routine
def mod_morse1(x, sigma, width):
    """threebody_mod.f90:mod_morse1."""
    mod_morse1 = f32(0.0)
    if width > 0.0:
        arg = (sigma - x) / width
        arg = 0.5 * arg * arg
        # arg=-abs(arg)
        if arg > 30.0:
            mod_morse1 = f32(0.0)
        else:
            # mod_morse1=(arg-1.0)*exp(arg)
            mod_morse1 = expf(-arg)
        # mod_morse1=1.0/(1.0+arg*arg)
    # write(6,*) x,sigma,width,mod_morse1
    return mod_morse1


@routine
def grad_mod_morse1(x, sigma, width):
    """threebody_mod.f90:grad_mod_morse1."""
    grad_mod_morse1 = f32(0.0)
    if width > 0.0:
        arg = (sigma - x) / width
        # arg=-abs(arg)
        if arg < -30.0:
            grad_mod_morse1 = f32(0.0)
        else:
            # if(x.gt.sigma) then
            grad_mod_morse1 = -(arg * expf(arg))  # Gradient is +ve if x > sigma
            # else
            #    grad_mod_morse1=arg*exp(arg) !Gradient is -ve if x < sigma
            # end if
            grad_mod_morse1 = grad_mod_morse1 / width
    return grad_mod_morse1


@routine
def save_threeb_energies(im):
    """threebody_mod.f90:save_threeb_energies -- accept the new triplet energies of im and
    of its triplet neighbours.
    """
    if td.nthreeb > 0:
        td.ergold_triplet[im] = td.ergnew_triplet[im]
        # IMPORTANT! Include the neighbours of im which may have been affected by the movement of im!
        if td.nthreebneighbmol[im] > 0:
            for it in range(1, td.nthreebneighbmol[im] + 1):
                imref = td.threebneighbmol[it, im]
                td.ergold_triplet[imref] = td.ergnew_triplet[imref]
    return
