"""Translation of fortran_src/coordconstraint_mod.f90: coordination-number constraints.

The file holds two Fortran modules; both sets of routines are here, in the Fortran order:

* module ``coordconstraint_mod`` (the type ``coordconstraint``, see
  ``common/coordconstraint_mod.py``): ``init_constraint``, ``find_constraints``,
  ``find_type``, ``en_coordcon``, ``en_coordcon_total``;
* module ``coordconstraints_mod`` (data in ``common/coordconstraints_mod.py``):
  ``reallocate1d_cc``, ``initialise_coords``, ``en_mole_coord``, ``en_coord_total``,
  ``en_change_mole_coord``, ``save_coord_energies``, ``check_coord_bonded``,
  ``writecoordcons``.

A ``type(coordconstraint)`` dummy is the Python object, modified in place (and also
returned where docs/SIGNATURES.md lists ``this``).

``call reallocate1d_r(this%ds, n)`` (an allocatable dummy, which Python cannot rebind)
is ``this.ds = reallocation_routines.reallocate1d_r(this.ds, n)``: that module returns
the (re)allocated array (see library/reallocation_routines.py).

Bug FORTRAN_BUGS.md 1.2 is fixed (``# FIX``); the "stale-slot family" (section 3) is
reproduced: ``ishortest``/``dshortest`` are not reset between central atoms, the search
for ``ia`` in the list of ``ja`` stops at slot 2, and slot 1 of an atom without bonds is
read (its stale value from an earlier call).
"""

from ..common import coordconstraint_mod as cc_type
from ..common import coordconstraints_mod
from ..fortran import (
    PI4ATAN,
    Ch,
    cosf,
    f32,
    fmt,
    fstr,
    int_array,
    lw,
    nint,
    r4,
    read_internal,
    sqrtf,
    substr,
    units,
)
from ..trace import routine
from . import parse, reallocation_routines


# --------------------------------------------------------------------------- helpers
def _reallocate_bonds(this, nbondccsize, natom):
    """The four ``call reallocate2d_*(this%iabondcc|wtcc|tbondcc|adistcc, nbondccsize, natom)``."""
    this.iabondcc = reallocation_routines.reallocate2d_i(this.iabondcc, nbondccsize, natom)
    this.wtcc = reallocation_routines.reallocate2d_i(this.wtcc, nbondccsize, natom)
    this.tbondcc = reallocation_routines.reallocate2d_r(this.tbondcc, nbondccsize, natom)
    this.adistcc = reallocation_routines.reallocate2d_r(this.adistcc, nbondccsize, natom)


def _check_same_molecules(this, natom, mref, npass):
    """The block "Check whether the same molecules appear as a neighbour to each molecule"
    that find_constraints repeats verbatim for passes 1-4 (only the pass number differs).
    """
    nt = 0
    sumsame = 0
    for ia in range(1, natom + 1):
        nt = nt + 1
        if this.nbondcc[ia] > 0:
            for j in range(1, this.nbondcc[ia] + 1):
                if mref[ia] == mref[this.iabondcc[j, ia]]:
                    sumsame = sumsame + 1
                    units.write(6, lw(f'find_constraints - pass {npass}> Neighbouring atoms in same molecule: ',
                                      ia, this.iabondcc[j, ia]))
            for j in range(1, this.nbondcc[ia] - 1 + 1):
                jm = mref[this.iabondcc[j, ia]]
                k = j
                while k < this.nbondcc[ia]:
                    k = k + 1
                    km = mref[this.iabondcc[k, ia]]
                    if km == jm:
                        sumsame = sumsame + 1
                        units.write(6, lw(f'find_constraints - pass {npass}> Neighbouring atoms in same molecule: ',
                                          ia, this.iabondcc[j, ia], this.iabondcc[k, ia]))
    units.write(6, lw(f'find_constraints - pass {npass}> ', Ch(this.atc, 3),
                      ' Molecules listed as neighbours to themselves (should be zero) ', sumsame))


# =========================================================================== module coordconstraint_mod
@routine
def init_constraint(this, sometext, len_sometext=256):
    """coordconstraint_mod.f90:init_constraint (module coordconstraint_mod) -- set up one
    coordination constraint from the text of a ``coordcon`` line.

    ``len_sometext`` is Fortran ``len(sometext)`` (the caller's ``text`` is character*256).
    Returns this (modified in place).
    """
    nworddim = 100
    nwordspc = 10
    ncf = int_array(nworddim)
    ncl = int_array(nworddim)
    nwords = 0

    nchartext = len_sometext
    nwords = parse.parse(sometext, nwords, ncf, ncl, nwordspc, nchartext, nworddim)
    units.write(6, lw('init_constraint> ', nwords, ' ', substr(sometext, ncf[1], ncl[nwords])))
    this.atc = fstr(substr(sometext, ncf[1], ncl[1]), 3)
    this.nc = read_internal(substr(sometext, ncf[2], ncl[2])).int(this.nc)
    this.cctemp = read_internal(substr(sometext, ncf[3], ncl[3])).real(this.cctemp)
    this.ccdistort = read_internal(substr(sometext, ncf[4], ncl[4])).real(this.ccdistort)
    this.ccgaussweight = read_internal(substr(sometext, ncf[5], ncl[5])).real(this.ccgaussweight)
    this.ccgausswidth = read_internal(substr(sometext, ncf[6], ncl[6])).real(this.ccgausswidth)
    iread = 6
    this.ns = 0
    this.ntbangles = 0
    this.ntbintraangles = 0
    this.tbweight = f32(0.0)
    this.tbdistort = f32(1.0)
    this.tbgaussweight = f32(0.0)
    this.tbgausswidth = f32(1.0)
    this.tbintraweight = f32(0.0)
    this.tbintradistort = f32(1.0)
    this.tbintragaussweight = f32(0.0)
    this.tbintragausswidth = f32(1.0)
    this.hastb = False
    this.hastbintra = False
    while iread < nwords:
        iread = iread + 1
        if not this.hastb and not this.hastbintra:
            this.hastb = substr(sometext, ncf[iread], ncl[iread]) == 'angles'
            this.hastbintra = substr(sometext, ncf[iread], ncl[iread]) == 'anglesintra'
            if not this.hastb and not this.hastbintra:
                this.ns = this.ns + 1
                this.ats = reallocation_routines.reallocate1d_c(this.ats, 3, this.ns)
                this.ds = reallocation_routines.reallocate1d_r(this.ds, this.ns)
                this.dsmax = reallocation_routines.reallocate1d_r(this.dsmax, this.ns)
                this.ats[this.ns] = fstr(substr(sometext, ncf[iread], ncl[iread]), 3)
                iread = iread + 1
                this.ds[this.ns] = read_internal(substr(sometext, ncf[iread], ncl[iread])).real(this.ds[this.ns])
                iread = iread + 1
                this.dsmax[this.ns] = read_internal(substr(sometext, ncf[iread], ncl[iread])).real(this.dsmax[this.ns])
        elif this.hastb and not this.hastbintra:
            this.tbweight = read_internal(substr(sometext, ncf[iread], ncl[iread])).real(this.tbweight)
            iread = iread + 1
            this.tbdistort = read_internal(substr(sometext, ncf[iread], ncl[iread])).real(this.tbdistort)
            iread = iread + 1
            this.tbgaussweight = read_internal(substr(sometext, ncf[iread], ncl[iread])).real(this.tbgaussweight)
            iread = iread + 1
            this.tbgausswidth = read_internal(substr(sometext, ncf[iread], ncl[iread])).real(this.tbgausswidth)
            while iread < nwords and not this.hastbintra:
                iread = iread + 1
                this.hastbintra = substr(sometext, ncf[iread], ncl[iread]) == 'anglesintra'
                if not this.hastbintra:
                    this.ntbangles = this.ntbangles + 1
                    this.tbangles = reallocation_routines.reallocate1d_r(this.tbangles, this.ntbangles)
                    this.tbangles[this.ntbangles] = read_internal(
                        substr(sometext, ncf[iread], ncl[iread])).real(this.tbangles[this.ntbangles])
        elif this.hastbintra:
            this.tbintraweight = read_internal(substr(sometext, ncf[iread], ncl[iread])).real(this.tbintraweight)
            iread = iread + 1
            this.tbintradistort = read_internal(substr(sometext, ncf[iread], ncl[iread])).real(this.tbintradistort)
            iread = iread + 1
            this.tbintragaussweight = read_internal(
                substr(sometext, ncf[iread], ncl[iread])).real(this.tbintragaussweight)
            iread = iread + 1
            this.tbintragausswidth = read_internal(
                substr(sometext, ncf[iread], ncl[iread])).real(this.tbintragausswidth)
            while iread < nwords - 1:
                iread = iread + 1
                this.ntbintraangles = this.ntbintraangles + 1
                this.tbintraats = reallocation_routines.reallocate1d_c(this.tbintraats, 3, this.ntbintraangles)
                this.tbintraangles = reallocation_routines.reallocate1d_r(this.tbintraangles, this.ntbintraangles)
                this.tbintraats[this.ntbintraangles] = fstr(substr(sometext, ncf[iread], ncl[iread]), 3)
                iread = iread + 1
                this.tbintraangles[this.ntbintraangles] = read_internal(
                    substr(sometext, ncf[iread], ncl[iread])).real(this.tbintraangles[this.ntbintraangles])
    items = []
    for i in range(1, this.ns + 1):
        items += [Ch(this.ats[i], 3), ' ', this.ds[i]]
    units.write(6, lw('Coordination constraint ', Ch(this.atc, 3), ' with coordination number ',
                      this.nc, ' created. ', *items))
    if this.hastb:
        units.write(6, lw(' 3-body angles: ', this.ntbangles,
                          *[this.tbangles[i] for i in range(1, this.ntbangles + 1)]))
        piconv = PI4ATAN / r4(180.0)
        this.tbcosines = reallocation_routines.reallocate1d_r(this.tbcosines, this.ntbangles)
        for i in range(1, this.ntbangles + 1):
            this.tbcosines[i] = cosf(piconv * this.tbangles[i])
    if this.hastbintra:
        items = []
        for i in range(1, this.ntbintraangles + 1):
            items += [Ch(this.tbintraats[i], 3), ' ', this.tbintraangles[i]]
        units.write(6, lw(' intra 3-body angles: ', this.ntbintraangles, ' ', *items))
        piconv = PI4ATAN / r4(180.0)
        this.tbintracosines = reallocation_routines.reallocate1d_r(this.tbintracosines, this.ntbintraangles)
        for i in range(1, this.ntbintraangles + 1):
            this.tbintracosines[i] = cosf(piconv * this.tbintraangles[i])
    return this


@routine
def find_constraints(this, nmol, natom, naf, nal, xadd, yadd, zadd, xa, ya, za, mref,
                     itype, ntype, atype, cfrac, dmax):
    """coordconstraint_mod.f90:find_constraints (module coordconstraint_mod) -- assign the
    coordination "bonds" of constraint ``this`` from the current configuration.

    Returns this (modified in place).
    """
    from . import minimagec

    ishortest = None
    dshortest = None
    ndonors = None
    nacceptors = None
    ndonorsmax = None
    donors = None
    acceptors = None
    units.write(6, lw('Coordination constraint ', Ch(this.atc, 3), ' with coordination number ',
                      this.nc, ' bonded to: ',
                      *[x for i in range(1, this.ns + 1) for x in (Ch(this.ats[i], 3), ' ', this.ds[i])]))
    # Get the atomic fraction of the central atom type
    cacfrac = f32(0.0)
    it = find_type(this.atc, ntype, atype)
    if it > 0:
        cacfrac = cfrac[it]
    # Sum the atomic fractions of the second atom types
    sacfrac = f32(0.0)
    for ia in range(1, this.ns + 1):
        it = find_type(this.ats[ia], ntype, atype)
        if it > 0:
            sacfrac = sacfrac + cfrac[it]
    if sacfrac > 0.0:
        this.alternatenc = nint(f32(this.nc) * cacfrac / sacfrac)
    # Initialise arrays to store the bonds for each atom
    nbondccsize = max(this.nc, this.alternatenc)
    this.nbondcc = reallocation_routines.reallocate1d_i(this.nbondcc, natom)
    this.nearnbondcc = reallocation_routines.reallocate1d_i(this.nearnbondcc, natom)
    this.expecnbondcc = reallocation_routines.reallocate1d_i(this.expecnbondcc, natom)
    _reallocate_bonds(this, nbondccsize, natom)
    #      call reallocate1d_r(this%distcc,nbondccsize)
    this.nbondcc[:] = 0  # Set the inital number of bonds for each atom to zero
    this.nearnbondcc[:] = 0
    this.expecnbondcc[:] = 0  # Expected number of bonds for this coordination constraint
    # First job is to prepare a list of atoms of type atc
    ic = 0
    for ia in range(1, natom + 1):
        if atype[itype[ia]] == this.atc:
            ic = ic + 1
    this.atclistsize = ic
    this.atclist = reallocation_routines.reallocate1d_i(this.atclist, this.atclistsize)
    ic = 0
    for ia in range(1, natom + 1):
        if atype[itype[ia]] == this.atc:
            ic = ic + 1
            this.atclist[ic] = ia
            this.expecnbondcc[ia] = this.nc
        else:  # See if this atom is one of the alternate atom types
            it1 = 1
            found = atype[itype[ia]] == this.ats[it1]
            while it1 < this.ns and not found:
                it1 = it1 + 1
                found = atype[itype[ia]] == this.ats[it1]
            if found:
                this.expecnbondcc[ia] = this.alternatenc
    # Next is to prepare the list of nearest neigbours of atc up to the limit ns.
    # NB: ishortest/dshortest are allocated once and NOT reset between central atoms
    # (stale-slot family, FORTRAN_BUGS.md section 3; Rust coordconstraint.rs:466).  Their
    # initial contents are uninitialised in the Fortran; 0 here.
    ishortest = reallocation_routines.reallocate1d_i(ishortest, this.nc)
    dshortest = reallocation_routines.reallocate1d_r(dshortest, this.nc)
    sumoutside = 0
    for ic in range(1, this.atclistsize + 1):
        ia = this.atclist[ic]
        im = mref[ia]
        nshortest = 0
        for ja in range(1, natom + 1):
            jm = mref[ja]
            # First check ja is not in the same molecule as ia, and also is on a different molecule to those on the list of shortest distances
            found = im == jm
            it = 0
            while it < nshortest and not found:
                it = it + 1
                found = jm == mref[ishortest[it]]
                #               found=ja.eq.ishortest(it)
            if not found:
                it = 0
                while it < this.ns and not found:
                    it = it + 1
                    found = this.ats[it] == atype[itype[ja]]
                if found:
                    # Get the separation of ia and ja
                    xdif = (xadd[im] + xa[ia]) - (xadd[jm] + xa[ja])
                    ydif = (yadd[im] + ya[ia]) - (yadd[jm] + ya[ja])
                    zdif = (zadd[im] + za[ia]) - (zadd[jm] + za[ja])
                    xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                    dsq = xdif * xdif + ydif * ydif + zdif * zdif
                    d12 = sqrtf(dsq)
                    # Step through the list of shortest distances and if d12 is smaller than one of these signal the fact
                    it = 0
                    found = False
                    while it < nshortest and not found:
                        it = it + 1
                        found = d12 < dshortest[it]
                    if found:
                        if nshortest < this.nc:
                            nshortest = nshortest + 1
                        # move all the existing values up 1 step
                        it2 = nshortest
                        while it2 > it:
                            it1 = it2 - 1
                            ishortest[it2] = ishortest[it1]
                            dshortest[it2] = dshortest[it1]
                            it2 = it1
                        # Finally insert the new values at it
                        ishortest[it] = ja
                        dshortest[it] = d12
                    elif nshortest < this.nc:
                        # Otherwise tag the value on the end of the list, but only if the number in the list is less than the maximum
                        nshortest = nshortest + 1
                        ishortest[nshortest] = ja
                        dshortest[nshortest] = d12
        # Finally save the list of nearest neighbours and sum the number of atoms outside the specified limit.
        for it in range(1, this.nc + 1):
            ja = ishortest[it]
            # Get the expected bond distance for atom ja
            # NB: if fewer than nc neighbours were found, ishortest(it) is the stale value of
            # the previous central atom (or 0 for the first: then atype(itype(0)) is '' and
            # not found; the Fortran reads out of bounds, Rust coordconstraint.rs:105 stops).
            it1 = 0
            found = False
            while it1 < this.ns and not found:
                it1 = it1 + 1
                found = atype[itype[ja]] == this.ats[it1]
            if found:  # Just to be safe, only set the atom bond numbers if the atom is correctly recognised.
                if dshortest[it] > this.dsmax[it1]:
                    sumoutside = sumoutside + 1  # Count the number of atoms outside the specified maximum distance
                dsfound = this.ds[it1]
                #               if(this%nbondcc(ia).lt.this%nc.and.dshortest(it).le.this%dsmax(it1)) then !Only set the bonds if they are within range
                this.nbondcc[ia] = this.nbondcc[ia] + 1
                if this.nbondcc[ia] > nbondccsize:
                    nbondccsize = int(this.nbondcc[ia])
                    _reallocate_bonds(this, nbondccsize, natom)
                this.iabondcc[this.nbondcc[ia], ia] = ja
                this.wtcc[this.nbondcc[ia], ia] = 1
                this.tbondcc[this.nbondcc[ia], ia] = dsfound
                this.adistcc[this.nbondcc[ia], ia] = dshortest[it]
                if itype[ia] != itype[ja]:
                    # Only set the alternate bonds if the atoms are of different types. A maximum of this%alternatenc such bonds are allowed
                    # Need to insert this bond in the order of shortest - longest
                    # First find the shortest
                    is_ = 0
                    found = False
                    while is_ < this.nbondcc[ja] and not found:
                        is_ = is_ + 1
                        found = dshortest[it] < this.adistcc[is_, ja]
                    if found:
                        # If a longer value was found, then insert before the longer value. First move all existing values up 1
                        this.nbondcc[ja] = this.nbondcc[ja] + 1
                        if this.nbondcc[ja] > nbondccsize:
                            nbondccsize = int(this.nbondcc[ja])
                            _reallocate_bonds(this, nbondccsize, natom)
                        is2 = int(this.nbondcc[ja])
                        while is2 > is_:
                            is1 = is2 - 1
                            this.iabondcc[is2, ja] = this.iabondcc[is1, ja]
                            this.wtcc[is2, ja] = this.wtcc[is1, ja]
                            this.tbondcc[is2, ja] = this.tbondcc[is1, ja]
                            this.adistcc[is2, ja] = this.adistcc[is1, ja]
                            is2 = is1
                        # Then save the new values at is
                        this.iabondcc[is_, ja] = ia
                        this.wtcc[is_, ja] = 1
                        this.tbondcc[is_, ja] = dsfound
                        this.adistcc[is_, ja] = dshortest[it]
                    else:  # Otherwise simply add new values at end of list if the list is not full
                        this.nbondcc[ja] = this.nbondcc[ja] + 1
                        if this.nbondcc[ja] > nbondccsize:
                            nbondccsize = int(this.nbondcc[ja])
                            _reallocate_bonds(this, nbondccsize, natom)
                        this.iabondcc[this.nbondcc[ja], ja] = ia
                        this.wtcc[this.nbondcc[ja], ja] = 1
                        this.tbondcc[this.nbondcc[ja], ja] = dsfound
                        this.adistcc[this.nbondcc[ja], ja] = dshortest[it]
    # Check whether the same molecules appear as a neighbour to each molecule
    _check_same_molecules(this, natom, mref, 1)
    # If the number of bonds for a given atom is greater than expected, we need to remove the extra atoms from the list and also from the lists of the atoms
    # to which they are bonded.
    for ia in range(1, natom + 1):
        # If nbondcc > expecnbondcc it means this atom occurs too often in the bonding lists of the other atoms, so it should be removed
        it = int(this.nbondcc[ia])
        while it > this.expecnbondcc[ia]:
            ja = this.iabondcc[it, ia]
            # Remove atom ia as neighbour of ja
            # NB: the search stops at it1 = 2, slot 1 is only tested when it is the last
            # one (stale-slot family, Rust coordconstraint.rs:166).  If nbondcc(ja) is 0
            # the Fortran reads iabondcc(0,ja); here the unused row 0 (value 0).
            it1 = int(this.nbondcc[ja])
            found = ia == this.iabondcc[it1, ja]
            while it1 > 2 and not found:
                it1 = it1 - 1
                found = ia == this.iabondcc[it1, ja]
            if found:
                # Save the values for this atom
                tbondsave = this.tbondcc[it1, ja]  # noqa: F841
                adistsave = this.adistcc[it1, ja]  # noqa: F841
                # Remove ia from the list of ja
                for it2 in range(it1, this.nbondcc[ja] - 1 + 1):
                    this.iabondcc[it2, ja] = this.iabondcc[it2 + 1, ja]
                    this.wtcc[it2, ja] = this.wtcc[it2 + 1, ja]
                    this.tbondcc[it2, ja] = this.tbondcc[it2 + 1, ja]
                    this.adistcc[it2, ja] = this.adistcc[it2 + 1, ja]
                this.nbondcc[ja] = this.nbondcc[ja] - 1
            it = it - 1
        this.nbondcc[ia] = it  # Finally reset the reduced number of bonds for ia.
    # Check whether the same molecules appear as a neighbour to each molecule
    _check_same_molecules(this, natom, mref, 2)
    # Need to consider the case where an atom is UNDERbonded - can suitable combinations be found??
    for ic in range(1, this.atclistsize + 1):
        # For simplicity we only look for atoms of type this%atc which are under-bonded - the opposite bonding is set automatically when a bond is found.
        ia = this.atclist[ic]
        im = mref[ia]
        it = int(this.nbondcc[ia])
        while it < this.expecnbondcc[ia]:
            #            write(6,*) 'Atom ',ia,' type ',atype(itype(ia)),' ',this%nbondcc(ia),this%expecnbondcc(ia)
            it = it + 1
            # Need to find an atom of one of the types this%ats which is ALSO underbonded.
            ja = 0
            found = False
            it1 = 0
            while ja < natom and not found:
                ja = ja + 1
                found = this.nbondcc[ja] < this.expecnbondcc[ja] and mref[ja] != mref[ia]  # Don't waste time on atoms which are in the same molecule as ia
                if found:
                    it1 = 1
                    found = atype[itype[ja]] == this.ats[it1]
                    while it1 < this.ns and not found:
                        it1 = it1 + 1
                        found = atype[itype[ja]] == this.ats[it1]
                    if found:
                        # Need to check molecule containing ja is not already listed as a neighbour of the molecule containing ia
                        # NB: slot 1 is read even when ia has no bonds (stale value of an
                        # earlier call; stale-slot family, Rust coordconstraint.rs:192).
                        it2 = 1
                        found = mref[ja] != mref[this.iabondcc[it2, ia]]
                        while it2 < this.nbondcc[ia] and found:
                            it2 = it2 + 1
                            found = mref[ja] != mref[this.iabondcc[it2, ia]]
            it1save = it1
            #            write(6,*) 'Atom ',ja,' type ',atype(itype(ja)),' ',this%ats(it1save),' ',this%nbondcc(ja),this%expecnbondcc(ja)
            # If one was found, now look for other possibilities, choosing the one that is nearest to ia
            if found:
                # measure the distance from ia to ja
                xdif = xadd[mref[ia]] + xa[ia] - xadd[mref[ja]] - xa[ja]
                ydif = yadd[mref[ia]] + ya[ia] - yadd[mref[ja]] - ya[ja]
                zdif = zadd[mref[ia]] + za[ia] - zadd[mref[ja]] - za[ja]
                xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                dsqmin = xdif * xdif + ydif * ydif + zdif * zdif
                jamin = ja
                # Step through the remaining atoms to see if any are closer
                while ja < natom:
                    ja = ja + 1
                    found = this.nbondcc[ja] < this.expecnbondcc[ja] and mref[ja] != mref[ia]
                    if found:  # Check that this atom is on the list possible second atoms
                        it1 = 1
                        found = atype[itype[ja]] == this.ats[it1]
                        while it1 < this.ns and not found:
                            it1 = it1 + 1
                            found = atype[itype[ja]] == this.ats[it1]
                        if found:
                            # Need to check molecule containing ja is not already listed as a neighbour of the molecule containing ia
                            it2 = 1
                            found = mref[ja] != mref[this.iabondcc[it2, ia]]
                            while it2 < this.nbondcc[ia] and found:
                                it2 = it2 + 1
                                found = mref[ja] != mref[this.iabondcc[it2, ia]]
                    if found:
                        xdif = xadd[mref[ia]] + xa[ia] - xadd[mref[ja]] - xa[ja]
                        ydif = yadd[mref[ia]] + ya[ia] - yadd[mref[ja]] - ya[ja]
                        zdif = zadd[mref[ia]] + za[ia] - zadd[mref[ja]] - za[ja]
                        xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                        dsq = xdif * xdif + ydif * ydif + zdif * zdif
                        if dsq < dsqmin:
                            jamin = ja
                            dsqmin = dsq
                            it1save = it1
                this.nbondcc[jamin] = this.nbondcc[jamin] + 1
                this.iabondcc[this.nbondcc[jamin], jamin] = ia
                this.wtcc[this.nbondcc[jamin], jamin] = 1
                this.tbondcc[this.nbondcc[jamin], jamin] = this.ds[it1save]
                this.adistcc[this.nbondcc[jamin], jamin] = sqrtf(dsqmin)
                this.nbondcc[ia] = this.nbondcc[ia] + 1
                this.iabondcc[this.nbondcc[ia], ia] = jamin
                this.wtcc[this.nbondcc[ia], ia] = 1
                this.tbondcc[this.nbondcc[ia], ia] = this.ds[it1save]
                this.adistcc[this.nbondcc[ia], ia] = this.adistcc[this.nbondcc[jamin], jamin]
                # Strictly the distance between ia and ja should be inserted in the correct order.
                # However from this point onwards adistcc is not used, so we won't bother for the time being.
                # In any case it is likely this atom (ja) is further away than other atoms already in the list.
    # Check whether the same molecules appear as a neighbour to each molecule
    _check_same_molecules(this, natom, mref, 3)
    # Save the number of nearest neighbours for each atom.
    this.nearnbondcc[:] = this.nbondcc
    # Assess how well the bonding criterion has been satisfied
    nt = 0
    sumratio = f32(0.0)
    for ia in range(1, natom + 1):
        nt = nt + 1
        if this.expecnbondcc[ia] == this.nbondcc[ia]:
            sumratio = sumratio + 1
    if nt > 0:
        sumratio = 100.0 * sumratio / f32(nt)
    ia = this.atclist[this.atclistsize]
    items = []
    for it in range(1, this.nc + 1):
        # NB: an unfilled ishortest(it) is 0 here: atype(itype(0)) is written as 3 blanks
        items += [ishortest[it], ' ', Ch(atype[itype[ishortest[it]]], 3), ' ', dshortest[it]]
    units.write(6, lw('find_constraints: ', ia, ' ', Ch(atype[itype[ia]], 3), ' ', *items))
    units.write(6, lw('find_constraints: ', Ch(this.atc, 3), ' Fraction outside specified distance(s) = ',
                      f32(sumoutside) / f32(this.nc * this.atclistsize)))
    units.write(6, lw('find_constraints: ', Ch(this.atc, 3), ' Percentage constraints correct = ', sumratio))
    # If these atoms are bonded via tb angles, need to define the corresponding bond distances
    if this.hastb:
        piconv = PI4ATAN / r4(180.0)  # noqa: F841
        for ic in range(1, this.atclistsize + 1):
            ia = this.atclist[ic]
            im = mref[ia]
            i = 0
            #            write(6,*) ia,this%nbondcc(ia),this%nearnbondcc(ia)
            while i < this.nearnbondcc[ia] - 1:
                i = i + 1
                ja1 = this.iabondcc[i, ia]
                jm1 = mref[ja1]
                # Get the separation of ia and ja1
                #               tb1=this%tbondcc(i,ia) !Use the expected separation distance to determine the three-body distance constraint.
                tb1 = this.adistcc[i, ia]  # Use the actual separation to determine the three-body distance constraint.
                j = i
                while j < this.nearnbondcc[ia]:
                    j = j + 1
                    ja2 = this.iabondcc[j, ia]
                    jm2 = mref[ja2]
                    # Get the separation of ia and ja2
                    #                  tb2=this%tbondcc(j,ia) !Use the expected separation to determine the three-body distance constraint
                    # FIX (FORTRAN_BUGS.md 1.2): the Fortran uses adistcc(i,ia), building the angle
                    # from the same arm twice; the separation of ia and ja2 is adistcc(j,ia)
                    # (as the Rust, coordconstraint.rs find_constraints angle pass).
                    tb2 = this.adistcc[j, ia]  # Use the actual separation to determine the three-body distance constraint.
                    # Get the separation of ja1 and ja2
                    xdif = (xadd[jm1] + xa[ja1]) - (xadd[jm2] + xa[ja2])
                    ydif = (yadd[jm1] + ya[ja1]) - (yadd[jm2] + ya[ja2])
                    zdif = (zadd[jm1] + za[ja1]) - (zadd[jm2] + za[ja2])
                    xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                    dist12 = xdif * xdif + ydif * ydif + zdif * zdif
                    dsq = tb1 * tb1 + tb2 * tb2 - 2.0 * tb1 * tb2 * this.tbcosines[1]  # If there are more than one angle we can check which one this particular pair is closer to.
                    tb = sqrtf(dsq)
                    dsave = abs(dist12 - dsq)
                    itb = 1
                    # If more than one angle specified, find the closest to dist12
                    while itb < this.ntbangles:
                        itb = itb + 1
                        dsq = tb1 * tb1 + tb2 * tb2 - 2.0 * tb1 * tb2 * this.tbcosines[itb]
                        dtest = abs(dist12 - dsq)
                        if dtest < dsave:
                            dsave = dtest
                            tb = sqrtf(dsq)
                    # It is possible atom ja2 already occurs in the coordination shell of ja1 (due to direct bonding to different centre atoms).
                    # To avoid this, we replace the distance between ja1 and ja2 if it is shorter, but don't increase the number of bonds
                    itb = int(this.nearnbondcc[ja1])
                    found = False
                    while itb < this.nbondcc[ja1] and not found:
                        itb = itb + 1
                        found = ja2 == this.iabondcc[itb, ja1]
                    if found:
                        tb = min(abs(this.tbondcc[itb, ja1]), tb)
                        this.tbondcc[itb, ja1] = -tb
                        # Need to repeat this for ja2 - first find which atom number ja1 is
                        itb = int(this.nearnbondcc[ja2])
                        found = False
                        while itb < this.nbondcc[ja2] and not found:
                            itb = itb + 1
                            found = ja1 == this.iabondcc[itb, ja2]
                        if found:
                            this.tbondcc[itb, ja2] = -tb
                    else:
                        # Otherwise we need to check the molecule containing ja2 does not already occur in the coordination shell of ja1
                        # NB: slot 1 read even if ja1/ja2 has no bonds (stale-slot family)
                        itb = 1
                        found = mref[this.iabondcc[itb, ja1]] != mref[ja2]
                        while itb < this.nbondcc[ja1] and found:
                            itb = itb + 1
                            found = mref[this.iabondcc[itb, ja1]] != mref[ja2]
                        # Repeat for ja2
                        itb = 1
                        found = found and mref[this.iabondcc[itb, ja2]] != mref[ja1]
                        while itb < this.nbondcc[ja2] and found:
                            itb = itb + 1
                            found = found and mref[this.iabondcc[itb, ja2]] != mref[ja1]
                        if found:
                            this.nbondcc[ja1] = this.nbondcc[ja1] + 1
                            this.nbondcc[ja2] = this.nbondcc[ja2] + 1
                            nbondccsize = int(max(nbondccsize, this.nbondcc[ja1], this.nbondcc[ja2]))
                            _reallocate_bonds(this, nbondccsize, natom)
                            this.iabondcc[this.nbondcc[ja1], ja1] = ja2
                            this.iabondcc[this.nbondcc[ja2], ja2] = ja1
                            this.wtcc[this.nbondcc[ja1], ja1] = 2
                            this.wtcc[this.nbondcc[ja2], ja2] = 2
                            this.tbondcc[this.nbondcc[ja1], ja1] = -tb  # Negative will indicate it has a different weighting to direct bonds.
                            this.tbondcc[this.nbondcc[ja2], ja2] = -tb
                            this.adistcc[this.nbondcc[ja1], ja1] = sqrtf(dist12)  # Negative will indicate it has a different weighting to direct bonds.
                            this.adistcc[this.nbondcc[ja2], ja2] = this.adistcc[this.nbondcc[ja1], ja1]
        # Check whether the same molecules appear as a neighbour to each molecule
        _check_same_molecules(this, natom, mref, 4)
    units.write(6, lw('Got to tbintra', this.hastbintra))
    if this.hastbintra:
        # For intra three-body angles, the object is to identify atoms of type tbintra ats bound (same molecule) to ia
        # and identify a neighbour of ia with which they are to form the specified angle. If more than one of these atoms is specified (nbtintraangles>1), then
        # a different neighbour must be chosen for each.
        # Keep track of donors and acceptors on each molecule
        ndonors = reallocation_routines.reallocate1d_i(ndonors, nmol)
        nacceptors = reallocation_routines.reallocate1d_i(nacceptors, nmol)
        ndonorsmax = reallocation_routines.reallocate1d_i(ndonorsmax, nmol)
        #         call reallocate1d_i(nacceptorsmax,nmol)
        # Step through the intra angles
        for itb in range(1, this.ntbintraangles + 1):
            ndonors[:] = 0
            nacceptors[:] = 0
            ndonaccsize = 0
            # The maximum number of donors for each molecule is set to the number of atoms of type this%tbintraats(1) on each molecule, i.e. the no. of bonds it can DONATE
            it1 = find_type(this.tbintraats[itb], ntype, atype)  # Type number of atom to be searched for
            for im in range(1, nmol + 1):
                ndonorsmax[im] = 0
                for ia in range(naf[im], nal[im] + 1):
                    if itype[ia] == it1:
                        ndonorsmax[im] = ndonorsmax[im] + 1
            # First pass (coordconstraint_mod.f90:759-861), then the
            # "Second pass to try to increase the number of donor molecules for those which are
            # underbonded" (862-965).  The two loops are statement-for-statement identical in
            # the Fortran apart from the form of the entry test (a logical variable in the
            # first, a direct if in the second, same expression), so the body is written once
            # and run twice.
            for _pass in (1, 2):
                for ic in range(1, this.atclistsize + 1):
                    ia = this.atclist[ic]
                    im = mref[ia]
                    # Search the atoms of molecule im to find which match tbintraats. We expect one constraint for each listed intra angle, but only up to the number of atoms of listed types in molecule im
                    for ja1 in range(naf[im], nal[im] + 1):
                        jm1 = mref[ja1]
                        # (Maximum number of acceptors assumed equal to maximum number of donors)
                        found = nacceptors[jm1] < ndonorsmax[jm1] and atype[itype[ja1]] == this.tbintraats[itb] and ja1 != ia
                        if found:
                            # Get the separation of ia and ja1
                            xdif = xa[ja1] - xa[ia]
                            ydif = ya[ja1] - ya[ia]
                            zdif = za[ja1] - za[ia]
                            xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                            dsq = xdif * xdif + ydif * ydif + zdif * zdif
                            tb1 = sqrtf(dsq)
                            # Search the near-neighbour bonds of ia and find that which is the closest to the expected separation from ja1
                            inearbond = 0
                            dsave = f32(0.0)
                            jamin = 0
                            dsqmin = f32(0.0)
                            dsqtbmin = f32(0.0)
                            found = False
                            while inearbond < this.nearnbondcc[ia]:
                                inearbond = inearbond + 1
                                ja2 = this.iabondcc[inearbond, ia]
                                jm2 = mref[ja2]
                                # First, check the maximum number of donors for jm2 is not already reached
                                # (second pass: Check whether jm2 has any unsaturated bonds)
                                found = ndonors[jm2] < ndonorsmax[jm2]
                                # Check that ja2 is not already occur in the coordination shell of ja1
                                it1 = 0
                                while it1 < this.nbondcc[ja1] and found:
                                    it1 = it1 + 1
                                    found = ja2 != this.iabondcc[it1, ja1]
                                # Check that ja2 is not already a near neighbour of ja1
                                it1 = 0
                                while it1 < this.nbondcc[ja2] and found:
                                    it1 = it1 + 1
                                    found = ja1 != this.iabondcc[it1, ja2]
                                # If accepted jm1 will become a donor for molecule jm2, and jm2 will become an acceptor for jm1
                                # Next, check jm1 is not already a donor for jm2
                                ic1 = 0
                                while ic1 < ndonors[jm2] and found:
                                    ic1 = ic1 + 1
                                    found = jm1 != donors[ic1, jm2]
                                # Two molecules cannot be both donor and acceptor for each other, so we need to ensure jm1 is not already an acceptor for jm2
                                ic2 = 0
                                while ic2 < nacceptors[jm2] and found:
                                    ic2 = ic2 + 1
                                    found = jm1 != acceptors[ic2, jm2]
                                if found:
                                    # Get the separation of ia and ja2
                                    tb2 = this.adistcc[inearbond, ia]
                                    # Get expected separation of ja1 and ja2
                                    dsqtb = tb1 * tb1 + tb2 * tb2 - 2.0 * tb1 * tb2 * this.tbintracosines[itb]
                                    # Get the separation of ja1 and ja2
                                    xdif = xadd[jm1] + xa[ja1] - (xadd[jm2] + xa[ja2])
                                    ydif = yadd[jm1] + ya[ja1] - (yadd[jm2] + ya[ja2])
                                    zdif = zadd[jm1] + za[ja1] - (zadd[jm2] + za[ja2])
                                    xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                                    dsq = xdif * xdif + ydif * ydif + zdif * zdif
                                    dtest = abs(dsq - dsqtb)
                                    if dtest < dsave or dsave == 0.0:
                                        dsave = dtest
                                        jamin = ja2
                                        dsqmin = dsq
                                        dsqtbmin = dsqtb
                            if jamin != 0:  # At least one valid bond was found
                                ja2 = jamin
                                jm2 = mref[ja2]
                                tb = sqrtf(dsqtbmin)
                                this.nbondcc[ja1] = this.nbondcc[ja1] + 1
                                this.nbondcc[ja2] = this.nbondcc[ja2] + 1
                                nbondccsize = int(max(nbondccsize, this.nbondcc[ja1], this.nbondcc[ja2]))
                                _reallocate_bonds(this, nbondccsize, natom)
                                this.iabondcc[this.nbondcc[ja1], ja1] = ja2
                                this.iabondcc[this.nbondcc[ja2], ja2] = ja1
                                this.wtcc[this.nbondcc[ja1], ja1] = 3
                                this.wtcc[this.nbondcc[ja2], ja2] = 3
                                this.tbondcc[this.nbondcc[ja1], ja1] = -tb
                                this.tbondcc[this.nbondcc[ja2], ja2] = -tb
                                this.adistcc[this.nbondcc[ja1], ja1] = sqrtf(dsqmin)  # Negative will indicate it has a different weighting to direct bonds.
                                this.adistcc[this.nbondcc[ja2], ja2] = this.adistcc[this.nbondcc[ja1], ja1]
                                nacceptors[jm1] = nacceptors[jm1] + 1
                                ndonors[jm2] = ndonors[jm2] + 1
                                ndonaccsize = int(max(ndonaccsize, nacceptors[jm1], ndonors[jm2]))
                                donors = reallocation_routines.reallocate2d_i(donors, ndonaccsize, nmol)
                                acceptors = reallocation_routines.reallocate2d_i(acceptors, ndonaccsize, nmol)
                                acceptors[nacceptors[jm1], jm1] = jm2
                                donors[ndonors[jm2], jm2] = jm1
            # Check whether any molecule pairs are both donor and acceptor to each other - this is not allowed
            sumsame = 0
            ndonorstotal = 0
            nacceptorstotal = 0
            for im in range(1, nmol + 1):
                nacceptorstotal = nacceptorstotal + nacceptors[im]
                ndonorstotal = ndonorstotal + ndonors[im]
                if ndonors[im] < ndonorsmax[im]:
                    ia = naf[im]
                # Check the acceptors of im
                for ic1 in range(1, nacceptors[im] + 1):
                    jm = acceptors[ic1, im]
                    # check if jm is a donor of im
                    for ic2 in range(1, ndonors[im] + 1):
                        if donors[ic2, im] == jm:
                            sumsame = sumsame + 1
            units.write(6, lw('find_constraints - pass 5> ', Ch(this.atc, 3), itb,
                              ' Number of donor-donor pairs (should be zero) ', sumsame))
            units.write(6, lw('find_constraints - pass 5> Total numbers of acceptors and donors: ',
                              int(nacceptorstotal), int(ndonorstotal)))
        # check for unphysically short distances between tbintraats type atoms - list the respective molecules and their bonds
        # (commented-out diagnostic in the Fortran, coordconstraint_mod.f90:996-1035)
    return this


@routine
def find_type(alabel, ntype, atype):
    """coordconstraint_mod.f90:find_type (module coordconstraint_mod) -- number of atom
    type ``alabel`` in ``atype(1:ntype)``, 0 if absent.
    """
    found = False
    find_type = 0
    i = 0
    while i < ntype and not found:
        i = i + 1
        found = alabel == atype[i]
    if found:
        find_type = i
    return find_type


@routine
def en_coordcon(this, ia, xadd, yadd, zadd, xa, ya, za, mref, itype, qtype, npairref, amassr, tethered):
    """coordconstraint_mod.f90:en_coordcon (module coordconstraint_mod) -- the coordination
    constraint energy of atom ia in coordination constraint this.

    Returns (en_coordcon, this).
    """
    from . import minimagec, molfit

    en_coordcon = f32(0.0)
    nbond = this.nbondcc[ia]
    if nbond > 0:
        im = mref[ia]
        x1 = xadd[im] + xa[ia]
        y1 = yadd[im] + ya[ia]
        z1 = zadd[im] + za[ia]
        for j in range(1, nbond + 1):
            ja = this.iabondcc[j, ia]
            jm = mref[ja]
            bothtethered = tethered[im] and tethered[jm] and qtype[ia] == 1 and qtype[ja] == 1  # Don't calculate pair energies if both molecules are tethered and both atoms are Q-type
            if not bothtethered:
                x2 = xadd[jm] + xa[ja]
                y2 = yadd[jm] + ya[ja]
                z2 = zadd[jm] + za[ja]
                xdif = x2 - x1
                ydif = y2 - y1
                zdif = z2 - z1
                # apply minimum image convention
                xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                dsq = xdif * xdif + ydif * ydif + zdif * zdif
                dtest = sqrtf(dsq)
                tb = this.tbondcc[j, ia]
                wtms = this.cctemp
                distort = this.ccdistort
                gweight = this.ccgaussweight
                gwidth = this.ccgausswidth
                if this.wtcc[j, ia] == 2:
                    wtms = this.tbweight
                    distort = this.tbdistort
                    gweight = this.tbgaussweight
                    gwidth = this.tbgausswidth
                elif this.wtcc[j, ia] == 3:
                    wtms = this.tbintraweight
                    distort = this.tbintradistort
                    gweight = this.tbintragaussweight
                    gwidth = this.tbintragausswidth
                wtm = amassr[npairref[itype[ia], itype[ja]]] / abs(tb)
                en = molfit.en_atom_inter_scale_coord(dtest, tb, wtm, wtms, distort, gweight, gwidth)
                #     write(6,*) 'en_coordcon> ',ia,ja,tb,dtest,wtm,wtms,width,en
                en_coordcon = en_coordcon + en
        # (commented-out three-body energy calculation in the Fortran, coordconstraint_mod.f90:1118-1170)
    return en_coordcon, this


@routine
def en_coordcon_total(this, natom, xadd, yadd, zadd, xa, ya, za, mref, itype, qtype, npairref, amassr, tethered):
    """coordconstraint_mod.f90:en_coordcon_total (module coordconstraint_mod) -- the total
    coordination constraint energy of all atoms in coordination constraint this.

    NB: unlike en_coordcon this divides by tb (not abs(tb)) and scales distort by
    sqrt(wtms), as the Fortran does.
    Returns en_coordcon_total.
    """
    from . import minimagec, molfit

    en_coordcon_total = f32(0.0)
    for ia in range(1, natom + 1):
        if this.nbondcc[ia] > 0:
            im = mref[ia]
            x1 = xadd[im] + xa[ia]
            y1 = yadd[im] + ya[ia]
            z1 = zadd[im] + za[ia]
            for j in range(1, this.nbondcc[ia] + 1):
                ja = this.iabondcc[j, ia]
                jm = mref[ja]
                #               if(ja.gt.ia) then ! Prevent double counting in sum of all coordination constraint energies
                tb = this.tbondcc[j, ia]
                wtms = this.cctemp
                distort = this.ccdistort
                gweight = this.ccgaussweight
                gwidth = this.ccgausswidth
                if this.wtcc[j, ia] == 2:
                    wtms = this.tbweight
                    distort = this.tbdistort
                    gweight = this.tbgaussweight
                    gwidth = this.tbgausswidth
                elif this.wtcc[j, ia] == 3:
                    wtms = this.tbintraweight
                    distort = this.tbintradistort
                    gweight = this.tbintragaussweight
                    gwidth = this.tbintragausswidth
                wtm = amassr[npairref[itype[ia], itype[ja]]] / tb
                distort = sqrtf(wtms) * distort
                x2 = xadd[jm] + xa[ja]
                y2 = yadd[jm] + ya[ja]
                z2 = zadd[jm] + za[ja]
                xdif = x2 - x1
                ydif = y2 - y1
                zdif = z2 - z1
                # apply minimum image convention
                xdif, ydif, zdif = minimagec.minimagec(xdif, ydif, zdif)
                dsq = xdif * xdif + ydif * ydif + zdif * zdif
                dtest = sqrtf(dsq)
                bothtethered = tethered[im] and tethered[jm] and qtype[ia] == 1 and qtype[ja] == 1  # Don't calculate pair energies if both molecules are tethered and both atoms are Q-type
                if not bothtethered:
                    en_coordcon_total = en_coordcon_total + molfit.en_atom_inter_scale_coord(
                        dtest, tb, wtm, wtms, distort, gweight, gwidth)
                #               end if
    return en_coordcon_total


# =========================================================================== module coordconstraints_mod
@routine
def reallocate1d_cc(array, newsize):
    """coordconstraint_mod.f90:reallocate1d_cc (module coordconstraints_mod) -- reallocate
    an allocatable array of type(coordconstraint), keeping the existing elements.

    The array is a Python list with index 0 unused (see common/coordconstraints_mod.py).
    An allocated list is resized *in place*; an unallocated one (None) cannot be, so when
    the actual argument is ``coordconstraints_mod.coordconstraints`` (the only caller,
    read_pcof) the module variable is rebound.  New elements (uninitialised in the
    Fortran) are fresh coordconstraint objects.
    Returns the (new) array as well (SIGNATURE_EXCEPTIONS.md); callers may ignore it.
    """
    is_module_array = array is coordconstraints_mod.coordconstraints
    if array is not None:
        oldsize = len(array) - 1
    else:
        oldsize = 0
    # Do nothing if array is already the requested size!
    if oldsize == newsize:
        return array
    if newsize > 0:
        if array is not None:
            # arraysave(i)=array(i) ... array(i)=arraysave(i): the elements keep their values
            del array[min(oldsize, newsize) + 1 :]
            new = array
        else:
            new = [None]
        while len(new) - 1 < newsize:
            new.append(cc_type.coordconstraint())
    else:
        new = None
    if is_module_array:
        coordconstraints_mod.coordconstraints = new
    return new


@routine
def initialise_coords(nmol):
    """coordconstraint_mod.f90:initialise_coords (module coordconstraints_mod)."""
    # Initialise energy storage arrays
    coordconstraints_mod.ergold_coord = reallocation_routines.reallocate1d_r(coordconstraints_mod.ergold_coord, nmol)
    coordconstraints_mod.ergnew_coord = reallocation_routines.reallocate1d_r(coordconstraints_mod.ergnew_coord, nmol)


@routine
def en_mole_coord(im, naf, nal, xadd, yadd, zadd, xa, ya, za, mref, itype, qtype, npairref, amassr, tethered):
    """coordconstraint_mod.f90:en_mole_coord (module coordconstraints_mod) -- total energy
    of molecule im due all coordination constraints.

    Returns en_mole_coord.
    """
    #      write(6,*) 'en_mole_coord> ',ncoordconstraints,im
    en_mole_coord = f32(0.0)
    for ic in range(1, coordconstraints_mod.ncoordconstraints + 1):
        c = coordconstraints_mod.coordconstraints[ic]
        for ia in range(naf[im], nal[im] + 1):
            if c.nbondcc[ia] > 0:
                en, _ = en_coordcon(c, ia, xadd, yadd, zadd, xa, ya, za, mref, itype, qtype, npairref, amassr, tethered)
                #         write(6,*) 'en_mole_coord>',ic,im,ia,en
                en_mole_coord = en_mole_coord + en
    return en_mole_coord


@routine
def en_coord_total(natom, xadd, yadd, zadd, xa, ya, za, mref, itype, qtype, npairref, amassr, tethered):
    """coordconstraint_mod.f90:en_coord_total (module coordconstraints_mod) -- total energy
    of all atoms due all coordination constraints.

    Returns en_coord_total.
    """
    en_coord_total = f32(0.0)
    for ic in range(1, coordconstraints_mod.ncoordconstraints + 1):
        en_coord_total = en_coord_total + en_coordcon_total(
            coordconstraints_mod.coordconstraints[ic], natom, xadd, yadd, zadd, xa, ya, za, mref, itype,
            qtype, npairref, amassr, tethered)
    return en_coord_total


@routine
def en_change_mole_coord(im, naf, nal, xadd, yadd, zadd, xa, ya, za, mref, itype, qtype, npairref, amassr, tethered):
    """coordconstraint_mod.f90:en_change_mole_coord (module coordconstraints_mod) -- change
    in total energy of molecule im due all coordination constraints.

    Returns en_change_mole_coord.
    """
    en_change_mole_coord = f32(0.0)
    if coordconstraints_mod.ncoordconstraints > 0:
        coordconstraints_mod.ergnew_coord[im] = en_mole_coord(im, naf, nal, xadd, yadd, zadd, xa, ya, za, mref,
                                                              itype, qtype, npairref, amassr, tethered)
        en_change_mole_coord = coordconstraints_mod.ergnew_coord[im] - coordconstraints_mod.ergold_coord[im]
    return en_change_mole_coord


@routine
def save_coord_energies(im):
    """coordconstraint_mod.f90:save_coord_energies (module coordconstraints_mod).

    NB: molfit calls it only after accepted internal-rotation moves, not after
    whole-molecule moves (stale-slot family, FORTRAN_BUGS.md section 3) -- the caller's
    business; reproduced there.
    """
    if coordconstraints_mod.ncoordconstraints > 0:
        coordconstraints_mod.ergold_coord[im] = coordconstraints_mod.ergnew_coord[im]


@routine
def check_coord_bonded(ia, ja, mref, itype):
    """coordconstraint_mod.f90:check_coord_bonded (module coordconstraints_mod) -- .true.
    if atoms ia and ja are coordination-bonded in any constraint.

    Returns check_coord_bonded.
    """
    check_coord_bonded = False
    ic = 0
    while ic < coordconstraints_mod.ncoordconstraints and not check_coord_bonded:
        ic = ic + 1
        c = coordconstraints_mod.coordconstraints[ic]
        jb = 0
        nbond = c.nbondcc[ia]
        #         nbond=min(coordconstraints(ic)%nbondcc(ia),coordconstraints(ic)%expecnbondcc(ia))
        while jb < nbond and not check_coord_bonded:
            jb = jb + 1
            check_coord_bonded = bool(ja == c.iabondcc[jb, ia])
            #            check_coord_bonded=mref(ja).eq.mref(coordconstraints(ic)%iabondcc(jb,ia))
        # Now check the reverse if atom types are different
        if itype[ia] != itype[ja] and not check_coord_bonded:
            jb = 0
            nbond = c.nbondcc[ja]
            #            nbond=min(coordconstraints(ic)%nbondcc(ja),coordconstraints(ic)%expecnbondcc(ja))
            while jb < nbond and not check_coord_bonded:
                jb = jb + 1
                check_coord_bonded = bool(ia == c.iabondcc[jb, ja])
                #               check_coord_bonded=mref(ia).eq.mref(coordconstraints(ic)%iabondcc(jb,ja))
    return check_coord_bonded


@routine
def writecoordcons(natom, itype, atype):
    """coordconstraint_mod.f90:writecoordcons (module coordconstraints_mod) -- write the
    coordination bonds of every atom to coordcon.txt.
    """
    units.open(29, 'coordcon.txt', status='unknown')
    for ia in range(1, natom + 1):
        for ic in range(1, coordconstraints_mod.ncoordconstraints + 1):
            c = coordconstraints_mod.coordconstraints[ic]
            nbond = int(c.nbondcc[ia])
            nearnbond = c.nearnbondcc[ia]
            #            write(6,*) 'writecoordcon> ',ia,atype(itype(ia)),ic,nbond
            if nbond > 0:
                formatline = fstr(fmt('(a,i2.2,a)', '(i6,1x,a,1x,i2,1x,i2,', nbond, '(i6,1x,a,1x,e10.3))'), 256)
                items = []
                for j in range(1, nbond + 1):
                    items += [c.iabondcc[j, ia], Ch(atype[itype[c.iabondcc[j, ia]]], 3), c.tbondcc[j, ia]]
                    #                     ,coordconstraints(ic)%adistcc(j,ia),j=1,nbond)
                units.write(29, fmt(formatline, ia, Ch(atype[itype[ia]], 3), nearnbond, nbond, *items))
            else:
                units.write(29, fmt('(i6,1x,a,1x,i2)', ia, Ch(atype[itype[ia]], 3), nbond))
    units.close(29)
