"""Translation of fortran_src/ergpairreference.f: modules ergpairstore and ergpairreference.

Module data: ``common/ergpairstore.py`` (imported as ``eps``) and
``common/ergpairreference.py`` (imported as ``erd``).  Both Fortran modules contain a
routine called ``save_energies``; the ergpairreference one (the one every program uses)
keeps the name ``save_energies``, the ergpairstore one is ``ergpairstore_save_energies``.
"""

from ..common import ergpairreference as erd  # module ergpairreference (data)
from ..common import ergpairstore as eps  # module ergpairstore (data)
from ..fortran import idiv, int_array, lw, real_array, units
from ..trace import routine

# ----------------------------------------------------------------------- module ergpairstore


@routine
def initialise_ergpairstore(nnbrefmax):
    """ergpairreference.f:initialise_ergpairstore (module ergpairstore) -- allocate the store.

    Returns nothing.
    """
    eps.nnbref = 0
    if eps.ergold is not None:
        eps.ergold = None
        eps.ergnew = None
        eps.nbref = None
    eps.ergold = real_array(nnbrefmax)
    eps.ergnew = real_array(nnbrefmax)
    eps.nbref = int_array(nnbrefmax)


@routine
def set_ergnew():
    """ergpairreference.f:set_ergnew (module ergpairstore) -- ergnew(1:nnbref) = ergold.

    Returns nothing.
    """
    for ibref in range(1, eps.nnbref + 1):
        eps.ergnew[ibref] = eps.ergold[ibref]
    return


@routine
def ergpairstore_save_energies():
    """ergpairreference.f:save_energies (module ergpairstore) -- ergold(1:nnbref) = ergnew.

    Renamed from ``save_energies`` (the name is also used by module ergpairreference).
    Returns nothing.
    """
    for ibref in range(1, eps.nnbref + 1):
        eps.ergold[ibref] = eps.ergnew[ibref]
    return


# end module ergpairstore

# ------------------------------------------------------------------ module ergpairreference
# This version uses a single dimension array to store pair energies


@routine
def initialise_ergpair(nmol):
    """ergpairreference.f:initialise_ergpair (module ergpairreference) -- allocate pair store.

    Allocates the arrays to store pair energies and their reference integers, sized
    from the cell neighbour list (get_nnbrefmax).  Returns nothing.
    """
    from . import cellneighbourlist  # use cellneighbourlist

    erd.nmolhere = nmol
    erd.nnbrefmax = cellneighbourlist.get_nnbrefmax(erd.nnbrefmax)
    # Allocate the arrays to store pair energies and reference integers
    if erd.ergpair is not None:
        erd.ergpair = None
    if erd.molergpair is not None:
        erd.molergpair = None
    if erd.ergpairref is not None:
        erd.ergpairref = None
    if erd.nergpairm is not None:
        erd.nergpairm = None
    erd.nergpair = 0
    erd.nergpairmax = 100 * nmol
    erd.ergpair = real_array(erd.nergpairmax)
    erd.molergpair = int_array(erd.nnbrefmax, nmol)
    erd.ergpairref = int_array(erd.nnbrefmax, nmol)
    erd.nergpairm = int_array(nmol)
    for im in range(1, nmol + 1):
        erd.nergpairm[im] = 0
    erd.nnbref = 0
    if erd.ergold is not None:
        erd.ergold = None
        erd.ergnew = None
        erd.nbref = None
    erd.ergold = real_array(erd.nnbrefmax)
    erd.ergnew = real_array(erd.nnbrefmax)
    erd.nbref = int_array(erd.nnbrefmax)
    units.write(6, lw("initialise_ergpair> nnbrefmax = ", erd.nnbrefmax, erd.ergnew.shape[0] - 1))
    return


@routine
def set_ergpair(im, jm, energy):
    """ergpairreference.f:set_ergpair (module ergpairreference) -- store energy of pair im,jm.

    Returns nothing.
    """
    if erd.nergpairm[im] < erd.nnbrefmax and erd.nergpairm[jm] < erd.nnbrefmax and energy != 0.0:
        # Store the energy
        erd.nergpair = erd.nergpair + 1
        if erd.nergpair > erd.nergpairmax:
            reallocate_ergpair()
        erd.ergpair[erd.nergpair] = energy
        # Update the energy reference arrays
        if jm == im:
            erd.nergpairm[im] = 1
            erd.molergpair[erd.nergpairm[im], im] = im
            erd.ergpairref[erd.nergpairm[im], im] = erd.nergpair
        else:
            erd.nergpairm[im] = erd.nergpairm[im] + 1
            erd.nergpairm[jm] = erd.nergpairm[jm] + 1
            erd.molergpair[erd.nergpairm[im], im] = jm
            erd.molergpair[erd.nergpairm[jm], jm] = im
            erd.ergpairref[erd.nergpairm[im], im] = erd.nergpair
            erd.ergpairref[erd.nergpairm[jm], jm] = erd.nergpair
    return


@routine
def reallocate_ergpair():
    """ergpairreference.f:reallocate_ergpair (module ergpairreference) -- grow ergpair by 20%.

    Returns nothing.
    """
    nergpairmaxsave = erd.nergpairmax
    #      write(6,*) 'reallocate_ergpair> ',nergpairmax
    ergpairsave = real_array(nergpairmaxsave)
    for i in range(1, nergpairmaxsave + 1):
        ergpairsave[i] = erd.ergpair[i]
    erd.ergpair = None
    erd.nergpairmax = idiv((erd.nergpairmax + 10) * 12, 10)
    #      write(6,*) 'reallocate_ergpair> ',nergpairmax
    erd.ergpair = real_array(erd.nergpairmax)
    for i in range(1, nergpairmaxsave + 1):
        erd.ergpair[i] = ergpairsave[i]
    ergpairsave = None
    return


@routine
def get_neighbsum(neighbsum):
    """ergpairreference.f:get_neighbsum (module ergpairreference) -- mean pairs per molecule.

    Returns neighbsum.
    """
    neighbsum = 0
    nmol = erd.nmolhere
    for im in range(1, nmol + 1):
        neighbsum = neighbsum + int(erd.nergpairm[im])
    if nmol > 0:
        neighbsum = idiv(neighbsum, nmol)
    return neighbsum


@routine
def make_neighbour_list_from_mols(im):
    """ergpairreference.f:make_neighbour_list_from_mols (module ergpairreference).

    Loads the stored pair energies of molecule im into nbref/ergold/ergnew.
    Returns nothing.
    """
    erd.nnbref = int(erd.nergpairm[im])
    if erd.ergold is None:
        erd.ergold = real_array(erd.nnbref)
        erd.ergnew = real_array(erd.nnbref)
        erd.nbref = int_array(erd.nnbref)
    elif erd.ergold.shape[0] - 1 < erd.nnbref:
        erd.ergold = None
        erd.ergnew = None
        erd.nbref = None
        erd.ergold = real_array(erd.nnbref)
        erd.ergnew = real_array(erd.nnbref)
        erd.nbref = int_array(erd.nnbref)
    for ibref in range(1, erd.nnbref + 1):
        erd.nbref[ibref] = erd.molergpair[ibref, im]
        erd.ergold[ibref] = erd.ergpair[erd.ergpairref[ibref, im]]
        erd.ergnew[ibref] = erd.ergold[ibref]
    return


@routine
def save_energies():
    """ergpairreference.f:save_energies (module ergpairreference) -- ergold(1:nnbref) = ergnew.

    Returns nothing.
    """
    for ibref in range(1, erd.nnbref + 1):
        erd.ergold[ibref] = erd.ergnew[ibref]
    return


@routine
def save_erg_pairs(im):
    """ergpairreference.f:save_erg_pairs (module ergpairreference) -- store ergold back.

    Returns nothing.
    """
    for ibref in range(1, erd.nnbref + 1):
        erd.ergpair[erd.ergpairref[ibref, im]] = erd.ergold[ibref]
    return


# end module ergpairreference
