"""Translation of fortran_src/reallocation_routines.f90: ``module reallocation_routines``.

Use for reallocating memory.  Each routine reallocates an allocatable array to a new
shape, keeping the overlapping elements (allocate a save copy, copy, deallocate,
allocate, copy back).  The new elements are "left uninitialised" in the Fortran; here
they are zero / ``.false.`` / blank (``fortran.arrays.resize``).

Python cannot reallocate an array in place, so -- unlike the SIGNATURES.md entry, which
lists no return value -- every routine **returns the (re)allocated array**, ``None``
when the Fortran leaves it deallocated.  Call as::

    molecom_mod.labrot = reallocation_routines.reallocate1d_i(molecom_mod.labrot, n)

When the size is unchanged the same object is returned (the Fortran returns early).
``size(array, d)`` is ``array.shape[d] - 1`` (arrays from ``*_array`` helpers, lower
bound 1).
"""

from ..fortran import char_array, int_array, logical_array, real_array, resize
from ..trace import routine


def _oldsize(array, rank):
    """oldsize(d) = size(array,d) if allocated(array), else 0."""
    if array is not None:
        return [array.shape[d] - 1 for d in range(rank)]
    return [0] * rank


def _all_gt0(sizes):
    for s in sizes:
        if not s > 0:
            return False
    return True


def _realloc(array, oldsize, newsize, make, always_allocate=False):
    """The common body after the 'already the requested size' test.

    if(all oldsize > 0) save a copy; deallocate(array);
    if(all newsize > 0) then allocate(array(newsize)); copy back the overlap.
    ``always_allocate`` reproduces reallocate3d_i, which allocates without testing
    newsize (a non-positive extent gives a zero-sized, allocated array).
    """
    if _all_gt0(oldsize):
        arraysave = array
    else:
        arraysave = None
    array = None  # if(allocated(array)) deallocate(array)
    if always_allocate or _all_gt0(newsize):
        dims = [max(n, 0) for n in newsize]
        if arraysave is not None:
            array = resize(arraysave, *dims)
        else:
            array = make(*dims)
    # if(allocated(arraysave)) deallocate(arraysave)
    return array


def _char_maker(ncharsize):
    def make(*dims):
        return char_array(*dims)

    return make


@routine
def reallocate1d_i(array, newsize):
    """reallocation_routines.f90:reallocate1d_i -- integer array(:).

    reallocates an array, new values are left uninitialised
    Returns the reallocated array (see the module docstring).
    """
    oldsize = _oldsize(array, 1)
    # Do nothing if array is already the requested size!
    if oldsize[0] == newsize:
        return array
    return _realloc(array, oldsize, [newsize], int_array)


@routine
def reallocate2d_i(array, newsize1, newsize2):
    """reallocation_routines.f90:reallocate2d_i -- integer array(:,:).

    reallocates an array, new values are left uninitialised
    Returns the reallocated array.
    """
    newsize = [newsize1, newsize2]
    oldsize = _oldsize(array, 2)
    # Do nothing if array is already the requested size!
    if oldsize[0] == newsize[0] and oldsize[1] == newsize[1]:
        return array
    return _realloc(array, oldsize, newsize, int_array)


@routine
def reallocate3d_i(array, newsize1, newsize2, newsize3):
    """reallocation_routines.f90:reallocate3d_i -- integer array(:,:,:).

    reallocates an array, new values are left uninitialised
    NB: unlike the other routines this one allocates without testing newsize > 0.
    Returns the reallocated array.
    """
    newsize = [newsize1, newsize2, newsize3]
    oldsize = _oldsize(array, 3)
    # Do nothing if array is already the requested size!
    if oldsize[0] == newsize[0] and oldsize[1] == newsize[1] and oldsize[2] == newsize[2]:
        return array
    array = _realloc(array, oldsize, newsize, int_array, always_allocate=True)
    #        write(6,*) 'reallocate3d_i> ',size(array,1),size(array,2),size(array,3)
    return array


@routine
def reallocate4d_i(array, newsize1, newsize2, newsize3, newsize4):
    """reallocation_routines.f90:reallocate4d_i -- integer array(:,:,:,:).

    reallocates an array, new values are left uninitialised
    Returns the reallocated array.
    """
    newsize = [newsize1, newsize2, newsize3, newsize4]
    oldsize = _oldsize(array, 4)
    # Do nothing if array is already the requested size!
    if (oldsize[0] == newsize[0] and oldsize[1] == newsize[1] and oldsize[2] == newsize[2]
            and oldsize[3] == newsize[3]):
        return array
    return _realloc(array, oldsize, newsize, int_array)


@routine
def reallocate1d_r(array, newsize):
    """reallocation_routines.f90:reallocate1d_r -- real array(:).

    reallocates an array, new values are left uninitialised
    Returns the reallocated array.
    """
    oldsize = _oldsize(array, 1)
    # Do nothing if array is already the requested size!
    if oldsize[0] == newsize:
        return array
    return _realloc(array, oldsize, [newsize], real_array)


@routine
def reallocate2d_r(array, newsize1, newsize2):
    """reallocation_routines.f90:reallocate2d_r -- real array(:,:).

    reallocates an array, new values are left uninitialised
    Returns the reallocated array.
    """
    newsize = [newsize1, newsize2]
    oldsize = _oldsize(array, 2)
    # Do nothing if array is already the requested size!
    if oldsize[0] == newsize[0] and oldsize[1] == newsize[1]:
        return array
    return _realloc(array, oldsize, newsize, real_array)


@routine
def reallocate3d_r(array, newsize1, newsize2, newsize3):
    """reallocation_routines.f90:reallocate3d_r -- real array(:,:,:).

    reallocates an array, new values are left uninitialised
    Returns the reallocated array.
    """
    newsize = [newsize1, newsize2, newsize3]
    oldsize = _oldsize(array, 3)
    # Do nothing if array is already the requested size!
    if oldsize[0] == newsize[0] and oldsize[1] == newsize[1] and oldsize[2] == newsize[2]:
        return array
    return _realloc(array, oldsize, newsize, real_array)


@routine
def reallocate4d_r(array, newsize1, newsize2, newsize3, newsize4):
    """reallocation_routines.f90:reallocate4d_r -- real array(:,:,:,:).

    reallocates an array, new values are left uninitialised
    Returns the reallocated array.
    """
    newsize = [newsize1, newsize2, newsize3, newsize4]
    oldsize = _oldsize(array, 4)
    # Do nothing if array is already the requested size!
    if (oldsize[0] == newsize[0] and oldsize[1] == newsize[1] and oldsize[2] == newsize[2]
            and oldsize[3] == newsize[3]):
        return array
    return _realloc(array, oldsize, newsize, real_array)


@routine
def reallocate1d_c(array, ncharsize, newsize):
    """reallocation_routines.f90:reallocate1d_c -- character(len=ncharsize) array(:).

    reallocates an array, new values are left uninitialised (here: blank)
    The elements are unpadded strings already no longer than ncharsize, so ncharsize
    only documents the length.
    Returns the reallocated array.
    """
    oldsize = _oldsize(array, 1)
    # Do nothing if array is already the requested size!
    if oldsize[0] == newsize:
        return array
    return _realloc(array, oldsize, [newsize], _char_maker(ncharsize))


@routine
def reallocate2d_c(array, ncharsize, newsize1, newsize2):
    """reallocation_routines.f90:reallocate2d_c -- character(len=ncharsize) array(:,:).

    reallocates an array, new values are left uninitialised (here: blank)
    Returns the reallocated array.
    """
    newsize = [newsize1, newsize2]
    oldsize = _oldsize(array, 2)
    # Do nothing if array is already the requested size!
    if oldsize[0] == newsize[0] and oldsize[1] == newsize[1]:
        return array
    return _realloc(array, oldsize, newsize, _char_maker(ncharsize))


@routine
def reallocate1d_l(array, newsize):
    """reallocation_routines.f90:reallocate1d_l -- logical array(:).

    reallocates an array, new values are left uninitialised (here: .false.)
    Returns the reallocated array.
    """
    oldsize = _oldsize(array, 1)
    # Do nothing if array is already the requested size!
    if oldsize[0] == newsize:
        return array
    return _realloc(array, oldsize, [newsize], logical_array)
