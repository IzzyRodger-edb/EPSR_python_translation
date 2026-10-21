"""Translation of fortran_src/minimagec.f: minimum image convention.

The arguments must be real*4 (``numpy.float32``) values, as in the Fortran; they are
not converted here because these routines sit in the innermost loops.
"""

from ..common import lattice
from ..fortran import f32, nint
from ..trace import routine


@routine
def minimagec(x, y, z):
    """minimagec.f:minimagec -- minimum image of the Cartesian point (x,y,z).

    subroutine to apply the minimum image convention of a point in Cartesian
    coordinates (x,y,z) relative to an arbitrary set of lattice vectors
    a,b,c. Note that if specified as non-cubic the supplied coordinates do not have to
    be normalised, since the routine normalises them in any case.
    Returns (x, y, z).
    """
    if lattice.cubic:

        x, y, z = minimage(x, y, z)
        return x, y, z

    else:
        #
        # fractional coordinates relative to (a,b,c)
        #
        fc = z / lattice.clatz
        fb = (y - fc * lattice.claty) / lattice.blaty
        fa = (x - fb * lattice.blatx - fc * lattice.clatx) / lattice.alatx
        #
        # minimum image convention
        #
        fa, fb, fc = minimage(fa, fb, fc)
        #
        # calculate corrected coordinates
        #
        z = fc * lattice.clatz
        y = fb * lattice.blaty + fc * lattice.claty
        x = fa * lattice.alatx + fb * lattice.blatx + fc * lattice.clatx

    return x, y, z


@routine
def minimage(fa, fb, fc):
    """minimagec.f:minimage -- minimum image in box-fraction units (cubic boundaries).

    subroutine to apply the minimum image convention of a point in Cartesian
    coordinates assuming cubic boundary conditions
    Returns (fa, fb, fc).
    """
    #
    # minimum image convention
    #
    fc = fc - f32(nint(fc))
    fb = fb - f32(nint(fb))
    fa = fa - f32(nint(fa))
    return fa, fb, fc
