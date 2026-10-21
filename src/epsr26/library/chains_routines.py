"""Translation of fortran_src/chains_routines.f: ``match`` and ``shortest_path``, shared
by the chains, rings and lsi programs."""

from ..fortran import fstr, int_array, substr
from ..trace import routine


@routine
def match(string1, string2):
    """chains_routines.f:match -- determines whether string1 is the same as string2 up any
    * that may be present in string2. If no * is present, the whole strings are matched.

    string1 is character*3, string2 character*4.  Returns match.
    """
    # (called millions of times by find_chains/find_rings/find_lsi: the character*3/*4
    # truncation is only applied when needed, and index() is str.find)
    if len(string1) > 3:
        string1 = fstr(string1, 3)
    if len(string2) > 4:
        string2 = fstr(string2, 4)

    nstar = string2.find("*") + 1  # index(string2,'*')

    if nstar > 0:
        nsearch = nstar - 1
        if nsearch > 0:
            match = substr(string1, 1, nsearch) == substr(string2, 1, nsearch)
        else:
            match = True
        return match

    else:
        match = string1 == string2
        return match


@routine
def shortest_path(im1, nz, nmol, nneighb, neighb, rneighb1, rneighb2, spvaluet, nspvaluet, mmol, mbond):
    """chains_routines.f:shortest_path -- calculates the shortest (topological) path to im1
    from all other molecules in the box, out to a 'diameter' of nz.

    nspvaluet(nmol,nz) receives the number of times each chain length occurs.
    Arrays are filled in place; nothing is returned.
    """
    # integer, dimension(:), allocatable :: sourceatom1,sourceatom2
    sourceatom1 = int_array(nmol)
    sourceatom2 = int_array(nmol)

    # Initialise the shortest topological path to other molecules larger than the number of
    # molecules in the box.

    #      do im=1,nmol
    #         spvaluet(im)=nmol+1
    #         do icount=1,nz
    #            nspvaluet(im,icount)=0
    # (filling with constants: done with slices)
    spvaluet[1:nmol + 1] = nmol + 1
    nspvaluet[1:nmol + 1, 1:nz + 1] = 0
    # Set the shortest topological distances to other molecules in the vicinity of im1
    nneighb1 = 1
    rneighb1[nneighb1] = im1
    sourceatom1[nneighb1] = im1
    # Distance to molecule at origin is zero
    spvaluet[im1] = 0
    icount = 0
    ineighb1 = 0
    # limit the maximum value of icount needed
    while icount < nz:
        icount = icount + 1
        # form a list of the neighbours of all the molecules listed in rneighb1
        nneighb2 = 0
        ineighb1save = ineighb1  # noqa: F841  Need to save the current molecule in list to ensure we do not refer back to it
        while ineighb1 < nneighb1:
            ineighb1 = ineighb1 + 1
            in_ = 0
            while in_ < nneighb[rneighb1[ineighb1]]:
                in_ = in_ + 1
                im2 = int(neighb[rneighb1[ineighb1], in_])
                # Need to ensure im2 does not refer back to the molecule that precedes rneighb1(ineighb1) in the chain
                im2precede = sourceatom1[ineighb1]
                if im2 != im2precede:
                    #               write(6,*) icount,im1,im2,sourceatom1(ineighb1)
                    # set the value of shortest path if it is smaller than the current value
                    if icount < spvaluet[im2]:
                        spvaluet[im2] = icount
                    # Count number of times this chain length occurs.
                    nspvaluet[im2, icount] = nspvaluet[im2, icount] + 1
                    # Check whether molecule im2 is already in the list
                    found = False
                    intest = 0
                    while not found and intest < nneighb1:
                        intest = intest + 1
                        found = im2 == rneighb1[intest]
                    # Also check this molecule does not occur in the new list either
                    intest = 0
                    while not found and intest < nneighb2:
                        intest = intest + 1
                        found = im2 == rneighb2[intest]
                    # If not found, set the new list
                    if not found and nneighb2 < mmol:
                        nneighb2 = nneighb2 + 1
                        rneighb2[nneighb2] = im2
                        sourceatom2[nneighb2] = rneighb1[ineighb1]  # Save the source atom for the current set of neighbours
        # Save the new values
        in_ = 0
        while in_ < nneighb2 and nneighb1 < mmol:
            in_ = in_ + 1
            nneighb1 = nneighb1 + 1
            rneighb1[nneighb1] = rneighb2[in_]
            sourceatom1[nneighb1] = sourceatom2[in_]
        #         write(6,*)
        #         write(6,*) 'shortest_path> ',im1,icount,nneighb1,(rneighb1(in),spvaluet(rneighb1(in)),in=1,nneighb1)
    return
