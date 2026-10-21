"""Translation of fortran_src/clusters.f: ``program standalone`` (the ``clusters``
program), ``run_clusters`` and ``find_clusters``.

Cluster-size distribution of molecules: two molecules are bonded when a pair of their
atoms of the specified types is within a distance range; also the fraction of
molecules in percolating clusters.  ``run_clusters``/``find_clusters`` have namesakes
in local_density.f (programs/local_density.py); each program module keeps its own.
The library modules are imported inside the routines.  See docs/translation_notes/O.md.
"""

from ..common import clustercom
from ..common import clusters
from ..common import molecom
from ..common import symbols
from ..common.epsr_dimension import mmol
from ..fortran import f32
from ..fortran import findex_padded
from ..fortran import fmt
from ..fortran import fstr
from ..fortran import ifix
from ..fortran import int_array
from ..fortran import lw
from ..fortran import r4
from ..fortran import read_internal
from ..fortran import sqrtf
from ..fortran import substr
from ..fortran import units
from ..trace import routine


def main():
    """clusters.f:standalone -- program to perform the EPSR clusters command.

    2009/3/30 Attempt to make this a standalone program
    """
    from ..library import file_and_directory_commands

    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_clusters()


_F101 = "(/'Cluster ',i3,': fraction containing 2 or more molecules',f10.5)"


@routine
def run_clusters():
    """clusters.f:run_clusters -- programme to estimate number of clusters in a
    distribution of molecules. Atoms are defined to be in a cluster when they are within
    a specified distance range of a neighbour.

    2007-01-25 Modified and simplified to correct for a small counting error when
    estimating which molecules belong to which cluster.
    2007-01-26 Modified to estimate the percentage of clusters which are percolating,
    i.e. spanning the box in one or more dimensions.
    """
    from ..library import (
        ato_read_and_save,
        file_and_directory_commands,
        read_and_save_files,
        setup_symbols_generic,
    )

    p = clusters
    ireff = file_and_directory_commands.ireff
    # integer*4 ifound1,ifound2,ierrout
    # character*20 text
    #
    # setup the clusters symbols
    #
    setup_symbols_generic.setup_symbols_clusters()
    #
    # now setup the input file
    #
    nin = 11
    nout = 6
    ierrout = 0
    nout, ierrout = file_and_directory_commands.setup_input_file(nin, nout, ierrout)
    #
    # if an error was encountered on setup, exit now
    #
    if ierrout != 0:
        return
    #
    # read the input file
    #
    nout, ierrout = read_and_save_files.read_input_file(nin, nout, ierrout)
    #
    # if an error was encountered on read, exit now
    #
    if ierrout != 0:
        return
    #
    # set up the number of cluster size bins
    #
    iref = 2
    p.nz = read_internal(symbols.symbol_value[iref]).int(p.nz)
    if p.nz > p.mz:
        p.nz = p.mz
    #
    # set up current accumulator
    #
    iref = iref + 1
    p.nsumt = read_internal(symbols.symbol_value[iref]).int(p.nsumt)
    #
    # check the number of distributions does not exceed the array size
    #
    p.ndist = symbols.nfiles
    if p.ndist > p.mdist:
        p.ndist = p.mdist
    #
    # setup the type number of the first and second atoms used to define
    # whether two molecules are bonded, and the maximum and minimum distances
    #
    for id_ in range(1, p.ndist + 1):
        it = 0
        ifound1 = 0
        ifound2 = 0
        iref1 = ireff(symbols.nsymboluptodata + 1, id_)
        iref2 = iref1 + 1
        iref3 = iref2 + 1
        iref4 = iref3 + 1
        while it < molecom.ntype and (ifound1 == 0 or ifound2 == 0):
            it = it + 1
            if molecom.type[it] == symbols.symbol_value[iref1]:
                ifound1 = it
            if molecom.type[it] == symbols.symbol_value[iref2]:
                ifound2 = it
        p.ltyp[id_] = ifound1
        p.ltyp1[id_] = ifound2
        p.atom1[id_] = fstr(molecom.type[ifound1], 4)
        p.atom2[id_] = fstr(molecom.type[ifound2], 4)
        p.damin[id_] = read_internal(symbols.symbol_value[iref3]).real(p.damin[id_])
        p.damax[id_] = read_internal(symbols.symbol_value[iref4]).real(p.damax[id_])
        units.write(6, lw(id_, p.ltyp[id_], p.ltyp1[id_], p.damin[id_], p.damax[id_]))
    units.write(6, lw(p.nz, p.ndist))
    #
    # setup step size and z values
    #
    p.zstep = r4(0.5) / f32(p.nz)
    p.zhalf = r4(0.5) * p.zstep
    for j in range(1, p.nz + 1):
        p.z[j] = f32(j - 1)

    read_and_save_files.initgofz(p.nz, p.ndist, p.gofz, p.mz, p.mdist)

    #
    # normalise atomic coordinates to box size
    #
    ato_read_and_save.normalise_coordinates()
    #
    # set up the table of molecular bonds
    #
    find_clusters()
    #
    # axis title
    #
    p.titler = "Cluster_size"
    for id_ in range(1, p.ndist + 1):
        #
        # set the column title for this pair of atoms
        #
        ilent1 = findex_padded(p.atom1[id_], " ", 4) - 1
        ilent2 = findex_padded(p.atom2[id_], " ", 4) - 1
        p.titledat[id_] = fstr("   " + substr(p.atom1[id_], 1, ilent1) + "-" + substr(p.atom2[id_], 1, ilent2), 29)

    # output the data

    p.ext = "n"

    p.nsumt = read_and_save_files.outgofz(symbols.fnameinp, p.ext, p.nz, p.ndist, p.nsumt, p.z, p.gofz, p.gzt,
                                          p.gzt2, p.titler, p.titledat, p.mz, p.mdist)

    #
    # For each distribution calculate the fraction of molecules in clusters of 2 or more
    #
    for id_ in range(1, p.ndist + 1):
        sum1 = f32(0.0)
        sum2 = f32(0.0)
        for i in range(1, p.nz + 1):
            sum1 = sum1 + p.z[i] * p.gzt[i, id_]
            if p.z[i] > 1.0:
                sum2 = sum2 + p.z[i] * p.gzt[i, id_]
        p.fract2[id_] = sum2 / sum1

    #
    # write input file
    #

    # update the current value of nsumt and ndist

    icount = 0
    text = "nsize"
    icount = file_and_directory_commands.search_for_symbol(text, icount)
    iref = ireff(icount, 1)
    symbols.symbol_value[iref] = fstr(lw(p.nz), 256)
    iref = iref + 1
    symbols.symbol_value[iref] = fstr(lw(p.nsumt), 256)
    iref = iref + 1
    symbols.symbol_value[iref] = fstr(lw(p.ndist), 256)

    read_and_save_files.save_files_others(14)

    for id_ in range(1, p.ndist + 1):
        #
        # for each cluster distribution indicate the mean number of molecules per cluster
        #
        units.write(14, fmt(_F101, id_, p.fract2[id_]))

    units.close(14)

    return


_F995 = "(/'For cluster ',i2,' there are ',i5,' clusters of average size = ',e12.5,' molecules')"
_F994 = "(/'Clusters> For cluster ',i2,' fraction of molecules which percolate is ',f10.5)"


@routine
def find_clusters():
    """clusters.f:find_clusters -- bonded-neighbour lists, cluster labelling,
    percolation test and cluster-size histogram for each distribution."""
    from ..library import minimagec

    p = clusters
    m = molecom
    c = clustercom
    neighb = int_array(mmol, p.mdist)
    mbond = int_array(mmol)
    labcluster = int_array(mmol)
    lencluster = int_array(mmol)
    ncllab = int_array(mmol)
    nfdim = p.mdist  # noqa: F841  (not used)
    p.zstep = r4(0.5) / f32(p.nz)
    pi = r4(3.141592653)  # noqa: F841  (not used)
    nmol = c.nmol
    naf, nal = c.naf, c.nal
    xadd, yadd, zadd = c.xadd, c.yadd, c.zadd
    xa, ya, za = m.xa, m.ya, m.za
    itype_l = m.itype.tolist()  # plain-int copy of itype for the inner loops
    dmax = c.dmax
    gofz = p.gofz
    nz = p.nz
    mi = minimagec.minimagec
    #
    # look at each requested pair of atoms
    #
    for id_ in range(1, p.ndist + 1):
        units.write(6, lw(p.ltyp[id_], p.ltyp1[id_], p.damin[id_], p.damax[id_]))
        lt = int(p.ltyp[id_])
        lt1 = int(p.ltyp1[id_])
        damin_id = p.damin[id_]
        damax_id = p.damax[id_]
        #
        # set up for each molecule a list of bonded neighbours.
        #
        # mbond keeps track of how many molecules each molecule is bonded to. -1
        # indicates that no specified atom pairs were found for this molecule, so that
        # it is not included in the cluster analysis.
        #
        for i in range(1, nmol + 1):
            mbond[i] = -1
        #
        # step through all atoms
        #
        for im1 in range(1, nmol - 1 + 1):
            for ia1 in range(int(naf[im1]), int(nal[im1]) + 1):
                it1 = itype_l[ia1]
                #
                # now look at all atoms in other (distinct) molecules and determine which
                # ones are bonded
                #
                for im2 in range(im1 + 1, nmol + 1):
                    for ia2 in range(int(naf[im2]), int(nal[im2]) + 1):
                        it2 = itype_l[ia2]
                        #
                        # find out if this pair of atoms match the specified types
                        #
                        if (it1 == lt and it2 == lt1) or (it2 == lt and it1 == lt1):
                            #
                            # O.K. so we have found a pair of molecules with atoms of the correct labels
                            #
                            # zero the bond number accumulator if this has not already been set
                            #
                            if mbond[im1] < 0:
                                mbond[im1] = 0
                            if mbond[im2] < 0:
                                mbond[im2] = 0
                            #
                            # calculate atom1 to atom2 distance, applying periodic boundary conditions
                            #
                            xc1 = xadd[im1]
                            yc1 = yadd[im1]
                            zc1 = zadd[im1]
                            xa1 = xa[ia1] + xc1
                            ya1 = ya[ia1] + yc1
                            za1 = za[ia1] + zc1
                            xc2 = xadd[im2]
                            yc2 = yadd[im2]
                            zc2 = zadd[im2]
                            xa2 = xa[ia2] + xc2
                            ya2 = ya[ia2] + yc2
                            za2 = za[ia2] + zc2
                            xdif = xa2 - xa1
                            ydif = ya2 - ya1
                            zdif = za2 - za1
                            #
                            # apply periodic boundary conditions
                            #
                            xdif, ydif, zdif = mi(xdif, ydif, zdif)
                            dsq12 = xdif * xdif + ydif * ydif + zdif * zdif
                            da12 = sqrtf(dsq12) * dmax
                            #
                            # check that distance is within bounds prescribed for this distribution
                            #
                            if da12 >= damin_id and da12 <= damax_id:
                                #
                                # O.K. found a bond - save molecule numbers
                                #
                                mbond[im1] = mbond[im1] + 1
                                mbond[im2] = mbond[im2] + 1
                                neighb[im1, mbond[im1]] = im2
                                neighb[im2, mbond[im2]] = im1
                                #     write(6,*) im1,im2,mbond(im1),mbond(im2)
            #        write(6,*) im1,mbond(im1),(neighb(im1,im2),im2=1,mbond(im1))
        #
        # Generate initial cluster labels - ncllab(im1). Each valid molecule is put in its own
        # unique cluster. If a molecule is not valid (i.e. does not satisfy the type requirements
        # so that mbond = -1) it is given a fictitious cluster label nmol+1. Such a cluster label
        # can never exist because there cannot be more clusters than the number of molecules in
        # the box.

        ncluster = 0
        for im1 in range(1, nmol + 1):
            if mbond[im1] >= 0:
                ncluster = ncluster + 1
                ncllab[im1] = ncluster
                labcluster[ncluster] = ncluster
            else:
                ncllab[im1] = nmol + 1
        nclusterold = 0
        #
        # begin loop to look for clusters
        #
        while nclusterold != ncluster:

            # save the previous number of clusters

            nclusterold = ncluster

            # step through molecules looking at the bonds to other molecules

            for im1 in range(1, nmol + 1):

                # Ignore molecules of the wrong type

                if mbond[im1] >= 0:
                    nseek = int(ncllab[im1])
                    in_ = 0
                    while in_ < mbond[im1]:

                        # This molecule is bonded to at least one neighbour

                        # Find out the lowest cluster label (if any) of bonded neighbours

                        in_ = in_ + 1
                        im2 = neighb[im1, in_]
                        if ncllab[im2] < nseek:
                            nseek = int(ncllab[im2])

                    # Set all these bonded atoms into the same cluster

                    ncllab[im1] = nseek
                    in_ = 0
                    while in_ < mbond[im1]:
                        in_ = in_ + 1
                        im2 = neighb[im1, in_]
                        ncllab[im2] = nseek

            # Now count how many distinct clusters there are

            ncluster = 0

            for im1 in range(1, nmol + 1):
                if mbond[im1] >= 0:

                    # Run through the existing clusters and check whether the cluster this molecule belongs to
                    # already exists.

                    iflag = 0
                    ic = 0
                    while ic < ncluster and iflag == 0:

                        ic = ic + 1
                        if ncllab[im1] == labcluster[ic]:
                            iflag = ic

                    # If this cluster does not exist, we start a new one.

                    if iflag == 0:

                        ncluster = ncluster + 1
                        labcluster[ncluster] = ncllab[im1]

            units.write(6, lw(ncluster, nclusterold))

        # now work out the length of each cluster

        units.write(6, lw(ncluster))
        for ic in range(1, ncluster + 1):
            sum_ = f32(0)
            for im in range(1, nmol + 1):
                if ncllab[im] == labcluster[ic]:
                    sum_ = sum_ + r4(1.0)
            lencluster[ic] = ifix(sum_)
            #         write(6,*) ic,labcluster(ic),lencluster(ic)

        # calculate the average size of each cluster

        lensum = 0
        for ic in range(1, ncluster + 1):
            lensum = lensum + int(lencluster[ic])
        if ncluster > 0:
            avlen = f32(lensum) / f32(ncluster)
        else:
            avlen = f32(0.0)
        units.write(6, fmt(_F995, id_, ncluster, avlen))

        # Look to see which clusters are percolating, i.e. span the simulation box in at least
        # one dimension.

        npercolate = 0
        for ic in range(1, ncluster + 1):

            # First find a molecule that is in this cluster

            im1 = 0
            iflag = 0
            while im1 < nmol and iflag == 0:

                im1 = im1 + 1
                if ncllab[im1] == labcluster[ic]:
                    iflag = im1

            # Record the maximum separation of all the other molecules in this cluster from this
            # starting molecule, applying periodic boundary conditions

            xmin = f32(0.0)
            ymin = f32(0.0)
            zmin = f32(0.0)
            xmax = f32(0.0)
            ymax = f32(0.0)
            zmax = f32(0.0)
            for im2 in range(1, nmol + 1):
                if ncllab[im2] == labcluster[ic]:

                    xdif = xadd[im2] - xadd[im1]
                    ydif = yadd[im2] - yadd[im1]
                    zdif = zadd[im2] - zadd[im1]
                    xdif, ydif, zdif = mi(xdif, ydif, zdif)

                    # Store the largest and smallest separations

                    if xdif > xmax:
                        xmax = xdif
                    if xdif < xmin:
                        xmin = xdif
                    if ydif > ymax:
                        ymax = ydif
                    if ydif < ymin:
                        ymin = ydif
                    if zdif > zmax:
                        zmax = zdif
                    if zdif < zmin:
                        zmin = zdif

            #        write(6,*) ic,xmin,xmax,ymin,ymax,zmin,zmax

            # Apply periodic boundary conditions to the difference between maximum and minimum
            # coordinates

            xdif = xmax - xmin
            ydif = ymax - ymin
            zdif = zmax - zmin
            diftest = r4(1.0) - p.damax[id_] / dmax

            # Test for overlap. If any one of these differences is less than the maximum separation
            # for this cluster, it is regarded as percolating.

            #        write(6,*) ic,xdif,ydif,zdif,diftest

            if xdif >= diftest or ydif >= diftest or zdif >= diftest:

                # Cluster is percolating.

                npercolate = npercolate + int(lencluster[ic])

        # Fraction of molecules in percolating clusters:-

        fracpercolate = f32(npercolate) / f32(lensum)

        units.write(6, fmt(_F994, id_, fracpercolate))

        # form histogram of number of bonded neighbours

        # form histogram of cluster sizes

        for ic in range(1, ncluster + 1):
            iref = int(lencluster[ic]) + 1
            if iref <= nz:
                gofz[iref, id_] = gofz[iref, id_] + r4(1.0)
            else:
                gofz[nz, id_] = gofz[nz, id_] + r4(1.0)

        # normalise the distributions to the number of clusters

        for ic in range(2, nz + 1):
            gofz[ic, id_] = gofz[ic, id_] / f32(ncluster)

        # save the percolation fraction

        gofz[1, id_] = fracpercolate
        #
        # repeat for other distributions
        #
    return
