"""Translation of fortran_src/local_density.f: ``program standalone`` (the
``local_density`` program), ``run_clusters`` and ``find_clusters``.

Despite its name this is the multiple-bond-definition variant of clusters.f (by RH):
each cluster distribution can list several (atom1, atom2, rmin, rmax) bond
definitions (``clusters_multibond.inc``).  It reads a ``.CLUSTERS.dat`` file
(setup_symbols_clusters).  ``run_clusters``/``find_clusters`` have namesakes in
clusters.f (programs/clusters.py); each program module keeps its own.
The library modules are imported inside the routines.  See docs/translation_notes/O.md.
"""

from ..common import clustercom
from ..common import clusters_multibond
from ..common import commands
from ..common import molecom
from ..common import symbols
from ..common.epsr_dimension import mmol
from ..fortran import Ch
from ..fortran import f32
from ..fortran import findex_padded
from ..fortran import fmt
from ..fortran import fstr
from ..fortran import ifix
from ..fortran import int_array
from ..fortran import lw
from ..fortran import r4
from ..fortran import read_internal
from ..fortran import setsub
from ..fortran import sqrtf
from ..fortran import substr
from ..fortran import units
from ..trace import routine


def main():
    """local_density.f:standalone -- program to perform the EPSR local_density command.

    2009/3/30 Attempt to make this a standalone program
    """
    from ..library import file_and_directory_commands

    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_clusters()


_F100 = "(a,'-',a,',')"
_F101 = "(/'Cluster ',i3,': fraction containing 2 or more molecules',f10.5)"


@routine
def run_clusters():
    """local_density.f:run_clusters -- programme to estimate number of clusters in a
    distribution of molecules. Atoms are defined to be in a cluster when they are within
    a specified distance range of a neighbour (several bond definitions per distribution).

    2007-01-25 Modified and simplified to correct for a small counting error when
    estimating which molecules belong to which cluster.
    2007-01-26 Modified to estimate the percentage of clusters which are percolating,
    i.e. spanning the box in one or more dimensions.
    """
    from ..library import (
        ato_read_and_save,
        file_and_directory_commands,
        parse,
        read_and_save_files,
        setup_symbols_generic,
    )

    p = clusters_multibond
    ireff = file_and_directory_commands.ireff
    # integer*4 ifound1,ifound2,ierrout,natom1s,natom2s,nrmins,nrmaxs
    # character*20 text
    # RH: text to allow processing of multiple bonds
    # character*3020 atom1stext,atom2stext,rminstext,rmaxstext
    natom1s = natom2s = nrmins = nrmaxs = 0
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

    # RH: Actually read the .dat parameters into symbols arrays
    units.write(6, lw("Ran setup_input_file ok."))

    nout, ierrout = read_and_save_files.read_input_file(nin, nout, ierrout)
    units.write(6, lw("Ran read_input_file ok."))
    #
    # if an error was encountered on read, exit now
    #
    if ierrout != 0:
        return

    # RH:
    # Using the variables from the input file (*CLUSTERS.dat), which
    # a restored in the char array symbol_value.
    #
    # Values read:
    #   nz -- the max cluster size
    #   nsumt -- number of accumulations
    # and the values for each cluster distribuition

    #
    # set up the number of cluster size bins
    #
    iref = 2
    p.nz = read_internal(symbols.symbol_value[iref]).int(p.nz)
    #
    # set up current accumulator
    #
    iref = iref + 1
    p.nsumt = read_internal(symbols.symbol_value[iref]).int(p.nsumt)
    #
    # check the number of distributions does not exceed the array size
    #

    # RH: nfiles is declared in symbols.inc, and set in
    # setup_symbols_clusters() in setup_symbols.f, and it is the number
    # of distributions
    p.ndist = symbols.nfiles
    if p.ndist > p.mdist:
        p.ndist = p.mdist
    #
    # setup the type number of the first and second atoms used to define
    # whether two molecules are bonded, and the maximum and minimum
    # distances
    #
    ncf, ncl = commands.ncf, commands.ncl
    for id_ in range(1, p.ndist + 1):
        # RH: nsymboluptodata is the iref value where the specification for
        # each cluster dist begins in symbol table
        iref1 = ireff(symbols.nsymboluptodata + 1, id_)  # atom1s
        iref2 = iref1 + 1  # atom2s
        iref3 = iref2 + 1  # rmins
        iref4 = iref3 + 1  # rmaxs
        # RH: Need to process the multiple atom pairs and their cutoffs:

        # RH: parse atom1s, atom2s, rmins, and rmaxs separately and store into arrays...
        atom1stext = fstr(symbols.symbol_value[iref1], 3020)
        atom2stext = fstr(symbols.symbol_value[iref2], 3020)
        rminstext = fstr(symbols.symbol_value[iref3], 3020)
        rmaxstext = fstr(symbols.symbol_value[iref4], 3020)

        natom1s = parse.parse(atom1stext, natom1s, ncf, ncl, 12, commands.nstsize, commands.nworddim)
        for iw in range(1, natom1s + 1):
            p.atom1s[id_, iw] = fstr(read_internal(substr(atom1stext, ncf[iw], ncl[iw])).word(p.atom1s[id_, iw]), 4)
            p.atom1s[id_, iw] = fstr(substr(atom1stext, ncf[iw], ncl[iw]), 4)

        natom2s = parse.parse(atom2stext, natom2s, ncf, ncl, 12, commands.nstsize, commands.nworddim)
        for iw in range(1, natom2s + 1):
            p.atom2s[id_, iw] = fstr(read_internal(substr(atom2stext, ncf[iw], ncl[iw])).word(p.atom2s[id_, iw]), 4)

        nrmins = parse.parse(rminstext, nrmins, ncf, ncl, 1, 256, commands.nworddim)
        for iw in range(1, nrmins + 1):
            p.rmins[id_, iw] = read_internal(substr(rminstext, ncf[iw], ncl[iw])).real(p.rmins[id_, iw])

        nrmaxs = parse.parse(rmaxstext, nrmaxs, ncf, ncl, 1, 256, commands.nworddim)
        for iw in range(1, nrmaxs + 1):
            p.rmaxs[id_, iw] = read_internal(substr(rmaxstext, ncf[iw], ncl[iw])).real(p.rmaxs[id_, iw])

        # Check that have same number of atom1s, atom2s, rmins, and rmaxs
        if (natom1s != natom2s) or (natom1s != nrmins) or (natom1s != nrmaxs):
            # The continued character constant is padded to column 72 by gfortran (two blanks
            # before "bonds"), as the Fortran binary prints.
            units.write(6, lw("Error: numbers of atoms and cutoffs in cluster  bonds definitions not equal."))
            return
        # RH: store the number of bond definitions for this cluster disribution
        p.numbdefs[id_] = natom1s

        # RH: Loop over cluster bonds for this distribution
        for ncbs in range(1, natom1s + 1):
            units.write(6, lw("cluster dist: ", id_))
            units.write(6, lw("bond def number: ", ncbs))
            units.write(6, lw("atom1 ", Ch(p.atom1s[id_, ncbs], 4)))
            units.write(6, lw("atom2 ", Ch(p.atom2s[id_, ncbs], 4)))
            units.write(6, lw("rmin ", p.rmins[id_, ncbs]))
            units.write(6, lw("rmax ", p.rmaxs[id_, ncbs]))

            # RH: Find atom types in cluster bond. Loop through the atom types
            # in type() array until find the atoms specified in .dat file
            it = 0
            ifound1 = 0
            ifound2 = 0
            while it < molecom.ntype and (ifound1 == 0 or ifound2 == 0):
                it = it + 1
                if molecom.type[it] == p.atom1s[id_, ncbs]:
                    ifound1 = it
                if molecom.type[it] == p.atom2s[id_, ncbs]:
                    ifound2 = it

            units.write(6, lw("atom1, atom2 types:", ifound1, ifound2))

            # RH: parse in multibond definitions into single bond definitions:
            p.ltyp[id_, ncbs] = ifound1
            p.ltyp1[id_, ncbs] = ifound2
            p.atom1[id_, ncbs] = fstr(molecom.type[ifound1], 4)
            p.atom2[id_, ncbs] = fstr(molecom.type[ifound2], 4)
            #            damin(id) = rmins(id,ncbs)
            #            damax(id) = rmaxs(id,ncbs)
            #            write(6,*) id, ncbs,ltyp(id),ltyp1(id),damin(id),damax(id)

        # end of loop over bond defs in this distribution
    # end of loop over distributions

    units.write(6, lw("nz,ndist", p.nz, p.ndist))
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
    # RH: For multi bond defs, the atom pairs are listed in the title
    # above the appropriate column
    p.titler = "Cluster_size"
    for id_ in range(1, p.ndist + 1):
        nchars = 0
        icharidx = 1
        p.titledat[id_] = ""  # zero char array
        for ib in range(1, p.numbdefs[id_] + 1):
            # set the column title for this cluster distribution
            ilent1 = findex_padded(p.atom1[id_, ib], " ", 4) - 1
            ilent2 = findex_padded(p.atom2[id_, ib], " ", 4) - 1
            nchars = ilent1 + ilent2 + 2 + 1
            p.titledat[id_] = setsub(p.titledat[id_], icharidx, icharidx + nchars,
                                     fmt(_F100, substr(p.atom1[id_, ib], 1, ilent1),
                                         substr(p.atom2[id_, ib], 1, ilent2)), 29)
            icharidx = icharidx + nchars
        p.titledat[id_] = setsub(p.titledat[id_], icharidx - 2, 29, " ", 29)  # rm last comma etc

    # output the data

    p.ext = "n"

    p.nsumt = read_and_save_files.outgofz(symbols.fnameinp, p.ext, p.nz, p.ndist, p.nsumt, p.z, p.gofz, p.gzt,
                                          p.gzt2, p.titler, p.titledat, p.mz, p.mdist)

    #
    # For each distribution calculate the fraction of molecules in
    # clusters of 2 or more
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
    text = "nsumt"
    icount = file_and_directory_commands.search_for_symbol(text, icount)
    iref = ireff(icount, 1)
    symbols.symbol_value[iref] = fstr(lw(p.nsumt), 256)
    iref = iref + 1
    symbols.symbol_value[iref] = fstr(lw(p.ndist), 256)

    read_and_save_files.save_files_others(14)

    for id_ in range(1, p.ndist + 1):
        #
        # for each cluster distribution indicate the mean number of
        # molecules per cluster
        #
        units.write(14, fmt(_F101, id_, p.fract2[id_]))

    units.close(14)

    return


_F995 = "(/'For cluster ',i2,' there are ',i5,' clusters of average size = ',e12.5,' molecules')"
_F994 = "(/'Clusters> For cluster ',i2,' fraction of molecules which percolate is ',f10.5)"


@routine
def find_clusters():
    """local_density.f:find_clusters -- bonded-neighbour lists (any of the bond
    definitions of the distribution, each molecule pair counted once), cluster
    labelling, percolation test and cluster-size histogram."""
    from ..library import minimagec

    p = clusters_multibond
    m = molecom
    c = clustercom
    neighb = int_array(mmol, p.mdist)
    mbond = int_array(mmol)
    labcluster = int_array(mmol)
    lencluster = int_array(mmol)
    ncllab = int_array(mmol)
    # integer*4 bonded
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

    # look at each requested pair of atoms
    for id_ in range(1, p.ndist + 1):
        nb = int(p.numbdefs[id_])
        lts = [0] + [int(p.ltyp[id_, ib]) for ib in range(1, nb + 1)]
        lt1s = [0] + [int(p.ltyp1[id_, ib]) for ib in range(1, nb + 1)]
        #
        # set up for each molecule a list of bonded neighbours.
        #
        # mbond keeps track of how many molecules each molecule is bonded
        # to. -1 indicates that no specified atom pairs were found for this
        # molecule, so that it is not included in the cluster analysis.
        #
        for i in range(1, nmol + 1):
            mbond[i] = -1

        # step through all atoms
        for im1 in range(1, nmol - 1 + 1):
            for ia1 in range(int(naf[im1]), int(nal[im1]) + 1):
                it1 = itype_l[ia1]
                # now look at all atoms in other (distinct) molecules and determine
                # which ones are bonded
                for im2 in range(im1 + 1, nmol + 1):
                    for ia2 in range(int(naf[im2]), int(nal[im2]) + 1):
                        it2 = itype_l[ia2]
                        #
                        # find out if this pair of atoms match the specified types
                        #

                        # RH: Need to change this bit to deal with multiple bond
                        # definitions:
                        for ib in range(1, nb + 1):
                            if (it1 == lts[ib] and it2 == lt1s[ib]) or (it2 == lts[ib] and it1 == lt1s[ib]):
                                # O.K. so we have found a pair of molecules with atoms of the
                                # correct labels
                                #
                                # zero the bond number accumulator if this has not already been set
                                if mbond[im1] < 0:
                                    mbond[im1] = 0
                                if mbond[im2] < 0:
                                    mbond[im2] = 0
                                # calculate atom1 to atom2 distance, applying periodic boundary
                                # conditions
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
                                if da12 >= p.rmins[id_, ib] and da12 <= p.rmaxs[id_, ib]:
                                    # O.K. found a bond - save molecule numbers

                                    # RH: Need to establish if these molecules are already bonded to
                                    # each other: loop through molecules already bonded to, and only add
                                    # to mbonds and neighb list if not already bonded:
                                    bonded = -1
                                    for im1i in range(1, int(mbond[im1]) + 1):
                                        if neighb[im1, im1i] == im2:
                                            bonded = 1
                                    if bonded == -1:
                                        mbond[im1] = mbond[im1] + 1
                                        mbond[im2] = mbond[im2] + 1
                                        neighb[im1, mbond[im1]] = im2
                                        neighb[im2, mbond[im2]] = im1
                                    #  write(6,*) im1,im2,mbond(im1),mbond(im2)
                        # e of loop thru bond definitions
                    # e of loop thru other atoms
                # e of loop thru other molecs
            # e of loop thru molec's atoms
        # e of loop thru molecs

        # RH: Have now worked out which molecules are bonded to which, now
        # need to calculate the clusters

        #
        # Generate initial cluster labels - ncllab(im1). Each valid molecule
        # is put in its own unique cluster. If a molecule is not valid
        # (i.e. does not satisfy the type requirements so that mbond = -1)
        # it is given a fictitious cluster label nmol+1. Such a cluster
        # label can never exist because there cannot be more clusters than
        # the number of molecules in the box.

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
                    # Run through the existing clusters and check whether the cluster
                    # this molecule belongs to already exists.
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
            # write(6,*) ic,labcluster(ic),lencluster(ic)

        # calculate the average size of each cluster

        lensum = 0
        for ic in range(1, ncluster + 1):
            lensum = lensum + int(lencluster[ic])
        if ncluster > 0:
            avlen = f32(lensum) / f32(ncluster)
        else:
            avlen = f32(0.0)
        units.write(6, fmt(_F995, id_, ncluster, avlen))

        # RH: When I want to calc RoG etc this is where I do it

        # Look to see which clusters are percolating, i.e. span the
        # simulation box in at least one dimension.

        npercolate = 0
        for ic in range(1, ncluster + 1):

            # First find a molecule that is in this cluster

            im1 = 0
            iflag = 0
            while im1 < nmol and iflag == 0:

                im1 = im1 + 1
                if ncllab[im1] == labcluster[ic]:
                    iflag = im1

            # Record the maximum separation of all the other molecules in this
            # cluster from this starting molecule, applying periodic boundary
            # conditions

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

            # write(6,*) ic,xmin,xmax,ymin,ymax,zmin,zmax

            # Apply periodic boundary conditions to the difference between
            # maximum and minimum coordinates

            xdif = xmax - xmin
            ydif = ymax - ymin
            zdif = zmax - zmin
            # NB: damax(id) of clusters_multibond.inc is never set in local_density: it is 0
            # (COMMON storage is zero-initialised), so diftest = 1.0 (Rust local_density.rs).
            diftest = r4(1.0) - p.damax[id_] / dmax

            # Test for overlap. If any one of these differences is less than the
            # maximum separation for this cluster, it is regarded as
            # percolating.

            # write(6,*) ic,xdif,ydif,zdif,diftest

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
        units.write(6, lw("number of clusters", ncluster))

        # save the percolation fraction

        gofz[1, id_] = fracpercolate
        #
        # repeat for other distributions
        #
    # e of loop thru cluster distributions
    return
