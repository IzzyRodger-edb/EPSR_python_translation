"""Translation of fortran_src/nwts.f: ``program standalone`` and ``subroutine nwts`` --
neutron weights for a DCS, 1st- or 2nd-order difference data set, driven by a
``.NWTS.dat`` keyword file (``setup_symbols_nwts``) and the .ato file it names.

The ``symbol_setup.eq.'nwts'`` branch of ``check_symbol`` is part of process_symbols.f.
"""

from ..common import clustercom, commands, directory_name, epsr_dimension, molecom, symbols
from ..fortran import (
    Ch,
    FortranEOF,
    FortranIOError,
    char_array,
    f32,
    fmt,
    fstr,
    int_array,
    iostat_of,
    lw,
    read_internal,
    real_array,
    substr,
    units,
)
from ..library import ato_manipulation_rev4
from ..library import ato_read_and_save
from ..library import file_and_directory_commands
from ..library import parse
from ..library import read_and_save_files
from ..trace import routine

misotope = 5  # parameter(misotope=5)

# FORMAT statements ('$' dropped: written with advance=False)
_F100 = ("(/'nwts> Program to calculate inter- and intra-molecular',' weightings ',"
         "/1x,'for DCS, 1st- or 2nd-order difference data.'/)")
_F202 = ("(/'nwts> The following components were found in this file'/1x,"
         "/1x,'Component no., label,  atomic fraction, chemical symbol'/1x,"
         "/((5x,i5,8x,a3,5x,e12.5,5x,a3)))")
_F219 = "(/'nwts> How many samples (1,2, or 3)? (0 to quit) ')"
_F235 = "(/'nwts> Number of components in input file different from number of atom types in .ato file')"
_F207 = "(/'nwts> For sample 1, get the scattering lengths for',' all components in the sample'/)"
_F223 = "(/'nwts> Get the scattering lengths for all',' components in the sample'/)"
_F209 = "(/'nwts> Did not recognise chemical symbol ',a3,' for atom type ',a3//'nwts> Try typing it again: ')"
_F212 = ("(/'nwts> Did not find massno ',i3,' for chemical',' symbol ',a3,' and abundance ',f10.5,"
         "//'Try typing again: ')")
_FISO = "(i3,1x,((a3,1x,f10.5,1x,f10.5)))"
_F108 = "(1x,a3,1x,a3)"
_F112 = "('tot',i1.1,'.wts')"
_F109 = "('tot.wts')"
_F113 = "(/'nwts> Writing TOTAL weights to file ',a)"
_F114 = "(a8)"
_F115 = "('dif.wts')"
_F116 = "(/'nwts> Writing 1st order difference weights to file ',a)"
_F117 = "('xx.wts')"
_F118 = "(1x,'Writing XX weights to file ',a)"
_F119 = "('xy.wts')"
_F120 = "(1x,'Writing XY (=XH) weights to file ',a)"
_F121 = "('yy.wts')"
_F122 = "(1x,'Writing YY (=HH) weights to file ',a)"


def main():
    """nwts.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    nwts()


@routine
def nwts():
    """nwts.f:nwts -- Program nwts.

    Calculates weights for a given input data (S(Q)) file or set of files, based on a
    specified .ATO file.  Currently set up for 70 components (number of atom types).
    """
    from ..library import setup_symbols_generic

    mcomp = epsr_dimension.mcomp
    mpairs = epsr_dimension.mpairs
    ncf = commands.ncf
    ncl = commands.ncl
    mol = molecom
    # character*256 wtsfile
    batom = real_array(mcomp, misotope, 2)
    abund = real_array(mcomp, misotope, 2)
    cprod = real_array(mpairs)
    bscatsum = real_array(2)
    binter = real_array(mcomp, 3)
    asum = f32(0.0)
    wintra = real_array(mpairs, 3)
    winter = real_array(mpairs, 3)
    wtsum = real_array(3)
    iflagsubs = int_array(mcomp)
    nisotope = int_array(mcomp, 2)
    # integer massno,normtot,normdif,ncomponent,ncharwts
    massno = normtot = normdif = ncomponent = 0
    # character*8 partlabel(mpairs),fnameb ; character*256 fname ; character*3 chsymtest
    partlabel = char_array(mpairs)
    # NB: these locals are only set by the reads below in the Fortran
    nrho = norder = 0
    fracmix = f32(0.0)
    ifind = 0
    ierr = 0

    # Get parameters from input file

    #
    # setup the nwts symbols
    #
    setup_symbols_generic.setup_symbols_nwts()
    #
    # now setup the input file
    #
    nin = 11
    nout = 6
    ierrout = 0  # NB: not initialised in the Fortran; setup_input_file sets it
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
    sv = symbols.symbol_value
    units.write(6, fmt(_F100))
    # 99    format(a80)
    #
    # Set up the base name for the weights file. This will be the same as the input file.

    iref = 1
    wtsfile = fstr(substr(symbols.fname_run, 1, symbols.ilenfname_run), 256)
    ncharwts = symbols.ilenfname_run
    #
    # set up whether to normalise to atoms or molecules.
    #
    iref = iref + 1
    nrho = read_internal(sv[iref]).int(nrho)
    #
    # get the no. of atoms per molecule ratio in case it is needed
    #
    atmolrat = f32(mol.natomr) / f32(clustercom.nmol)
    items = []
    for it in range(1, mol.ntype + 1):
        items += [it, Ch(mol.type[it], 3), mol.cfrac[it], Ch(mol.atsymbol[it], 3)]
    units.write(6, fmt(_F202, *items))

    # How many samples

    #      write(6,219)
    iref = iref + 1
    norder = read_internal(sv[iref]).int(norder)
    #
    # the program assumes
    #
    # norder = 1 means this is single dataset
    # norder = 2 means this is a first-order difference experiment
    # norder = 3 means this is a second-order difference experiment

    # Read fraction of sample 1 in mixture (only used if norder=3)

    iref = iref + 1
    fracmix = read_internal(sv[iref]).real(fracmix)

    # Normalisation on totals

    iref = iref + 1
    normtot = read_internal(sv[iref]).int(normtot)

    # Normalisation on first order differences

    iref = iref + 1
    normdif = read_internal(sv[iref]).int(normdif)

    # Check number of components in input file is same as atom types

    iref = iref + 1
    ncomponent = read_internal(sv[iref]).int(ncomponent)
    if ncomponent != mol.ntype:

        units.write(6, fmt(_F235))
        return

    #
    # set flags to say which (if any) atoms have been substituted
    #
    for it in range(1, mol.ntype + 1):
        iflagsubs[it] = 0
    nsubs = 0  # noqa: F841
    #
    # Get the scattering lengths for the first sample
    #
    if norder > 1:
        units.write(6, fmt(_F207))
    else:
        units.write(6, fmt(_F223))
    for it in range(1, mol.ntype + 1):

        # Skip over next two entries in input file as they will be determined from the ato file

        iref = iref + 2
        #
        # cannot see "Q" type atoms
        #
        if mol.qtypea[it] == 0:

            chsymtest = fstr(mol.atsymbol[it], 3)
            massno = 0
            #
            # test to see if this chemical symbol exists in the table
            #
            btest, ifind = ato_manipulation_rev4.scatt_len(chsymtest, massno, ifind)
            while ifind == 0:
                units.write(6, fmt(_F209, Ch(chsymtest, 3), Ch(mol.type[it], 3)), advance=False)
                try:
                    chsymtest = fstr(units.list_read(5).word(chsymtest), 3)
                    ierr = 0
                except (FortranEOF, FortranIOError) as e:
                    ierr = iostat_of(e)  # noqa: F841
                btest, ifind = ato_manipulation_rev4.scatt_len(chsymtest, massno, ifind)
            mol.atsymbol[it] = chsymtest
            #
            # for hydrogen atoms we need to check whether they exchange with other hydrogens
            #
            if mol.atsymbol[it] == "H":
                iref = iref + 1
                mol.iexchange[it] = read_internal(sv[iref]).int(mol.iexchange[it])
            else:
                iref = iref + 1
                mol.iexchange[it] = 0
            sv[iref] = fstr(fmt("(i1)", mol.iexchange[it]), 256)
            #
            # Read the atomic mass numbers and corresponding abundances
            #
            iref = iref + 1
            commands.sometext = fstr(sv[iref], commands.nstsize)
            commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 5, commands.nstsize,
                                          commands.nworddim)
            #      write(6,'(1x,i5,/((2(i5,1x))))') nwords,(ncf(j),ncl(j),j=1,nwords)
            #
            # the list of mass numbers and abundances must go in pairs  - any odd values at the end will be ignored
            #
            asum = f32(0)
            nsum = 0
            while 2 * nsum < commands.nwords:

                nsum = nsum + 1
                j2 = 2 * nsum
                j1 = j2 - 1
                massno = read_internal(substr(commands.sometext, ncf[j1], ncl[j1])).int(massno)
                abund[it, nsum, 1] = read_internal(substr(commands.sometext, ncf[j2], ncl[j2])).real(
                    abund[it, nsum, 1])
                asum = asum + abund[it, nsum, 1]
                if massno == 0:
                    batom[it, nsum, 1] = btest
                else:
                    batom[it, nsum, 1], ifind = ato_manipulation_rev4.scatt_len(chsymtest, massno, ifind)
                    while ifind == 0:
                        units.write(6, fmt(_F212, massno, Ch(mol.atsymbol[it], 3), abund[it, nsum, 1]),
                                    advance=False)
                        massno = units.list_read(5).int(massno)
                        batom[it, nsum, 1], ifind = ato_manipulation_rev4.scatt_len(chsymtest, massno, ifind)
            #
            # Renormalise abundances to total abundance
            #
            for isum in range(1, nsum + 1):
                abund[it, isum, 1] = abund[it, isum, 1] / asum
            nisotope[it, 1] = nsum
            items = [nisotope[it, 1]]
            for j in range(1, nisotope[it, 1] + 1):
                items += [Ch(mol.type[it], 3), batom[it, j, 1], abund[it, j, 1]]
            units.write(6, fmt(_FISO, *items))

            # Determine if this atom will be substituted - this value is ignored if there is only 1 sample.

            iref = iref + 1
            iflagsubs[it] = read_internal(sv[iref]).int(iflagsubs[it])

            # If atom is substituted, get the mass numbers and abundances of substituted atom

            if iflagsubs[it] == 1:

                # Get the abundances of the substituted atom

                iref = iref + 1
                commands.sometext = fstr(sv[iref], commands.nstsize)
                commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 5, commands.nstsize,
                                              commands.nworddim)
                #
                # the list of mass numbers and abundances must go in pairs  - any odd values at the end will be ignored
                #
                asum = f32(0)
                nsum = 0
                while 2 * nsum < commands.nwords:

                    nsum = nsum + 1
                    j2 = 2 * nsum
                    j1 = j2 - 1
                    massno = read_internal(substr(commands.sometext, ncf[j1], ncl[j1])).int(massno)
                    abund[it, nsum, 2] = read_internal(substr(commands.sometext, ncf[j2], ncl[j2])).real(
                        abund[it, nsum, 2])
                    asum = asum + abund[it, nsum, 2]
                    if massno == 0:
                        batom[it, nsum, 2] = btest
                    else:
                        batom[it, nsum, 2], ifind = ato_manipulation_rev4.scatt_len(chsymtest, massno, ifind)
                        while ifind == 0:
                            units.write(6, fmt(_F212, massno, Ch(mol.atsymbol[it], 3), abund[it, nsum, 2]),
                                        advance=False)
                            massno = units.list_read(5).int(massno)
                            batom[it, nsum, 2], ifind = ato_manipulation_rev4.scatt_len(chsymtest, massno, ifind)

                # Renormalise abundances to total abundance
                #
                for isum in range(1, nsum + 1):
                    abund[it, isum, 2] = abund[it, isum, 2] / asum
                nisotope[it, 2] = nsum

            else:

                # Otherwise set the mass numbers and abundances to the same as the first

                iref = iref + 1
                nisotope[it, 2] = nisotope[it, 1]
                for isum in range(1, nisotope[it, 2] + 1):

                    batom[it, isum, 2] = batom[it, isum, 1]
                    abund[it, isum, 2] = abund[it, isum, 1]

            items = [nisotope[it, 2]]
            for j in range(1, nisotope[it, 2] + 1):
                items += [Ch(mol.type[it], 3), batom[it, j, 2], abund[it, j, 2]]
            units.write(6, fmt(_FISO, *items))

        else:

            # Otherwise skip over the entries for this atom type

            iref = iref + 4

    # we now need to check that components which
    # are stated to be exchangeable through the sample do indeed all have the
    # the same scattering length.
    #
    # Step through the components and if iexchange is 1 form the average scattering
    # for all the atoms with the same chemical symbol.
    #
    for ii in range(1, 3):
        bscatsum[ii] = 0.0
    cfracsum = f32(0.0)
    for it in range(1, mol.ntype + 1):
        if mol.iexchange[it] == 1:
            cfracsum = cfracsum + mol.cfrac[it]
            for ii in range(1, nisotope[it, 1] + 1):
                term = mol.cfrac[it] * abund[it, ii, 1] * batom[it, ii, 1]
                bscatsum[1] = bscatsum[1] + term
            for ii in range(1, nisotope[it, 2] + 1):
                term = mol.cfrac[it] * abund[it, ii, 2] * batom[it, ii, 2]
                bscatsum[2] = bscatsum[2] + term
    if cfracsum > 0.0:
        for ii in range(1, 3):
            bscatsum[ii] = bscatsum[ii] / cfracsum
    else:
        for ii in range(1, 3):
            bscatsum[ii] = 0.0
    #
    # go through each atomic type and if iexchange is 1, then replace the individual
    # scattering lengths with a single average scattering length
    #
    for it in range(1, mol.ntype + 1):
        if mol.iexchange[it] == 1:
            nisotope[it, 1] = 1
            abund[it, 1, 1] = 1.0
            batom[it, 1, 1] = bscatsum[1]
            nisotope[it, 2] = 1
            abund[it, 1, 2] = 1.0
            batom[it, 1, 2] = bscatsum[2]
    #
    # if this is a 2nd-order difference experiment, then third sample is
    # treated as a mixture of the first two. If this is not the case then
    # you can either use this programme to set up the weights for each DCS file
    # separately, and so fit each total DCS, or else write another version!
    #
    if norder == 3:
        frac1 = fracmix
        frac2 = (f32(1.0) - fracmix)
    else:
        frac1 = f32(1.0)
        frac2 = f32(0.0)
    #
    # for each sample calculate both intra-molecular and inter-molecular weights
    # for the total differential cross-sections
    #
    # total number of unique pairs expected
    #
    npart = mol.ntype * (mol.ntype + 1) // 2
    #
    # calculate the weights for the first one or two samples
    #
    # weights for the third sample if present will be calculated
    # as an average of the first two samples
    #
    ndata = min(2, norder)
    #
    # go through all the components and calculate the isotope averaged scattering
    # length for datasets 1 and 2
    #
    for idata in range(1, ndata + 1):
        for it1 in range(1, mol.ntype + 1):
            #
            # average scattering length of this component
            #
            binter[it1, idata] = 0.0
            if mol.qtypea[it1] == 0:
                for ib in range(1, nisotope[it1, idata] + 1):
                    term1 = abund[it1, ib, idata] * batom[it1, ib, idata]
                    binter[it1, idata] = binter[it1, idata] + term1
        #
        # step through all component pairs and write down the intra- and inter-
        # molecular weights
        #
        # check the number of partials generated
        #
        ipart = 0
        #
        # sum of weights
        #
        wtsum[idata] = 0.0
        #
        # first components
        #
        for it1 in range(1, mol.ntype + 1):
            binter1 = binter[it1, idata]

            # Get the normalisation in the event of normalising to single atom scattering

            if normtot == 2 and mol.qtypea[it1] == 0:

                wtsum[idata] = wtsum[idata] + mol.cfrac[it1] * binter1 * binter1

            #
            # step through second components
            #
            for it2 in range(it1, mol.ntype + 1):
                ipart = ipart + 1
                if mol.qtypea[it1] == 0 and mol.qtypea[it2] == 0:
                    #
                    # if idata=1 write a label for this pair of components
                    #
                    if idata == 1:
                        partlabel[ipart] = fstr(fmt(_F108, Ch(mol.type[it1], 3), Ch(mol.type[it2], 3)), 8)
                        #
                        # atomic fraction product for this pair
                        #
                        if it1 == it2:
                            cprod[ipart] = mol.cfrac[it1] * mol.cfrac[it1]
                        else:
                            cprod[ipart] = f32(2.0) * mol.cfrac[it1] * mol.cfrac[it2]
                    #
                    # average scattering length of second component
                    #
                    binter2 = binter[it2, idata]
                    #
                    # scattering length product for intermolecular terms
                    #
                    winter[ipart, idata] = binter1 * binter2
                    #
                    # for intra-molecular terms involving H if both components are non-exchanging
                    # and no. of isotopes is the same for both components, need to sum
                    # the scattering length product for each isotope.
                    #
                    # Note that no checking has been done to see if this component pair
                    # occurs on the same molecule. Strictly if it does not the intramolecular
                    # weight should be set to zero. However EPSR will never find
                    # such pairs.
                    #
                    # !  (the per-isotope sum for non-exchanging H-H pairs is commented out)
                    wintra[ipart, idata] = winter[ipart, idata]
                    if normtot < 2:
                        term = cprod[ipart] * winter[ipart, idata]
                        wtsum[idata] = wtsum[idata] + term
    #
    # now do third sample if present
    #
    if norder == 3:
        wtsum[3] = 0
        for it1 in range(1, mol.ntype + 1):
            if mol.qtypea[it1] == 0:
                #
                # average scattering length for this component
                #
                binter[it1, 3] = frac1 * binter[it1, 1] + frac2 * binter[it1, 2]
                # !                 write(6,*) norder,it1,binter(it1,3),frac1,frac2
        ipart = 0
        for it1 in range(1, mol.ntype + 1):

            # Set up normalisation when normalising to single atom scattering

            if normtot == 2:

                # NB: binter(it1,3) of a Q atom is never set (uninitialised local); 0 here
                wtsum[3] = wtsum[3] + mol.cfrac[it1] * binter[it1, 3] * binter[it1, 3]

            for it2 in range(it1, mol.ntype + 1):
                ipart = ipart + 1
                if mol.qtypea[it1] == 0 and mol.qtypea[it2] == 0:
                    #
                    # for the intermolecular terms, we just have the random mixing result
                    #
                    winter[ipart, 3] = binter[it1, 3] * binter[it2, 3]
                    #
                    # for intramolecular terms, we have to mix the intra-molecular terms in proportion
                    # unless the atoms are exchangeable, in which case you get just the random mixing result
                    #
                    if mol.iexchange[it1] == 0 and mol.iexchange[it2] == 0:
                        wintra[ipart, 3] = frac1 * wintra[ipart, 1] + frac2 * wintra[ipart, 2]
                    else:
                        #
                        # otherwise we just have the random mixing result
                        #
                        wintra[ipart, 3] = winter[ipart, 3]
                    if normtot < 2:
                        term = cprod[ipart] * winter[ipart, 3]
                        wtsum[3] = wtsum[3] + term
    #
    # write out the total weights for the dataset(s)
    #
    fname = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                 + substr(wtsfile, 1, ncharwts), 256)
    #
    # find out how many characters there are in this name
    #
    nchar = directory_name.nchardir + ncharwts
    for idata in range(1, norder + 1):

        nstruct = normtot
        #
        # generate file name to write to
        #
        if norder > 1:
            fnameb = fstr(fmt(_F112, idata), 8)
            nchart = nchar + 8
        else:
            fnameb = fstr(fmt(_F109), 8)
            nchart = nchar + 7
        #
        # concatenate these two filenames
        #
        fname = fstr(substr(fname, 1, nchar) + fnameb, 256)
        units.write(6, fmt(_F113, substr(fname, 1, nchart)))
        units.open(10, fname, status="unknown")
        #
        # signal whether these weights are per atom or per molecule
        #
        if nrho == 2:
            if nstruct == 0:
                rhofac = atmolrat
            elif nstruct == 1:
                rhofac = f32(1.0) / atmolrat
            else:
                rhofac = f32(1.0)
        else:
            rhofac = f32(1.0)
        units.write(10, lw(f32(1.0)))
        for ipart in range(1, npart + 1):
            it1 = mol.n1extr[ipart]
            it2 = mol.n2extr[ipart]
            if mol.qtypea[it1] == 0 and mol.qtypea[it2] == 0:
                winterwr = winter[ipart, idata]
                wintrawr = wintra[ipart, idata]
                # !                  write(6,*) idata,ipart,wintra(ipart,idata),winter(ipart,idata)
                #
                # relative contribution to total
                #
                relativewt = winter[ipart, idata] * cprod[ipart]
                #
                # If required normalise the weights to the total sum of weights
                #
                if nstruct > 0:
                    winterwr = winterwr / wtsum[idata]
                    wintrawr = wintrawr / wtsum[idata]
                    relativewt = relativewt / wtsum[idata]
                units.write(10, fmt(_F114, Ch(partlabel[ipart], 8)))
                units.write(10, lw(wintrawr * rhofac, winterwr * rhofac, relativewt))
        units.close(10)
    #
    # now if norder .eq. 2, take the first order difference between samples 1 and 2
    # and write out the weights (but not if both intra and inter weights are zero
    #
    if norder == 2:

        nstruct = normdif
        #
        # generate file name to write to
        #
        fnameb = fstr(fmt(_F115), 8)
        #
        # concatenate these two filenames
        #
        fname = fstr(substr(fname, 1, nchar) + fnameb, 256)
        units.write(6, fmt(_F116, substr(fname, 1, nchar + 7)))
        units.open(10, fname, status="unknown")
        #
        # signal whether these weights are per atom or per molecule
        #
        if nrho == 2:
            if nstruct > 0:
                rhofac = f32(1.0) / atmolrat
            else:
                rhofac = atmolrat
        else:
            rhofac = f32(1.0)
        units.write(10, lw(f32(1.0)))
        wtsumdif = wtsum[1] - wtsum[2]
        for ipart in range(1, npart + 1):
            it1 = mol.n1extr[ipart]
            it2 = mol.n2extr[ipart]
            if mol.qtypea[it1] == 0 and mol.qtypea[it2] == 0:
                #
                # only write out those partials which involve a substituted atom
                #
                if iflagsubs[it1] == 1 or iflagsubs[it2] == 1:
                    winterwr = winter[ipart, 1] - winter[ipart, 2]
                    wintrawr = wintra[ipart, 1] - wintra[ipart, 2]
                    #
                    # relative contribution to total
                    #
                    relativewt = cprod[ipart] * winterwr
                    #
                    # If required normalise the weights to the total sum of weights
                    #
                    if nstruct == 1 and wtsumdif != 0.0:
                        winterwr = winterwr / wtsumdif
                        wintrawr = wintrawr / wtsumdif
                        relativewt = relativewt / wtsumdif
                    units.write(10, fmt(_F114, Ch(partlabel[ipart], 8)))
                    units.write(10, lw(wintrawr * rhofac, winterwr * rhofac, relativewt))
        units.close(10)
    #
    # Do second order difference if appropriate
    #
    if norder == 3:
        #
        # it will be assumed here that partial structure factors have
        # been normalised
        #
        nstruct = 1
        #
        # generate 3 filenames to write to
        #
        fnameb = fstr(fmt(_F117), 8)
        #
        # concatenate these two filenames
        #
        fname = fstr(substr(fname, 1, nchar) + fnameb, 256)
        units.write(6, lw())
        units.write(6, fmt(_F118, substr(fname, 1, nchar + 6)))
        units.open(10, fname, status="unknown")
        #
        # signal whether these weights are per atom or per molecule
        #
        if nrho == 2:
            rhofac = f32(1.0) / atmolrat
        else:
            rhofac = f32(1.0)
        units.write(10, lw(f32(1.0)))
        fnameb = fstr(fmt(_F119), 8)
        #
        # concatenate these two filenames
        #
        fname = fstr(substr(fname, 1, nchar) + fnameb, 256)
        units.write(6, lw())
        units.write(6, fmt(_F120, substr(fname, 1, nchar + 6)))
        units.open(11, fname, status="unknown")
        #
        # signal whether these weights are per atom or per molecule
        #
        units.write(11, lw(f32(1.0)))
        fnameb = fstr(fmt(_F121), 8)
        #
        # concatenate these two filenames
        #
        fname = fstr(substr(fname, 1, nchar) + fnameb, 256)
        units.write(6, lw())
        units.write(6, fmt(_F122, substr(fname, 1, nchar + 6)))
        units.open(12, fname, status="unknown")
        #
        # signal whether these weights are per atom or per molecule
        #
        units.write(12, lw(f32(1.0)))
        #
        # set up the total atomic fractions for X (non-substituted)
        # and Y (substituted) atoms and their respective average
        # scattering lengths
        #
        cx = f32(0)
        cy = f32(0)
        cxbx = f32(0)
        cyby1 = f32(0)
        cyby2 = f32(0)
        for it in range(1, mol.ntype + 1):
            if mol.qtypea[it] == 0:
                if iflagsubs[it] == 1:
                    cy = cy + mol.cfrac[it]
                    cyby1 = cyby1 + mol.cfrac[it] * binter[it, 1]
                    cyby2 = cyby2 + mol.cfrac[it] * binter[it, 2]
                else:
                    cx = cx + mol.cfrac[it]
                    cxbx = cxbx + mol.cfrac[it] * binter[it, 1]
        #
        # write down the sum of weights for each of the YY, XY
        # and XX terms
        #
        wtsumyy = frac1 * wtsum[1] + frac2 * wtsum[2] - wtsum[3]
        wtsumxy = f32(2.0) * cxbx * (cyby1 - cyby2)
        wtsumxx = cxbx * cxbx
        #
        # also write down the overall weighting of the YY term in the
        # first order difference sample 1 - sample 2
        #
        wtsumyydif = cyby1 * cyby1 - cyby2 * cyby2
        #
        # write down the weights of the YY and XY terms in the 1st sample
        #
        wtsumyy1 = cyby1 * cyby1
        wtsumxy1 = f32(2.0) * cxbx * cyby1
        #
        # now step through all the partials and write down their
        # weights. It is assumed the data have been normalised.
        #
        for ipart in range(1, npart + 1):
            #
            # get atom type numbers for this pair
            #
            it1 = mol.n1extr[ipart]
            it2 = mol.n2extr[ipart]
            if mol.qtypea[it1] == 0 and mol.qtypea[it2] == 0:
                #
                # YY terms
                #
                if iflagsubs[it1] == 1 and iflagsubs[it2] == 1:
                    winteryy = (frac1 * winter[ipart, 1] + frac2 * winter[ipart, 2] - winter[ipart, 3]) / wtsumyy
                    if mol.iexchange[it1] == 1 or mol.iexchange[it2] == 1:
                        wintrayy = (frac1 * wintra[ipart, 1] + frac2 * wintra[ipart, 2]
                                    - wintra[ipart, 3]) / wtsumyy
                    else:
                        wintrayy = f32(0.0)
                else:
                    winteryy = f32(0.0)
                    wintrayy = f32(0.0)
                relativeyy = cprod[ipart] * winteryy
                #
                # XY terms
                #
                if iflagsubs[it1] != iflagsubs[it2]:
                    winterxy = (winter[ipart, 1] - winter[ipart, 2] - wtsumyydif * winteryy) / wtsumxy
                    wintraxy = (wintra[ipart, 1] - wintra[ipart, 2] - wtsumyydif * wintrayy) / wtsumxy
                elif (iflagsubs[it1] == 1 and iflagsubs[it2] == 1
                      and mol.iexchange[it1] == 0 and mol.iexchange[it2] == 0):
                    winterxy = f32(0.0)
                    wintraxy = (wintra[ipart, 1] - wintra[ipart, 2] - wtsumyydif * wintrayy) / wtsumxy
                else:
                    winterxy = f32(0.0)
                    wintraxy = f32(0.0)
                relativexy = cprod[ipart] * winterxy
                #
                # XX terms
                #
                if iflagsubs[it1] == 0 and iflagsubs[it2] == 0:
                    winterxx = (winter[ipart, 1] - wtsumxy1 * winterxy - wtsumyy1 * winteryy) / wtsumxx
                    wintraxx = (wintra[ipart, 1] - wtsumxy1 * wintraxy - wtsumyy1 * wintrayy) / wtsumxx
                elif (iflagsubs[it1] == 1 and iflagsubs[it2] == 1
                      and mol.iexchange[it1] == 0 and mol.iexchange[it2] == 0):
                    winterxx = f32(0.0)
                    wintraxx = (wintra[ipart, 1] - wtsumxy1 * wintraxy - wtsumyy1 * wintrayy) / wtsumxx
                else:
                    winterxx = f32(0.0)
                    wintraxx = f32(0.0)
                relativexx = cprod[ipart] * winterxx
                #
                # ignore small contributions
                #
                if winteryy != 0.0 or wintrayy != 0.0:
                    units.write(12, fmt(_F114, Ch(partlabel[ipart], 8)))
                    units.write(12, lw(wintrayy * rhofac, winteryy * rhofac, relativeyy))
                if winterxy != 0.0 or wintraxy != 0.0:
                    units.write(11, fmt(_F114, Ch(partlabel[ipart], 8)))
                    units.write(11, lw(wintraxy * rhofac, winterxy * rhofac, relativexy))
                if winterxx != 0.0 or wintraxx != 0.0:
                    units.write(10, fmt(_F114, Ch(partlabel[ipart], 8)))
                    units.write(10, lw(wintraxx * rhofac, winterxx * rhofac, relativexx))
        units.close(10)
        units.close(11)
        units.close(12)
    #
    # write the .ato file in case anything was changed
    #
    ato_read_and_save.moloutc(symbols.fnameato)
    #
    # Finally rewrite the input file, in case any values changed
    #
    read_and_save_files.save_files_others(14)

    return
