"""Translation of fortran_src/epsrwts.f: ``program standalone`` and ``subroutine epsrwts``
-- neutron weights for a DCS, 1st- or 2nd-order difference data set, from a .ato file
and answers on standard input.
"""

from ..common import clustercom, commands, directory_name, epsr_dimension, molecom, symbols
from ..fortran import (
    Ch,
    FortranEOF,
    FortranIOError,
    char_array,
    f32,
    findex,
    fmt,
    fstr,
    int_array,
    iostat_of,
    lw,
    real_array,
    substr,
    units,
)
from ..library import ato_manipulation_rev4
from ..library import ato_read_and_save
from ..library import file_and_directory_commands
from ..library import parse
from ..trace import routine

misotope = 5  # parameter(misotope=5)

# FORMAT statements ('$' dropped: written with advance=False)
_F499 = "(/'epsrwts> No .ato file specified! - routine will exit'/)"
_F599 = "(/'epsrwts> Problems with .ato file ',a/'Routine will exit now'/)"
_F100 = ("(/'epsrwts> Program to calculate inter- and intra-molecular',' weightings ',"
         "/1x,'for DCS, 1st- or 2nd-order difference data.'/)")
_F200 = "(/'epsrwts> Is the output to be per atom (1), or per molecule (2)? ')"
_F202 = ("(/'epsrwts> The following components were found in this file'/1x,"
         "/1x,'Component no., label,  atomic fraction, chemical symbol'/1x,"
         "/((5x,i5,8x,a3,5x,e12.5,5x,a3)))")
_F219 = "(/'epsrwts> How many samples (1,2, or 3)? (0 to quit) ')"
_F207 = "(/'epsrwts> For sample 1, get the scattering lengths for',' all components in the sample'/)"
_F223 = "(/'epsrwts> Get the scattering lengths for all',' components in the sample'/)"
_F209 = ("(/'epsrwts> Did not recognise chemical symbol ',a3,' for atom type ',a3"
         "//'epsrwts> Try typing it again: ')")
_F103 = ("(/'epsrwts> For hydrogen component ',a3,':-'//' If it exchanges with atoms on other molecules type 1,'"
         ",' if not type 0: ')")
_F211 = "(/'epsrwts> For component ',a3//' Type 0 for a natural isotope or mass number',' for a specific isotope: ')"
_F224 = "(/' and its abundance (0.0-1.0): ')"
_F212 = "(/'epsrwts> Did not find this massno for chemical',' symbol ',a3,' - try typing again'/)"
_F213 = "(/'epsrwts> For sample 2:-'//' Give the component number of ','a substituted atom, (0 to exit) ')"
_F155 = "(/'epsrwts> Cannot substitute on a Q atom!')"
_F214 = ("(/'epsrwts> For sample 3, specify the mixing fraction'//' (This should be the same as that used',"
         "' in the .hdb_dat file): ')")
_F108 = "(1x,a3,1x,a3)"
_F220 = "(/'epsrwts> Type basename of file to output weights to: ')"
_F111 = "(a)"
_F201 = "(/'epsrwts> For total data ',i2,' has data been normalised',' (1) or not (0)? ')"
_F112 = "('tot',i1.1,'.wts')"
_F109 = "('tot.wts')"
_F113 = "(/'epsrwts> Writing TOTAL weights to file ',a)"
_F114 = "(a8)"
_F301 = "(/'epsrwts> Has first order difference data data been normalised',' (1) or not (0)? ')"
_F115 = "('dif.wts')"
_F116 = "(/'epsrwts> Writing 1st order difference weights to file ',a)"
_F117 = "('xx.wts')"
_F118 = "(1x,'Writing XX weights to file ',a)"
_F119 = "('xy.wts')"
_F120 = "(1x,'Writing XY (=XH) weights to file ',a)"
_F121 = "('yy.wts')"
_F122 = "(1x,'Writing YY (=HH) weights to file ',a)"


def main():
    """epsrwts.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    epsrwts()


def _read5_err(kind, old):
    """``read(5,*,err=label) x``: (True, value), or (False, old) on a conversion error
    (the caller jumps to the label); end of file is not handled (fatal, as in the Fortran)."""
    r = units.list_read(5)
    try:
        if kind == "int":
            return True, r.int(old)
        return True, r.real(old)
    except FortranIOError:
        return False, old


@routine
def epsrwts():
    """epsrwts.f:epsrwts -- Program EPSRWTS.

    Calculates weights for a given input data (S(Q)) file or set of files, based on a
    specified .ATO file.  Currently set up for 70 components (number of atom types).
    """
    mcomp = epsr_dimension.mcomp
    mpairs = epsr_dimension.mpairs
    ncf = commands.ncf
    ncl = commands.ncl
    mol = molecom
    # character*256 readfile ; character*4 ext ; character*1 iaccept
    iaccept = ""  # NB: not initialised in the Fortran; search_for_files sets it
    batom = real_array(mcomp, misotope, 2)
    abund = real_array(mcomp, misotope, 2)
    cprod = real_array(mpairs)
    binter = real_array(mcomp, 3)
    wintra = real_array(mpairs, 3)
    winter = real_array(mpairs, 3)
    wtsum = real_array(3)
    iflagsubs = int_array(mcomp)
    nisotope = int_array(mcomp, 2)
    # integer*4 massno ; character*8 partlabel(mpairs),fnameb ; character*80 fname
    partlabel = char_array(mpairs)
    # character*3 chsymtest
    # implicitly typed locals passed to minitc (epsr16.inc is not included)
    rmaxcell = cellmax = upol = rho = amubar = amubar2 = f32(0.0)
    ndist = 0
    ierr = 0
    # NB: these locals are uninitialised in the Fortran before the reads that set them;
    # a null list-directed value keeps the old value (0 as in the Rust port)
    nrho = norder = massno = isubs = nstruct = 0
    fracmix = f32(0.0)
    ifind = 0
    #
    # generate the name of the .ato file
    #
    ext = ".ato"
    ilenext = 4
    if commands.nwords < 2:
        commands.nwords = 0
        while commands.nwords <= 0:
            iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)
            commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 5, commands.nstsize,
                                          commands.nworddim)
            if commands.nwords <= 0:
                units.write(6, fmt(_F499))
                return
        readfile = fstr(substr(commands.sometext, ncf[1], ncl[1]), 256)
        ilenreadfile = ncl[1] - ncf[1] + 1
    else:
        readfile = fstr(substr(commands.sometext, ncf[2], ncl[2]), 256)
        ilenreadfile = ncl[2] - ncf[2] + 1
    #
    # strip out any extension that was supplied - extension must be .ato, otherwise
    # nothing is done.
    #
    ilendot = findex(readfile, ".") - 1
    if ilendot > 0:
        readfile = fstr(substr(readfile, 1, ilendot), 256)
        ilenreadfile = ilendot
    #
    # create the full name of the input file by adding on the directory name if not present
    #
    symbols.fnameato, symbols.ilenfnameato = file_and_directory_commands.get_full_filename(
        readfile, ilenreadfile, symbols.fnameato, symbols.ilenfnameato)
    symbols.fnameato = fstr(substr(symbols.fnameato, 1, symbols.ilenfnameato) + ".ato", 256)
    symbols.ilenfnameato = symbols.ilenfnameato + 4
    #
    # initialise number of atom and molecule types to zero
    #
    mol.ntype = 0
    clustercom.nmoltype = 0
    #
    # read .ato file
    #
    mol.sizefactor = f32(1.0)
    (rmaxcell, cellmax, upol, rho, amubar, amubar2, ndist,
     ierr) = ato_read_and_save.minitc(symbols.fnameato, rmaxcell, cellmax, upol, rho, amubar, amubar2,
                                      ndist, ierr)
    if rho == 0.0 or ierr != 0:
        units.write(6, fmt(_F599, Ch(symbols.fnameato, 256)))
        return
    units.write(6, lw())
    units.write(6, fmt(_F100))
    # 99    format(a80)
    units.write(6, fmt(_F200), advance=False)
    nrho = units.list_read(5).int(nrho)
    #
    # get the no. of atoms per molecule ratio in case it is needed
    #
    atmolrat = f32(mol.natomr) / f32(clustercom.nmol)
    items = []
    for it in range(1, mol.ntype + 1):
        items += [it, Ch(mol.type[it], 3), mol.cfrac[it], Ch(mol.atsymbol[it], 3)]
    units.write(6, fmt(_F202, *items))
    while True:
        # 203 continue
        units.write(6, fmt(_F219), advance=False)
        norder = units.list_read(5).int(norder)
        #
        # the program assumes
        # norder = 1 means this is single dataset
        # norder = 2 means this is a first-order difference experiment
        # norder = 3 means this is a second-order difference experiment
        #
        if norder < 1:
            return
        if norder > 3:
            continue  # go to 203
        break
    #
    # set flags to say which (if any) atoms have been substituted
    #
    for it in range(1, mol.ntype + 1):
        iflagsubs[it] = 0
    nsubs = 0
    #
    # Get the scattering lengths for the first sample
    #
    if norder > 1:
        units.write(6, fmt(_F207))
    else:
        units.write(6, fmt(_F223))
    for it in range(1, mol.ntype + 1):
        #
        # cannot see "Q" type atoms
        #
        if mol.qtypea[it] == 0:
            while True:
                # 206 continue
                chsymtest = fstr(mol.atsymbol[it], 3)
                massno = 0
                #
                # test to see if this chemical symbol exists in the table
                #
                btest, ifind = ato_manipulation_rev4.scatt_len(chsymtest, massno, ifind)
                if ifind == 0:
                    units.write(6, fmt(_F209, Ch(chsymtest, 3), Ch(mol.type[it], 3)), advance=False)
                    try:
                        mol.atsymbol[it] = fstr(units.list_read(5).word(mol.atsymbol[it]), 3)
                        ierr = 0
                    except (FortranEOF, FortranIOError) as e:
                        ierr = iostat_of(e)
                    continue  # go to 206
                break
            #
            # for hydrogen atoms we need to check whether they exchange with other hydrogens
            #
            if mol.atsymbol[it] == "H":
                mol.iexchange[it] = -1
                while mol.iexchange[it] != 0 and mol.iexchange[it] != 1:
                    units.write(6, fmt(_F103, Ch(mol.type[it], 3)), advance=False)
                    try:
                        mol.iexchange[it] = units.list_read(5).int(mol.iexchange[it])
                        ierr = 0
                    except (FortranEOF, FortranIOError) as e:
                        ierr = iostat_of(e)
                if mol.iexchange[it] != 1:
                    mol.iexchange[it] = 0
            else:
                mol.iexchange[it] = 0
            asum = f32(0)
            nsum = 1
            while True:
                # 205 continue
                units.write(6, fmt(_F211, Ch(mol.type[it], 3)), advance=False)
                ok, massno = _read5_err("int", massno)
                if not ok:
                    continue  # err=205
                units.write(6, fmt(_F224), advance=False)
                ok, abund[it, nsum, 1] = _read5_err("real", abund[it, nsum, 1])
                if not ok:
                    continue  # err=205
                if massno == 0:
                    batom[it, nsum, 1] = btest
                else:
                    batom[it, nsum, 1], ifind = ato_manipulation_rev4.scatt_len(chsymtest, massno, ifind)
                    if ifind == 0:
                        units.write(6, fmt(_F212, Ch(mol.atsymbol[it], 3)))
                        continue  # go to 205
                #     write(6,*) it,nsum,batom(it,nsum,1)
                asum = asum + abund[it, nsum, 1]
                if asum < 1.0:
                    nsum = nsum + 1
                    continue  # go to 205
                else:
                    for isum in range(1, nsum + 1):
                        abund[it, nsum, 1] = abund[it, nsum, 1] / asum
                break
            nisotope[it, 1] = nsum
    #
    # we now need to check that components which
    # are stated to be exchangeable through the sample do indeed all have the
    # the same scattering length.
    #
    # Step through the components and if iexchange is 1 form the average scattering
    # for all the atoms with the same chemical symbol.
    #
    bscatsum = f32(0.0)
    cfracsum = f32(0.0)
    for it in range(1, mol.ntype + 1):
        if mol.iexchange[it] == 1:
            cfracsum = cfracsum + mol.cfrac[it]
            for ii in range(1, nisotope[it, 1] + 1):
                term = mol.cfrac[it] * abund[it, ii, 1] * batom[it, ii, 1]
                bscatsum = bscatsum + term
    if cfracsum > 0.0:
        bscatsum = bscatsum / cfracsum
    else:
        bscatsum = f32(0.0)
    #
    # go through each atomic type and if iexchange is 1, then replace the individual
    # scattering lengths with a single average scattering length
    #
    for it in range(1, mol.ntype + 1):
        if mol.iexchange[it] == 1:
            nisotope[it, 1] = 1
            abund[it, 1, 1] = 1.0
            batom[it, 1, 1] = bscatsum
    #
    # now for the second sample, if appropriate identify those atoms
    # which are substituted and given the relevant scattering lengths
    # for the substituted atoms.
    #
    if norder > 1:
        while True:
            # 204 continue
            units.write(6, fmt(_F213), advance=False)
            #
            # obviously it is important to ensure you do not specify an atom
            # which previously was treated as if only 1 isotope exists. The
            # programme does not do any checking on this.
            #
            ok, isubs = _read5_err("int", isubs)
            if not ok:
                continue  # err=204
            if isubs >= 1 and isubs <= mol.ntype:
                if mol.qtypea[isubs] != 0:
                    units.write(6, fmt(_F155))
                    continue  # goto 204
                nsubs = nsubs + 1
                iflagsubs[isubs] = 1
                #
                # read in mass number and abundance
                #
                asum = f32(0)
                nsum = 1
                while True:
                    # 208 continue
                    units.write(6, fmt(_F211, Ch(mol.type[isubs], 3)), advance=False)
                    ok, massno = _read5_err("int", massno)
                    if not ok:
                        continue  # err=208
                    units.write(6, fmt(_F224), advance=False)
                    ok, abund[isubs, nsum, 2] = _read5_err("real", abund[isubs, nsum, 2])
                    if not ok:
                        continue  # err=208
                    chsymtest = fstr(mol.atsymbol[isubs], 3)
                    batom[isubs, nsum, 2], ifind = ato_manipulation_rev4.scatt_len(chsymtest, massno, ifind)
                    if ifind == 0:
                        units.write(6, fmt(_F212, Ch(mol.atsymbol[isubs], 3)))
                        continue  # go to 208
                    asum = asum + abund[isubs, nsum, 2]
                    if asum < 1.0:
                        nsum = nsum + 1
                        continue  # go to 208
                    else:
                        for isum in range(1, nsum + 1):
                            abund[isubs, nsum, 2] = abund[isubs, nsum, 2] / asum
                    break
                nisotope[isubs, 2] = nsum
                continue  # go to 204
            break
        #
        # for all the NON-substituted atoms, set the relevant scattering
        # lengths and parameters to those of 1st sample
        #
        for it in range(1, mol.ntype + 1):
            if iflagsubs[it] == 0:
                nisotope[it, 2] = nisotope[it, 1]
                for ib in range(1, nisotope[it, 2] + 1):
                    batom[it, ib, 2] = batom[it, ib, 1]
                    abund[it, ib, 2] = abund[it, ib, 1]
        #
        # we now need to check that components which have the same chemical symbol and which
        # are stated to be exchangeable through the sample do indeed all have the
        # the same scattering length.
        #
        # Step through the components and if iexchange is 1 form the average scattering
        # for all the atoms with the same chemical symbol.
        #
        bscatsum = f32(0.0)
        cfracsum = f32(0.0)
        for it in range(1, mol.ntype + 1):
            if mol.iexchange[it] == 1:
                cfracsum = cfracsum + mol.cfrac[it]
                for ii in range(1, nisotope[it, 2] + 1):
                    term = mol.cfrac[it] * abund[it, ii, 2] * batom[it, ii, 2]
                    bscatsum = bscatsum + term
        if cfracsum > 0.0:
            bscatsum = bscatsum / cfracsum
        else:
            bscatsum = f32(0.0)
        #
        # go through each atomic type and if iexchange is 1, then replace the individual
        # scattering lengths with a single average scattering length
        #
        for it in range(1, mol.ntype + 1):
            if mol.iexchange[it] == 1:
                nisotope[it, 2] = 1
                abund[it, 1, 2] = 1.0
                batom[it, 1, 2] = bscatsum
    #
    # if this is a 2nd-order difference experiment, then third sample is
    # treated as a mixture of the first two. If this is not the case then
    # you can either use this programme to set up the weights for each DCS file
    # separately, and so fit each total DCS, or else write another version!
    #
    if norder == 3:
        units.write(6, fmt(_F214), advance=False)
        fracmix = units.list_read(5).real(fracmix)
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
    if norder > 2:
        ndata = 2
    else:
        ndata = norder
    #
    # go through all the components and calculate the isotope averaged scattering
    # length for datasets 1 and 2
    #
    for idata in range(1, ndata + 1):
        for it1 in range(1, mol.ntype + 1):
            if mol.qtypea[it1] == 0:
                #
                # average scattering length of this component
                #
                binter[it1, idata] = 0.0
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
                    if (mol.iexchange[it1] == 0 and mol.iexchange[it2] == 0
                            and nisotope[it1, idata] == nisotope[it2, idata]
                            and mol.atsymbol[it1] == "H" and mol.atsymbol[it2] == "H"):
                        bintra = f32(0.0)
                        for ib in range(1, nisotope[it1, idata] + 1):
                            # also the abundances have to be identical, but we won't check
                            # since if they are not, then this is probably not a valid
                            # correlation in any case
                            term = abund[it1, ib, idata] * batom[it1, ib, idata] * batom[it2, ib, idata]
                            bintra = bintra + term
                        wintra[ipart, idata] = bintra
                    else:
                        wintra[ipart, idata] = winter[ipart, idata]
                    term = cprod[ipart] * winter[ipart, idata]
                    wtsum[idata] = wtsum[idata] + term
    #
    # now do third sample if present
    #
    if norder == 3:
        ipart = 0
        wtsum[3] = 0
        for it1 in range(1, mol.ntype + 1):
            if mol.qtypea[it1] == 0:
                #
                # average scattering length for this component
                #
                binter[it1, 3] = frac1 * binter[it1, 1] + frac2 * binter[it1, 2]
        for it1 in range(1, mol.ntype + 1):
            if mol.qtypea[it1] == 0:
                for it2 in range(it1, mol.ntype + 1):
                    if mol.qtypea[it2] == 0:
                        ipart = ipart + 1
                        #
                        # for the intermolecular terms, we just have the random mixing result
                        #
                        winter[ipart, 3] = binter[it1, 3] * binter[it2, 3]
                        #
                        # for intramolecular terms, we have to mix the intra-molecular terms in proportion
                        # unless the atoms are exchangeable, in which you get just the random mixing result
                        #
                        if mol.iexchange[it1] == 0 and mol.iexchange[it2] == 0:
                            wintra[ipart, 3] = frac1 * wintra[ipart, 1] + frac2 * wintra[ipart, 2]
                        else:
                            #
                            # otherwise we just have the random mixing result
                            #
                            wintra[ipart, 3] = winter[ipart, 3]
                        term = cprod[ipart] * winter[ipart, 3]
                        wtsum[3] = wtsum[3] + term
    #
    # write out the total weights for the dataset(s)
    #
    commands.nwords = 0
    fname = ""
    while commands.nwords == 0:
        units.write(6, fmt(_F220), advance=False)
        commands.sometext = fstr(units.read_line(5), commands.nstsize)
        commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 5, commands.nstsize,
                                      commands.nworddim)
        fname = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                     + substr(commands.sometext, ncf[1], ncl[1]), 80)
    #
    # find out how many characters there are in this name
    #
    nchar = directory_name.nchardir + ncl[1] - ncf[1] + 1
    for idata in range(1, norder + 1):
        units.write(6, fmt(_F201, idata), advance=False)
        #
        # ("Normalised" means divided by the sum of the weight factors)
        #
        nstruct = units.list_read(5).int(nstruct)
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
        fname = fstr(substr(fname, 1, nchar) + fnameb, 80)
        units.write(6, fmt(_F113, substr(fname, 1, nchart)))
        units.open(10, fname, status="unknown")
        #
        # signal whether these weights are per atom or per molecule
        #
        if nrho == 2:
            if nstruct == 1:
                rhofac = f32(1.0) / atmolrat
            else:
                rhofac = atmolrat
        else:
            rhofac = f32(1.0)
        units.write(10, lw(f32(1.0)))
        for ipart in range(1, npart + 1):
            it1 = mol.n1extr[ipart]
            it2 = mol.n2extr[ipart]
            if mol.qtypea[it1] == 0 and mol.qtypea[it2] == 0:
                winterwr = winter[ipart, idata]
                wintrawr = wintra[ipart, idata]
                #
                # relative contribution to total
                #
                relativewt = winter[ipart, idata] * cprod[ipart]
                #
                # If required normalise the weights to the total sum of weights
                #
                if nstruct == 1:
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
        units.write(6, fmt(_F301), advance=False)
        #
        # ("Normalised" means divided by the sum of the weight factors)
        #
        nstruct = units.list_read(5).int(nstruct)
        #
        # generate file name to write to
        #
        fnameb = fstr(fmt(_F115), 8)
        #
        # concatenate these two filenames
        #
        fname = fstr(substr(fname, 1, nchar) + fnameb, 80)
        units.write(6, fmt(_F116, substr(fname, 1, nchar + 7)))
        units.open(10, fname, status="unknown")
        #
        # signal whether these weights are per atom or per molecule
        #
        if nrho == 2:
            if nstruct == 1:
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
        fname = fstr(substr(fname, 1, nchar) + fnameb, 80)
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
        fname = fstr(substr(fname, 1, nchar) + fnameb, 80)
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
        fname = fstr(substr(fname, 1, nchar) + fnameb, 80)
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
    return

