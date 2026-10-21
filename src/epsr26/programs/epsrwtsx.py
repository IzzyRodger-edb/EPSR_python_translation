"""Translation of fortran_src/epsrwtsx.f: ``program standalone``, ``subroutine epsrwtsx``
and ``subroutine form_fact`` -- X-ray form-factor "weights" for the atom types of a .ato
file, from answers on standard input.
"""

from ..common import clustercom, commands, directory_name, epsr_dimension, molecom, symbols
from ..fortran import (
    Ch,
    FortranEOF,
    FortranIOError,
    char_array,
    dble_array,
    f32,
    findex,
    fmt,
    fstr,
    iostat_of,
    isign,
    lw,
    read_internal,
    substr,
    units,
)
from ..library import ato_read_and_save
from ..library import file_and_directory_commands
from ..library import parse
from ..trace import routine

# FORMAT statements ('$' dropped: written with advance=False)
_F499 = "(/'epsrwtsx> No .ato file specified! - routine will exit'/)"
_F599 = "(/'epsrwtsx> Problems with .ato file ',a/'Routine will exit now'/)"
_F100 = "(/'epsrwtsx> Program to setup the X-ray form factors for','for the specified .ato file.'/)"
_F200 = "(/'epsrwtsx> Is the output to be per atom (1)',', or per molecule (2)? ')"
_F202 = ("(/'epsrwtsx> The following components were found in this file:-'/1x,"
         "/1x,'Component no., label,  atomic fraction, chemical symbol'/1x,"
         "/((5x,i5,8x,a3,5x,e12.5,5x,a3)))")
_F223 = "(/'epsrwtsx> Get the form factor coefficients for all',' components in the sample'/)"
_F300 = "(/'epsrwtsx> For component ',a3,' Type form factor symbol: ')"
_F209 = ("(/'epsrwtsx> Did not recognise form factor symbol ',a4,' for atom type ',a3"
         "//'epsrwtsx> Try typing it again: ')")
_F210 = "(/'epsrwtsx> Specify q and delta for MAFF ','(0 0 for IAFF): ')"
_F224 = ("(/'epsrwtsx> Are data normalised to single atom scattering (1)',"
         "' or mean form factor squared (2)? ')")
_F225 = "(/'epsrwtsx> Use the IAFFs (1) or MAFFs (2) for normalisation? ')"
_F220 = "(/'epsrwtsx> Type basename of file to output weights to: ')"
_F111 = "(a)"
_F113 = "(/'epsrwtsx> Writing X-ray weights to file ',a)"
_F114 = "(1x,a3,1x,a4)"
_F115 = "(11(1x,f11.6))"
# form_fact
_FF_F301 = "(/'form_factor> X-ray file ',a15,' does not exist')"


def main():
    """epsrwtsx.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    epsrwtsx()


@routine
def epsrwtsx():
    """epsrwtsx.f:epsrwtsx -- Program EPSRWTSX.

    Sets up the X-ray form factors for the elements in the specified .ato file.
    Form factor data are read from the Waasmaier & Kirfel compilation (Acta Cryst.
    (1995). A51, 416-413): f0[k] = c + [SUM a_i*EXP(-b_i*(k^2))], i=1,5, where
    k = sin(theta)/lambda; columns a1 a2 a3 a4 a5 c b1 b2 b3 b4 b5.
    Currently set up for 70 components (number of atom types).
    """
    mcomp = epsr_dimension.mcomp
    ncf = commands.ncf
    ncl = commands.ncl
    mol = molecom
    # character*256 readfile ; character*4 ext ; character*1 iaccept
    iaccept = ""  # NB: not initialised in the Fortran; search_for_files sets it
    # real*8 ax(5,mcomp),bx(5,mcomp),cx(mcomp),alphax(mcomp),deltax(mcomp)
    ax = dble_array(5, mcomp)
    bx = dble_array(5, mcomp)
    # NB: cx(it) is never set for a Q atom in the Fortran (uninitialised); 0 here
    cx = dble_array(mcomp)
    alphax = dble_array(mcomp)
    deltax = dble_array(mcomp)
    # real*8 axt(5),bxt(5),cxt
    axt = dble_array(5)
    bxt = dble_array(5)
    cxt = 0.0
    # character*9 fnameb ; character*80 fname ; character*4 chsymtest,ffsym(mcomp)
    ffsym = char_array(mcomp)
    chsymtest = ""  # NB: uninitialised in the Fortran until read
    # integer*4 nrho,normx
    nrho = normx = 0
    iafformaff = 0  # NB: uninitialised in the Fortran; read with iostat=
    ifind = 0
    # implicitly typed locals passed to minitc (epsr16.inc is not included)
    rmaxcell = cellmax = upol = rho = amubar = amubar2 = f32(0.0)
    ndist = 0
    ierr = 0
    #
    # generate the name of the .ato file
    #
    ext = ".ato"
    fnameb = "X-ray.wts"
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

    # Check that a sensible value of nrho is specified. If not it is assumed to be
    # per molecule.

    if nrho < 1 or nrho > 2:
        nrho = 2

    # If nrho .eq. 2, then signal per molecule by making it -1

    if nrho == 2:
        nrho = -1

    #
    # get the no. of atoms per molecule ratio in case it is needed
    #
    atmolrat = f32(mol.natomr) / f32(clustercom.nmol)  # noqa: F841
    items = []
    for it in range(1, mol.ntype + 1):
        items += [it, Ch(mol.type[it], 3), mol.cfrac[it], Ch(mol.atsymbol[it], 3)]
    units.write(6, fmt(_F202, *items))
    # 203 continue
    norder = 1  # noqa: F841
    #
    # norder = 1 means this is single dataset (X-ray form factors have no isotope
    # dependence).
    #
    units.write(6, fmt(_F223))
    for it in range(1, mol.ntype + 1):
        #
        # cannot see "Q" type atoms
        #
        if mol.qtypea[it] == 0:
            # 206 continue
            units.write(6, fmt(_F300, Ch(mol.type[it], 3)), advance=False)
            try:
                chsymtest = fstr(units.list_read(5).word(chsymtest), 4)
                ierr = 0
            except (FortranEOF, FortranIOError) as e:
                ierr = iostat_of(e)
            #
            # test to see if this form factor symbol exists in the table
            #
            ifind = 0
            while ifind == 0:
                cxt, ifind = form_fact(chsymtest, axt, bxt, cxt, ifind, 5)
                if ifind == 0:
                    units.write(6, fmt(_F209, Ch(chsymtest, 4), Ch(mol.type[it], 3)), advance=False)
                    try:
                        chsymtest = fstr(units.list_read(5).word(chsymtest), 4)
                        ierr = 0
                    except (FortranEOF, FortranIOError) as e:
                        ierr = iostat_of(e)
            for i in range(1, 6):
                ax[i, it] = axt[i]
                bx[i, it] = bxt[i]
            cx[it] = cxt
            ffsym[it] = chsymtest

            # Read the alpha and delta values to define the modified atomic form factors (MAFF)
            # 0 0 means use the independent atomic form factors (IAFF)

            units.write(6, fmt(_F210), advance=False)
            try:
                r = units.list_read(5)
                alphax[it] = r.dble(alphax[it])
                deltax[it] = r.dble(deltax[it])
                ierr = 0
            except (FortranEOF, FortranIOError) as e:
                ierr = iostat_of(e)

        else:

            for i in range(1, 6):
                ax[i, it] = 0.0
                bx[i, it] = 0.0
            alphax[it] = 0.0
            deltax[it] = 0.0
            ffsym[it] = "nil"

    # Decide whether data are normalised to single atom scattering (1) or mean scattering
    # length squared.

    units.write(6, fmt(_F224), advance=False)
    normx = units.list_read(5).int(normx)

    # Decide whether to use the IAFFs or the MAFFs for the normalisation

    units.write(6, fmt(_F225), advance=False)
    try:
        iafformaff = units.list_read(5).int(iafformaff)
        ierr = 0
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)

    # Based on the values of nrho, normx and iafformaff, calculate the data type.

    ndatatype = isign(1, nrho) * (normx + 2 * (iafformaff - 1))  # noqa: F841

    # write out the total weights for the dataset(s)

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
    nchart = nchar + 9
    #
    # concatenate these two filenames
    #
    fname = fstr(substr(fname, 1, nchar) + fnameb, 80)
    units.write(6, fmt(_F113, substr(fname, 1, nchart)))
    units.open(10, fname, status="unknown")

    # signal whether these weights are per atom or per molecule (not used for
    # X-ray weights)

    units.write(10, lw(nrho, normx, iafformaff))
    for it in range(1, mol.ntype + 1):
        units.write(10, fmt(_F114, Ch(mol.type[it], 3), Ch(ffsym[it], 4)))

        # 28/1/2006 we add alpha and delta for modified atomic form factors if needed.
        # See equation 5, JCP, vol113, p 9149-9161 (2000), which however is incorrect!

        vals = [cx[it]]
        for i in range(1, 6):
            vals += [ax[i, it], bx[i, it]]
        vals += [alphax[it], deltax[it]]
        units.write(10, fmt(_F115, *vals))

    units.close(10)
    #
    # write the .ato file in case anything was changed
    #
    #     call moloutc(fnameato)
    return


@routine
def form_fact(chsym, axt, bxt, cxt, ifind, ndim):
    """epsrwtsx.f:form_fact -- read the X-ray form factor coefficients of symbol chsym
    (character*4) from f0_WaasKirf.dat.  Returns (cxt, ifind)."""
    ncf = commands.ncf
    ncl = commands.ncl
    # character*4 chsym ; character*15 fnameX ; real*8 axt(ndim),bxt(ndim),cxt
    chsym = fstr(chsym, 4)
    fnameX = "f0_WaasKirf.dat"

    ifind = 0

    # check the number of non-space characters in chsym

    ilenchsym = len(chsym)
    if ilenchsym < 1:
        return cxt, ifind

    # Open form factor file and skip some lines

    try:
        units.open(10, fnameX, status="old")
        ierr = 0
    except (FortranEOF, FortranIOError) as e:
        ierr = iostat_of(e)

    if ierr != 0:
        units.write(6, fmt(_FF_F301, Ch(fnameX, 15)))
        return cxt, ifind

    while ierr == 0 and ifind == 0:

        try:
            commands.sometext = fstr(units.read_line(10), commands.nstsize)
            ierr = 0
        except (FortranEOF, FortranIOError) as e:
            ierr = iostat_of(e)

        # Search sometext for the character string showing a valid data line

        if findex(commands.sometext, "#S") > 0:

            commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 5, commands.nstsize,
                                          commands.nworddim)

            # The third word should contain the text being searched for

            if substr(commands.sometext, ncf[3], ncl[3]) == substr(chsym, 1, ilenchsym):
                for _k in range(2):
                    try:
                        units.skip_record(10)
                        ierr = 0
                    except (FortranEOF, FortranIOError) as e:
                        ierr = iostat_of(e)
                try:
                    commands.sometext = fstr(units.read_line(10), commands.nstsize)
                    ierr = 0
                except (FortranEOF, FortranIOError) as e:
                    ierr = iostat_of(e)

                # Check that the number of words is correct

                commands.nwords = parse.parse(commands.sometext, commands.nwords, ncf, ncl, 5, commands.nstsize,
                                              commands.nworddim)
                if ierr == 0 and commands.nwords == 11:
                    r = read_internal(commands.sometext)
                    for i in range(1, 6):
                        axt[i] = r.dble(axt[i])
                    cxt = r.dble(cxt)
                    for i in range(1, 6):
                        bxt[i] = r.dble(bxt[i])
                    ifind = 1
    units.close(10)
    return cxt, ifind
