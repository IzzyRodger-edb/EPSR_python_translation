"""Translation of fortran_src/xwts.f: ``program standalone``, ``subroutine xwts`` and
``subroutine form_fact`` -- X-ray form-factor weights for the atom types of a .ato file,
driven by a ``.XWTS.dat`` keyword file (``setup_symbols_xwts``).

``form_fact`` is the same routine as in epsrwtsx.f (each Fortran file has its own copy).
"""

from ..common import commands, directory_name, epsr_dimension, molecom, symbols
from ..fortran import (
    Ch,
    FortranEOF,
    FortranIOError,
    char_array,
    dble_array,
    findex,
    fmt,
    fstr,
    iostat_of,
    read_internal,
    substr,
    units,
)
from ..library import file_and_directory_commands
from ..library import parse
from ..library import read_and_save_files
from ..trace import routine

# FORMAT statements
_F100 = "(/'xwts> Program to setup the X-ray form factors for','for the specified .ato file.'/)"
_F202 = ("(/'xwts> The following components were found in this file:-'/1x,"
         "/1x,'Component no., label,  atomic fraction, chemical symbol'/1x,"
         "/((5x,i5,8x,a3,5x,e12.5,5x,a3)))")
_F235 = "(/'xwts> Number of components in input file different from number of atom types in .ato file')"
_F223 = "(/'xwts> Get the form factor coefficients for all',' components in the sample'/)"
_F209 = ("(/'xwts> Did not recognise form factor symbol ',a4,' for atom type ',a3"
         "//'xwts> Try typing it again: ')")
_F114 = "(1x,a3,1x,a4)"
_F115 = "(11(1x,f11.6))"
# form_fact
_FF_F301 = "(/'form_factor> X-ray file ',a15,' does not exist')"


def main():
    """xwts.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    xwts()


@routine
def xwts():
    """xwts.f:xwts -- Program xwts.

    Sets up the X-ray form factors for the elements in the specified .ato file.
    Form factor data are read from the Waasmaier & Kirfel compilation (Acta Cryst.
    (1995). A51, 416-413): f0[k] = c + [SUM a_i*EXP(-b_i*(k^2))], i=1,5, where
    k = sin(theta)/lambda; columns a1 a2 a3 a4 a5 c b1 b2 b3 b4 b5.
    Currently set up for 70 components (number of atom types).
    """
    from ..library import setup_symbols_generic

    mcomp = epsr_dimension.mcomp
    mol = molecom
    # character*256 wtsfile
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
    # character*256 fname ; character*4 chsymtest,ffsym(mcomp)
    ffsym = char_array(mcomp)
    # integer*4 nrho,normx,ncharwts
    nrho = normx = iafformaff = ncomponent = 0
    ifind = 0
    ierr = 0
    units.write(6, fmt(_F100))

    # Get parameters from input file
    #
    # setup the xwts symbols
    #
    setup_symbols_generic.setup_symbols_xwts()
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

    units.write(6, fmt(_F100))
    # 99    format(a80)
    #
    # Set up the base name for the weights file. This will be the same as the input file.

    sv = symbols.symbol_value
    iref = 1
    wtsfile = fstr(substr(symbols.fname_run, 1, symbols.ilenfname_run), 256)
    ncharwts = symbols.ilenfname_run
    #
    # set up whether to normalise to atoms or molecules.
    #
    iref = iref + 1
    nrho = read_internal(sv[iref]).int(nrho)

    # Get type of normalisation employed

    iref = iref + 1
    normx = read_internal(sv[iref]).int(normx)

    # Decide on whether to use IAFFs or MAFFs for normalisation

    iref = iref + 1
    iafformaff = read_internal(sv[iref]).int(iafformaff)
    #
    # get the no. of atoms per molecule ratio in case it is needed
    #
    items = []
    for it in range(1, mol.ntype + 1):
        items += [it, Ch(mol.type[it], 3), mol.cfrac[it], Ch(mol.atsymbol[it], 3)]
    units.write(6, fmt(_F202, *items))
    # 203 continue
    #
    # norder = 1 means this is single dataset (X-ray form factors have no isotope
    # dependence).
    #

    # Check number of components in input file is same as atom types

    iref = iref + 1
    ncomponent = read_internal(sv[iref]).int(ncomponent)
    if ncomponent != mol.ntype:

        units.write(6, fmt(_F235))
        return

    units.write(6, fmt(_F223))
    for it in range(1, mol.ntype + 1):

        # Skip over next entry in input file as it is be determined from the ato file

        iref = iref + 1
        #
        # cannot see "Q" type atoms
        #
        if mol.qtypea[it] == 0:

            iref = iref + 1
            chsymtest = fstr(sv[iref], 4)
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
                        ierr = iostat_of(e)  # noqa: F841
            sv[iref] = chsymtest
            for i in range(1, 6):
                ax[i, it] = axt[i]
                bx[i, it] = bxt[i]
            cx[it] = cxt
            ffsym[it] = chsymtest

            # Read the alpha and delta values to define the modified atomic form factors (MAFF)
            # 0 0 means use the independent atomic form factors (IAFF)

            iref = iref + 1
            alphax[it] = read_internal(sv[iref]).dble(alphax[it])
            iref = iref + 1
            deltax[it] = read_internal(sv[iref]).dble(deltax[it])

        else:

            iref = iref + 3
            for i in range(1, 6):
                ax[i, it] = 0.0
                bx[i, it] = 0.0
            alphax[it] = 0.0
            deltax[it] = 0.0
            ffsym[it] = "nil"

    # write out the total weights for the dataset(s)

    fname = fstr(substr(directory_name.directory_name, 1, directory_name.nchardir)
                 + substr(wtsfile, 1, ncharwts) + ".wts", 256)
    units.open(10, fname, status="unknown")

    # signal whether these weights are per atom or per molecule (not used for
    # X-ray weights)

    units.write(10, fmt("(3(1x,i3))", nrho, normx, iafformaff))
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
    # Rewrite the input file in case anything was changed
    #
    read_and_save_files.save_files_others(14)

    return


@routine
def form_fact(chsym, axt, bxt, cxt, ifind, ndim):
    """xwts.f:form_fact -- read the X-ray form factor coefficients of symbol chsym
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
