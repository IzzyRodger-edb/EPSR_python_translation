"""Translation of fortran_src/introtcluster.f: ``program standalone`` and
``subroutine randomise`` (the ``randomise`` program) -- random rotation and
translation of all molecules in a cluster.
"""

from ..common import clustercom, commands, lattice, molecom, molecom_mod, symbols
from ..fortran import Ch, FortranStop, f32, findex, fmt, fstr, lw, substr, units
from ..library import ato_read_and_save
from ..library import auxiliary_routines
from ..library import file_and_directory_commands
from ..library import minimagec
from ..library import molfit
from ..library import parse
from ..trace import routine

_F499 = "(/'randomise> No .ato file specified! - routine will exit'/)"
_F599 = "(/'Problems with .ato file ',a/'Routine will exit now'/)"
_F193 = "('randomise> ',i5,' rotational groups in ',i5,' molecules ')"


def main():
    """introtcluster.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    randomise()
    raise FortranStop(0)


@routine
def randomise():
    """introtcluster.f:randomise -- perform a random rotation and translation of all
    molecules in a cluster."""
    # use molecom_mod
    # character*256 readfile
    # character*4 ext
    # character*1 iaccept
    iaccept = ""  # NB: not initialised in the Fortran; search_for_files sets it
    #
    # generate the name of the .ato file
    #
    ext = ".ato"
    ilenext = 4
    if commands.nwords < 2:
        commands.nwords = 0
        while commands.nwords <= 0:
            iaccept = file_and_directory_commands.search_for_files(ext, ilenext, iaccept)
            commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                          commands.nstsize, commands.nworddim)
            if commands.nwords <= 0:
                units.write(6, fmt(_F499))
                return
        readfile = fstr(substr(commands.sometext, commands.ncf[1], commands.ncl[1]), 256)
        ilenreadfile = commands.ncl[1] - commands.ncf[1] + 1
    else:
        readfile = fstr(substr(commands.sometext, commands.ncf[2], commands.ncl[2]), 256)
        ilenreadfile = commands.ncl[2] - commands.ncf[2] + 1
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
    molecom.ntype = 0
    clustercom.nmoltype = 0
    #
    # read .ato file
    #
    molecom.sizefactor = f32(1.0)
    # rmaxcell, cellmax, upol, rho, amubar, amubar2 (real) and ndist (integer) are
    # implicitly typed locals here (epsr16.inc is not included)
    rmaxcell = cellmax = upol = rho = amubar = amubar2 = f32(0.0)
    ndist = 0
    ierr = 0
    (rmaxcell, cellmax, upol, rho, amubar, amubar2, ndist,
     ierr) = ato_read_and_save.minitc(symbols.fnameato, rmaxcell, cellmax, upol, rho, amubar, amubar2,
                                      ndist, ierr)
    if rho == 0.0 or ierr != 0:
        units.write(6, fmt(_F599, Ch(symbols.fnameato, 256)))
        return
    units.write(6, lw(clustercom.dmax, molecom.dmaxmol))
    ato_read_and_save.normalise_coordinates()
    #
    # save the current rotational step sizes, and substitute 2*pi rotations
    # for the present purpose
    #
    stepris = molecom.stepri
    steprs = molecom.stepr
    # For rotational moves, only perform the randomisation if the relevant step size is .gt. 0
    if stepris > 0.0:
        molecom.stepri = f32(1.0)
    if steprs > 0.0:
        molecom.stepr = f32(1.0)
    for im in range(1, clustercom.nmol + 1):

        if not clustercom.tethered[im]:

            # perform random translation of whole molecule

            if molecom.stepm > 0.0:

                fa = (auxiliary_routines.ran1() - f32(0.5))
                fb = (auxiliary_routines.ran1() - f32(0.5))
                fc = (auxiliary_routines.ran1() - f32(0.5))

                # Convert to Cartesian coordinates

                #
                # calculate corrected coordinates
                #
                clustercom.zadd[im] = fc * lattice.clatz
                clustercom.yadd[im] = fb * lattice.blaty + fc * lattice.claty
                clustercom.xadd[im] = fa * lattice.alatx + fb * lattice.blatx + fc * lattice.clatx

                # ensure all molecules are in the central box ( + or - 0.5)

                clustercom.xadd[im], clustercom.yadd[im], clustercom.zadd[im] = minimagec.minimagec(
                    clustercom.xadd[im], clustercom.yadd[im], clustercom.zadd[im])
                #               call minimagec(phix(im),phiy(im),phiz(im))

            # perform random whole molecule rotation if requested
            if molecom.stepr > 0.0:
                molfit.rotatm(im, molecom.stepr)
    # perform random internal molecule rotation
    # !      if(ngrrot.gt.0.and.stepri.gt.0.0) then
    # !            do ig=1,ngrrot
    # !               if(rotgrp(ig).eq.'ROT') call rotati(ig) !Only rotate free rotation groups at this stage
    # !            end do
    # !      endif
    units.write(6, fmt(_F193, molecom_mod.ngrrot, clustercom.nmol))
    #
    # set the rotational step sizes back to what they were on input
    #
    molecom.stepri = stepris
    molecom.stepr = steprs
    #
    # write out the new molecule coordinates
    #
    units.write(6, lw(clustercom.dmax, molecom.dmaxmol))
    ato_read_and_save.unnormalise_coordinates()
    ato_read_and_save.moloutc(symbols.fnameato)
    return
