"""Translation of fortran_src/system_commands.f."""

from ..common import commands, system_commands
from ..fortran import FortranEOF, FortranIOError, f32, fstr, r4, read_internal, substr, units
from . import parse
from ..trace import routine

@routine
def get_system_commands():
    """system_commands.f:get_system_commands -- read system_commands.txt (current folder).

    initialise the system commands, then set each one named in the file.
    """
    #     write(6,*) 'Got to get_system_commands'
    #
    # initialise the system commands
    #
    system_commands.system_ls = ""
    system_commands.ilensystem_ls = 0
    system_commands.system_cd = ""
    system_commands.ilensystem_cd = 0
    system_commands.system_md = ""
    system_commands.ilensystem_md = 0
    system_commands.system_edit = ""
    system_commands.ilensystem_edit = 0
    system_commands.system_del = ""
    system_commands.ilensystem_del = 0
    #
    # termination character(s) for the end of a directory name
    # This will be added to the end of directory names if it is not present
    #
    system_commands.system_termination = ""
    system_commands.ilensystem_termination = 0
    #
    # characters to join commands in one line
    #
    system_commands.system_join = ""
    system_commands.ilensystem_join = 0
    #
    # command to invoke gnuplot
    #
    system_commands.system_gnu = ""
    system_commands.ilensystem_gnu = 0

    # Command to invoke jmol

    system_commands.system_jmol = ""
    system_commands.ilensystem_jmol = 0

    # Command to invoke Mopac

    system_commands.system_mopac = ""
    system_commands.ilensystem_mopac = 0

    # Stores the allowed number of orbitals in Mopac. 300 is the default.

    system_commands.system_mopac_norbitals = "300"
    system_commands.ilensystem_mopac_norbitals = 3

    # System binaries location

    system_commands.system_binaries = ""
    system_commands.ilensystem_binaries = 0

    # System binaries extension

    system_commands.system_bin_ext = ""
    system_commands.ilensystem_bin_ext = 0
    #
    # commands to do 2d and 3d plotting
    #
    system_commands.system_plot3d = ""
    system_commands.ilensystem_plot3d = 0
    system_commands.system_plot3djmol = ""
    system_commands.ilensystem_plot3djmol = 0
    system_commands.system_plot2d = ""
    system_commands.ilensystem_plot2d = 0
    system_commands.system_plotato = ""
    system_commands.ilensystem_plotato = 0
    system_commands.system_pgout = ""
    system_commands.ilensystem_pgout = 0
    system_commands.system_quotes = ""
    system_commands.ilensystem_quotes = 0
    system_commands.gnu_pointtype = 0
    system_commands.gnu_pointsize = f32(0.0)
    #
    # open the system_command file in the current directory
    # If a fault occurs, the Windows Command Prompt commands are assumed
    #
    try:
        units.open(10, "system_commands.txt", status="old")
        ierr = 0
    except FortranIOError as e:
        ierr = e.iostat
    while ierr == 0:
        commands.sometext = ""
        try:
            commands.sometext = fstr(units.read_line(10), commands.nstsize)
        except FortranEOF:
            ierr = -1
        except FortranIOError as e:
            ierr = e.iostat
        #         write(6,*) sometext(1:len_trim(sometext))
        #
        # parse this text to see if contains any words
        #
        # len(sometext) = 3020 = nstsize
        commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 11,
                                      commands.nstsize, commands.nworddim)
        #     write(6,*) nwords,(' ',sometext(ncf(i):ncl(i)),i=1,nwords)
        if commands.nwords > 1:
            word1 = substr(commands.sometext, commands.ncf[1], commands.ncl[1])
            if word1 == "system_ls":
                nlast = commands.nwords
                _, nlast, system_commands.system_ls, system_commands.ilensystem_ls = parse.concatenate(
                    2, nlast, system_commands.system_ls, system_commands.ilensystem_ls,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
                #              write(6,102) ilensystem_ls,system_ls(1:ilensystem_ls)
                # 102   format(1x,i3,1x,a)
            elif word1 == "system_cd":
                nlast = commands.nwords
                _, nlast, system_commands.system_cd, system_commands.ilensystem_cd = parse.concatenate(
                    2, nlast, system_commands.system_cd, system_commands.ilensystem_cd,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "system_md":
                nlast = commands.nwords
                _, nlast, system_commands.system_md, system_commands.ilensystem_md = parse.concatenate(
                    2, nlast, system_commands.system_md, system_commands.ilensystem_md,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "system_edit":
                nlast = commands.nwords
                _, nlast, system_commands.system_edit, system_commands.ilensystem_edit = parse.concatenate(
                    2, nlast, system_commands.system_edit, system_commands.ilensystem_edit,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "system_del":
                nlast = commands.nwords
                _, nlast, system_commands.system_del, system_commands.ilensystem_del = parse.concatenate(
                    2, nlast, system_commands.system_del, system_commands.ilensystem_del,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "system_termination":
                nlast = commands.nwords
                _, nlast, system_commands.system_termination, system_commands.ilensystem_termination = parse.concatenate(
                    2, nlast, system_commands.system_termination, system_commands.ilensystem_termination,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "system_join":
                nlast = commands.nwords
                _, nlast, system_commands.system_join, system_commands.ilensystem_join = parse.concatenate(
                    2, nlast, system_commands.system_join, system_commands.ilensystem_join,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "system_gnu":
                nlast = commands.nwords
                _, nlast, system_commands.system_gnu, system_commands.ilensystem_gnu = parse.concatenate(
                    2, nlast, system_commands.system_gnu, system_commands.ilensystem_gnu,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "system_jmol":
                nlast = commands.nwords
                _, nlast, system_commands.system_jmol, system_commands.ilensystem_jmol = parse.concatenate(
                    2, nlast, system_commands.system_jmol, system_commands.ilensystem_jmol,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "system_binaries":
                nlast = commands.nwords
                _, nlast, system_commands.system_binaries, system_commands.ilensystem_binaries = parse.concatenate(
                    2, nlast, system_commands.system_binaries, system_commands.ilensystem_binaries,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "system_bin_ext":
                nlast = commands.nwords
                _, nlast, system_commands.system_bin_ext, system_commands.ilensystem_bin_ext = parse.concatenate(
                    2, nlast, system_commands.system_bin_ext, system_commands.ilensystem_bin_ext,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "system_plot3d":
                nlast = commands.nwords
                _, nlast, system_commands.system_plot3d, system_commands.ilensystem_plot3d = parse.concatenate(
                    2, nlast, system_commands.system_plot3d, system_commands.ilensystem_plot3d,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "system_plot3djmol":
                nlast = commands.nwords
                _, nlast, system_commands.system_plot3djmol, system_commands.ilensystem_plot3djmol = parse.concatenate(
                    2, nlast, system_commands.system_plot3djmol, system_commands.ilensystem_plot3djmol,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "system_plot2d":
                nlast = commands.nwords
                _, nlast, system_commands.system_plot2d, system_commands.ilensystem_plot2d = parse.concatenate(
                    2, nlast, system_commands.system_plot2d, system_commands.ilensystem_plot2d,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "system_plotato":
                nlast = commands.nwords
                _, nlast, system_commands.system_plotato, system_commands.ilensystem_plotato = parse.concatenate(
                    2, nlast, system_commands.system_plotato, system_commands.ilensystem_plotato,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "system_pgout":
                nlast = commands.nwords
                _, nlast, system_commands.system_pgout, system_commands.ilensystem_pgout = parse.concatenate(
                    2, nlast, system_commands.system_pgout, system_commands.ilensystem_pgout,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "system_quotes":
                nlast = commands.nwords
                _, nlast, system_commands.system_quotes, system_commands.ilensystem_quotes = parse.concatenate(
                    2, nlast, system_commands.system_quotes, system_commands.ilensystem_quotes,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "system_mopac":
                nlast = commands.nwords
                _, nlast, system_commands.system_mopac, system_commands.ilensystem_mopac = parse.concatenate(
                    2, nlast, system_commands.system_mopac, system_commands.ilensystem_mopac,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "system_mopac_norbitals":
                nlast = commands.nwords
                _, nlast, system_commands.system_mopac_norbitals, system_commands.ilensystem_mopac_norbitals = parse.concatenate(
                    2, nlast, system_commands.system_mopac_norbitals, system_commands.ilensystem_mopac_norbitals,
                    commands.nwords, commands.ncf, commands.ncl, commands.sometext, 256)
            elif word1 == "gnu_pointtype":
                try:
                    system_commands.gnu_pointtype = read_internal(
                        substr(commands.sometext, commands.ncf[2], commands.ncl[commands.nwords])
                    ).int(system_commands.gnu_pointtype)
                except FortranEOF:
                    ierr = -1
                except FortranIOError as e:
                    ierr = e.iostat
            elif word1 == "gnu_pointsize":
                try:
                    system_commands.gnu_pointsize = read_internal(
                        substr(commands.sometext, commands.ncf[2], commands.ncl[commands.nwords])
                    ).real(system_commands.gnu_pointsize)
                except FortranEOF:
                    ierr = -1
                except FortranIOError as e:
                    ierr = e.iostat
    units.close(10)
    if system_commands.gnu_pointtype == 0:
        system_commands.gnu_pointtype = 8
    if system_commands.gnu_pointsize == 0.0:
        system_commands.gnu_pointsize = r4(0.3)
    return
