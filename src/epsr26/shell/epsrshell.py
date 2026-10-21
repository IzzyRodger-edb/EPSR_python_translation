"""Translation of fortran_src/epsrshell_standalone.f: the EPSR command shell.

It keeps the Fortran's prompt, banner, help text, command names and aliases, script mode
(``runflag.txt``, ``ss``/``ps``/``es``, ``#`` comments, rewinding the script at its end,
pausing and resuming across shell restarts) and ``FChistory.txt``.  What changes:

* programs run in this process (``programs.run_program``) instead of through
  ``call system(...)``; the shell's own COMMON state and open units are set aside while a
  program runs, as if it ran in a separate process;
* ``ls``, ``pwd``, ``cd``, ``md``, ``del``, ``edit``, ``system``/``!cmd`` are built-ins
  (``builtins.py``) that work on Linux;
* Ctrl+C and shutdown are handled (``signals.py``): no residual processes.

New commands: ``exit``/``quit`` (as ``bye``), ``jobs`` (live child processes),
``trace on [pattern] | off | summary`` (call tracing).
"""

import os
import sys
import traceback

from .. import common, programs
from ..common import commands, directory_name, version
from ..fortran import fmt, fstr, lw, substr, units
from ..fortran.stop import FortranEOF, FortranIOError, FortranStop
from ..library import parse
from ..trace import routine
from . import builtins, history, processes, signals

# commands that run a program (process_command in the Fortran), in its order
PROGRAM_COMMANDS = [
    "setup", "makeato", "readgpr", "readjmol", "readmole", "readxyz", "makemole", "fmole", "mixato",
    "addato", "makelattice", "makelatticeato", "dockato", "changeato", "randomise", "epsrwts", "nwts",
    "epsrwtsx", "xwts", "bonds", "epsr", "readpdb", "sharm", "sdf", "sdfcube", "clusters", "partials",
    "coord", "chains", "rings", "lsi", "voids", "plot", "plot2d", "splot2d", "plot3d", "plot3djmol",
    "plotato", "triangles", "torangles", "fluctuations", "writexyz", "runjmol", "mapgrentropy", "mapgr",
]

# Fortran binary name used by run_program -> key in programs.PROGRAMS
_PROGRAM_KEYS = {"upset": "upset", "mapgrquad": "mapgr"}

_F101 = "(/'EPSRshell> Welcome to EPSR version ',a,//'Type \"help\" or \"?\" for a list of commands')"
_F600 = "(/'Binaries folder: ',a)"


def _write(text=""):
    units.write(6, text)


@routine
def main(argv=None):
    """epsrshell_standalone.f:epsrshell -- the main command loop."""
    home = os.getcwd()
    history.set_home(home)
    history.setup()
    signals.install(home)
    signals.state.script_pauser = _pause_on_exit
    other = signals.take_lock()
    common.reset()
    units.reset()

    history.get_history_list()
    version.version = version.date
    version.ilenversion = version.ilendate
    _write(fmt(_F101, substr(version.version, 1, version.ilenversion)))
    if other is not None:
        _write(f"\nEPSRshell> Warning: another EPSRshell is running in this folder: PID {other}")
    # get the system commands for the particular operating system in use
    _get_system_commands()
    signals.state.keep_viewers = _keep_viewers()
    # The programs run inside this shell (Python translation), not from system_binaries
    _write(fmt(_F600, f"{os.path.dirname(os.path.dirname(__file__))} (Python, in-process)"))
    # set initial directory name - this is called the home directory
    # (file_and_directory_commands.f:get_directory_name; the separator is always '/')
    builtins.set_current(home)
    _write(fmt("(/'Home folder is: ',a)", substr(directory_name.directory_name, 1, directory_name.nchardir)))
    directory_name.ncharhomedir = directory_name.nchardir
    directory_name.home_directory = directory_name.directory_name

    # initialize the text stream
    commands.sometext = ""
    # get the global script status, and start the script if it is not running
    commands.ncharscript = 0
    start_script(1)

    while commands.runflag in ("interactive", "script"):
        commands.sometext = ""
        # read a line of text from the appropriate input stream
        if commands.runflag == "interactive":
            _write(lw())
            try:
                commands.sometext = fstr(history.readline("EPSRshell>"), commands.nstsize)
            except KeyboardInterrupt:
                # Ctrl+C at the prompt: discard the line (as bash does)
                sys.stdout.write("\n")
                continue
            except EOFError:
                # Ctrl+D: as 'bye'
                sys.stdout.write("\n")
                commands.sometext = "bye"
        elif commands.runflag == "script":
            try:
                commands.sometext = fstr(units.read_line(commands.ncommand_stream), commands.nstsize)
            except (FortranEOF, FortranIOError):
                commands.sometext = ""
            # parse the text to see if it contains any words
            commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                          commands.nstsize, commands.nworddim)
            # if no. of words is zero rewind the input file and try again
            # if this fails, abandon the run
            if commands.nwords <= 0:
                units.rewind(commands.ncommand_stream)
                try:
                    commands.sometext = fstr(units.read_line(commands.ncommand_stream), commands.nstsize)
                except (FortranEOF, FortranIOError) as e:
                    ierr1 = -1 if isinstance(e, FortranEOF) else e.iostat
                    _write(lw("Error condition ", ierr1, " after rewind"))
                    raise SystemExit(0)
        # now parse it into separate words
        commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 5,
                                      commands.nstsize, commands.nworddim)
        # see if the first character of the first word contains a # - if so we will ignore this line
        if commands.nwords > 0:
            ilenhash = substr(commands.sometext, commands.ncf[1], commands.ncl[commands.nwords]).find("#") + 1
        else:
            ilenhash = 1
        # if running a script, reconstruct the original text and write it out
        if commands.runflag == "script":
            nlast = commands.nwords
            _, nlast, output_string, nchartext = parse.concatenate(
                1, nlast, "", 0, commands.nwords, commands.ncf, commands.ncl, commands.sometext)
            _write(lw())
            _write(lw(substr(output_string, 1, nchartext)))
        if commands.nwords > 0 and ilenhash == 0:
            # find which command was called, if any, and act on it
            process_command()
        check_script_status()
        # if the global runflag has been set to interactive, this means the script is to be abandoned.
        if commands.runflagglobal == "interactive":
            if commands.runflag == "script":
                units.close(commands.ncommand_stream)
                commands.ncharscript = 0
                commands.runflag = "interactive"
        else:
            # if the current program is running the script and script_status has been set to
            # pausing, stop running the script and go to interactive mode, but leaving open
            # the possibility that the script can be restarted again
            if commands.runflag == "script":
                if commands.script_status in ("pausing", "ps"):
                    commands.script_status = "notrunning"
                    _write_runflag("script", commands.script_status, substr(commands.script_fname, 1, commands.ncharscript))
                    # close the script file
                    units.close(commands.ncommand_stream)
                    commands.runflag = "interactive"
                    if commands.runflagglobal == "exit":
                        commands.runflag = "exit"
                elif commands.runflagglobal == "exit":
                    _write_runflag("interactive")
                    commands.runflag = "exit"
    signals.shutdown()
    return 0


def _write_runflag(*words):
    """``open(10,file='runflag.txt')`` + ``write(10,*) words`` + ``close(10)``."""
    try:
        units.open(10, "runflag.txt", status="unknown")
        if len(words) == 1:
            units.write(10, lw(words[0]))
        else:
            units.write(10, lw(words[0], " ", words[1], " ", words[2]))
        units.close(10)
    except FortranIOError:
        pass


def _pause_on_exit():
    """Shutdown: a script this shell is running is left 'notrunning' so it resumes next time."""
    if commands.runflag == "script" and commands.ncharscript > 0:
        _write_runflag("script", "notrunning", substr(commands.script_fname, 1, commands.ncharscript))


def _get_system_commands():
    """system_commands.f:get_system_commands (the translated routine, when available)."""
    try:
        from ..library import system_commands as scmod
    except ImportError:
        return
    try:
        scmod.get_system_commands()
    except Exception:
        pass


def _keep_viewers():
    """``system_keep_viewers yes`` in system_commands.txt keeps gnuplot/Jmol open at exit."""
    try:
        with open("system_commands.txt") as fh:
            for line in fh:
                w = line.split()
                if len(w) > 1 and w[0] == "system_keep_viewers":
                    return w[1].lower() in ("yes", "y", "true", "1")
    except OSError:
        pass
    return False


def _words():
    return [substr(commands.sometext, commands.ncf[i], commands.ncl[i]) for i in range(1, commands.nwords + 1)]


@routine
def process_command():
    """epsrshell_standalone.f:process_command -- dispatch one command line."""
    first = substr(commands.sometext, commands.ncf[1], commands.ncl[1]).lower()
    commands.sometext = fstr(substr(commands.sometext, 1, commands.ncf[1] - 1) + first
                             + commands.sometext[commands.ncl[1]:], commands.nstsize)
    commands.command_text = fstr(first, 20)
    commands.ilencommand_text = commands.ncl[1] - commands.ncf[1] + 1
    # save the name of the command being executed
    from ..common import symbols

    symbols.symbol_setup = commands.command_text
    symbols.ilensetup = len(commands.command_text)
    ct = commands.command_text
    words = _words()
    if ct in ("help", "?"):
        help()
    elif ct == "system":
        builtins.system(words)
    elif first.startswith("!") and len(first) > 1:
        builtins.system(words, text=commands.sometext[commands.ncf[1]:].strip())
    elif ct == "pwd":
        builtins.pwd(words)
    elif ct == "cd":
        builtins.cd(words)
    elif ct in ("ls", "dir"):
        builtins.ls(words)
    elif ct == "md":
        builtins.md(words)
    elif ct in ("ed", "edit"):
        builtins.edit(words)
    elif ct == "del":
        builtins.delete(words)
    elif ct == "ss":
        start_script(2)
    elif ct == "ps":
        pause_script()
    elif ct == "es":
        end_script()
    elif ct == "jobs":
        builtins.jobs(words)
    elif ct == "trace":
        builtins.trace_cmd(words)
    elif ct in PROGRAM_COMMANDS:
        run_program()
    elif ct in ("bye", "e", "q", "exit", "quit"):
        # if the current program is running the script, set the script status
        # to 'notrunning'. This means when EPSR starts up next time it will automatically
        # start the script running again.
        if commands.runflag == "script":
            commands.script_status = "notrunning"
            _write_runflag("script", commands.script_status, substr(commands.script_fname, 1, commands.ncharscript))
            # close the script file
            units.close(commands.ncommand_stream)
            commands.runflag = "interactive"
        history.save_history_list()
        quit()
    else:
        _write(fmt("(/'process_command> Command \"',a,'\" not recognised.'//'  -- Will attempt to run command anyway')",
                   substr(commands.sometext, commands.ncf[1], commands.ncl[1])))
        # the Fortran tried system_binaries/<command>; here it runs through /bin/sh
        builtins.system(words, text=commands.sometext.strip())


# FORMAT statements of subroutine help, verbatim from the Fortran (command names -> format)
_HELP_FORMATS = {
    ('help', '?'): "(/'Lists the valid commands. Type: help or ?'/)",
    ('cd',): "(/'Sets the working folder. Type: cd <FOLDER_NAME>'/)",
    ('pwd',): "(/'Shows the current working folder. Type: pwd'/)",
    ('ls', 'dir'): "(/'Lists the contents of the working folder.',' Type: ls (or dir) <FILENAME>'/)",
    ('ss',): "(/'Starts an EPSR script file. Type: ss FILENAME'/)",
    ('ps',): "(/'Pauses a script if one is running. Type: ps'/)",
    ('es',): "(/'Ends a script if one is running. Type: es'/)",
    ('epsr',): "(/'Runs EPSR for the specified input files. Type: epsr FILENAME'/)",
    ('plot',): "(/'Sets up the plot file for GNUPlot. Type PLOT <PLOTTYPE>'/)",
    ('bye', 'e'): "(/'Exits from EPSR. Type: bye or e'/)",
    ('randomise',): "(/'Randomizes rotations and positions of molecules.'//'Type randomise <.ato FILENAME>'/)",
    ('fmole',): '(/\'Molecular "shake" to generate molecules. \'//\'Type fmole <.ATO FILENAME> <NTIMES> \'/)',
    ('mixato',): "(/'Mix .ato files. '//'Type mixato <No. of .ATO files to be mixed> '/)",
    ('changeato',): "(/'Change parameters in a .ato file. '//'Type changeato <.ato FILENAME> '/)",
    ('nwts', 'epsrwts'): "(/'Sets up the neutron .wts files from an input filename'//'Type nwts <FILENAME> '/)",
    ('makeato',): "(/'Makes a .ato file (a single molecule or headgroup). '//'Type makeato <.ato FILENAME> '/)",
    ('readgpr',): "(/'Converts a ghemical .gpr file to .mol file. '//'Type readgpr <.gpr FILENAME> '/)",
    ('readjmol',): "(/'Converts a Jmol .jmol file to .mol file. '//'Type readjmol <.jmol FILENAME> '/)",
    ('readmole',): "(/'Converts a ghemical .gpr or Jmol .jmol file to .mol file. '//'Type readmole <FILENAME>, readmole <.EXT> <FILENAME>, where .EXT=.gpr or .jmol '/)",
    ('readxyz',): "(/'Converts a .xyz file to a .ato file. '//'Type readxyz <FILENAME>, readxyz <FILENAME> <vanderWaalsFactor> '/)",
    ('makemole',): "(/'Makes an .ato file for a single molecule from a template .mol file'//'Type makemole <.mol FILENAME> ')",
    ('setup',): "(/'Setup an input file for EPSR or one of the auxiliary programs'//'Type setup <PROGRAM NAME> ')",
    ('bonds',): "(/'Calculate the bond lengths for each molecule in a .ato file.'//'Type bonds <.ATO FILENAME> '/)",
    ('sharm',): "(/'Runs SHARM for the specified input file.'//'Type: sharm <FILENAME>'/)",
    ('sdf',): "(/'Runs SDF for the specified input file.'//'Type: sdf <FILENAME>'/)",
    ('sdfcube',): "(/'Runs SDFCUBE for the specified input file.'//'Type: sdfcube <FILENAME>'/)",
    ('clusters',): "(/'Runs clusters for the specified input file.'//' Type: clusters <FILENAME>'/)",
    ('voids',): "(/'Runs voids for the specified input file.'//' Type: voids <FILENAME>'/)",
    ('plot3d',): "(/'Runs plot3d for the specified input file.'//' Type: plot3d <FILENAME>'/)",
    ('plot3djmol',): "(/'Runs plot3djmol for the specified input file.'//' Type: plot3djmol <FILENAME>'/)",
    ('triangles',): "(/'Runs tri-angles for the specified input file.'//' Type: triangles <FILENAME>'/)",
    ('coord',): "(/'Runs coord for the specified input file.'//' Type: coord <FILENAME>'/)",
    ('chains',): "(/'Runs chains for the specified input file.'//' Type: chains <FILENAME>'/)",
    ('torangles',): "(/'Runs torangles for the specified input file.'//' Type: torangles <FILENAME>'/)",
    ('partials',): "(/'Runs partials for the specified input file.'//' Type: partials <FILENAME>'/)",
    ('md',): "(/'Creates a new directory or folder with specified name.'//' Type: md <FOLDERNAME>'/)",
    ('ed', 'edit'): '(/\'Invokes the editor specified in "system_commands.txt".\'//\' Type: edit <FILENAME>\'/)',
    ('system',): "(/'Invokes the system command in the subsequent words.'//' Type: system <COMMAND> '/)",
    ('del',): "(/'Delete the specified files.'//' Type: del <FILENAME>'/)",
    ('plot2d',): "(/'Runs plot2d for the specified input file.'//' Type: plot2d <FILENAME>'/)",
    ('splot2d',): "(/'Runs splot2d for the specified input file.'//' Type: splot2d <FILENAME>'/)",
    ('plotato',): "(/'Plots the specified .ato file using GNUplot, PGplot or ATOMS.'//' Type: plotato <FILENAME>'/)",
    ('xwts', 'epsrwtsx'): "(/'Sets up X-ray form factor coefficients for a .ato file. '//'Type xwts <FILENAME> '/)",
    ('rings',): "(/'Runs rings for specified input file. '//'Type rings <FILENAME> '/)",
    ('dockato',): "(/'Docks the first .ato file onto the second. '//'Type dockato '/)",
    ('runjmol',): "(/'Runs Jmol.'//'Type runjmol '/)",
    ('readpdb',): "(/'Limited read of .pdb file. Type: readpdb FILENAME'/)",
    ('writexyz',): "(/'Writes xyz for specified ato file. Type: writexyz FILENAME'/)",
    ('mapgrentropy',): "(/'Calculates 6D entropy from .SHARM or .SDF files. Type: mapgrentropy FILENAME'/)",
    ('mapgr',): "(/'Shows g(r) in specified directions (views) from .SHARM or .SDF files. Type: mapgr FILENAME'/)",
}
_F100 = "(/' setup'/' makeato',/' runjmol',/' readgpr',/' readjmol',/' readmole',/' readxyz',/' readpdb',/' makemole',/' fmole',/' bonds',/' mixato',/' addato',/' makelattice',/' makelatticeato',/' mapgr',/' mapgrentropy',/' dockato',/' changeato',/' randomise',/' nwts',/' xwts',//' setup',//' epsr',/' partials',/' triangles',/' torangles',/' chains',/' rings',/' clusters',/' voids',/' coord'/' sharm',/' sdf',/' sdfcube',/' writexyz',//' plot  -  Starts the plot menu',/' plotato  -  Plots the specified .ato file using GNUplot, PGplot, or ATOMS',//' plot2d - plots a 2D density surface using PGPLOT/PGXTAL',/' splot2d - plots a 2D density surface using GNUplot',/' plot3d - plots a 3D density surface using PGPLOT/PGXTAL',/' plot3djmol - Plots a density surface from SHARM/SDF using Jmol',//' system  -  Invoke the specified system command.',/' cd  -  Sets the working folder.',/' pwd  -  Shows the current working folder.'/' ls or dir  -  Lists the contents of the working folder.',/' md  -  Makes a new folder.',/' ed or edit  -  Invokes the specified editor.',/' del  -  Delete a file or files.',/' ss  -  Starts an EPSR script file.',/' ps  -  Pauses a script if one is running.',/' es  -  Ends a script if one is running.',//' bye or e',/)"
_F99 = "(' For more help on a particular command ','type: help command_name.',//' <> means optional parameter')"

_HELP_EXTRA = """
 Also (Python version): exit or quit, !<command> (as system), jobs (running
 child processes), trace on <pattern> | trace off | trace summary (call tracing)."""


@routine
def help():
    """epsrshell_standalone.f:help -- the command list, or help on one command."""
    if commands.nwords < 2:
        _write(lw())
        _write(lw("The list of valid commands is:- "))
        _write(fmt(_F100))
        _write(fmt(_F99))
        _write(_HELP_EXTRA)
        return
    # evaluate which command is being sought
    commands.command_text = fstr(substr(commands.sometext, commands.ncf[2], commands.ncl[2]), 20)
    ct = commands.command_text
    for names, f in _HELP_FORMATS.items():
        if ct in names:
            _write(fmt(f))
            return
    _write(lw("This command is not implemented."))
    _write(lw())
    _write(lw("For a list of valid commands type: help"))
    _write(lw())


@routine
def start_script(nfirst):
    """epsrshell_standalone.f:start_script -- start the script named on the command line
    (nfirst=2) or resume a paused one at startup (nfirst=1)."""
    # first find out whether a script is already running
    check_script_status()
    if commands.script_status != "notrunning":
        _write(fmt("(/'start script> A script may already be running - ','will run in interactive mode')"))
        commands.runflag = "interactive"
        return
    # generate the script file name if one is not defined so far
    if nfirst > 1 and commands.ncharscript == 0:
        get_script_fname(nfirst)
    # if no file was generated exit now
    if commands.ncharscript == 0:
        if nfirst > 1:
            _write(fmt("(/'start script> Script file name not specified')"))
        commands.runflag = "interactive"
        return
    # open the current script file
    commands.ncommand_stream = 19
    try:
        units.open(commands.ncommand_stream, commands.script_fname, status="old")
    except FortranIOError:
        _write(fmt("(/'start_script> Problems with script file ',a)", substr(commands.script_fname, 1, commands.ncharscript)))
        commands.ncharscript = 0
        commands.runflag = "interactive"
        return
    # set the script running
    _write(fmt("(/'start script> Running script file: ',a)", substr(commands.script_fname, 1, commands.ncharscript)))
    commands.script_status = "running"
    # set the runflag, and write it to the runflag.txt file
    commands.runflag = "script"
    commands.runflagglobal = "script"
    _write_runflag("script", commands.script_status, substr(commands.script_fname, 1, commands.ncharscript))


@routine
def check_script_status():
    """epsrshell_standalone.f:check_script_status -- read runflag.txt (the global run status
    of EPSR in this folder)."""
    ncft = [0] * 101
    nclt = [0] * 101
    # first make sure runflag.txt exists - if not create it
    if not os.path.isfile("runflag.txt"):
        _write_runflag("interactive")
        commands.ncharscript = 0
    # now read some text. We expect to find either 1, 2 or 3 words, depending on the value
    # of the first. Any problems with this, then the runflagglobal is set to interactive
    testtext = ""
    try:
        with open("runflag.txt", encoding="latin-1") as fh:
            testtext = fstr(fh.readline().rstrip("\r\n"), 256)
    except OSError:
        pass
    nwordst = parse.parse(testtext, 0, ncft, nclt, 5, 256, 100)
    if nwordst > 0:
        commands.runflagglobal = fstr(substr(testtext, ncft[1], nclt[1]), 20)
        # If runflagglobal is 'exit' then script_status='ps' or 'pausing' will pause the script if runflag='script'
        if commands.runflagglobal == "exit":
            if nwordst > 1:
                commands.script_status = fstr(substr(testtext, ncft[2], nclt[2]), 20)
        elif commands.runflagglobal == "script" and nwordst > 2:
            commands.script_status = fstr(substr(testtext, ncft[2], nclt[2]), 20)
            _, _, commands.script_fname, commands.ncharscript = parse.concatenate(
                3, nwordst, commands.script_fname, commands.ncharscript, nwordst, ncft, nclt, testtext)
        else:
            commands.runflagglobal = "interactive"
            commands.ncharscript = 0
    else:
        commands.runflagglobal = "interactive"
        commands.ncharscript = 0
    # if no script file was set up - for whatever reason, the local run status is set to 'interactive'
    if commands.runflagglobal == "interactive":
        commands.runflag = commands.runflagglobal
        commands.script_status = "notrunning"
        _write_runflag("interactive")


@routine
def get_script_fname(nfirst):
    """epsrshell_standalone.f:get_script_fname -- the script name from the command words."""
    if commands.nwords < nfirst:
        commands.ncharscript = 0
        return
    nlast = commands.nwords
    _, _, commands.script_fname, commands.ncharscript = parse.concatenate(
        nfirst, nlast, commands.script_fname, commands.ncharscript, commands.nwords, commands.ncf, commands.ncl,
        commands.sometext)


@routine
def pause_script():
    """epsrshell_standalone.f:pause_script -- ask a running script to pause."""
    # get the current status of any runs in progress
    check_script_status()
    if commands.script_status != "running":
        _write(lw())
        _write(lw("No need to pause script"))
        return
    commands.script_status = "pausing"
    _write_runflag("script", commands.script_status, substr(commands.script_fname, 1, commands.ncharscript))


@routine
def end_script():
    """epsrshell_standalone.f:end_script -- end the script (it cannot resume)."""
    _write_runflag("interactive")


@routine
def quit():
    """epsrshell_standalone.f:quit -- leave the shell."""
    # if a script is running in the current program, set it to not running before quiting
    if commands.runflag == "script" and commands.script_status == "running":
        pause_script()
    signals.shutdown()
    raise SystemExit(0)


@routine
def run_program():
    """epsrshell_standalone.f:run_program -- run the program named by the command.

    The Fortran built ``<system_binaries>/<command> <directory_name> <command> <args>`` and
    passed it to ``call system``; here the same arguments go to the translated program,
    which runs in this process.
    """
    history.save_history_list()  # Save the current history list for access by other programs
    if commands.nwords <= 0:
        return
    # Get the program to execute. This should always be lower case
    command = commands.command_text.lower()
    commandsuffix = ""
    # For the command 'setup' Vista gets upset
    if command == "setup":
        command = "upset"
    # Translate other commands
    if command == "readgpr":
        command, commandsuffix = "readmole", ".gpr"
    elif command == "readjmol":
        command, commandsuffix = "readmole", ".jmol"
    elif command == "readms":
        command, commandsuffix = "readmole", ".ms"
    elif command == "readmsi":
        command, commandsuffix = "readmole", ".msi"
    elif command == "mapgr":
        command = "mapgrquad"
    # concatenate the all the words into the command line
    text = ""
    nchartext = 0
    if commands.nwords > 1:
        _, _, text, nchartext = parse.concatenate(2, commands.nwords, "", 0, commands.nwords, commands.ncf,
                                                  commands.ncl, commands.sometext)
    dname = substr(directory_name.directory_name, 1, directory_name.nchardir)
    args = [dname, command]
    if commandsuffix:
        args.append(commandsuffix)
    if nchartext > 0:
        args += substr(text, 1, nchartext).split()
    _write(fmt("(/a)", " ".join([command] + args)))
    key = _PROGRAM_KEYS.get(command, command)
    in_script = commands.runflag == "script"
    snap = common.snapshot()
    saved_units = units.save_state()
    interrupted = False
    try:
        with signals.running(key):
            programs.run_program(key, args)
            interrupted = signals.stop_requested()
    except programs.NotTranslated as e:
        _write(f"\nrun_program> {e}")
    except KeyboardInterrupt:
        interrupted = True
        _write(f"\n\nEPSRshell> {key} interrupted. The files it was writing may be incomplete.")
    except FortranStop:
        pass
    except SystemExit:
        pass
    except Exception:
        # the Fortran program would have crashed; the shell carries on
        _write(f"\nrun_program> {key} stopped with an error:")
        traceback.print_exc()
    finally:
        common.restore(snap)
        units.restore_state(saved_units)
    if interrupted and in_script:
        _write("\nEPSRshell> Pausing the script (type 'ss' to resume it).")
        pause_script()
    history.get_history_list()  # Get the current list - possibly modified by program.


def entry():
    """Console-script entry point ``epsrshell``."""
    try:
        status = main()
    except SystemExit as e:
        status = e.code if isinstance(e.code, int) else 0
    finally:
        signals.shutdown()
        processes.terminate_all()
    sys.exit(status or 0)
