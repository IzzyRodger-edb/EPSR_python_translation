"""Translation of fortran_src/FCreadline_windows.f90: ``module jmk_readline_module``.

The Windows stand-in for GNU readline: ``readline`` writes the prompt without advancing
and reads one line from standard input.  The history routines are empty, as in the
Fortran.  (The shell will later replace ``readline`` with GNU readline.)
"""

from ..fortran import fmt, fstr, units
from ..trace import routine


@routine
def readline(line, prompt, len_line=256):
    """FCreadline_windows.f90:readline -- prompt and read one line from unit 5.

    ``len_line`` is the declared length of the caller's ``line`` (Fortran
    ``character(len=*)``).
    Returns line.
    """
    # write(*,'(a,1x)',advance='no') prompt(1:len_trim(prompt))
    # NB: fformat.fmt drops a trailing 1x at the end of the record, but gfortran writes the
    # blank on a non-advancing WRITE (checked with gfortran 'EPSR>' -> 'EPSR> '), so the
    # blank is added here explicitly.
    units.write(6, fmt("(a)", prompt.rstrip(" ")) + " ", advance=False)
    # read(*,'(a)') line
    line = fstr(units.read_line(5), len_line)
    return line


@routine
def show_history_list():
    """FCreadline_windows.f90:show_history_list -- empty."""
    return


@routine
def save_history_list():
    """FCreadline_windows.f90:save_history_list -- empty."""
    return


@routine
def get_history_list():
    """FCreadline_windows.f90:get_history_list -- empty."""
    return

# the test program (commented out in the Fortran)
# program test_program
#   use jmk_readline_module
#   character(len=256):: line
#   call get_history_list()
#   write(*,*) ... ' Enter text and then edit it. "q" quits; "h" displays history:'
#   do
#     call readline(line,'EPSRshell>') ! read editable input line
#     if (line=='q') exit
#   enddo
#   call save_history_list()
# end program test_program
