"""Translation of fortran_src/plot.f: ``program standalone`` (the EPSR ``plot`` program).

The plot command itself is ``plot_routines.plot_new``.  With the default plot symbols
the Fortran aborts in a list-directed read (see docs/translation_notes/H.md); the
Python does the same.
"""

from ..library import file_and_directory_commands, plot_routines


def main():
    """plot.f:standalone -- program to perform the EPSR plot command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    plot_routines.plot_new()
