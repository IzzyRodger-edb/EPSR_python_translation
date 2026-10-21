"""Translation of fortran_src/makemole.f: ``program standalone`` running ``makemole``.

The builder itself (``subroutine makemole``) is in ato_manipulation_rev4.f and is
translated in library/ato_manipulation_rev4.py.
"""

from ..library import ato_manipulation_rev4
from ..library import file_and_directory_commands


def main():
    """makemole.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    ato_manipulation_rev4.makemole()
