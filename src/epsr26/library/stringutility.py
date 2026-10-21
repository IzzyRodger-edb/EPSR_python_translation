"""Translation of fortran_src/stringutility.f90: ``module String_Utility``.

``StrUpCase``/``StrLowCase`` are ``strupcase``/``strlowcase`` (lower case, as in
SIGNATURES.md).  The result has ``len(Input_String)``; character values are stored
without trailing blanks, which case conversion never touches, so the Python string
returned is the Fortran result without its trailing blanks.
"""

from ..trace import routine

# CHARACTER( * ), PRIVATE, PARAMETER
LOWER_CASE = "abcdefghijklmnopqrstuvwxyz"
UPPER_CASE = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


@routine
def strupcase(input_string):
    """stringutility.f90:StrUpCase -- copy of input_string with a-z made upper case.

    Returns output_string.
    """
    # -- Copy input string
    output_string = list(input_string)
    # -- Loop over string elements
    for i in range(1, len(output_string) + 1):
        # -- Find location of letter in lower case constant string
        n = LOWER_CASE.find(output_string[i - 1]) + 1
        # -- If current substring is a lower case letter, make it upper case
        if n != 0:
            output_string[i - 1] = UPPER_CASE[n - 1]
    return "".join(output_string)


@routine
def strlowcase(input_string):
    """stringutility.f90:StrLowCase -- copy of input_string with A-Z made lower case.

    Returns output_string.
    """
    # -- Copy input string
    output_string = list(input_string)
    # -- Loop over string elements
    for i in range(1, len(output_string) + 1):
        # -- Find location of letter in upper case constant string
        n = UPPER_CASE.find(output_string[i - 1]) + 1
        # -- If current substring is an upper case letter, make it lower case
        if n != 0:
            output_string[i - 1] = LOWER_CASE[n - 1]
    return "".join(output_string)
