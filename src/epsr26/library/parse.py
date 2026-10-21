"""Translation of fortran_src/parse.f: EPSR's word splitter and word joiner."""

from ..fortran import fstr, substr
from ..trace import routine


@routine
def parse(sometext, nwords, ncf, ncl, nspcword, nchartext, nworddim):
    """parse.f:parse -- split sometext into blank-separated words.

    ncf(i), ncl(i) receive the first and last character numbers of word i.  The search
    stops once more than nspcword consecutive spaces are found, which is how EPSR skips
    the free text that follows each value in its input files.

    Returns nwords (ncf, ncl are filled in place).
    """
    lt = len(sometext)

    def ch(i):  # sometext(i:i); blank beyond the stored (unpadded) text
        return sometext[i - 1] if i <= lt else " "

    # First bypass any blank space before the first non-space character
    nspace = 0
    nwords = 0
    i = 1
    while ch(i) == " " and i < nchartext:
        i = i + 1
    nfirst = i
    nlast = i
    while i < nchartext and nspace <= nspcword:
        i = i + 1
        # Check whether the current letter is a space
        if ch(i) == " ":
            # A space is found. Check whether this is the end of a word
            if nlast >= nfirst:
                # New word found
                if nwords < nworddim:
                    nwords = nwords + 1
                    ncf[nwords] = nfirst
                    ncl[nwords] = nlast
                nspace = 1
            else:
                nspace = nspace + 1
            nfirst = i + 1
        else:
            nlast = i
    # Get the final word in case the text does not end with a space
    # First check whether the last letter is a space
    if ch(i) != " ":
        if nlast >= nfirst:
            # New word found
            if nwords < nworddim:
                nwords = nwords + 1
                ncf[nwords] = nfirst
                ncl[nwords] = nlast
    return nwords


@routine
def concatenate(nfirst, nlast, return_string, ncharstring, nwords, ncf, ncl, input_string, len_return_string=256):
    """parse.f:concatenate -- join words nfirst..nlast with single spaces (max 256 chars).

    ``len_return_string`` is the declared length of the caller's return_string
    (Fortran ``len(return_string)``).

    Returns (nfirst, nlast, return_string, ncharstring): the routine also clips the
    caller's nfirst and nlast.
    """
    output_string = ""  # character(len=256)
    # maximum number of characters in a word
    ncharword = 256
    ncharstring = 0
    if nlast > nwords:
        nlast = nwords
    if nfirst < 1:
        nfirst = 1
    # set a flag to say whether the character limit has been reached
    iflag = 0  # noqa: F841
    iw = nfirst
    while iw <= nlast and ncharstring < ncharword:
        ncharadd = ncl[iw] - ncf[iw] + 1
        nchartest = ncharstring + ncharadd
        # Need to allow for the space between words if there is more than one word in the string
        if iw > nfirst:
            nchartest = nchartest + 1
        if nchartest > ncharword:
            ncharadd = ncharword - ncharstring
            # Again need to allow for space between words if there is more than one word already
            if iw > nfirst:
                ncharadd = ncharadd - 1
            iflag = 1  # noqa: F841
        if iw == nfirst:
            output_string = fstr(substr(input_string, ncf[iw], min(ncl[iw], ncf[iw] + ncharadd - 1)), 256)
            ncharstring = ncharadd
        else:
            # test the final number of characters in the text
            if ncharadd > 0:
                output_string = fstr(
                    substr(output_string, 1, ncharstring) + " "
                    + substr(input_string, ncf[iw], min(ncl[iw], ncf[iw] + ncharadd - 1)), 256)
            ncharstring = ncharstring + ncharadd + 1
        iw = iw + 1
    # Make we don't exceed the length of return_string
    if ncharstring > len_return_string:
        ncharstring = len_return_string
    return_string = fstr(substr(output_string, 1, ncharstring), len_return_string)
    return nfirst, nlast, return_string, int(ncharstring)
