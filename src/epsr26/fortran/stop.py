"""Fortran ``stop`` and I/O error conditions as Python exceptions."""


class FortranStop(Exception):
    """``stop [code]`` / ``call exit(code)``: ends the current program.

    The shell (and ``programs.run_program``) catches it; ``code`` is the exit status.
    """

    def __init__(self, code=0, message=""):
        super().__init__(message or f"STOP {code}")
        self.code = code
        self.message = message


class FortranEOF(Exception):
    """End of file in a READ (``end=`` branch, ``iostat < 0``)."""

    iostat = -1


class FortranIOError(Exception):
    """Error in an OPEN/READ/WRITE (``err=`` branch, ``iostat > 0``)."""

    def __init__(self, iostat=5010, message=""):
        super().__init__(message)
        self.iostat = iostat


def iostat_of(exc):
    """The ``iostat=`` value for an I/O exception: -1 at end of file, >0 on error.

        ios = 0
        try:
            x = units.list_read(10).real(x)
        except (FortranEOF, FortranIOError) as e:
            ios = iostat_of(e)
    """
    if isinstance(exc, FortranEOF):
        return -1
    return getattr(exc, "iostat", 5010) or 5010
