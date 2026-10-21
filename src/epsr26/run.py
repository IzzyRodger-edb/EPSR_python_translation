"""Run one EPSR program outside the shell, exactly as the Fortran binary would be run.

    python -m epsr26.run <program> <arguments...>

is the counterpart of ``fortran_build/<program> <arguments...>``; e.g.
``python -m epsr26.run epsr ./ epsr water.EPSR.inp``.  Used by the validation harness.
"""

import sys

from . import trace
from .programs import run_program


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        sys.stderr.write(__doc__)
        return 2
    trace.enable_from_environment()
    try:
        status = run_program(argv[0], argv[1:])
    finally:
        if trace.is_enabled():
            import os

            if os.environ.get("EPSR_TRACE") == "summary":
                trace.summary(sys.stderr)
            trace.disable()
    return status


if __name__ == "__main__":
    sys.exit(main())
