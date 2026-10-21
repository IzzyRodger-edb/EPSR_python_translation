# Package A: small shared libraries (translation notes)

## Files

| Fortran | Python |
|---|---|
| lngamma_mod.f90 | `library/lngamma_mod.py` (routines), `common/lngamma_mod.py` (data: `mgamma`, `gamln`) |
| reallocation_routines.f90 | `library/reallocation_routines.py` (no module data) |
| molecom_mod.f90 | `library/molecom_mod.py` (routines), `common/molecom_mod.py` (data) |
| stringutility.f90 (module String_Utility) | `library/stringutility.py` (`strupcase`, `strlowcase`; LOWER_CASE/UPPER_CASE as module constants) |
| minimagec.f | `library/minimagec.py` |
| bondnumbers.f | `library/bondnumbers.py` |
| system_commands.f | `library/system_commands.py` |
| FCreadline_windows.f90 (module jmk_readline_module) | `library/FCreadline_windows.py` |
| pfroutines_mod.f90 (module pfroutines) | `library/pfroutines_mod.py` (routines), `common/pfroutines.py` (data: `rratio`, `intenprobfunction`) |
| get_molecules.f | `library/get_molecules.py` |
| file_and_directory_commands.f | `library/file_and_directory_commands.py` |

43 routines in all (2 + 11 + 2 + 2 + 2 + 1 + 1 + 4 + 3 + 1 + 14), every one translated completely.

## What integrators must know

1. **reallocate\* return the array.** NumPy cannot grow an array in place, so every
   `reallocateNd_x` returns the reallocated array (`None` when the Fortran leaves it
   deallocated: a non-positive new size). Write
   `molecom_mod.labrot = reallocation_routines.reallocate1d_i(molecom_mod.labrot, n)`.
   SIGNATURES.md lists no return value (logged in SIGNATURE_EXCEPTIONS.md). The same
   object comes back when the size is unchanged (the Fortran returns early). New elements
   are 0 / blank / `.false.` (the Fortran leaves them uninitialised). `size(array,d)` is
   `shape[d]-1`, so do not reallocate arrays that were allocated with a lower bound of 0.
   `reallocate3d_i` allocates even if a new extent is <= 0 (zero-sized but allocated, as in
   the Fortran); the others leave the array deallocated in that case.
2. **Name clash library vs data module.** `library/molecom_mod.py` and
   `common/molecom_mod.py` share a name (likewise `lngamma_mod`). A file that needs both
   must alias one, e.g. `from . import molecom_mod as molecom_mod_routines` next to
   `from ..common import molecom_mod`. `pfroutines` has no clash (library file
   `pfroutines_mod.py`, data `common/pfroutines.py`).
3. **Two `lngamma`s.** `lngamma_mod.lngamma(ic)` (real*8 table lookup of log(ic!)) and
   `pfroutines_mod.lngamma(argin)` (real*4, uses `gamma`). Check the link list for which
   one a program uses; inside pfroutines the module's own one is called.
4. **pfroutines arrays are (0:nr)**: created with `real_array((0, nr))`; Fortran
   `size(rratio)` is `rratio.shape[0]` (not `-1`).
5. **Imports of other packages' modules are inside the functions** (`auxiliary_routines`
   in pfroutines_mod.getrratio, get_nterms, setup_input_file; `ato_read_and_save` and
   `ato_manipulation_rev4` in get_molecules), so these modules import before those files
   exist. They are still called through the module (`auxiliary_routines.ran1()`), so
   tracing works.
6. `minimagec`/`minimage` do not convert their arguments: pass `numpy.float32` values
   (array elements of real arrays, or f32 expressions). Returns `(x, y, z)`.
7. `testsize_rotgrps(newsizein=None)`, `testsize_labrot(newsizein=None)`: `optional`
   argument as a keyword defaulting to None.
8. `search_for_files(ext, ilenext, iaccept, len_iaccept=1)` and
   `readline(line, prompt, len_line=256)`: extra keyword = declared length of the
   `character(len=*)` output dummy. `search_for_files` and `setup_input_file` leave the
   selected name in `commands.sometext`, as in the Fortran.
9. `setup_input_file` returns `(nout, ierrout)` as SIGNATURES says, although nout is
   never assigned.

## Uninitialised variables / quirks reproduced

- get_molecules: `vanderWaalsRadius(melements)` is an uninitialised local; elements not
  read from vanderWaalsRadii.txt are 0 here. `ir` (read with list-directed input before
  any assignment) starts at 0 (only matters for a null value such as `,`).
- get_molecules reads `vanderWaalsFactor` and the type pairs from unit 5 without `end=`:
  end of input raises `FortranEOF` (the Fortran stops with a runtime error). Same for
  `read(5,101) iaccept` in search_for_files and `read(5,101) sometext` in make_directory.
- get_system_commands: a failed internal read of `gnu_pointtype`/`gnu_pointsize` sets
  `ierr` (iostat) and so ends the reading of system_commands.txt, as in the Fortran.
- get_command_from_args parses the character*256 argument with `nchartext = nstsize`
  (3020). Fine unless an argument is longer than ~250 characters (the Fortran would read
  past `text`; parse.py sees blanks).
- `read(u,'(a)',iostat=ierr) sometext` at end of file leaves sometext unchanged (helper
  `_read_a`), which the loops in get_nterms/search_for_files rely on.
- set_directory_name: `'.'` inside a path (`deep/./x`) resets to the home folder, as in
  the Fortran (checked).

## Compatibility-layer gaps (helpers written locally)

- `fformat.fmt` does not support the `$` edit descriptor: the formats with `,$` are
  written without it and output with `units.write(..., advance=False)`; the original
  format is kept as a comment.
- `fmt('(a,1x)', p)` drops the trailing `1x`, but gfortran writes the blank on a
  non-advancing write (checked: `'EPSR>'` -> `'EPSR> '`). `readline` appends the blank
  by hand. It would be better if fmt kept trailing X on non-advancing output.
- No libm `tgamma` in fmath: `pfroutines_mod._tgamma` calls glibc through ctypes
  (Python's `math.gamma` is CPython's own implementation and need not be bit-identical).
- No `getcwd`: `file_and_directory_commands._getcwd` (os.getcwd, status 0).
- `_feq`/`_fne` for Fortran character comparison of blank-padded substrings, and
  `_read_a` for `read(u,'(a)',iostat=)`; candidates for `fortran/text.py`/`units.py`.

## Checks done (gfortran -O2, byte comparison of stdout)

- pfroutines: getpfprobability/lngamma/getrratio for 5 values of rfract (stub ran1): identical.
- get_system_commands on EPSR-test/system_commands.txt (+ a gnu_pointsize line) and
  minimagec on a non-cubic box (20 points): identical.
- get_command_from_args (with arguments), set_directory_name (sub, a/../a/./x, .., home),
  get_full_filename, get_nterms (existing and missing file), search_for_symbol,
  setup_input_file (file found / run-name mismatch / missing): identical.
- setup_input_file -> search_for_files interactive: accept, `u` backspacing past the
  start of the list, `e`, typed name, and end of input: identical.
- get_molecules on an 8-atom system with stubbed normalise/centreofmass: identical
  (merging, bond lists, renumbering, moltype list, all output).
- bondnumbers, testsize_rotgrps/labrot, lngamma_mod.lngamma, StrUpCase/StrLowCase: identical.
- reallocation routines: Python unit checks (grow/shrink/unchanged/None/zero size).
