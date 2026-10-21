# Package D: .ato read/write, EPSR symbol table, cell neighbour list, pair energies

## Files

| Fortran | Python |
|---|---|
| ato_read_and_save.f | `library/ato_read_and_save.py`: `minitc`, `normalise_coordinates`, `unnormalise_coordinates`, `moloutc` |
| setup_epsr_rev1.f | `library/setup_epsr_rev1.py`: `read_epsr16`, `read_epsr17`, `setup_symbols_epsr` |
| cellneighbourlist.f, module routines | `library/cellneighbourlist.py`: `initialise_cellneighbourlist`, `makecellneighbourlist`, `make_neighbour_list_from_cells`, `get_nnbrefmax`, `findiaibic`, `getcellreference`, `findperiodiccell` |
| cellneighbourlist.f, free routines | same file: `cellcoord`, `atomcellno`, `get_date_and_time`, `get_daysthisyearsofar`, `get_currentmonth` |
| cellneighbourlist.f, module data | `common/cellneighbourlist.py` |
| cellneighbourlist.f, `common /time/` in get_date_and_time | `common/time.py` (startyear, startabsolutetime, startmusecs) |
| ergpairreference.f, modules ergpairstore + ergpairreference | `library/ergpairreference.py` (all 10 routines) |
| ergpairreference.f, module data | `common/ergpairstore.py`, `common/ergpairreference.py` |

26 routines, all translated in full: 4 + 3 + 12 + 10, counting the second `save_energies`.

## What integrators must know

1. **Library and data modules have the same name** (`cellneighbourlist`,
   `ergpairreference`). Inside the library files the data is imported as `cnl` / `erd` /
   `eps`. A caller that needs both should alias one, e.g.
   `from ..common import cellneighbourlist as cnl` next to `from . import cellneighbourlist`.
   The same goes for `common/time.py`: import it as `from ..common import time as time_common`.
   It does not shadow the stdlib `time` because imports are absolute.
2. **Two `save_energies`.** `ergpairreference.save_energies()` is the one in module
   ergpairreference, which every program uses. The one in module ergpairstore is
   `ergpairreference.ergpairstore_save_energies()`. Logged in SIGNATURE_EXCEPTIONS.md.
3. **Imports of other packages' modules are made inside the routines**: auxiliary_routines,
   ato_manipulation_rev4, bondnumbers, molecom_mod (data and routines), reallocation_routines,
   threebody_mod, coordconstraints_mod, file_and_directory_commands, read_and_save_files.
   This avoids the ato_manipulation_rev4 <-> ato_read_and_save cycle. Calls still go through
   the module attribute, so tracing works. `initialise_ergpair` imports `cellneighbourlist`
   inside the function as well.
4. **reallocate1d_r** follows package A's convention: it returns the array, so I write
   `molecom_mod.scalefacmol = reallocation_routines.reallocate1d_r(molecom_mod.scalefacmol, n)`.
   The Fortran `scalefacmol=1.0` after it is `[:] = 1.0`. That fails if the new size is <= 0,
   because A then returns None. Fortran assignment to an unallocated array is undefined too.
5. **coordconstraints_mod** (package H): I use `coordconstraints_mod.coordconstraints[ir]`
   (a 1-based list) with attributes `.ns`, `.ds[j]`, `.dsmax[j]`, as H wrote them.
6. **minitc** returns `(rmax, cell, upol, rho, amubar, amubar2, ndist, ierr)`. `amubar2` is
   only accumulated (`amubar2 = amubar2 + amud*amud`, never zeroed), exactly as in the
   Fortran, so the caller's incoming value matters. Pass a float32. The Rust zeroes it in the
   `nbar == 0` branch, but the Fortran doesn't, and the Python follows the Fortran.
7. **cosd/sind**: minitc declares `external cosd,sind`, so it uses `auxiliary_routines.cosd/sind`.
   moloutc does not declare them, so it uses gfortran's intrinsics (`cosd_intrinsic`,
   `sind_intrinsic`), as in the Rust port.
8. **moloutc with natomsphere unallocated**: the Fortran reads `natomsphere(itype(ia))` in the
   atom loop *before* its "allocate if not allocated" block. gfortran crashes with SIGSEGV. The
   Python raises TypeError at the same point. Programs that never call ctype hit this, as the
   Rust notes.
9. **read_epsr16 / read_epsr17 have no callers** in any Fortran file (dead code). They are
   translated in full anyway. `nout` is returned as SIGNATURES says. `ifnameato`, `ifnamepcof`,
   `npairspcof` and `ntype` are dead locals there, because the include files that hold them
   are not included.
10. `get_date_and_time` uses `runtime.date_and_time()` values; `get_daysthisyearsofar`
    returns `(days, latestmonth)` per SIGNATURES (the dummy may be clipped to 13).
11. Explicit imports from `..fortran` are used instead of `import *`, so that
    `ruff --select F` is clean.

## Uninitialised variables and choices

- minitc: the implicit real locals that are read with `iostat=` start at f32(0) once per call
  and keep their value across the atom-type loop, as Fortran stack locals do. These are
  `epsit sigt amass qt rqt natomsph radminsph radmaxsph pseudofacsph`, plus `imolread`.
  The choice only matters when a read fails.
- minitc `iexch`, read without iostat, and read_epsr `nq`: a null value keeps 0.
- `get_currentmonth`: `month` is blank when values is outside 1..12 (it is uninitialised in
  the Fortran).
- `make_neighbour_list_from_cells` / `get_nnbrefmax`: the uninitialised `icellref` passed
  to getcellreference is 0. It is overwritten anyway.

## Verified

I built a Fortran driver (`drv.f`, scratch only) against the reference objects in
`fortran_build/`, from the `epsr` link list without `epsr_standalone_rev3.o`, plus a
matching Python driver. Each run does minitc -> normalise_coordinates ->
initialise_cellneighbourlist -> get_nnbrefmax -> make_neighbour_list_from_cells ->
initialise_ergpair / set_ergpair / get_neighbsum / make_neighbour_list_from_mols ->
unnormalise_coordinates -> moloutc -> setup_symbols_epsr, and dumps the full symbol table.

- **All 137 .ato files in the repository** (EPSR-test, test_problems, validation; this includes
  5 non-cubic boxes, a DIH group and tethered molecules): the written .ato file and the whole
  stdout are **byte-identical** to the Fortran. That covers the saved ran1 state, the minitc
  and cell-list messages, and the `write(6,*)` values.
- Plain minitc+moloutc round trip of EPSR-test/200DMSO400Water.ato reproduces the input
  file byte for byte.
- read_epsr17 and read_epsr16 with synthetic old-format .inp files, and the "file not found"
  branch: stdout and all symbol values are identical to the Fortran. The only difference is
  my driver printing the never-set `fname_pcof`: NUL bytes in Fortran static storage, blank in
  Python. That is a convention of the common modules, not of this translation.

Not verified: the coordination-constraint and three-body branches of
(un)normalise_coordinates (no test file defines them), and get_date_and_time against the
Fortran (it prints wall-clock time).
