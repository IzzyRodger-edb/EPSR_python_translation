# Package G: molfit.f, containerforces_mod.f90, threebody_mod.f90

## Files

| Fortran | Python |
|---|---|
| `molfit.f` (32 routines) | `src/epsr26/library/molfit.py` |
| `containerforces_mod.f90` routines (3) | `src/epsr26/library/containerforces_mod.py` |
| `containerforces_mod.f90` module data | `src/epsr26/common/containerforces_mod.py` |
| `threebody_mod.f90` routines (10) | `src/epsr26/library/threebody_mod.py` |
| `threebody_mod.f90` module data | `src/epsr26/common/threebody_mod.py` |

That is 45 routines, all translated in full.

## Verification (done; scripts in the session scratchpad, not checked in)

Each check builds the same input in Fortran and Python and compares the output as hex
bit patterns (`z8`):

1. `rotatm`, `movem` (including the doubling loop, from a 1e-12 step), `en_atom_inter_scale`,
   `en_dihedral`, `en_atom_inter_scale_coord`, `en_gauss`/`en_harmonic`, `en_morse*`, all
   extracted from molfit.f and compiled with gfortran using the reference flags (non-cubic
   box): **bit-identical**.
2. `rotati`, `en_atom_intra` (bonds, angles, image atoms, dihedral exclusion,
   refine_intra), `en_mole_intra`, `en_mole_dihedral`, `en_change_tethered`,
   `check_tethered`: **bit-identical**.
3. `find_triplets` (stdout too), `en_mole_triplets`, `pr_mole_triplets`,
   `truncate_threeb`, `mod_morse*`, `grad_mod_morse1`, `get_triplet_pot_cons` (stdout),
   `get_containerforce_list`, `en_containerforce` (all 4 branches), `dist_containerforce`,
   against the real module sources: **bit-identical, stdout identical**.
4. **8 full `molfit` cycles** on 27 three-atom molecules with bonds, angles and a ROT
   group each (shake branch, internal rotations every `rotfreq`, whole-molecule moves,
   step adaptation, ergpair bookkeeping), linking the reference `fortran_build/*.o`
   against the Python molfit plus packages D/E's ergpairreference, cellneighbourlist and
   centreofmassi: every coordinate, step size, move counter, fadds/faddbs and every stored
   pair energy is **bit-identical**.

## Decisions

- **Reallocation**: I use package A's `reallocation_routines.reallocate*` (each returns the
  new array, see SIGNATURE_EXCEPTIONS.md). I wrote no private reallocation helper.
- **Module-data imports**: where a data module and a library module have the same name I
  alias the data module: `from ..common import threebody_mod as td` (and `tb_d`, `cf`/`cf_d`,
  `erg_d`, `cc_d`, `cnl`). Data is always read through the module at call time.
- **Speed**: COMMON arrays (never reallocated) and module arrays are aliased into locals at
  the start of each routine. `ergpairreference.ergold/ergnew/nbref` are aliased only
  inside `en_change`, because `make_neighbour_list_from_mols` can reallocate them.
  `molfit` re-reads them through the module after every call to that routine.
- **`check_coord_bonded` short-circuit** (`en_change`, `en_mole_pair`): it is only called
  when `ncoordconstraints > 0`. With 0 constraints the Fortran function returns `.false.`
  with no side effects, so the result is the same; the Rust does the same.
- **Constant folding**: `exp(15.0)` (the `limit` parameter) and `exp(30.0)` are module
  constants `f32(math.exp(x))`. gfortran folds them at compile time, and these values
  equal glibc `expf`.
- **Timers**: `call system_clock(time1)` with the integer(8) `timer.time1`. gfortran's
  count is CLOCK_MONOTONIC in ns (rate 1e9, checked with gfortran 15), so the private helper
  `_system_clock8()` is `time.monotonic_ns()`. `sumtime1/2` are wall-clock values, so they are
  never byte-identical (stdout line 1999 and unit 35 in runepsr17).
- `xs, ys, zs, tol` are COMMON (`molecom`), not locals: `movem` and `molfit` set
  `molecom.xs` etc., and `check_tethered` reads `molecom.tol`.

## Bugs and quirks reproduced

- **FORTRAN_BUGS.md section 2**: whole-molecule translation step doubled.
  `wholetrans[mt] = movem(im, wholetrans[mt])` stores the doubled `stepsz` back,
  marked `NB (FORTRAN_BUGS.md section 2, reproduced)`.
- **FORTRAN_BUGS.md 1.1**: none of the fix belongs in the threebody module routines. The
  Rust fix is entirely in `write_pcof_threebody` (the save path of `read_and_save_files`,
  package B), so there is no `# FIX` site here.
- **`get_triplet_pot_cons`** subtracts cosines from radians (FORTRAN_BUGS.md section 3). It
  only affects a printed value; reproduced.
- **`find_triplets`** prints `ltyp1` on all three `size(atomNlist)` lines; reproduced.
  `size()` of an unallocated list (a zero count) is written as 0, which is undefined in
  the Fortran.
- **`en_mole_triplets_old`** assigns a *local* `en_mole_triplets` and never sets its own
  result; it returns 0.0 (NB). It is never called.
- **`get_containerforce_list`**: `typecentral/typesecondary` are uninitialised when no
  type matches; 0 is used (NB, as the Rust does).
- **`find_triplets`**: `ergold_triplet/ergnew_triplet` are freshly allocated
  (uninitialised) and set to zero here.
- **`molfit`**: the locals `sumevib, sumeang, sumedih, actualstep` are uninitialised in the
  Fortran, but they are always set before use.

## Contract issue (already in SIGNATURE_EXCEPTIONS.md, by package F)

`en_mole_intra` returns `(en_mole_intra, sumedih)` as SIGNATURES.md says. The
Fortran also accumulates `sumevib`/`sumeang` (through `en_atom_intra`), and with this contract
those values are lost. Inside this package that only affects `en_change_intra(_rot)`'s
returned `sumevib/sumeang`, which only commented-out writes use. `adjust_intra_wts`
(package F) prints them, so that line will show zeros until both sides change the contract
to `(en_mole_intra, sumevib, sumeang, sumedih)`. To fix it, change `en_mole_intra`'s return and its
three call sites (`en_change_intra`, `en_change_intra_rot`, and the callers in
ato_manipulation_rev4).

## For the integrator

- `library/molfit.py` imports `ato_manipulation_rev4` (for `centreofmassi`), which imports
  `ato_read_and_save` (package D). Until that file exists, `import epsr26.library.molfit`
  fails, but only because of that missing sibling. With a stub module it imports and runs.
  The other four modules import cleanly now.
- Callers pass `nisum` as the COMMON variable (`epsrinput_par.nisum` in runepsr17) and must
  assign the returned value back.
- `ruff check --select F` is clean on all five files. Unused locals that are kept for
  fidelity (`distmax12`, `d12star`, ...) carry `# noqa: F841`.
