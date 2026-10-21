# Package F: ato_manipulation_rev4.f lines 2344-end

Output: `src/epsr26/library/ato_manipulation_rev4_part2.py`, to be merged into
`library/ato_manipulation_rev4.py` (package E). 23 routines, all translated in full, in
Fortran order: fmolec1, get_direct_bonds, get_rotgrp, get_whole_rotgrp, rotatecoords,
dihedral, generate_mp, crystfrac, crystarg, bragginten, radd, sharm_diffuse, jhalf,
truncconv, change_boxsize, atogen, press, stog, lorch, abs_energy_ep, abs_energy_ep_old,
adjust_intra_wts, get_maxmolsize.

## Merging

- Part-1 names used by bare name (ruff F821 until the merge, each marked
  `# noqa: F821 (part 1)`): `en_tethered` (fmolec1), `generate_fractional_coordinates`,
  `generate_cartesian_coordinates` (change_boxsize), `initialise_diffuse_values`
  (sharm_diffuse). With those noqa markers `ruff check --select F` is clean.
- Imports to carry into the merged module. Library modules whose name is also a
  common data module are imported with a `_lib` suffix:
  `from . import cellneighbourlist as cellneighbourlist_lib`, `containerforces_mod as
  containerforces_mod_lib`, `coordconstraint_mod as coordconstraint_mod_lib`,
  `ergpairreference as ergpairreference_lib`, `threebody_mod as threebody_mod_lib`; plus
  `auxiliary_routines, minimagec, molfit, reallocation_routines`. `read_and_save_files`
  (outgofzs_sharm) is imported inside sharm_diffuse (it did not exist yet and may form a
  cycle).
- Import check: `molfit` imports part 1, which imports `ato_read_and_save` (package D,
  missing when this was written). With an in-memory stub for that module the file imports
  and a synthetic smoke test of atogen, press, generate_mp, stog, lorch, truncconv,
  get_whole_rotgrp runs with float32 results throughout.
- Private helpers (the compat layer lacks them): `_omp_get_max_threads()`,
  `_omp_get_num_threads()` (both 1), `_atan2d_intrinsic(y, x)` = gfortran's ATAN2D for
  real*4, i.e. `atan2f(y,x) * (r4(180)/PI4ATAN)` (= 57.2957764). I checked this against
  gfortran 15 on 2000 points: it is NOT the correctly rounded 180/pi (57.29578). It would
  be worth adding `atan2d_intrinsic` next to `sind_intrinsic`/`cosd_intrinsic` in fmath.
  `_intenhkl_ref` (below).

## Quirks reproduced / decisions

- **fmolec1 store loop starts at inn=2** (FORTRAN_BUGS.md section 2): reproduced, marked NB.
- fmolec1 OpenMP: `num_threads_input <= 0` gives `num_threads = 1/2 = 0` with
  OMP_NUM_THREADS=1, and `num_threads_input` is returned as 0. That is what the Fortran
  prints (`0 1 1`). OMP_SET_NUM_THREADS is a no-op. The parallel loops run in serial order.
- fmolec1: sumevib/sumeang/sumedih are uninitialised locals, only accumulated, never
  used: they start at 0.0.
- `en_mole_intra` follows SIGNATURES.md (returns energy, sumedih), so
  **adjust_intra_wts prints 0.0 for sumevib and sumeang**, and its ratios are
  1.0/1.0. Stdout only. See SIGNATURE_EXCEPTIONS.md; fix both sides together.
- get_rotgrp/get_whole_rotgrp return `(labrott, nrott)`. labrott is
  `allocatable, intent(out)`, so it is set to None on entry and grown with
  `reallocation_routines.reallocate1d_i` (which returns the array, package A's
  convention). This matches E's makemole calls. Logged in SIGNATURE_EXCEPTIONS.md.
- generate_mp: `it`, `jt` passed to get_containerforce_list are never set (uninitialised).
  0 is used, so the container-force lists come out empty (FORTRAN_BUGS.md section 2, Rust
  potential.rs:561). The unused locals ratt, summax, nterms, decayp are still computed.
- bragginten: more than nhklmax reflections would write out of bounds in the Fortran.
  Here NumPy raises IndexError (Rust: `unsupported`). The `write(6,*) ib,hi,ki,li` debug
  line per reflection is kept. At the end the allocatables are set to None as the
  Fortran deallocates them; if diffusedefined, hihkl/kihkl/lihkl/qhkl/intenhkl are kept
  for sharm_diffuse.
- **sharm_diffuse**: the run-time FORMAT is the reference build's patched one,
  `'(1h(,'//il1fmt//...` (tools/build_fortran_reference.sh). The unpatched
  `'(1h('//il1fmt` is invalid. I checked the title strings against gfortran
  (`(12,0,12;3,0)`).
- sharm_diffuse indexes `intenhkl(nhkl,npairs)` as `intenhkl(ip,ib)`. `_intenhkl_ref`
  maps it to the element gfortran addresses, `(ip-1)+(ib-1)*nhkl` in column-major order.
  When nhkl > npairs that runs past the end of the array (undefined behaviour), so the
  routine raises RuntimeError up front, as the Rust does. In practice it only works for
  nhkl <= npairs.
- sharm_diffuse: cosd, sind and atan2d are gfortran intrinsics there (not declared
  external), and so are cosd/sind in change_boxsize. `intenfile` sums with weightx/cwt
  (real*8) are done in double and stored as real. Not reproducible: sharm_diffuse
  deallocates qhkl but not hihkl, so a later bragginten call would fail its
  `deallocate(hihkl,kihkl,lihkl,qhkl)` in the Fortran. Python does not fail there.
- press: prpot is real*8. It is the real*4 product converted on assignment, and
  `ergnew = ergnew + prwt + prpot` is (real4 sum) + double, stored as real. presst and
  pressthreeb are real*8. pressc and prthreebout follow the Rust casts. prthreebout is
  returned unchanged when nthreeb = 0. Speed shortcut: `check_coord_bonded` is skipped when
  ncoordconstraints == 0, where it returns .false. without looking at anything (same
  result).
- press `do if=1,nfile` uses the variable name `if_`.
- abs_energy_ep_old never sets its result (it assigns `abs_energy_ep`, a local there):
  returns 0.0. Unused by the engine.
- abs_energy_ep writes `type(...)` with `a3`: passed as `Ch(type, 3)`.
- Hot loops (atogen, press, fmolec1) alias module arrays to locals. Neither routine
  reallocates them, except ergold/ergnew/nbref, which are re-fetched after each
  make_neighbour_list_from_mols/make_neighbour_list_from_cells call. `minimagec.minimagec` is
  bound to a local at routine entry, so tracing still works when it is enabled before
  the call.
- coordconstraints_mod.coordconstraints is used 1-based (`[None] + ...`, package H).
