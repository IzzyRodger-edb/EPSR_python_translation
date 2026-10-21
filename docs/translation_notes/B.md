# Package B: auxiliary_routines.f, read_and_save_files.f

Files: `src/epsr26/library/auxiliary_routines.py` (30 routines),
`src/epsr26/library/read_and_save_files.py` (11 routines), test
`tests/library/test_ran1.py` with `tests/library/data/ran1_driver.f90` and
`tests/library/data/ran1_reference.txt`.

## Verification (all against gfortran, same flags as tools/build_fortran_reference.sh)

- `ran1`/`initran1`/`saveran1`: bit-exact over 6 seeds x 1000 deviates plus a
  save/restore/continue block, including the saveran1 state line (pytest).
- All other auxiliary routines (w_int, w_intr, r_int, r_int_titles, indexlastdot,
  get_atomicnumber, cosd/sind/acosd/asind/atan2d, convert*, perpd, normvec, rotateeuler,
  smoof/smood/smoop, dlm/maxdlm, cgcoef/maxcgcoef): a scratch driver comparing bit
  patterns and written files (including 105 columns split over .f01/.f02, w_intr, and
  stdout) -- byte-identical.
- outgofz, outgofz_sharm, outgofzs_sharm: 4 accumulation passes (write, read back with
  r_int, average) -- byte-identical files and stdout.
- read_input_file -> process_commands -> read_ato/read_pcof, then save_files_others:
  a driver doing get_command_from_args, setup_symbols_epsr, setup_input_file,
  read_input_file, save_files_others(14) (no simulation), run on all 49 EPSR-test ex
  cases with packages A/C/D/H code: .EPSR.inp, .pcof, .ato and stdout byte-identical in
  48 cases; ex108 differs only by FIX 1.3 (see below), the divergence documented for
  the Rust.  Covers tb lines (ex110-114) and coordcon lines (ex115-119).
- Not exercised: save_files / force_save_files (interactive prompts), the error branch
  of convertalphabetagamma.

## Bug fixes (FORTRAN_BUGS.md section 1), marked `# FIX (FORTRAN_BUGS.md 1.x)`

- 1.1 `save_files_others`: the re-parsed tb values are kept in locals (as Rust
  `write_pcof_threebody`), so memory keeps read_pcof's cosines / width 1.0 / abs()
  amplitudes; the written line is unchanged.  The `nthreebangsmax` update (wrong index,
  FORTRAN_BUGS 3) and the reallocate calls are kept on the module data as written
  (no-ops in practice).  Checked: after save on ex110 memory still holds cos(180)=-1.
- 1.4 `save_files_others`: `igauss = igauss - 1` on a Gaussian read error.
- 1.3 `read_pcof` (new pair created for a pair missing from the .pcof): coulfac=1,
  efacp=1, coordex=0.  FORTRAN_BUGS.md names "save_pcof" but there is no such Fortran
  routine; the Rust fix is in pcof.rs read_pcof, i.e. this branch of read_pcof.  Not in
  my brief explicitly, but the guide says all five are fixed the Rust way; it is the only
  difference seen in the 49-case comparison (ex108, coulfac column of the new pairs).

## Uninitialised variables / quirks (reproduced)

- `indexlastdot` returns 0 when there is no dot (the comment says len+1).
- `read_pcof`: `ireset, iinit, roverlap, rminfac` are implicit locals; roverlap/rminfac
  start at 0.0 when their symbol cannot be read.  The lengths returned by concatenate
  (`ilenvalue`, `ilensymbol`, `ilensymbol_value`) start at 0 (outputs only).
- `read_pcof`: the q/Q exemption from roverlap tests the numbers line, not the labels
  (FORTRAN_BUGS 2, Rust pcof.rs:258) -- reproduced.
- `read_pcof`/`save_files_others`: `ierr` is the iostat of the *last* read, as in the
  Fortran (Rust makes some errors sticky; only matters for malformed files).  Words past
  nwords are read from the stale ncf/ncl entries of COMMON /commands_i/, as in Fortran.
- `read_input_file`, `read_pcof`: a formatted read at end of file leaves `sometext`
  unchanged, so the last line is processed once more (as Fortran and Rust do).
- `outgofz*`: `lptin` uninitialised when r_int fails before reading; passed as 0 (unused).
- `w_intr` builds its title format with nblock repeats instead of iwr; same output.
- `w_int`: `recl=` of the open is not modelled (records never exceed it).
- `convertalphabetagamma`: on the error returns, the intent(out) values not yet set keep
  the caller's values (the Python takes them as arguments and returns them).
- `ran1`: `-idum` is wrap32'ed (idum = -2**31 case) and mcall wraps at 32 bits.

## Conventions an integrator needs

- `ran1()` takes no arguments; state is `common/ran1com`.
- `cosd/sind/acosd/asind/atan2d` here are EPSR's own (degradconv = 3.141592653/180
  folded in single precision); not gfortran's intrinsics (fmath.*_intrinsic).
- `dlm`, `cgcoef` read `common/lngamma_mod.gamln`/`mgamma` (call fctrln first, as Fortran).
- `outgofz` and `outgofz_sharm` are identical in the Fortran and share a private body.
- Reallocations: `tb.x = reallocation_routines.reallocate..(tb.x, ...)` and
  `coordconstraints_mod.coordconstraints = coordconstraint_mod.reallocate1d_cc(...)`
  (package A/H convention, see SIGNATURE_EXCEPTIONS.md).
- `readline` is called with `len_line=commands.nstsize` (sometext is character*3020);
  `init_constraint` with `len_sometext=256`.
- Other-package modules are imported inside the routines (cycles: process_symbols
  imports this module).

## Private helpers (compat-layer wishes)

- `_star_repeat`: fformat has no `*(...)` unlimited-repeat group; formats built with
  `i1.1` for 2*nangs >= 10 contain `*(`; replaced by `999(` (same output).
- fformat treats a repeat count of 0 (`0(...)`) as 1; harmless where it occurs here
  (no data left), but the compat layer should honour 0 and `*`.
- `_read_real`/`_read_int`: `read(text,*,iostat=ierr) x` -> (x, ierr).
- `_dexp`/`_dlog`: real*8 exp/log with Fortran's inf/-inf/NaN instead of Python errors.
