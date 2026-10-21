# Package I: epsr_standalone_rev3.f lines 1-2522 (`program epsr_standalone`, `runepsr17`)

Files written:
- `src/epsr26/programs/epsr_standalone_rev3.py`: `main()` (the program) and `runepsr17()`, one
  long function in Fortran statement order, plus two private helpers (`_ios`, `_iteration_finished`).
- `src/epsr26/common/sumvalues.py`: hand-written COMMON `/sumvalues/` (fc1, fc2 real*8; add1, add2,
  sumtot1, sumtot2, errort real*8(2)), with `reset()` and `common.register`. **Also declared in
  `sumgr` and `sumsq`** (same file, packages J/K): they must use `from ..common import sumvalues` and
  read `sumvalues.fc1` / `sumvalues.fc2` (set by runepsr17 just before calling them).

## For the integrator

- **Bare-name calls to routines from part2/part3** (F821, each line carries `# noqa: F821`, which can
  be removed after the merge): `inita, get_wts, initg, initpt, interp, lsparabola1, openf, gensq,
  sumsq, sumgr, entrpb, entrpc, ptgnsq, repulsivept, sumpe, norma, ftback, uniform, multiply_pcoefex`.
  Return conventions used (from SIGNATURES.md): `iinit, nz = inita(...)`; `nsum = initg(...)`
  (assigned back to `epsrinput_par.nsum`, or `.nesum` in the `nsums < 0` branch);
  `nqmind[in_], nqmaxd[in_] = interp(...)`; `oiseb = entrpc(oiseb, wav)`;
  `sentb = entrpb(sentb, nsum)`; `(nparab, aparab, bparab, cparab, rfacwrite, aav, bav, lastqf,
  nextqf) = lsparabola1(9, fname, ...)`; the others return nothing.
- **Library and Fortran-module data modules are imported inside `runepsr17`/`main`** (they did not
  exist yet; a top-level import would make the module fail to import). Names assumed:
  library `file_and_directory_commands` (get_command_from_args, setup_input_file, search_for_symbol,
  ireff), `setup_epsr_rev1` (setup_symbols_epsr), `read_and_save_files` (read_input_file,
  save_files_others), `ato_read_and_save` (normalise/unnormalise_coordinates),
  `ato_manipulation_rev4` (get_maxmolsize, fmolec1, adjust_intra_wts, atogen, bragginten,
  change_boxsize, press, sharm_diffuse, stog, and `initialise_timer` of module `timer`),
  `minimagec`, `molfit`, `auxiliary_routines` (convertphibthetacphic, convertalphabetagamma, ran1,
  w_int, w_intr), `cellneighbourlist` (initialise_cellneighbourlist, get_date_and_time -- the latter
  is outside the Fortran module but in the same file), `parse`.
  Data: `common.timer` (sumtime1, sumtime2), `common.diffuse_mod`, `common.threebody_mod`
  (nthreeb, sume_triplet_dist), `common.molecom_mod` (qscaleinputfactor, selfbroadpwr,
  coarsegraining, rhocg). They can be moved to the top of the module once they exist.
- `w_int`/`w_intr` are given the titles as plain strings (`titleq`/`titler` are character*14,
  `titlepsf`/`titledat` char_array of character*29 values); if those routines need the declared
  lengths for `(a)` output they must pad themselves (or we add keywords).
- The `do 3000` loop: `trace.context["iter"] = iter` at the top of every iteration; the end of an
  iteration is marked `# ---- 3000 continue` and followed by `_iteration_finished(iter, niter)`
  (no-op). NB for the graceful Ctrl+C: all result files (`.f .s .d .q .t .u .v .g .l .y .r .x .w .o .p
  .z`, norma) are written only when `iter == niter`, tested **before** the end of the iteration, so
  to "make the current iteration the last" the shell will have to lower `niter` before that test
  (e.g. a second hook, or a flag checked where `if iter == niter:` is) -- the hook at the loop end
  can only stop the loop after a non-final iteration (then no result files for that run, as with a
  Fortran kill). Nothing was changed in the physics/order for this.

## Deliberate fix
- `# FIX (FORTRAN_BUGS.md 1.5)`: pressure-control keyword: the Fortran accepts `ca` but its branch
  tests `ac`; now both `ca` and `ac` are accepted and select the c+a axes (as Rust pressure.rs).

## Uninitialised variables (all set to 0, which is what the Rust reference does)
- `aparab, bparab, cparab, nparab, aav, bav, lastqf, nextqf` when `iinit == 1` (lsparabola1 not
  called): 0 (Rust `parab = None`).
- `epcup, epccoeff` on non-final iterations (written to `.erg`/`.out`): 0, or the value from the
  previous iteration of the same run (they are function locals that persist across iterations,
  same as Rust `self.summary`).
- `ctrllevel, ctrlbaseline` (only set in the `bparab != 0` branch of the automatic ereq): 0, kept
  across iterations (Rust `self.ctrllevel`).
- `rfacd` (only set on the last iteration, written every iteration in `.out`): zero array.
- `rfacraw` on non-final iterations = `rfacwrite` exactly as the Fortran says (lsparabola1's value,
  or 0 with iinit = 1) -- FORTRAN_BUGS.md section 2 / Rust `rfacraw_last`.
- `eminmaxdir` (a typo for eminmaxdif in a `write(6,*)`, never set): 0.
- `overallfac` (passed to ptgnsq, never set, unused by ptgnsq): 0.
- `rmaxgrsave` if its `read(...,iostat=ierr)` fails, `navmaxparab/nparabstartfactor` on a null read:
  0 (Rust `unwrap_or(0)`). `hashsearch` starts as ''.

## Quirks reproduced
- `text='rmaxgr'; call search_for_symbol(text,icount)` does not reset `icount` (keeps the `inter`
  value); passed through as in the Fortran.
- Lattice refinement: `clatnew = clatnew*blatfac` (not clatfac) kept, marked NB.
- `add` used as a scratch variable in the Bragg convolution is the COMMON `epsr16.add`; `twopi` is
  `lattice.twopi`; both are written to the COMMON as in the Fortran.
- `datagr2(mpz,mpairs)` is passed to `stog` whose dummy is `gofr(mpz,mfile)`: same leading dimension,
  so passing the array unchanged is equivalent.
- `open(8,...,status='OLD')` and the list-directed reads of the symbol table have no `err=`/`iostat=`
  in the Fortran: a missing data file / bad value is a Fortran runtime error; here the
  `FortranIOError`/`FortranEOF` propagates out of `runepsr17` (run_program does not catch it).
- Unit 35 is never opened: writes go to `fort.35` in the current directory (units default), unit 9
  (`.erg`) is opened by `openf`, read to its end by lsparabola1 (so the write appends) unless
  `iinit == 1` (then it is overwritten), closed after the loop.
- In the data-file read loop the list-directed values are assigned one by one (an error part way
  leaves the earlier items assigned, as gfortran); the `(i7,3g14.7)` GENIE read assigns only on success.

## Arithmetic notes
- weightx/backx (X-rayf0) and cwt/ciwt/cwttq (weights) are real*8: every `sum=sum+w*x` is
  `f32(sum + w*x)` (double add, then rounded), `cwttemp = f32(float(ratinter1)*cwttq + ...)`.
- `fc1 = float(r4(1.0)/f32(nsum))`, `fc2 = float(f32(nsum-1))*fc1` (single then real*8, as Rust).
- `twopi**1.5 -> powf(twopi, 1.5)`, `vol**(1.0/3.0) -> powf(vol, r4(1)/r4(3))`,
  `(1+rat)**(1.0/real(npcontrols)) -> powf`, `(porodq/q)**porodpwr -> powf`, `nisumm = ifix(0.3*intra)`.

## Helper I would like in the compat layer
- An `iostat` helper: the pattern `try: ...; ierr = 0 / except (FortranEOF, FortranIOError) as e:
  ierr = -1 or e.iostat` appears ~15 times; I used a private `_ios(e)`.

## Checks done
- `python -c "import epsr26.programs.epsr_standalone_rev3"` OK; `ruff check --select F` clean.
- Smoke test with all dependencies mocked (3 iterations, lattice refinement, Bragg block, automatic
  ereq, pressure control `abc`): runs to the end and writes `.out`, `.erg`, `fort.35`, `.inten` in
  the expected formats. Not yet compared with the Fortran (needs the other packages).
