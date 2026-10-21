# Package M: setup and weights programs (wave 2)

Files (all in `src/epsr26/programs/`, each with `def main()` translating the `program`):

| Fortran | Python | routines |
|---|---|---|
| makemole.f | makemole.py | program only (builder is `library/ato_manipulation_rev4.makemole`) |
| mixato.f | mixato.py | mixato |
| introtcluster.f | introtcluster.py | randomise (program `randomise`; ends with `stop`) |
| fmole.f | fmole.py | fmole, initpt_fmole, update_ato |
| nwts.f | nwts.py | nwts |
| xwts.f | xwts.py | xwts, form_fact |
| epsrwts.f | epsrwts.py | epsrwts |
| epsrwtsx.f | epsrwtsx.py | epsrwtsx, form_fact |
| changeato.f | changeato.py | changeato |
| readpdb.f | readpdb.py | run_readpdb, readpdb, upcase |
| makeato.f | makeato.py | makeato |

Every routine is translated completely; no stubs.

## Decisions and quirks

* **fmole runs makemole** for each molecule template (`sometext='makemole '//template`,
  then `call makemole()`), as the Fortran does; the Rust port refuses this path. Verified
  with the validation/water steps *without* `delete Water.mol` (makemole regenerates
  Water.ato inside fmole): IDENTICAL.
* **readpdb** follows the patched reference (explicit `ctype` interface): `ctype` gets the
  11 positional arguments, i.e. `natomsph=0.5*sigma, rhosph=1.0, radminsph=0.0`, with
  radmaxsph/pseudofacsph absent. `dihedral(im, ...)` with absolute atom numbers and the
  never-set `nal` of one-atom molecules are reproduced. DO-variable values after loops are
  set explicitly (`_do_end`) where the main read loop's `ia<matom .and. im<mmol` test sees them.
* **makeato**: the reference dies with SIGSEGV in `moloutc` (natomsphere/pseudofacmol never
  allocated). The Python catches the resulting TypeError only when those arrays are
  unallocated, discards the buffered unit output (`units.reset()`, so the .ato file stays
  empty as in the reference), prints gfortran's SIGSEGV line on stderr and ends with
  `FortranStop(139)` (exit 139, same as the Fortran). The Fortran's stdout is lost in the
  crash too (buffered pipe); the Python's prompts are printed -- stdout only.
* Uninitialised locals (`ierrout` before `setup_input_file`, `iaccept` before
  `search_for_files`, read targets such as `nrho`, `massno`, `fracmix`, template counters
  in `update_ato`, `cx(it)` of Q atoms in epsrwtsx/xwts, `binter(it,3)` of Q atoms in nwts)
  start at 0, marked `# NB:`.
* `epsrwts`: the abundance renormalisation loop divides `abund(it,nsum,k)` (not
  `abund(it,isum,k)`) nsum times -- reproduced.
* `read(5,*,err=label)` (epsrwts) jumps back only on a conversion error; end of file is not
  caught (fatal, as in the Fortran).
* `a3` output of character*3 variables uses `Ch(s,3)` (first run of mixato showed the
  right-justified difference in stdout).
* Programs that do not include `epsr16.inc` (randomise, mixato, epsrwts, epsrwtsx) pass
  implicitly typed *locals* rmaxcell, cellmax, upol, rho, amubar, amubar2, ndist to
  `minitc`; changeato and fmole pass the epsr16 COMMON variables.

## Verification (tools/run_ex_py.py, resume mode, one case at a time)

All output files of every step of these programs are IDENTICAL. stdout logs were also
compared (after masking folder names, timestamps and `/bin/sh:` vs `sh:`): identical, except
gfortran's runtime note "IEEE_DENORMAL" after randomise (ex35, water) and the lost Fortran
stdout of makeato (ex30).

| case | step | program | Fortran s | Python s |
|---|---|---|---:|---:|
| ex20 | 1-3 | nwts | 0.27 / 0.13 / 0.12 | 1.22 / 1.38 / 1.49 |
| ex21 | 1-4 | nwts | 0.23 / 0.12 / 0.17 / 0.01 | 1.93 / 1.38 / 1.64 / 0.48 |
| ex22 | 1-3 | nwts | 0.17 / 0.12 / 0.17 | 1.43 / 1.29 / 1.64 |
| ex23 | 1-3 | xwts | 0.17 / 0.07 / 0.12 | 1.02 / 0.87 / 0.92 |
| ex24 | 1-4 | epsrwts | 0.17 / 0.12 / 0.11 / 0.12 | 1.43 / 1.39 / 1.18 / 1.43 |
| ex25 | 1-3 | epsrwtsx | 0.18 / 0.07 / 0.07 | 2.00 / 1.24 / 1.53 |
| ex26 | 1-3 | changeato | 0.12 / 0.02 / 0.17 | 1.33 / 0.58 / 6.48 |
| ex27 | 1-2 | changeato | 0.12 / 0.01 | 6.06 / 0.42 |
| ex28 | 1-2 | makemole | 0.02 / 0.04 | 0.42 / 0.37 |
| ex29 | 1-3 | makemole | 0.01 / 0.01 / 0.05 | 0.48 / 0.42 / 0.37 |
| ex30 | 1-2 | makeato (exit 139 both) | 0.02 / 0.01 | 0.98 / 0.42 |
| ex31 | 1 | mixato | 0.03 | 0.48 |
| ex32 | 1 | mixato | 0.07 | 0.57 |
| ex33 | 1 / 2-3 | mixato / fmole | 0.04 / 0.07 / 0.01 | 0.43 / 1.48 / 1.03 |
| ex34 | 1-2 | fmole | 0.07 / 0.02 | 0.98 / 0.43 |
| ex35 | 1 / 2,4,5 / 3 | mixato / randomise / fmole | 0.07 / 0.04, 0.04, 0.12 / 0.03 | 0.57 / 0.63, 0.73, 1.54 / 1.33 |
| ex39 | 1 / 3 / 4 | mixato / randomise / fmole | 0.08 / 0.05 / 0.12 | 0.78 / 0.67 / 2.03 |
| ex44 | 1 | readpdb | 0.04 | 0.67 |
| ex45 | 1 | readpdb | 0.02 | 0.58 |
| ex49 | 1 | readpdb | 0.01 | 0.47 |
| validation/water | 1-7 | makemole, mixato, randomise, fmole 100, nwts x3 | 0.02, 0.07, 0.12, 0.21, 0.17, 0.17, 0.12 | 0.48, 1.36, 1.59, 19.31, 1.46, 1.38, 1.38 |
| water with Water.mol kept | 4 | fmole 100 (+ internal makemole) | 0.22 | 17.64 |

ex39 step 2 (makelatticeato, another package) DIFFERS / python exit 1 at the time of the run;
the steps after it start from the Fortran folder and are IDENTICAL.

## Notes for other packages

* None of the library modules needed changes.
