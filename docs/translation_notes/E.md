# Package E: ato_manipulation_rev4.f lines 1-2343

## Files

- `src/epsr26/library/ato_manipulation_rev4.py`: 21 routines, all translated in full.
  - Module procedures: `initialise_timer` (module timer) and `initialise_diffuse_values` (module diffuse_mod).
  - The rest: `ftype`, `ctype`, `scatt_len`, `fractcoord`, `wholecoord`,
    `generate_fractional_coordinates`, `generate_cartesian_coordinates`, `centreofmassi`, `centreofmass`,
    `en_tethered`, `gauss6pot`, `chsph`, `sum1`, `sum2`, `resetq`, `truncate`, `truncatec`, `truncatecweeks`,
    `makemole`.
- Module data (TRANSLATION_GUIDE 5.2, each with `reset()` and registered):
  - `src/epsr26/common/timer.py`
  - `src/epsr26/common/bragginten_mod.py`
  - `src/epsr26/common/reciprocallattice_mod.py`
  - `src/epsr26/common/diffuse_mod.py`

## Merging with part 2 (package F)

- `get_rotgrp` / `get_whole_rotgrp` appear inside makemole (lines 606-616) only as `INTERFACE` blocks. They are
  defined in part 2 (lines 2579, 2682). makemole calls them by bare name:
  `labrott1, nrott1 = get_rotgrp(ir1, ir2, labrott1, nrott1)`, with `# noqa: F821 (part 2)`. This is the only
  F821 in the file. Part 2 returns `(labrott, nrott)` (logged by F in SIGNATURE_EXCEPTIONS.md).
- Part 2 calls `en_tethered`, `initialise_diffuse_values` and `generate_fractional/cartesian_coordinates` from
  this file by bare name. The signatures agree.
- Private names are `_do`, `_iread`, `_read_a`, `_read5_real` and the FORMAT constants `_F<label>`
  (`_F100_ftype`, `_F100_ctype`, `_F100_makemole` where labels repeat). Part 2's private names are
  `_omp_*`, `_atan2d_intrinsic`, `_intenhkl_ref` and `_F<label>_<routine>`, so there is no clash as of this
  writing. The two import blocks need merging.

## Conventions used

- Imports:
  - The library module procedures of `molecom_mod.f90` and `lngamma_mod.f90` are imported as `molecom_mod_lib` and
    `lngamma_mod_lib`.
  - The data modules keep their own names (`common.molecom_mod`).
- `reallocate1d_r` and the other reallocation routines are called as `arr = reallocation_routines.reallocate1d_r(arr, n)`,
  following package A's convention.
- `ctype` takes its five Fortran `optional` arguments as keywords defaulting to `None`.
- `gauss6pot` returns `(gauss6pot, sigg)`, because it clips `sigg` to 0.29.
- `scatt_len` returns `(scatt_len, ifind)`.
- `fractcoord`, `wholecoord` and `generate_*` return the coordinates they assign.

## Arithmetic decisions

- `coul**6` in gauss6pot is `c3 = c*c*c; c3*c3`. gfortran expands the constant power this way, and the Rust
  potential.rs does the same.
- `sum1`/`sum2`:
  - Integer literals that meet reals are written `f32(3)` and `f32(4)`.
  - The constant `1.0/6.0` is folded in single precision (`r4(1.0)/r4(6.0)`).
- makemole:
  - The bond-angle distance is `b13 = b21*b21 + b23*b23 - f32(2)*b21*b23*cosd(bangle)`.
  - `external cosd,sind` means EPSR's own `auxiliary_routines.cosd`/`sind`, also used for the unit-cell volume.
  - `vol = f32(natomr)/rho` and `dmax = powf(vol, r4(1.0)/r4(3.0))`.
- `truncate` uses `PI4ATAN`. `truncatecweeks` uses `erfcf`.

## Quirks reproduced (all marked `# NB` in the code)

- **makemole, wholerot groups.** The Fortran assigns `rotgrp(ir)='ROT'` *before*
  `testsize_rotgrps(ir)` enlarges the array.
  - The store goes past the end of the old array and is lost.
  - The new element is uninitialised heap memory, which is zero bytes in the reference.
  - The Python stores 'ROT' only if `ir` is already in range. After the resize it sets the new element to
    `'\x00\x00\x00'` (Rust makemole.rs:517).
  - Checked: the labels match the Fortran `.ato` (`DIH`, `ROT`, `^@^@^@`, `ROT`) on a custom input.
  - `moloutc` (package D) must write `rotgrp` as is, so that the NUL bytes reach the file.
- **makemole, uninitialised locals.** These are set to 0:
  - `nspace`: only read for `.mol`-type files. A column-map line in a `.gmol` file then divides by zero; `idiv`
    raises ZeroDivisionError, where the Fortran gets SIGFPE.
  - `nrow`, `nrowmax`, `ncolumnmax`, `nstar`.
  - `nmap`: incremented but never used.
  - `nrowoffset`, `ncoloffset`, `irmap`, `icmap`, `nread`: these keep their previous value when an
    `iostat=` internal read fails.
- **makemole, read loops at end of standard input.** `do while(ierr.ne.0) read(5,*,iostat=ierr)` loops for ever, as
  the Fortran does. This applies to the temperature, density and potential-parameter prompts. The bond and angle
  prompt loops run once whatever happens.
- **makemole, the `atom` keyword.** The test is `naread.lt.mstarlabel` (angles, not atoms), as written.
- **makemole, `pseudofac` with no value.** This reads `ncf(2):ncl(2)`, which still holds the previous line's word
  positions, as in the Fortran.
- **initialise_diffuse_values.**
  - The counting loop is `do m=0,mdiffmax` with no step, but the filling loop uses step `mdiffstep`. So
    `l1ex`/`n1ex` are allocated larger than the final `ndiffuse`, and `ndiffuse` ends at the smaller count.
    This is reproduced.
  - `mdiffstep=0` (or `ldiffstep=0`) is a DO loop with a zero step. `_do()` raises ValueError; the Rust reports
    "unsupported".
- **ftype.** A bare `stop` raises `FortranStop(0)`.

## Format notes

- The `$` edit descriptor (formats 101 and 102) is dropped from the format string and written with
  `advance=False`. The compat FORMAT parser does not accept `$`; a `$` descriptor in fformat would be a
  welcome addition.
- Character variables written with `a3` are passed as `Ch(x, 3)`.
- Internal writes to `sometext` are `fstr(fmt(...), nstsize)`. Writes into substrings use `setsub`.

## Verification (scratch harness)

The harness stubs `moloutc` and the program header, and uses a scratch copy of `get_rotgrp`. The `.ato` files
themselves need package D's `moloutc` and were not compared.

- **ex28** (DMSO, WaterMolecule): the rewritten `.mol` files and all of makemole's stdout are identical to the
  Fortran.
- **ex29** (Chain with stdin, Tip4p, Grid map file with stdin): `Chain.mol`, `Tip4p.mol`, `Grid.mol`, `Grid.atm`
  and stdout are identical to the reference and to `fortran_build/makemole`.
- **Custom inputs**:
  - `.gmol` with `dihedral` x2 on the same axis, `wholerot`, `introt`, `extra`, `qradius`, `unitcell` and a bad
    `vibtemp`.
  - `.mol` grid with a missing potential (prompt), `temperature -5` and a density prompt.
  - Both are identical to the Fortran in the `.mol`/`.atm` rewrite and stdout.

## Helpers written privately (in the module)

- `_do(start, end, step)`: the iteration range of a Fortran DO loop for any sign of step.
- `_iread(text, kind, old)`: an internal list-directed read of one item with `iostat`.
- `_read_a(unit, old, n)`: `read(u,'(a)',iostat=)`.
- `_read5_real(old)`: `read(5,*,iostat=)` of one real.

These, and the `$` descriptor, would be welcome in the compat layer.
