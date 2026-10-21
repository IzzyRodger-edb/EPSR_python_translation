# Package K: epsr_standalone_rev3.f lines 4737-6584

Output: `src/epsr26/programs/epsr_standalone_rev3_part3.py` (to be merged into
`programs/epsr_standalone_rev3.py`).

Routines (19, all complete): gensq, fit_expansion, expand, calc_gaussq, calc_gaussr,
calc_gaussr_grad, calc_pwrexpq, calc_pwrexpr, calc_pwrexpr_grad, dnext, ftback, uniform,
invert_cwt, generate_rp, generate_ep, lsparabola, lsparabola1, multiply_pot_coeffs,
multiply_pcoefex.

## Verification (bit-exact against gfortran -O2 -fno-tree-vectorize)

Scratch drivers built the Fortran routines in isolation and compared every float as hex
with the Python output, plus stdout and the written files:

* calc_pwrexpr_grad, calc_pwrexpr (with rmin > 0, ncfirst), calc_gaussr_grad, calc_gaussr,
  calc_gaussq, calc_pwrexpq, dnext, expand: identical.
* fit_expansion (with the real smoop and a deterministic stand-in ran1, same in both):
  wt, wtold, yout, stepr, acceptr, rejrat, chisqnew and the `fit_expansion>` line identical.
* invert_cwt (nfile=2, npairs=3, same stand-in ran1): cwtt, cwt and all stdout identical.
* uniform: identical where the Fortran writes; beyond the last point it fills, the Fortran
  array holds stack garbage (Python: zeros).
* gensq, ftback (real epsr16/molecom/molecom_mod data): identical.
* lsparabola, lsparabola1 on a synthetic .erg: returns, stdout and the .qdr files identical.
* generate_rp (modmorse + Gaussians, exponential and harmonic cores, gradlimit > 0),
  generate_ep (power-exponential and Gaussian bases), multiply_pot_coeffs, multiply_pcoefex
  (symbol rewrite through parse + list-directed internal write): all arrays, pcoefp,
  rminex, symbol_value and stdout identical. The ato_manipulation_rev4 routines
  (generate_mp, chsph, gauss6pot, abs_energy_ep) were replaced by the same simple stubs on
  both sides, so their own translations still need testing.

## For the integrator

* **Imports.** At module level: common (clustercom, commands, epsr16, epsrinput_par,
  molecom, potential_par, symbols, molecom_mod, and `containerforces_mod` /
  `lngamma_mod` imported as `containerforces_mod_data` / `lngamma_mod_data`), library
  (containerforces_mod, lngamma_mod, parse), fortran, trace, and `math`.
  `auxiliary_routines` (ran1, smoop) and `ato_manipulation_rev4` (generate_mp, chsph,
  gauss6pot, abs_energy_ep) are imported **inside** fit_expansion / invert_cwt /
  generate_rp / generate_ep, because those modules did not exist when I wrote this. They
  can be hoisted to module level once they exist.
* **No F821 names from the sibling parts.** My routines only call each other. Outside
  this part, callers are runepsr17 (part I: gensq, ftback, uniform, lsparabola1,
  multiply_pcoefex) and part J (initpt: generate_rp, generate_ep; ptgnsq: calc_gaussq,
  calc_pwrexpq, fit_expansion, generate_ep; repulsivept: generate_rp; get_wts: invert_cwt;
  initg: uniform).
* **Module-level names** that could clash when merging: `_EXP25`, `_F200_INVERT_CWT`,
  `_F100_LSPARABOLA`, `_mpairs` (alias of epsr_dimension.mpairs), `_lsp_read_record`
  (private helper for lsparabola/lsparabola1's 10-item `read(nchan,*,iostat=ierr)`, which
  assigns items in order, so a failed read keeps the items already read, as gfortran does).
* **ftback's first dummy `in` is `in_`** (`in` is a Python keyword). Pass it by position.
* **Signature notes.** No exceptions: SIGNATURES.md matches every routine. Callers must
  use the returns: `ncfirst = calc_*(...)`, `d, rootd = dnext(...)`,
  `stepr, acceptr, rejrat, chisqnew = fit_expansion(...)`, and the 9-tuples of
  lsparabola / lsparabola1.
* **Types callers must pass.** Real dummies must be numpy.float32 (e.g. sigmab, rmin,
  fbroad of the calc_* routines; rmin of ftback; epfact). Plain Python literals would give
  real*8 or weak-typed arithmetic. generate_ep passes `r4(0.0)` and `r4(1.0)`.

## Fortran quirks reproduced / uninitialised variables

* generate_rp: `add` is a COMMON variable (epsr16.inc), so every `add=` is written to
  `epsr16.add`.
* generate_rp sets `xav(npt,ic)=0.0`, not `xav1(npt,ic)`. Reproduced as written.
* generate_ep: `ncfirst` is an uninitialised local passed to calc_*_grad with rmin = 0, so
  it is never assigned or used. 0 is passed.
* lsparabola1: `navref` is uninitialised if no point after the first is read, and is then
  written to stdout. 0 is used (as in the Rust).
* lsparabola: `x`, `y`, `sum*` are real*8 but `xm` and `xmax` are implicit real*4, and
  `(1.0-epcfrac)*rfacwrite` is a real*4 product added to real*8. All the casts follow the
  Fortran exactly. In lsparabola1, `x` and `y` are real*4.
* Both lsparabola routines inhibit the quadratic fit (`delta=0`). The dead branch is
  translated anyway. The `delta` written to stdout is the straight-line determinant.
* invert_cwt: `add=(2.0*ran1()-1)*amp` is a real*4 expression stored in real*8, and
  `amp`/`fac`/`feedback1` are real*4. The `write(6,200)` lines and blank lines are kept
  (the Rust drops them).
* multiply_pcoefex: the internal list-directed write into the character*256
  `symbol_value(icount)` uses `fstr(lw(...), 256)`. If that record overflowed 256
  characters, gfortran would return an error in `iostat` (ignored by the Fortran); here
  the text is truncated instead. This should never happen in practice.
* gensq: `sumwavself`, `ycohfactor` and the coarse-graining `prodssff` product are
  computed but unused, as in the Fortran.
* exp(25.0): gfortran folds this constant (correctly rounded). It equals expf(25.0)
  (checked), held in `_EXP25`.

## Performance (speed-only changes; arithmetic and order unchanged)

* gensq and expand iterate over Python lists / slices of the arrays instead of
  indexing numpy elements one by one.
* invert_cwt works on `tolist()` copies of cwt/cwtt/carray/carrayold (Python floats
  are real*8) and copies them back at the end. It can run up to 1000 x 50000 MC steps, so
  this matters. ran1 is bound to a local at entry.
* Literals in the hot loops are plain Python floats meeting float32 operands
  (TRANSLATION_GUIDE 2.3), not `r4(...)` calls.
* Possible further speed-up (not done, against guide rule 5): the expand rows in
  fit_expansion's chi-squared and the gensq Q points are independent. Vectorising
  *across* rows with elementwise float32 numpy (keeping the per-row term order) would give
  identical bits; the Rust does this with rayon. fit_expansion is the main remaining cost
  of ptgnsq.

## Compat-layer helpers I would have liked

* A list-directed READ helper that assigns a fixed list of typed items and returns
  `iostat`, keeping partially read items (my `_lsp_read_record`).
