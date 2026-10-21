# Package Q: sharm, sdf, sdfcube (sharm_routines_rev1.f, sharm_rev2.f, sdf_rev2.f, sdfcube.f)

## Files

| Fortran | Python |
|---|---|
| sharm_routines_rev1.f, `module sharm_mod` (lines 1-30) | `common/sharm_mod.py` (module data, guide 5.2) |
| sharm_routines_rev1.f, routines | `library/sharm_routines_rev1.py` |
| sharm_rev2.f (`program standalone`, `run_sharm`, `sharm`) | `programs/sharm_rev2.py` |
| sdf_rev2.f (`program standalone`, `run_sdf`, `sdf`) | `programs/sdf_rev2.py` |
| sdfcube.f, `module sdfcube_mod` (lines 1-13) | `common/sdfcube_mod.py` |
| sdfcube.f (`program standalone`, `run_sdfcube`, `initcube`, `outcube`, `sdfcube`) | `programs/sdfcube.py` |

Routines translated in full (no stubs): get_l_and_n_values, form_cgcoeffs, getrefint,
getrefintnon, form_eulerangles, get_max_size_dlmn, form_dlmn, form_sum_prod_dlmn_cgcoeff,
populate_disk (never called: its calls are commented out in the Fortran); main/run_sharm/sharm;
main/run_sdf/sdf; main/run_sdfcube/initcube/outcube/sdfcube.

## Names

* The include files `sharm.inc` and `sdfcube.inc` are the generated modules `common.sharm`
  and `common.sdfcube`, whose names clash with the routines `sharm` (sharm_rev2.f) and
  `sdfcube` (sdfcube.f).  The programs import them as `sharm_inc` / `sdfcube_inc` so that
  the routines keep their Fortran names.
* `library/sharm_routines_rev1.py` imports `library/lngamma_mod.py` as `lngamma_mod_lib`
  (guide 11, same-name modules); dlm/maxdlm/cgcoef/maxcgcoef/scapro/vecpro/ran1 are
  auxiliary_routines.f's (`auxiliary_routines.*`).
* The three averaged-molecule xyz blocks of sharm and sdfcube (identical in the Fortran) are
  written by one helper `_write_xyz` per program, and the "is this atom part of the central
  molecule definition" test by `_central_definition`; the two "FINALLY we sum the products
  of coefficients" blocks of form_sum_prod_dlmn_cgcoeff (identical except for which molecule
  is at the origin and which dlmnc row is used) by `_sum_products`.  Each keeps the
  Fortran's statements and order.

## Patched FORMAT (reference build)

sharm_rev2.f and sdf_rev2.f build the column-title format
`'(1h('//il1fmt//',1h,,'//...` with no comma after the Hollerith `1h(`; gfortran rejects it
at run time.  The reference build (`tools/build_fortran_reference.sh`) inserts the comma,
`'(1h(,'//il1fmt//...`, and the Python translates that patched behaviour (titles
`(l1,l2,l;n1,n2)`, e.g. `(2,2,0;2,-2)`), marked in the code.

## Quirks reproduced

* sdf_rev2.f `run_sdf` does not include sharm_sdf.inc: its `pi`, `piconv` (declared) and the
  test angles `phim1 ... phillab` (implicit) are locals.  `sdf`/`form_dlmn` use the
  sharm_sdf.inc variables, which stay 0, so with `test` the "known" orientations are all
  zero (Rust sdf.rs says the same).  In sharm_rev2.f they are the COMMON ones.
* sdf: the check of the remaining intra-group distances uses `zdiff=ytemp(i2)-ztemp(i1)`.
* sdf: when a second molecule sharing atoms with a first one is removed, only nm, phim,
  cthetm, chim and iasave are shifted (not xcent.. or dlmnmc/chiphi).
* sdfcube `outcube`: `fc1=1.0-fc1` reads the uninitialised `fc1`; gfortran -O2 folds it to
  `fc1=1.0` (Rust sdfcube.rs, checked there in the object code).  Marked `NB`.
* sdfcube: for identical molecules there is no `iref=iref+2`, so `averagemol` is read from
  the `axiss1` symbol.
* sdfcube `sdfcube`: `vol` is a local (no include declares it); `xstepd/ystepd/zstepd` are
  declared `real` locally but are the sdfcube.inc COMMON variables.
* `(-1)**m` in form_sum_prod_dlmn_cgcoeff is an integer power (+-1 also for m < 0); written
  `1 if m % 2 == 0 else -1` (Python's `(-1)**m` is a float for m < 0).
* form_eulerangles: `at`, `bt` are uninitialised when the secondary axis is undefined; they
  are then only used with zero norms, so the result does not depend on them (zeros here).
* Uninitialised `ierrout` before setup_input_file, `ierr` before the iostat reads: 0.
* `itimestotal*100/ntimestotal` in the progress line is a division by zero when there are no
  pairs (ntimestotal = 0): the Fortran crashes (SIGFPE), the Python raises ZeroDivisionError.
  No test case gets there.

## Performance

The innermost loop of form_sum_prod_dlmn_cgcoeff (`_sum_products`) is the whole cost: per
term two libm calls (`cosf`, `sinf`) and ~10 float32 operations.  It is a scalar loop in the
Fortran order.  The rows it reads are copied once per call into Python lists of the same
numpy.float32 elements (list indexing is ~3x faster than ndarray indexing); gofz(iref,id) is
accumulated in a local and stored after each id (the same sequence of float32 additions).
The ctypes call of `cosf`/`sinf` costs ~1-1.7 us each, more than half of the time; a faster
exact cosf/sinf in `fortran/fmath.py` would speed up sharm/sdf almost 2x (not done: not my
package).

## Verification

`tools/run_ex_py.py` (resume mode), every test case using sharm/sdf/sdfcube (work folders
`work/Q/exNN`).  **Every sharm, sdf and sdfcube step is IDENTICAL** (all files), including
the two `test` cases (ex72 orient.SHARM.dat, ex78 dmsotest.SDF.dat) and the accumulation
steps (nsumt = 5: sdfcube reads back its cube file, sharm/sdf average with the .h01 file).
The stdout of the first ex70 sharm step was also compared by hand: identical apart from the
folder name.  Times (machine shared with a 10-core sweep; Fortran / Python, s):

| case | step: program (input) | result | Fortran | Python |
|---|---|---|---:|---:|
| ex70 | 1, 16: sharm (water) | IDENTICAL | 0.47 | 71.5, 69.5 |
| ex70 | 4, 7, 10, 13: sharm (nsumt 1-4, no calculation) | IDENTICAL | 0.07-0.13 | 1.1-1.3 |
| ex71 | 1: sharm (dmsowat) | IDENTICAL | 1.22 | 215.7 |
| ex72 | 1: sharm (watdmso) | IDENTICAL | 2.94 | 631.8 |
| ex72 | 2, 6: sharm (orient, `test`) | IDENTICAL | 0.47, 0.57 | 79.2, 82.8 |
| ex72 | 5: sharm (watdmso, no calculation) | IDENTICAL | 0.07 | 0.8 |
| ex73 | 1, 11: sdf (water) | IDENTICAL | 0.22, 0.27 | 31.8, 34.8 |
| ex73 | 3, 5, 7, 9: sdf (no calculation) | IDENTICAL | 0.07-0.12 | 0.9-1.0 |
| ex74 | 1: sdf (sc) | IDENTICAL | 0.12 | 4.7 |
| ex75 | 1, 6: sdfcube (water) | IDENTICAL | 0.17, 0.12 | 5.2, 5.2 |
| ex75 | 2-5: sdfcube (no calculation) | IDENTICAL | 0.06-0.07 | 0.7-0.8 |
| ex76 | 1: sdfcube (dmsowat) | IDENTICAL | 0.12 | 4.1 |
| ex77 | 1: sharm (dmso) | IDENTICAL | 3.03 | 663.6 |
| ex78 | 1, 4: sdf (dmsotest, `test`) | IDENTICAL | 0.12, 0.17 | 3.4, 3.7 |
| ex78 | 2: sdf (watso) | IDENTICAL | 0.42 | 51.1 |
| ex78 | 5: sdf (watso, no calculation) | IDENTICAL | 0.13 | 1.2 |
| ex79 | 1: sharm (axes) | IDENTICAL | 0.22 | 31.3 |
| ex79 | 2: sdfcube (axes) | IDENTICAL | 0.07 | 3.8 |

Other programs in these cases: mapgr and mapgrentropy (other packages) were IDENTICAL in
ex70, ex71, ex72, ex77 (run later); in the earlier ex73, ex74, ex78 runs they were not yet
translated (python exit 1, NotTranslated), so those steps were not compared.

## Open issues

* Speed: sharm with DMSO (ex72 step 1, ex77) takes ~11 min against 3 s (x220); see
  Performance.  A faster exact `cosf`/`sinf` in the compat layer is the main lever.
* No signature exceptions; no bugs found in other packages' modules.
