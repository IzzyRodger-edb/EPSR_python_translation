# Package J: epsr_standalone_rev3.f lines 2523-4736

Output: `src/epsr26/programs/epsr_standalone_rev3_part2.py`, to be merged into
`programs/epsr_standalone_rev3.py` between runepsr17 (package I) and gensq... (package K).

All 24 routines are translated in full: interp, inita (including the coarse-graining block), j1qr,
j1qr_mod, gaussianqr2, lorentzianqr2, pofr, initg, entrpa, entrpb, entrpc, sumgr, sumsq, sumpe,
initpt, ptgnsq, repulsivept, norma, openf, get_wts, wtype, wxinit, wxtype, wxcreate.
Return conventions follow SIGNATURES.md exactly; no exceptions were needed.

## Checked against the Fortran reference (byte-identical)
- `get_wts` X-ray branch (wxinit/wxtype/wxcreate): `X-rayf0x.dat` and `X-raywts.dat` of
  EPSR-test-ex9 reproduced byte for byte (state set by hand: types, cfraccg, pair order).
- `inita` r grid + `norma`: the `.uni` file of ex9 is byte-identical; `ncoeffp = 116` as in the `.pcof`.
- `inita` coarse-graining block: `formfactors.txt`, `intraformfactors.txt`,
  `intraformfactorsself.txt` of ex100 are byte-identical (natomsphere/rhosphere/radmin/radmax set
  as `ato_read_and_save.f` does, scalefacmol = 1, selfbroadpwr = 0).
- Everything else was smoke-tested only (runs, plausible values); the end-to-end test is still needed.

## Names from the sibling part-files (bare-name calls, each line has `# noqa: F821`)
`uniform`, `generate_rp`, `generate_ep`, `calc_gaussq`, `calc_pwrexpq`, `fit_expansion`,
`invert_cwt` (all package K, `_part3.py`). The `noqa` comments can be dropped after the merge.

## Imports
- The compat layer is imported by explicit names (`from ..fortran import f32, r4, ...`) instead of
  `import *`, so that `ruff --select F` reports real undefined names (with `*` every compat name is an
  F405 and sibling names are not reported). Merging with part1's `import *` is harmless.
- Uses `common.sumvalues` (package I's hand-written COMMON /sumvalues/; fc1, fc2 set by runepsr17),
  `common.molecom_mod`, `library.reallocation_routines` (reallocate1d_r returns the array),
  `library.auxiliary_routines.r_int`, `library.file_and_directory_commands.ireff`,
  `library.ato_manipulation_rev4.truncate/truncatec`, `library.parse.parse`.

## Decisions / uninitialised variables
- `inita`: `ireff(iqminref, ifilecount)` reads the undeclared, uninitialised local `ifilecount`;
  it is not used by ireff for iqminref <= nsymboluptodata, so 1 is passed (the Rust passes 1 too).
- `inita`: `volume` is the clustercom COMMON variable (inita includes clustercom.inc), so inita sets
  `clustercom.volume = alatx*blaty*clatz` as a side effect.
- `inita`: `sqrt(2.*alog(2.0))` is computed as `sqrtf(2*logf(2))` (gfortran folds it in single
  precision; identical result).
- `inita` coarse graining: `fracntermssff` is allocated but missing from the Fortran's
  `deallocate(...)` list, so a second inita call with coarse graining in the same program is a
  gfortran runtime error. Reproduced as `FortranStop(2, "Fortran runtime error: Attempting to allocate
  already allocated variable 'fracntermssff'")` (never happens in practice: inita runs once per program).
- `inita` coarse graining: unit 10 is reopened on `intraformfactorsself.txt` without closing
  `intraformfactors.txt`; Fortran (and `units.open`) closes the first file implicitly.
- `initg`: `lptin` is uninitialised; 0 is passed (r_int only assigns it).
- `ptgnsq`: `ncfirst`, `rejrat`, `chisqnew` uninitialised locals start at 0 (callees assign them).
- `ptgnsq` delfit.dat: the `if(abs(qtemp(i,id)).gt.0.0)` test uses `id` after the loop, i.e.
  `nfile+1`, an untouched (zero) column of the local `qtemp`, so the second branch is always taken
  (as in the Rust). `qtemp` is allocated with `mfile+1` columns so that nfile = mfile cannot raise.
  Local arrays (`qtemp`, `fittemp`) are zero-initialised; the rfactorep loop reads `fittemp` beyond
  `nqmax(id)` (0, matching the Rust).
- `ptgnsq` uses `weights.cwtpot` (COMMON, (mpairs,mpairs)) exactly as the Fortran;
  `get_wts` fills it as `cwtpot(in,ic)` for invert_cwt and then overwrites `cwtpot(ic,in)`.
- `entrpb`: FORTRAN_BUGS.md section 2 bug reproduced (loops `id=1,nfile`, not npairs).
- `get_wts` X-ray header: `read(10,'(1x,a3,1x,a4)') a1,a2` with a2 character*3 keeps the rightmost
  3 characters of the 4-character field (private helper `_read_a1a2`). Two more private helpers
  (`_read_wread`, `_read_xray_coeffs`) implement the `iostat=` list-directed reads; items read
  before an error keep their new values, as in Fortran.
- `get_wts`: `datanorm` is an integer array, so `datanorm(nxfile)=rhofac(in)` truncates (ifix).
- The dummy argument `in` (Python keyword) is spelled `in_` in interp, wtype, wxtype and in loops.
- Commented-out Fortran blocks of more than a few lines (inita .uni read, ptgnsq X-ray
  reconstruction, repulsivept pcoefexmax update) are referred to by line number rather than copied.

## Helper wished for in the compat layer
- An `iostat` helper for reads (`ierr = units.iostat(lambda: ...)`) would remove the
  try/except boilerplate repeated in many routines.
