# Package P: chains, rings, lsi, mapgr, mapgrentropy

| Fortran | Python | routines |
|---|---|---|
| chains_routines.f | library/chains_routines.py | match, shortest_path |
| chains_sp_criterion.f | programs/chains_sp_criterion.py | main (program standalone), run_chains, find_chains |
| rings_sp_criterion.f | programs/rings_sp_criterion.py | main, run_rings, find_rings |
| lsi.f | programs/lsi.py | main, run_lsi, find_lsi |
| plotsharm_routines.f | library/plotsharm_routines.py | smoosharm, quadinterp2d, quadinterp, sumoverm2, suminterp |
| mapgrquad.f | programs/mapgrquad.py | main, run_mapgr |
| mapgrentropy.f | programs/mapgrentropy.py | main, run_mapgrentropy |

Every routine is translated completely; signatures follow SIGNATURES.md (no exceptions
needed). `setup_symbols_*` (package L) are imported inside the run_* routines.

## Verification

End to end with `tools/run_ex_py.py` (resume mode), one case at a time. Every step of my
programs is IDENTICAL (all output files), and the stdout logs match the Fortran's apart
from the "Home folder" line (and gfortran's IEEE_UNDERFLOW note on stderr for mapgr).

| case | step(s) | program | Fortran | Python |
|---|---|---|---|---|
| ex62 | 1, 2 | chains | 1.02, 0.87 s | 54.6, 53.8 s (machine busy with a parallel sweep) |
| ex63 | 1, 2 | chains | 0.20, 0.22 s | 9.5, 9.1 s |
| ex64 | 1, 2 | rings | 0.69, 0.63 s | 40.6, 39.7 s |
| ex65 | 1, 2 | rings | 0.32, 0.32 s | 27.9, 29.9 s |
| ex66 | 1, 2, 3 | lsi | 0.38, 0.22, 0.22 s | 12.3, 12.1, 11.1 s |
| ex67 | 1, 2 | lsi | 0.47, 0.42 s | 29.6, 28.8 s |
| ex70 | 2,5,8,11,14,17 | mapgr | 0.02-0.03 s | 0.6-0.8 s |
| ex70 | 3,6,9,12,15,18 | mapgrentropy | 0.03-0.12 s | 8.1-9.5 s |
| ex71 | 2, 3 / 4 | mapgr / mapgrentropy | 0.07, 0.03 / 0.12 s | 1.3, 1.2 / 20.1 s |
| ex72 | 3, 7 / 4, 8 | mapgr / mapgrentropy | 0.04, 0.07 / 0.07, 0.07 s | 0.9, 0.8 / 6.5, 6.8 s |
| ex73 | 2,4,...,12 | mapgr | 0.02-0.12 s | 0.6-0.9 s |
| ex74 | 2 | mapgr | 0.02 s | 1.24 s |
| ex77 | 2 / 3 | mapgr / mapgrentropy | 0.07 / 0.12 s | 1.6 / 22.1 s |
| ex78 | 3, 6 | mapgrentropy | 0.04, 0.04 s | 2.0, 2.2 s |

(The mapgr/mapgrentropy cases were run with `--steps` = my steps only, the sharm/sdf
steps being package Q's.)

plotsharm_routines.f was also checked in isolation: a gfortran driver (reference flags,
linked with the mapgrquad object list) calling quadinterp2d, quadinterp, suminterp (both
branches of the phi-range reset), smoosharm and sumoverm2 (positive/negative/zero n1)
gives bit-identical results to the Python (hex comparison). quadinterp, quadinterp2d and
suminterp are not used by any program in the table.

## Quirks reproduced (all as in the Rust port)

* chains/rings: the "Not found enough items" message says `rings>` in both programs.
  `icount`, `sptest`, `spmin`, `spmax` are the COMMON (chains.inc) variables, so in
  find_rings the ring-closure check overwrites `spmin` that the stack walk then reuses.
* chains: the mean/std records are written to unit 14 after `save_files_others(14)` has
  closed it, so they go to `fort.14` (units.py opens `fort.N` for an unopened unit).
* lsi: `gofz=0.0` zeroes the whole array at the start of every distribution, so only the
  last distribution keeps non-zero histograms; `fc = 1/nsumt` uses the input nsumt.
* lsi: a neighbour for which no atom pair matched has `iwmin = 0` and `ia2min = 0`; the
  Fortran then reads `dfs1(id,0)`, `dfs2(id,0)` (= `dff2(id,mdist)`, `dfs1(id,mdist)` in
  COMMON /lsi_r/) and `xa(0)`, `ya(0)`, `za(0)` (= `stepmi`, `xa(matom)`, `ya(matom)` in
  COMMON /molecom_r/). Emulated by `_dfs_range` and `_atom_xyz` in lsi.py.
* lsi: `max(anglefound,-1.0)` / `min(...,+1.0)` with gfortran's NaN handling (`_amax`,
  `_amin`).
* lsi: `titledat(idref7..10)` use `atom3(id,1)(1:ilent1)` (atom1's length), as written.
* mapgr: `grwerr(:,1)` is never set (stack garbage in the Fortran, 0 here); masked by
  compare_outputs.py.
* mapgr/mapgrentropy: `cosd` in run_mapgr is EPSR's own (`external cosd`); `cosd`/`sind`
  in sumoverm2 are not declared external, so gfortran >= 14 uses its intrinsics
  (`cosd_intrinsic`/`sind_intrinsic`).
* mapgrentropy: the min/max radius and angle variables are uninitialised unless a
  minimum/maximum is found; 0 (as in the reference).
* `(-1)**n` (integer) is computed by `_m1pow` (correct also for negative n).

## Differences from the Rust port (not exercised by the tests)

* mapgr: the view angles are read with `read(...,*,iostat=itest)`; when the value holds
  fewer than 5 numbers gfortran gives an end-of-file (iostat -1) and the view is skipped.
  The Python does that; the Rust treats missing values as "unchanged" with itest = 0.

## Open issues / notes

* lsi with more than 3 distributions (`ndist3 = 10*ndist > mdist = 30`) indexes gofz out
  of bounds in the Fortran; the Python raises IndexError there (the Rust reports
  "unsupported").
* Speed: the cost is the per-pair atom loops (find_chains/find_rings/find_lsi call
  `match` ~15 million times for ex62). `match` skips the `fstr` truncation when the
  strings are already short enough and uses `str.find` for `index`; otherwise the code
  follows the Fortran loops as written.
* Large unused local arrays of mapgrquad/mapgrentropy (`out(0:mpixel3)`, `lvl`, `th`,
  `rgbbal`, ...) are not created (listed in comments).
