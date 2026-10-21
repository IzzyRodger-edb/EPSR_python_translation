# Package O: analysis programs (wave 2)

Fortran: `partials.f`, `coord.f`, `triangles.f`, `torangles.f`, `fluctuations.f`,
`local_density.f`, `clusters.f`, `voids.f`.

Output (one program module each, `def main()` = the Fortran `program standalone`):

| file | routines |
|---|---|
| `src/epsr26/programs/partials.py` | main, run_partials, find_partials |
| `src/epsr26/programs/coord.py` | main, run_coord, find_coord |
| `src/epsr26/programs/triangles.py` | main, run_triangles, find_triangles (+ helper `_add_flat`) |
| `src/epsr26/programs/torangles.py` | main, run_torangles, find_torangles (+ `_add_flat`) |
| `src/epsr26/programs/fluctuations.py` | main, run_fluctuations, find_fluctuations |
| `src/epsr26/programs/local_density.py` | main, run_clusters, find_clusters |
| `src/epsr26/programs/clusters.py` | main, run_clusters, find_clusters |
| `src/epsr26/programs/voids.py` | main, run_voids, find_voids, ranmine, occupypix, markpixel, search_voids, get_void_list, spheregr (+ `_pair_loop`, `_wrap_table`) |

The duplicated names (`run_clusters`/`find_clusters` in clusters.f and local_density.f,
`program standalone` everywhere) each live in their own program module. The `interface`
block for `find_voids` in `run_voids` (SIGNATURES lists it as a second `find_voids` at
voids.f:24) is only a declaration; there is one `find_voids`.

No signature exceptions: every routine here has no scalar dummies except the voids helpers,
which return what SIGNATURES lists (`occupypix`/`markpixel` -> npixfound, `search_voids`
-> nvoids, `get_void_list` -> nvoidlist, `ranmine`, `spheregr`).

## Verification (tools/run_ex_py.py, resume mode; every output file compared)

| case | program | steps | result | Fortran s/step | Python s/step |
|---|---|---|---|---|---|
| ex50 | partials (all pairs) | 2 | IDENTICAL | 0.4 | 63-68 |
| ex51 | partials (atom-c, atom-z, atom-r) | 6 | IDENTICAL | 0.12-0.18 | 7-16 |
| ex52 | coord (4 inputs incl. slab/cylinder/centre) | 8 | IDENTICAL | 0.07-0.23 | 1.9-5.8 |
| ex53 | triangles (4 inputs, ltyp 1 and bonded-average) | 8 | IDENTICAL | 0.12-0.47 | 7-31 |
| ex54 | torangles | 2 | IDENTICAL | 0.10-0.13 | 1.2-2.0 |
| ex55 | voids (2 distributions, fresh ran1) | 1 | IDENTICAL (run by hand: all 4 files + stdout) | 8 | 2266 |
| ex56 | voids (include list, accumulating) | 2 | IDENTICAL | 5.9-6.0 | 1398-1472 |
| ex57 | fluctuations | 3 | IDENTICAL | 0.12-0.26 | 13-14 |
| ex58 | local_density (multi-bond clusters) | 2 | IDENTICAL | 0.32 | 19-22 |
| ex59 | fluctuations | 2 | IDENTICAL | 0.12 | 7 |
| ex60 | clusters | 2 | IDENTICAL | 0.12-0.20 | 9.5-9.8 |
| ex61 | clusters | 2 | IDENTICAL | 0.12 | 5.5-5.7 |

The standard output (logs `f<n>.log` vs `p<n>.log`) is also identical in every step
(apart from the work-folder name in the "Home folder" line).

## Decisions and quirks

* `read(symbol_value(iref),'(a)') atomc` with `character*3 atomc` is `fstr(value, 3)`.
* `index(atom1(id),' ')-1` on the `character*4` atom labels is `findex_padded(.., ' ', 4) - 1`.
* Records written with `a3` from `character*4` labels pass `Ch(label, 4)` (the leftmost 3
  characters are written; a plain shorter string would be right-justified).
* Records written to unit 14 after `save_files_others(14)` has closed it go to `fort.14`
  (units handles the implicit open), as in gfortran.
* **triangles** uses gfortran's intrinsic `sind`/`cosd` (not declared external; Rust
  anal_common.rs): `sind_intrinsic`/`cosd_intrinsic`. **torangles** uses
  auxiliary_routines' `vecpro`/`scapro` (the copies in torangles.f are commented out).
* triangles/torangles: an angle of exactly 180 deg gives bin `nz+1`; `_add_flat` adds with
  the Fortran column-major address (Rust `add_flat`), so `nz = mz` would spill into the
  next column as in the Fortran instead of raising an IndexError.
* torangles: `ltypfmt`/`ltyp1fmt` are uninitialised in the Fortran when ltyp >= 1000; blank
  is used (unreachable in the tests). `toranglesvalues.txt` (unit 29) is never closed in
  the Fortran; it is flushed at program end (units.close_all), as gfortran does.
* fluctuations: `frac=1.0/real(nsumt)` is single precision converted to real*8.
* local_density: `damax(id)` of `clusters_multibond.inc` is never set (zero COMMON), so the
  percolation test uses `diftest = 1.0` (as the Rust). The error message of the
  "numbers of atoms and cutoffs" check is a character constant continued over two fixed-form
  lines; gfortran pads the first line to column 72, so the message has **two** blanks
  before "bonds" (checked with gfortran 15; the Rust port prints one -- untested branch).
* voids: `npixfound` (uninitialised if no atom occupies pixels) and `iscount` (incremented
  before it is first set) start at 0. If there is no unoccupied pixel the Fortran loops
  forever choosing a random origin pixel; the Python raises FortranStop(1) with a message
  (as the Rust reports it unsupported).
* voids: `voidpix` (integer*1) and `voidnum` (integer*2) are Python lists (values 0/1 and
  <= maxvoid cannot overflow); `get_void_list` evaluates the periodic wrap
  `it-npixtot*nint(real(it)/real(npixtot))` through a table built with the same float32
  expression (`_wrap_table`), once per index instead of 27 times per pixel. The inner
  pixel-pair loop of `find_voids` is in `_pair_loop` (same statements, same order) only to
  keep its variables local.

## Performance

Python is ~25-180x slower than the Fortran for these programs (numpy.float32 scalar
arithmetic, ~5 us per `ran1` call, ~4 us per `minimagec` call). voids is the extreme case:
its pair-correlation loop draws 3 random numbers per void pixel for each of 100 random
origins (~80 million iterations for ex55), ~31 us per iteration in Python, so a voids step takes 23-38 minutes against the Fortran's
6-8 s. A bit-exact speed-up is possible (bulk ran1 in plain ints plus elementwise float32
NumPy for the pair loop, with float32 count saturation at 2**24 for the histogram), but it
would no longer be a line-by-line translation, so it is not done.

## Bugs found in other modules

None that blocked this package.
