# Validation results (real measured data)

Generated from `work/val/<case>/summary.json` (`tools/run_ex_py.py validation/<case> --max-niter 2`),
2026-10-21, gfortran 15.2 `-O2 -fno-tree-vectorize`, glibc 2.43. Resume mode: every step starts from
the Fortran folder before that step. EPSR steps longer than 2 iterations are compared on a copy
with `niter 2` (both Fortran and Python); the Fortran chain itself runs the full step. A step is
IDENTICAL when every file in the folder matches byte for byte (date in `.out` and the timing
columns of `fort.35` masked). Times are wall-clock on a shared machine (other runs in parallel).

## Water (H2O, D2O, HDO 50:50; SANDALS)

| step | program | args | compared | result | Fortran [s] | Python [s] | ratio |
|---:|---|---|---|---|---:|---:|---:|
| 1 | makemole | Water | full step | IDENTICAL | 0.01 | 0.9 | 62x |
| 2 | mixato |  | full step | IDENTICAL | 0.07 | 1.3 | 18x |
| 3 | randomise | waterbox | full step | IDENTICAL | 0.12 | 1.3 | 11x |
| 4 | fmole | waterbox 100 | full step | IDENTICAL | 0.32 | 19.7 | 61x |
| 5 | nwts | H2O.NWTS.dat | full step | IDENTICAL | 0.12 | 2.1 | 18x |
| 6 | nwts | D2O.NWTS.dat | full step | IDENTICAL | 0.12 | 1.5 | 13x |
| 7 | nwts | HDO.NWTS.dat | full step | IDENTICAL | 0.17 | 1.2 | 7x |
| 8 | epsr | waterbox.EPSR.inp | full step | IDENTICAL | 2.69 | 420.5 | 156x |
| 9 | epsr | waterbox.EPSR.inp | niter 2 of 99 | IDENTICAL | 4.40 | 768.8 | 175x |
| 10 | epsr | waterbox.EPSR.inp | niter 2 of 100 | IDENTICAL | 6.84 | 1208.4 | 177x |
| 11 | epsr | waterbox.EPSR.inp | niter 2 of 100 | IDENTICAL | 5.46 | 852.1 | 156x |

11/11 steps identical.

## Amorphous SiO2 (NIMROD)

| step | program | args | compared | result | Fortran [s] | Python [s] | ratio |
|---:|---|---|---|---|---:|---:|---:|
| 1 | makemole | Si | full step | IDENTICAL | 0.01 | 0.7 | 49x |
| 2 | makemole | O | full step | IDENTICAL | 0.01 | 0.7 | 104x |
| 3 | mixato |  | full step | IDENTICAL | 0.17 | 1.3 | 7x |
| 4 | randomise | silica | full step | IDENTICAL | 0.23 | 2.1 | 9x |
| 5 | nwts | SiO2.NWTS.dat | full step | IDENTICAL | 0.27 | 2.1 | 8x |
| 6 | epsr | silica.EPSR.inp | full step | IDENTICAL | 5.28 | 657.7 | 125x |
| 7 | epsr | silica.EPSR.inp | niter 2 of 99 | IDENTICAL | 7.09 | 1043.5 | 147x |
| 8 | epsr | silica.EPSR.inp | niter 2 of 100 | IDENTICAL | 9.61 | 1261.8 | 131x |
| 9 | epsr | silica.EPSR.inp | niter 2 of 100 | IDENTICAL | 5.76 | 757.8 | 132x |

9/9 steps identical.

## Cu-BTC framework (NIMROD, Bragg peaks)

| step | program | args | compared | result | Fortran [s] | Python [s] | ratio |
|---:|---|---|---|---|---:|---:|---:|
| 1 | readpdb | CuBTC.pdb | full step | IDENTICAL | 0.04 | 1.1 | 28x |
| 2 | changeato | CuBTC.ato | full step | IDENTICAL | 0.07 | 1.1 | 15x |
| 3 | nwts | CuBTC.NWTS.dat | full step | IDENTICAL | 0.07 | 0.7 | 10x |
| 4 | epsr | CuBTC.EPSR.inp | full step | IDENTICAL | 2.94 | 718.7 | 244x |
| 5 | epsr | CuBTC.EPSR.inp | niter 2 of 99 | IDENTICAL | 4.88 | 1277.4 | 262x |
| 6 | epsr | CuBTC.EPSR.inp | niter 2 of 100 | IDENTICAL | 6.96 | 1289.0 | 185x |
| 7 | epsr | CuBTC.EPSR.inp | niter 2 of 100 | IDENTICAL | 4.73 | 841.4 | 178x |

7/7 steps identical.
