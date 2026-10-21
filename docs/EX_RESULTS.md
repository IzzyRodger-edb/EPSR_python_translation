# Test-problem results

All 115 cases of `test_problems/200DMSO400Water`, Fortran (gfortran 15.2 `-O2 -fno-tree-vectorize`,
glibc 2.43) against the Python, `tools/run_ex_py.py` resume mode (each step starts from the Fortran
folder before it), 2026-10-21. A step is identical when every file matches byte for byte (date in
`.out` and timing columns of `fort.35` masked). Times are wall-clock with up to 12 runs in parallel.

| case | programs | identical steps | Fortran [s] | Python [s] | slowdown | note |
|---|---|---:|---:|---:|---:|---|
| EPSR-test-ex1 | epsr | 3/3 | 11.5 | 1836.9 | 159x |  |
| EPSR-test-ex2 | epsr | 6/6 | 22.6 | 3648.7 | 161x |  |
| EPSR-test-ex3 | epsr | 6/6 | 22.7 | 3671.6 | 162x |  |
| EPSR-test-ex4 | epsr | 2/2 | 7.3 | 1255.1 | 172x |  |
| EPSR-test-ex5 | epsr | 3/3 | 9.1 | 1314.8 | 145x |  |
| EPSR-test-ex6 | epsr | 3/3 | 12.3 | 1849.0 | 151x |  |
| EPSR-test-ex7 | epsr | 3/3 | 11.1 | 1834.1 | 165x |  |
| EPSR-test-ex8 | epsr | 4/4 | 11.3 | 1696.4 | 150x |  |
| EPSR-test-ex9 | epsr | 2/2 | 12.3 | 2046.8 | 167x |  |
| EPSR-test-ex10 | epsr | 2/2 | 12.0 | 2043.6 | 170x |  |
| EPSR-test-ex11 | epsr | 2/2 | 7.9 | 1222.5 | 155x |  |
| EPSR-test-ex12 | epsr | 2/2 | 7.3 | 1222.8 | 167x |  |
| EPSR-test-ex13 | epsr | 3/3 | 11.2 | 1818.2 | 162x |  |
| EPSR-test-ex14 | epsr | 2/2 | 9.4 | 1413.8 | 151x |  |
| EPSR-test-ex15 | epsr | 2/2 | 8.8 | 1420.2 | 162x |  |
| EPSR-test-ex20 | nwts | 3/3 | 0.5 | 5.3 | 11x |  |
| EPSR-test-ex21 | nwts | 4/4 | 0.4 | 5.6 | 13x |  |
| EPSR-test-ex22 | nwts | 3/3 | 0.5 | 4.7 | 10x |  |
| EPSR-test-ex23 | xwts | 3/3 | 0.4 | 3.2 | 9x |  |
| EPSR-test-ex24 | epsrwts | 4/4 | 0.8 | 4.8 | 6x |  |
| EPSR-test-ex25 | epsrwtsx | 3/3 | 0.4 | 3.5 | 9x |  |
| EPSR-test-ex26 | changeato | 3/3 | 0.4 | 7.0 | 19x |  |
| EPSR-test-ex27 | changeato | 2/2 | 0.2 | 6.6 | 44x |  |
| EPSR-test-ex28 | makemole | 2/2 | 0.1 | 1.1 | 13x |  |
| EPSR-test-ex29 | makemole | 3/3 | 0.0 | 1.0 | 22x |  |
| EPSR-test-ex30 | makeato | 2/2 | 0.2 | 2.9 | 15x |  |
| EPSR-test-ex31 | mixato | 1/1 | 0.1 | 0.7 | 9x |  |
| EPSR-test-ex32 | mixato | 1/1 | 0.1 | 0.6 | 5x |  |
| EPSR-test-ex33 | fmole,mixato | 3/3 | 0.1 | 3.0 | 25x |  |
| EPSR-test-ex34 | fmole | 2/2 | 0.1 | 1.3 | 25x |  |
| EPSR-test-ex35 | fmole,mixato,randomise | 5/5 | 0.3 | 3.8 | 15x |  |
| EPSR-test-ex36 | addato | 2/2 | 0.4 | 3.4 | 9x |  |
| EPSR-test-ex37 | makelattice | 3/3 | 0.0 | 1.4 | 31x |  |
| EPSR-test-ex38 | makelatticeato | 4/4 | 0.1 | 2.6 | 26x |  |
| EPSR-test-ex39 | fmole,makelatticeato,mixato,randomise | 4/4 | 0.2 | 3.4 | 18x |  |
| EPSR-test-ex40 | readxyz | 2/2 | 0.0 | 0.8 | 27x |  |
| EPSR-test-ex41 | readxyz | 4/4 | 0.0 | 2.2 | 47x |  |
| EPSR-test-ex42 | writexyz | 5/5 | 0.4 | 7.1 | 18x |  |
| EPSR-test-ex43 | writexyz | 3/3 | 0.1 | 3.6 | 39x |  |
| EPSR-test-ex44 | readpdb | 1/1 | 0.1 | 0.8 | 10x |  |
| EPSR-test-ex45 | readpdb | 1/1 | 0.0 | 0.9 | 26x |  |
| EPSR-test-ex46 | bonds | 4/4 | 0.2 | 4.3 | 19x |  |
| EPSR-test-ex47 | alignmole | 2/2 | 0.1 | 1.2 | 23x |  |
| EPSR-test-ex48 | alignmole | 3/3 | 0.1 | 2.1 | 33x |  |
| EPSR-test-ex49 | readpdb | 1/1 | 0.0 | 0.4 | 31x |  |
| EPSR-test-ex50 | partials | 2/2 | 0.8 | 145.7 | 177x |  |
| EPSR-test-ex51 | partials | 6/6 | 0.9 | 73.4 | 83x |  |
| EPSR-test-ex52 | coord | 8/8 | 1.1 | 50.4 | 47x |  |
| EPSR-test-ex53 | triangles | 8/8 | 1.4 | 143.4 | 105x |  |
| EPSR-test-ex54 | torangles | 2/2 | 0.2 | 3.2 | 17x |  |
| EPSR-test-ex55 | voids | 1/1 | 11.0 | 1523.6 | 139x |  |
| EPSR-test-ex56 | voids | 2/2 | 11.1 | 1616.9 | 146x |  |
| EPSR-test-ex57 | fluctuations | 3/3 | 0.4 | 43.2 | 108x |  |
| EPSR-test-ex58 | local_density | 2/2 | 0.5 | 40.7 | 75x |  |
| EPSR-test-ex59 | fluctuations | 2/2 | 0.2 | 19.8 | 80x |  |
| EPSR-test-ex60 | clusters | 2/2 | 0.3 | 24.0 | 69x |  |
| EPSR-test-ex61 | clusters | 2/2 | 0.2 | 15.6 | 63x |  |
| EPSR-test-ex62 | chains | 2/2 | 1.9 | 89.8 | 48x |  |
| EPSR-test-ex63 | chains | 2/2 | 0.6 | 23.4 | 39x |  |
| EPSR-test-ex64 | rings | 2/2 | 1.0 | 68.1 | 65x |  |
| EPSR-test-ex65 | rings | 2/2 | 0.6 | 50.1 | 85x |  |
| EPSR-test-ex66 | lsi | 3/3 | 0.6 | 28.2 | 44x |  |
| EPSR-test-ex67 | lsi | 2/2 | 0.8 | 51.2 | 65x |  |
| EPSR-test-ex70 | mapgr,mapgrentropy,sharm | 18/18 | 1.6 | 166.7 | 102x |  |
| EPSR-test-ex71 | mapgr,mapgrentropy,sharm | 4/4 | 1.3 | 209.3 | 163x |  |
| EPSR-test-ex72 | mapgr,mapgrentropy,sharm | 8/8 | 3.8 | 619.9 | 164x |  |
| EPSR-test-ex73 | mapgr,sdf | 12/12 | 0.8 | 62.0 | 74x |  |
| EPSR-test-ex74 | mapgr,sdf | 2/2 | 0.1 | 4.1 | 33x |  |
| EPSR-test-ex75 | sdfcube | 6/6 | 0.5 | 11.7 | 25x |  |
| EPSR-test-ex76 | sdfcube | 1/1 | 0.1 | 3.3 | 28x |  |
| EPSR-test-ex77 | mapgr,mapgrentropy,sharm | 3/3 | 3.0 | 554.7 | 186x |  |
| EPSR-test-ex78 | mapgrentropy,sdf | 6/6 | 0.5 | 54.2 | 114x |  |
| EPSR-test-ex79 | sdfcube,sharm | 2/2 | 0.3 | 29.2 | 102x |  |
| EPSR-test-ex80 | dockato | 2/2 | 0.1 | 0.8 | 15x |  |
| EPSR-test-ex81 | dockato | 2/2 | 0.0 | 0.7 | 33x |  |
| EPSR-test-ex82 | dockato | 2/2 | 0.0 | 1.5 | 97x |  |
| EPSR-test-ex83 | dockato | 3/3 | 0.0 | 1.1 | 46x |  |
| EPSR-test-ex84 | dockato | 3/3 | 1.1 | 26.6 | 25x |  |
| EPSR-test-ex85 | dockato | 2/2 | 0.5 | 8.7 | 19x |  |
| EPSR-test-ex86 | dockato | 2/2 | 0.0 | 1.1 | 31x |  |
| EPSR-test-ex87 | dockato | 2/2 | 0.0 | 0.7 | 36x |  |
| EPSR-test-ex90 | epsr | 3/3 | 11.6 | 1804.6 | 156x |  |
| EPSR-test-ex91 | epsr | 3/3 | 13.6 | 2201.1 | 162x |  |
| EPSR-test-ex92 | epsr | 3/3 | 13.7 | 2259.2 | 165x |  |
| EPSR-test-ex93 | epsr | 3/3 | 11.1 | 1934.9 | 174x |  |
| EPSR-test-ex94 | epsr | 3/3 | 12.0 | 1989.8 | 165x |  |
| EPSR-test-ex95 | epsr | 3/3 | 13.1 | 2160.1 | 165x |  |
| EPSR-test-ex96 | epsr | 3/3 | 12.1 | 2167.0 | 180x |  |
| EPSR-test-ex97 | epsr | 3/3 | 12.6 | 2157.5 | 172x |  |
| EPSR-test-ex98 | epsr | 3/3 | 12.4 | 2208.9 | 177x |  |
| EPSR-test-ex99 | epsr | 3/3 | 13.0 | 2241.5 | 172x |  |
| EPSR-test-ex100 | epsr | 2/2 | 8.6 | 1433.2 | 166x |  |
| EPSR-test-ex101 | epsr | 2/2 | 8.3 | 1430.8 | 173x |  |
| EPSR-test-ex102 | epsr | 2/2 | 6.8 | 1152.4 | 170x |  |
| EPSR-test-ex103 | epsr | 2/2 | 7.1 | 1384.3 | 196x |  |
| EPSR-test-ex105 | epsr | 4/4 | 14.9 | 2314.1 | 155x |  |
| EPSR-test-ex106 | epsr | 2/2 | 7.6 | 1124.1 | 148x |  |
| EPSR-test-ex107 | epsr | 2/2 | 8.2 | 1255.2 | 154x |  |
| EPSR-test-ex108 | epsr | 1/2 | 8.0 | 1108.3 | 139x | FORTRAN_BUGS.md 1.3 fix: differs from the Fortran by design; every pass identical to the Rust port (tools/check_fixed_case.py) |
| EPSR-test-ex109 | epsr | 6/6 | 13.7 | 2078.6 | 152x |  |
| EPSR-test-ex110 | epsr | 0/2 | 9.1 | 1418.8 | 156x | FORTRAN_BUGS.md 1.1 fix: differs from the Fortran by design; every pass identical to the Rust port (tools/check_fixed_case.py) |
| EPSR-test-ex111 | epsr | 0/3 | 13.6 | 2189.5 | 161x | FORTRAN_BUGS.md 1.1 fix: differs from the Fortran by design; every pass identical to the Rust port (tools/check_fixed_case.py) |
| EPSR-test-ex112 | epsr | 2/2 | 5.8 | 880.1 | 151x |  |
| EPSR-test-ex113 | epsr | 2/2 | 7.2 | 1055.7 | 147x |  |
| EPSR-test-ex114 | epsr | 2/2 | 5.8 | 902.5 | 157x |  |
| EPSR-test-ex115 | epsr | 3/3 | 7.9 | 1215.3 | 154x |  |
| EPSR-test-ex116 | epsr | 3/3 | 8.1 | 1334.6 | 165x |  |
| EPSR-test-ex117 | epsr | 0/3 | 7.6 | 1300.9 | 171x | FORTRAN_BUGS.md 1.2 fix: differs from the Fortran by design; every pass identical to the Rust port (tools/check_fixed_case.py) |
| EPSR-test-ex118 | epsr | 0/3 | 10.6 | 1296.0 | 122x | FORTRAN_BUGS.md 1.2 fix: differs from the Fortran by design; every pass identical to the Rust port (tools/check_fixed_case.py) |
| EPSR-test-ex119 | epsr | 3/3 | 8.2 | 1227.0 | 149x |  |
| EPSR-test-ex120 | epsr | 2/2 | 4.3 | 778.6 | 179x |  |
| EPSR-test-ex121 | epsr | 2/2 | 1.4 | 289.3 | 208x |  |
| EPSR-test-ex122 | epsr | 2/2 | 4.3 | 586.1 | 136x |  |
| EPSR-test-ex123 | epsr | 2/2 | 1.6 | 321.4 | 195x |  |
| EPSR-test-ex124 | epsr | 2/2 | 2.5 | 371.3 | 149x |  |

**Total: 115 cases, 339/351 steps identical.** Every differing step is in one of the five
cases that exercise a bug fixed on purpose (docs/FORTRAN_BUGS.md section 1, fixed as in the Rust port); there
the Python is compared with the Rust port instead (tools/check_fixed_case.py, tools/compare_with_rust.py).
