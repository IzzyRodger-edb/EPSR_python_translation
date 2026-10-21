# Progress log

Read this first if you are picking up the work. Plan: [PLAN.md](../PLAN.md). Rules:
[TRANSLATION_GUIDE.md](TRANSLATION_GUIDE.md). Routine contract: [SIGNATURES.md](SIGNATURES.md) and
[SIGNATURE_EXCEPTIONS.md](SIGNATURE_EXCEPTIONS.md). Where each Fortran routine went:
[ROUTINE_INDEX.md](ROUTINE_INDEX.md).

## Summary (2026-10-21)

Everything `epsrshell` can run is translated (477 of 480 routines of the active sources; the other 3
are interface blocks), about 63,000 lines of Python. Verified byte for byte against the Fortran:

| test | result |
|---|---|
| test problems: 115 cases, 351 steps ([EX_RESULTS.md](EX_RESULTS.md)) | **339/351 steps identical**; the 12 others are all in the 5 cases that exercise the bugs fixed on purpose (ex108, ex110, ex111, ex117, ex118), and there every pass is identical to the Rust port, which has the same fixes |
| validation, water / SiO2 / Cu-BTC on measured data ([VALIDATION_RESULTS.md](VALIDATION_RESULTS.md)) | **27/27 steps identical** (11 + 9 + 7; EPSR steps of 99/100 iterations compared at 2 iterations) |
| shell (`tests/shell`, pexpect) | 6/6: built-ins, Ctrl+C at the prompt, Ctrl+D, no residual processes after `bye` and SIGHUP, lock file, script mode |
| graceful Ctrl+C during `epsr` (`tools/test_graceful_stop.py`) | pass: niter 3 interrupted in iteration 1 gives files identical to a Fortran niter-1 run (the saved `.EPSR.inp` keeps niter 3, as intended), clean exit, no residual process |
| compatibility layer, ran1, keyword tables (`pytest`) | 14/14 |

"Identical" = every file of the run folder, byte for byte, with only the date in `.out` and the
timing columns of `fort.35` masked (`tools/compare_outputs.py`).

## Environment (this machine)

- gfortran 15.2, glibc 2.43 (Ubuntu, WSL2), Python 3.14, NumPy 2.3.5, 12 cores.
- Fortran reference build: `tools/build_fortran_reference.sh $PWD/fortran_build` (gitignored, ~5
  min). It reproduces `EPSR-test-ex1/fortran_reference.tar.gz` byte for byte on this glibc.
  `fortran_build/epsrshell` was also built (for the shell comparison) by compiling
  `epsrshell_standalone.f` against the epsr objects.
- Python env: `python3 -m venv .venv && . .venv/bin/activate && pip install -e '.[dev]'`.
- Byte identity depends on glibc's libm (`sinf`, `expf`, `sinhf`, ...): regenerate the
  Fortran references on each machine and compare there; do not compare across machines.

## How to re-run the verification

```bash
. .venv/bin/activate
pytest -q                                                      # unit + shell tests
python tools/run_ex_py.py test_problems/200DMSO400Water/EPSR-test-ex1 --work work/ex1
python tools/run_all_ex_py.py --jobs 10                        # all 115 cases (~6 h on 12 cores)
python tools/run_ex_py.py validation/water --max-niter 2 --work work/val/water   # ~1.5 h
python tools/check_fixed_case.py test_problems/200DMSO400Water/EPSR-test-ex110 work/fixed/ex110  # vs Rust
python tools/test_graceful_stop.py work/graceful
```

`run_ex_py.py` is resume mode by default (each step starts from the Fortran folder before that step);
`--mode chain` runs the Python on its own from the raw inputs. The Rust binary used for the fixed
cases is `~/EPSRpy/rust/target/release/epsr-rs`.

## Run time: Python against Fortran

Across all 351 test-problem steps (median of Python time / Fortran time per program, machine shared
with other runs; short steps are dominated by Python start-up, ~0.4 s):

| program | steps | ratio (median) | range |
|---|---:|---:|---|
| epsr | 135 | 162x | 80x-211x |
| mapgrentropy | 12 | 144x | 46x-239x |
| voids | 3 | 139x | 125x-187x |
| sharm | 13 | 119x | 9x-190x |
| fluctuations, triangles, partials, local_density, rings, clusters | 31 | 68x-92x | |
| lsi, coord, chains, sdf, sdfcube, mapgr, torangles | 56 | 20x-56x | |
| setup and weights tools (makemole, mixato, randomise, fmole, nwts, readpdb, changeato, dockato, ...) | 98 | 6x-52x | mostly start-up |

Absolute examples: one EPSR-test iteration (600 molecules, 5 MC cycles) 2.3 s Fortran, ~4-5 min
Python; water validation step of 2 iterations (1000 molecules) 4.4 s vs 13 min; `voids` 7 s vs 15-38 min;
`sharm` on DMSO 3 s vs 11 min.

Stages of one `epsr` pass (EPSR-test-ex1; Python from `EPSR_TRACE=summary`, Fortran from `perf`
self-time shares of the 2.3 s run, so approximate; `tools/stage_times.py` compares the time-stamped
stages of the two logs directly):

| stage | Fortran | Python | ratio |
|---|---:|---:|---:|
| Monte Carlo moves and energies (`molfit` x5, `fmolec1` x2) | ~1.1 s | 132 s | ~120x |
| EP refinement (`ptgnsq`, almost all of it `fit_expansion`) | ~0.3 s | 56 s | ~190x |
| histograms, g(r) transforms, pressure (`atogen`, `stog`, `press`) | ~0.2 s | 59 s | ~300x |
| S(Q) from histograms (`gensq`) | ~0.02 s | 3.4 s | ~200x |
| reading/writing every file, setup (`inita`, `initpt`, `w_int`, `moloutc`, ...) | ~0.5 s | ~3 s | ~6x |

Where the Python time goes: every float32 operation is a NumPy scalar operation (~50-100 ns), and
every single-precision transcendental (`expf`, `sinf`, `cosf`) is a ctypes call into glibc (~1 µs),
which the rules require for byte identity (NumPy's float32 ufuncs do not use glibc).
Ideas for a faster path that keeps the bits (not done, outside the plan's first pass): a small C
extension for the libm calls; numba on the energy and histogram loops checked against this
translation; vectorising loops whose iterations are independent and contain no reductions
(`fit_expansion`'s chi-square rows, `voids` pixel pairs), which gives identical bits but departs from
the line-by-line rule.

## Infrastructure

| file | what |
|---|---|
| `src/epsr26/fortran/` | fmath (float32, libm via ctypes, powi, nint, atan2d...), fformat (`fmt`, `lw`, `Ch`), listio (list-directed/formatted input), units (I/O units, gfortran reopen rules), text (character helpers), arrays, runtime (getarg, system), stop (`iostat_of`) |
| `src/epsr26/common/` | generated include-file modules (`tools/gen/gen_common.py`) + hand-written Fortran-module data; `SaveBlock`; `reset()`, `snapshot()`/`restore()` |
| `src/epsr26/trace.py` | `@routine`; `EPSR_TRACE=calls|summary`, `EPSR_TRACE_PATTERN=a,b,c`; shell `trace on/off/summary` |
| `src/epsr26/programs/__init__.py` | program table and `run_program(name, args)` |
| `src/epsr26/run.py` | `python -m epsr26.run <program> <fortran args>` |
| `src/epsr26/shell/` | `epsrshell.py` (translation), `builtins.py`, `signals.py`, `processes.py`, `history.py` |
| `tools/run_ex_py.py`, `run_all_ex_py.py` | Fortran-vs-Python harness (commands.txt and steps.txt; resume/chain; `--max-niter`) |
| `tools/compare_with_rust.py`, `check_fixed_case.py` | Python vs Rust for the deliberately fixed cases |
| `tools/test_graceful_stop.py`, `stage_times.py` | Ctrl+C check; stage timing from the logs |
| `tools/gen/` | `gen_common.py`, `gen_signatures.py` (SIGNATURES.md, call_graph.md), `gen_routine_index.py`, `gen_setup_symbols_generic.py`, `merge_parts.py` |

## How the translation was done

Two waves of parallel translators, each given `TRANSLATION_GUIDE.md` and the signature contract;
per-package notes (uninitialised variables, quirks reproduced, decisions) are in
`docs/translation_notes/<pkg>.md`.

- Wave 1, the `epsr` link list: A small libraries, B auxiliary_routines + read_and_save_files,
  C process_symbols, D ato_read_and_save + setup_epsr_rev1 + cell lists + pair store,
  E/F ato_manipulation_rev4 (two halves, merged), G molfit + three-body + container forces,
  H coordination constraints + plot_routines, I/J/K epsr_standalone_rev3 (three parts, merged).
- Wave 2, all other programs: L setup_symbols_generic, M setup and weights programs, N1 dockato/addato,
  N2 small tools + setup/plot/plot3djmol, O analysis I, P chains/rings/lsi/mapgr, Q sharm/sdf/sdfcube.
  `runjmol` was translated by the integrator.

Integration fixes worth knowing: `en_mole_intra` returns all its sums (the signature generator
missed arguments modified through function references; fixed and re-checked); reallocating routines
return the new array; `units.open` follows gfortran when a connected unit is reopened; FORMAT accepts
`0x`, `0(...)`, `*(...)`; shell commands run as `sh -c` like `system(3)`.

## Shell: what was wrong with the Fortran shell on Linux (verified 2026-10-21)

Built `fortran_build/epsrshell` and drove it with scripts.

With **no** `system_commands.txt` (`system_termination` blank):

- `pwd` prints "No such folder!": `directory_name` has no trailing separator, so the existence check
  on `<dir>.` fails;
- `cd sub` **hangs forever**: `set_directory_name` loops on `index(temp_name, '')`, which is always 1
  for an empty separator. The only way out is Ctrl+C, which kills the shell but not what it started.

With the **Windows** `system_commands.txt` shipped in EPSR-test (`\`, `dir/b`, `cd/d`):

- every folder name gets a `\` appended, so `pwd` and every `cd` report "No such folder!";
- `ls` runs `dir/b <dir>`: "sh: dir/b: not found".

In both cases Ctrl+D crashes the shell ("Fortran runtime error: End of file" in
`FCreadline_windows.f90`), and programs started with `call system(...)` (and viewers started with
`... &`) are not tracked, so they outlive the shell.

The Python shell keeps the Fortran's model (the shell process stays in its start folder, where
`runflag.txt` and `FChistory.txt` live; `cd` changes `directory_name`, which every program receives
as its input folder) but implements the commands directly, runs programs in-process, and tracks every
child process group. Visible differences: the "Binaries folder" line names the Python package; `help`
has extra lines for `exit`/`quit`, `!cmd`, `jobs`, `trace`; an unknown command runs through `/bin/sh`
in the working folder.

## Open items

1. **Speed** (above). The literal translation is the reference; any faster path must be checked
   against it.
2. **Fortran defects found and reproduced, not fixed** (fixing changes outputs; candidates for
   docs/FORTRAN_BUGS.md):
   - `plot` aborts with default settings: `setup_plot_filenames`/`gnuplot` read symbols by fixed
     offsets that skip the later-added `f` and `b` symbols (translation_notes/H.md, N2.md);
   - `setup_symbols_generic` parses 3020 characters of a 256-character line and reads the next
     variable (the `.symbols` file name) on blank lines (L.md);
   - `add_image_atom` and `remove_molecule_types` loop forever on some inputs (C.md);
   - `makeato` segfaults in `moloutc` (two arrays never allocated; reproduced as exit 139, M.md);
   - `epsrwts` renormalises only the last abundance (M.md); `lsi` overruns `gofz` with more than 3
     distributions (P.md).
3. **Not exercised by any test**: `save_files`/`force_save_files` interactive saves; three-body and
   coordination branches of `normalise_coordinates`; `readmole` (`readgpr`/`readjmol`) has no Fortran
   source in `fortran_src/`, so the shell reports "no program 'readmole'", as the Fortran shell had no
   binary for it either.
4. **Chain mode** (PLAN.md section 6, final acceptance: the Python alone from raw inputs over a whole
   validation case) has not been run end to end; resume mode covered every step. With 99/100-iteration
   EPSR steps a full Python chain takes days per case at the current speed.
5. CI (GitHub Actions) is not set up.
