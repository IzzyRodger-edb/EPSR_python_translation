# EPSR_python_translation

A line-by-line Python translation of EPSR26 (Empirical Potential Structure Refinement,
A. K. Soper), written to be read and modified rather than to be fast. It keeps the
`epsrshell` interface and is validated byte for byte against the Fortran on the test
problems and on measured neutron data for water, amorphous silica and Cu-BTC.

Status and results: [docs/PROGRESS.md](docs/PROGRESS.md). Plan: [PLAN.md](PLAN.md).

## Summary of the work (2026-10-21)

Everything `epsrshell` can run is translated: the `epsr` engine, the 36 setup and analysis
programs, `setup`, `plot`, `plot3djmol`, `runjmol` and the shell itself. That is 477 of the 480
routines in the active Fortran sources (the other 3 are interface blocks), about 63,000 lines of
Python.

| check | result |
|---|---|
| test problems, 115 cases, 351 steps ([EX_RESULTS](docs/EX_RESULTS.md)) | **339 identical**; the other 12 are in the 5 cases that exercise the documented bug fixes (ex108, ex110, ex111, ex117, ex118), where every pass is identical to the Rust port, which has the same fixes |
| measured data ([VALIDATION_RESULTS](docs/VALIDATION_RESULTS.md)) | **water 11/11, SiO₂ 9/9, Cu-BTC 7/7 steps identical**: setup, scratch start, equilibration, empirical-potential refinement, restart (99/100-iteration steps compared at 2 iterations) |
| shell | built-ins, Ctrl+C, Ctrl+D, lock file, script mode, and no processes left after `bye` or a closed terminal |
| Ctrl+C during `epsr` | the current iteration is finished and saved; the files equal a Fortran run that stops at that iteration |

"Identical" means every file in the run folder, byte for byte; only the date in `.out` and the
timing columns of `fort.35` are masked.

**The two shell fixes.** Run on Linux, the Fortran shell fails in two ways:

- **Without `system_commands.txt`:** `cd` hangs forever (it searches for an empty path separator).
- **With the Windows `system_commands.txt` that ships with EPSR26:** every folder name gets a `\`,
  so `cd` and `pwd` fail, and `ls` runs `dir/b`.

In both cases Ctrl+D crashes the shell, and the programs it started outlive it. The Python shell
implements these commands directly and tracks every child process. It keeps the Fortran's banner,
prompt, help, commands, script mode and program questions.

**Time cost.** One `epsr` iteration of the 600-molecule test case takes 4–5 minutes in Python
against 2.3 s in Fortran. `epsr` is a median of ~160× slower across 135 steps. By stage:

| stage | slowdown |
|---|---:|
| Monte Carlo moves and energies | ~120× |
| empirical-potential fit | ~190× |
| histograms, g(r) and pressure | ~300× |
| file reading and writing | ~6× |

Analysis programs are 20–150× slower; setup tools are 6–50× slower, mostly Python start-up. The
cost comes from float32 scalar arithmetic and per-call access to glibc's libm. Both are needed for
byte identity.

**Open items:**

- The plan's final chain-mode run (the Python alone over a whole validation case) has not been
  done; it would take days per case at this speed.
- CI is not set up.
- Some Fortran bugs were found and reproduced rather than fixed: `plot` aborts with default
  settings, two editing routines in `changeato` can loop forever, and `makeato` segfaults. See
  docs/PROGRESS.md.

## Is a pure Python port worth it? (verdict)

**Yes, as the reference and working copy of EPSR. Not as the program people run for production
refinements.**

What it is good for:

- **Understanding and changing the code.** Every routine is one Python function with the Fortran
  name and docstring (`file.f:routine`). COMMON blocks are visible modules, and
  `trace on` / `EPSR_TRACE` shows which routines are called, when, and how long they take. A
  change can be checked against the Fortran with one command
  (`tools/run_ex_py.py`). This was the main goal, and it is met.
- **Trust.** Byte identity on 115 feature cases and three measured-data cases means any
  difference after a modification is caused by the modification.
- **Usability.** The shell now works properly on Linux and cleans up after itself.

Where it falls short:

- **Speed.** At ~160× slower, a typical refinement of a few hundred iterations goes from minutes
  to a day or more. A 100-iteration water run takes about 11 hours.

  This cannot be fixed within the pure-Python rules. The rules that make it a faithful reference
  (float32 scalars, Fortran loop order, glibc libm calls) are what make it slow, and vectorising
  NumPy changes the last bits.

Recommendation:

- Keep this translation as the readable, modifiable and validated specification.
- For production, use the Fortran or the Rust port, or a fast path built from this code: numba or
  a small C extension on the hot loops (the energies, `fit_expansion`, `atogen`/`stog`). Accept it
  only after the same byte-identity harness passes against this translation.
- Prototype physics changes here, on short runs, then port them to the fast path.

## Install

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e '.[dev]'
```

Requires Linux with glibc (single-precision maths calls glibc's `libm`, as the
gfortran build does), Python >= 3.10 and NumPy >= 2.

## Use

```bash
cd my_epsr_folder
epsrshell                 # or: python -m epsr26
```

The shell looks and behaves like the Fortran `epsrshell`: same banner, prompt, `help`,
command names, script mode (`ss`/`ps`/`es`, `runflag.txt`) and the questions each program
asks. Differences:

- `ls`, `pwd`, `cd` (also `cd ..`, `cd -`, `cd ~`), `md`, `del`, `edit` work on Linux
  without `system_commands.txt` (in the Fortran `cd` could hang and `ls` ran `dir/b`);
  `!cmd` or `system cmd` runs a shell command in the working folder.
- `bye`, `e`, `exit`, `quit` or Ctrl+D leave cleanly; every process the shell started
  (gnuplot, Jmol, `system` commands) is stopped, so nothing is left running.
- Ctrl+C at the prompt discards the line. Ctrl+C during `epsr` finishes the current
  iteration and saves every file consistently; a second Ctrl+C within 3 s (or Ctrl+C in
  any other program) aborts it. A running script is paused.
- `jobs` lists running child processes; `trace on [pattern]`, `trace off`,
  `trace summary` show which routines are called, when, and how long they take.

A single program can also be run like the Fortran binary:

```bash
python -m epsr26.run epsr ./ epsr water.EPSR.inp
EPSR_TRACE=summary python -m epsr26.run epsr ./ epsr water.EPSR.inp   # with timing per routine
```

It is slow: one EPSR iteration of the 600-molecule test case takes ~5 min against ~2 s
for the Fortran (see docs/PROGRESS.md for timings of every stage).

## Layout

| path | contents |
|---|---|
| `src/epsr26/fortran/` | Fortran compatibility layer (float32 arithmetic, glibc libm, FORMAT, list-directed I/O, units) |
| `src/epsr26/common/` | COMMON blocks and Fortran module data (generated from the `.inc` files) |
| `src/epsr26/library/` | shared routine files, one Python module per Fortran file |
| `src/epsr26/programs/` | one module per Fortran program (`epsr_standalone_rev3.py` is the engine) |
| `src/epsr26/shell/` | the shell, built-ins, signals, child processes, history |
| `docs/` | translation guide, routine index, signatures, call graph, progress, notes per package |
| `tools/` | Fortran reference build, comparison and test harness, generators |
| `fortran_src/` | the original Fortran, unchanged |

Finding the Python for a Fortran routine: [docs/ROUTINE_INDEX.md](docs/ROUTINE_INDEX.md).
Every Python function's docstring starts with the Fortran `file:routine` it translates.

## Testing

```bash
pytest -q                                              # compatibility layer, library, shell
tools/build_fortran_reference.sh $PWD/fortran_build   # once, ~5 min (gfortran)
python tools/run_ex_py.py test_problems/200DMSO400Water/EPSR-test-ex1   # Fortran vs Python, every file
python tools/run_ex_py.py validation/water --max-niter 2                # real data, shortened EPSR steps
```
