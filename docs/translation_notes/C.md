# Package C: process_symbols.f

Output: `src/epsr26/library/process_symbols.py`. All 22 routines are translated, in Fortran order,
with every branch, including the interactive ones the Rust skipped: `d`/`u` navigation, `help`, `ls`/`cd`,
`save`/`forcesave`, the "not saved" prompt, datafile/wtsfile `search`, the mapgr/sharm/sdf/nwts/xwts/bragginten
branches, and all changeato editing routines.

## Calls to other packages: lazy module proxies

`process_symbols.py` and `read_and_save_files.py` call each other, and most callees belong to other
packages. So the other library modules are reached through a private `_LazyModule` proxy defined at the top of the file.
The proxy looks up the attribute on every call, so tracing and `common.reset()` still work. Call
sites read like ordinary module calls (`auxiliary_routines.scapro(a, b)`). The proxy also accepts several
candidate module names:

| proxy | candidates (first that defines the name wins) |
|---|---|
| `ato_manipulation_rev4` | `library.ato_manipulation_rev4`, then `library.ato_manipulation_rev4_part2`. `change_boxsize` (Fortran line 3591) sits in part 2 until the files are merged |
| `fcreadline_windows` | `library.fcreadline_windows`, `library.FCreadline_windows`, `library.jmk_readline_module`. The Python name for FCreadline_windows.f90 was not fixed, so package A should use one of these |
| `lib_molecom_mod` | `library.molecom_mod` (`testsize_rotgrps`, `testsize_labrot`) |
| `molecom_mod` | `common.molecom_mod` (data; assignments are forwarded, e.g. `molecom_mod.ngrrot = ih1`) |
| `pfroutines_mod` | `library.pfroutines_mod` (`getpfprobability`, `getrratio`) |
| others | `auxiliary_routines`, `ato_read_and_save`, `file_and_directory_commands`, `get_molecules`, `plot_routines`, `read_and_save_files` |

Compat-layer wish: a shared lazy-import helper of this kind (`fortran.lazy_module(...)`) would let other
files with cycles use it as well.

Other private helpers:

| helper | what it does |
|---|---|
| `_read_int` / `_read_real` / `_read_word` | `read(text,*,iostat=ierr) x`. Returns `(x, ierr)` and keeps the old value on error |
| `_read_line5` | `read(5,'(a)',iostat=ierr) s` |
| `_parse_sometext` | the many `call parse(...)` calls on the COMMON `nwords`/`ncf`/`ncl` |
| `_word(k)` | `sometext(ncf(k):ncl(k))` |
| `_pr()` | `prompt(1:ilenprompt)` |
| `_elem_or0` | see uninitialised variables below |

A READ without `iostat` is translated as a plain `read_internal(...)`/`units.list_read(5)`. A bad value
or end of file then raises, as the Fortran runtime error does.

## Signatures

- The calls follow SIGNATURES.md, e.g. `nout, command, icount, ifilecount, ierrout = process_commands(...)`.
- `concatenate(2,nwords,symbol_value(iref),...)` passes the COMMON `nwords` as `nlast`. The returned
  `nlast` is assigned back to `commands.nwords`, which is the same aliasing as in the Fortran.
- `force_save_files(nin,nout,.true.)` gets a third actual argument that the routine does not declare.
  The Python passes `(nin, nout)`. This is logged in SIGNATURE_EXCEPTIONS.md.
- `scatt_len` is called as `btest, ifind = scatt_len(ctemp, massno, ifind)`.
- `change_boxsize`, `get_molecules`, `moloutc`, `ftype` and `ctype` return nothing.

## Uninitialised variables (choice made, marked `NB:` in the code)

| where | variable | Fortran | Python |
|---|---|---|---|
| `change_bond`, `add_bond`, `remove_bond`, `change_label` | READ targets (`ncha1`..., `it`, `iamol`) | old value kept on a null value | 0 |
| `align_molecules` | `labaxis` for an axis letter other than x,y,z,a,b,c | stack garbage | zeros |
| `set_molecule_scale_factors` | `rfract` when its read fails | uninitialised but still used | 0.0 |
| `change_tether`, `shift_molecules`, `align_molecules`, `set_molecule_scale_factors` | `moltyperead` / `mtyperead` | uninitialised | only used when ierr = 0 |
| `write_symbols_changeato` | `natomsphere`, `radminsphere`, `radmaxsphere`, `pseudofacmol` element when the array is unallocated or too short | reads garbage | `_elem_or0` gives 0.0 (as the Rust `get_or0` in changeato.rs) |
| `change_symbol_changeato` | local `ntf`/`ntl` read via `tokens(ntf(1):ntl(1))` when the value is blank | garbage | zeros, so the token is empty |
| `check_symbol`, abundances with nwords = 0 | `ncl(0)` | the word before the array (`ncf(201)`) | element 0 = 0, same empty result |
| `check_symbol`, ATOM CHECK | `type(it)(1:ncsearch)` with ncsearch > 3 | reads past the character*3 element | blank padding (as the Rust) |

## Quirks reproduced (not fixed)

- `change_tether` with the values given on the keyword line (`ntokens > 4`) does not parse the new
  `sometext`. It uses the stale `ncf`/`ncl` from the previous parse, so the words are offset.
- `rotate_molecules` given a molecule *type* (a non-integer first word) never rotates, because the test
  is `ierr.eq.0.and.nwords.gt.2`.
- `add_image_atom(.true.)` with the atom numbers on the keyword line (`ntokens > 3`) loops forever. The
  same tokens are re-read each pass and the loop only ends on a 0.
- `add_image_atom(.false.)` tests `mtype(im)` instead of `mtype(im1)` in the inner loop.
- `remove_molecule_types` loops forever at end of file on unit 5 (`readagain` stays `.true.`).
- `remove_molecule` moves rotational groups but not their dihedral data. `cleanup_types` does not
  shift `pseudofacmol`.
- check_symbol writes its "not an integer" message to unit 6 and its "not a real number" message to
  nout, as the Fortran does.
- The `process_commands` "Command not recognised" FORMAT has 2 data descriptors and 3 items, so
  format reversion produces 4 records: blank, the line with the command, blank, then the line with
  the text. This was checked against gfortran 15. The Rust (`symbols.rs`) writes it as a single line,
  so the Rust is wrong here.
- `write(6,300) mcoeff` with i3 and mcoeff = 1000 prints `***`.

## Intrinsics

The cubic-box branch of `change_symbol_changeato` calls `cosd`/`sind` without declaring them
EXTERNAL, so gfortran uses its intrinsics. The Python uses `cosd_intrinsic`/`sind_intrinsic` (as
`ato.rs` moloutc). `x**(1.0/3.0)` becomes `powf(x, r4(1.0)/r4(3.0))`.

## Testing done

- `python -c "import epsr26.library.process_symbols"` works, and `ruff check --select F` is clean.
- Smoke tests with stub siblings: the interactive loop (help, `<symbol> <value>`, type-check
  rejection, nfiles change with resequencing, `data 2`, d/u navigation, the not-saved prompt), and
  the changeato routines (cleanup_types, remove_molecule, rotate_mol, write_symbols_changeato,
  atlabel/rho/box changes). The list-directed internal writes give gfortran's text
  (`'  0.100000001'`, `' cubic'`).
- Not yet compared end to end with the Fortran: changeato needs packages A, B and D (EPSR-test ex26/ex27).
