# Package L: keyword tables of the auxiliary programs (setup_symbols_generic.f)

## Files

| Fortran | Python |
|---|---|
| setup_symbols_generic.f, `setup_symbols_generic` | `library/setup_symbols_generic.py`, hand-written (top of the file) |
| setup_symbols_generic.f, the 23 `setup_symbols_<program>` tables | same file, below the `# --- GENERATED BELOW` marker, written by `tools/gen/gen_setup_symbols_generic.py` |
| (test) | `tests/library/data/setup_symbols_driver.f`, `tests/library/data/setup_symbols/*.symbols`, `tests/library/data/setup_symbols_reference.txt`, `tests/library/test_setup_symbols_generic.py` |

All 24 routines are translated in full: setup_symbols_generic, changeato, partials, coord,
clusters, rings, chains, lsi, sharm, sdf, sdfcube, mapgrentropy, plot2d, splot2d, plot3d,
plot3djmol, mapgr, triangles, torangles, voids, fluctuations, nwts, xwts, bragginten. Their
signatures match SIGNATURES.md (`setup_symbols_generic(auxprogname, ssg_ierr)` returns
`ssg_ierr`; the others take nothing and return nothing). No signature exceptions.

## Generator

`python tools/gen/gen_setup_symbols_generic.py` keeps the file up to the marker and
rewrites everything after it. It translates the Fortran one statement at a time and keeps
the comments. It only knows the statement forms that occur in these tables (symbol/type/
text/default/iflag assignments, `nsymbol*=iref`, the defaults loop with its `write(6,100)`,
the internal read of `nfiles`, prompt/extension assignments, `index`/`len_trim`, and the
`write(symbol_text(iref),101) mring` / `write(symbol_default(iref),*) mchain` lines in rings,
chains and lsi). Any other statement stops the script, so nothing can be dropped silently.
FORMATs go to module constants: `_F100` (the same in every table) and `_F101_rings`,
`_F101_chains`, `_F101_lsi`. `mring` comes from `common.chains`; `mchain` comes from
`common.chains` (chains) or `common.lsi` (lsi), the include each routine uses. Both are 25
and 200.

## Verification

`tests/library/test_setup_symbols_generic.py` (2 tests, about 13 s, most of it the gfortran
build) passes:

1. The Python runs the same call sequence as the gfortran driver and its unit-6 output
   matches `setup_symbols_reference.txt` byte for byte (1271 lines). The sequence is:
   - all 23 tables, then sharm again (its stale read of `symbol_value(nsymbol+1)`);
   - setup_symbols_generic on 16 cases. One is a good file with ndist, a truncated
     symbol name, and a text cut off by 13+ blanks. Then: no ndist, no final newline,
     missing prompt, prompt without `>`, missing extension, ndist without dataprompt or
     dataprompttext, a blank line inside a symbol block, a bad key, end of file inside a
     block, missing default/type/text, a missing file, and leading blanks in the program
     name.

   After each call the test dumps every symbol (name, type, value, default, text, iflag),
   the counters, nfiles, prompt, dataprompt, extinp with their lengths, irefrho, and
   commands.nwords/ncf(1)/ncl(1). It also checks everything the routines write themselves.
2. When gfortran is present, the driver is rebuilt and the stored reference is checked
   against a fresh run. I also built it with the reference flags
   (`-O2 -fno-tree-vectorize -fopenmp`) and the output was identical.

## Quirks reproduced

- **parse reads past the end of sometext2.** setup_symbols_generic calls
  `parse(sometext2,...,nstsize,...)` with nstsize = 3020, but sometext2 is character*256.
  On a blank line (or at end of file), gfortran's parse scans on into the next local,
  `symbols_fname`, and finds its first word. Effects:
  - a blank line or end of file inside a symbol's default/type/text block (the
    `nwords.ge.0` test is always true) gives
    `Unrecognised value for symbol X : <name>.symbols in line N`;
  - after a successful call, `commands.nwords` = 1 and `ncf(1)` = 257 (the last read is
    EOF).

  The Python passes `_mem(sometext2, symbols_fname)` (the two padded variables joined) to
  parse and takes cmd/val from it. On lines with text the result is unchanged. The driver
  checks both effects.
- The "ndist not found" message continues a character constant onto the next source line,
  so it prints `therefore nodistributions specified`.
- `write(6,102) "Complete ...", "-- symbols() ..."`: format reversion gives two
  `/`-preceded records (fmt does this).
- `write(6,102) " command passed in: "//auxprogname` writes the full 256-character
  padded auxprogname.
- On an error return unit 14 stays open, as in the Fortran.
- sharm, sdf and sdfcube read `nfiles` from `symbol_value(nsymbol+1)` before setting the
  values, then set `nfiles=1`. mapgrentropy sets `nfiles=0`. These are reproduced
  as written.

## Not reproduced / open

- If the first `.symbols` line holding a keyword has a single word, and no earlier parse
  set `ncf(2)`, the Fortran takes `sometext2(0:n)`, which reads one byte before the
  variable. The Python's `substr` with i=0 does not emulate that byte. This needs a
  malformed file on the first line and never happens with valid files.
- The first `parse(auxprogname,...)` would also read past auxprogname if it were entirely
  blank (the memory beyond it belongs to the caller). This is not emulated.
- The callers (setup_generic.f, ...) belong to other packages. The tables only fill
  `symbols.*` (and `irefrho` in changeato).
