# Package H: coordconstraint_mod.f90, plot_routines.f

## Files

- `src/epsr26/common/coordconstraint_mod.py`: the derived type `type(coordconstraint)` as the
  class `coordconstraint` (all components; allocatables start as `None`, the others 0 / '' /
  False since gfortran leaves them uninitialised). The module has no variables (`reset()` does
  nothing). It is registered with `common.register` all the same.
- `src/epsr26/common/coordconstraints_mod.py`: `coordconstraints` (None, or a **Python list**
  with index 0 unused, so `coordconstraints[ic]` is `coordconstraints(ic)`), `ncoordconstraints`
  (0), `ergold_coord`, `ergnew_coord` (None until `initialise_coords`).
- `src/epsr26/library/coordconstraint_mod.py`: the 13 routines of both modules, in Fortran
  order. Each docstring names its module.
- `src/epsr26/library/plot_routines.py`: the 18 routines of plot_routines.f.

## Calling conventions an integrator needs

- `type(coordconstraint)` dummies are objects modified in place. As SIGNATURES lists,
  `init_constraint` and `find_constraints` return `this`, and `en_coordcon` returns
  `(en, this)`. Callers usually ignore `this`. Example:
  `coordconstraint_mod.find_constraints(coordconstraints_mod.coordconstraints[id], nmol, natom, ...)`.
- `reallocate1d_cc`: see SIGNATURE_EXCEPTIONS.md. read_pcof can call
  `coordconstraint_mod.reallocate1d_cc(coordconstraints_mod.coordconstraints, coordconstraints_mod.ncoordconstraints)`
  and ignore the result, then
  `coordconstraint_mod.init_constraint(coordconstraints_mod.coordconstraints[n], text)`.
  New elements are fresh `coordconstraint()` objects.
- `init_constraint(this, sometext, len_sometext=256)`: the keyword is Fortran `len(sometext)`.
- `ato_read_and_save.f` (normalise/unnormalise_coordinates) changes `ds`/`dsmax` element by
  element: `c.ds[j] = fractcoord(c.ds[j], dmax)`.
- The component arrays use the usual layout: `this.iabondcc[j, ia]`, `this.nbondcc[ia]`.
- Allocatable components are reallocated through `library/reallocation_routines.py`, which
  returns the new array. For example: `this.ds = reallocation_routines.reallocate1d_r(this.ds, n)`.
- Two routines from other files are imported inside the calling function. This avoids import
  cycles and lets the module import before those files exist:
  - `molfit.en_atom_inter_scale_coord`, called by `en_coordcon` and `en_coordcon_total`;
  - `minimagec.minimagec`.
  plot_routines likewise imports `file_and_directory_commands` and `auxiliary_routines`
  inside its functions.

## Bugs and quirks

- FIX (FORTRAN_BUGS.md 1.2): in the angle (`hastb`) pass of `find_constraints`,
  `tb2 = adistcc(j,ia)`, not `(i,ia)`. This is marked `# FIX`.
- The stale-slot family (FORTRAN_BUGS.md section 3) is reproduced:
  - `ishortest`/`dshortest` are not reset between central atoms.
  - Removing an over-bonded atom searches only down to slot 2.
  - Slot 1 of an atom with no bonds is read. Its value is stale from an earlier call (the
    arrays persist in `this` between calls, as the allocatables do), or 0 on the first call.
    In that case `mref[0]` (the unused element, 0) is read where the Fortran reads out of
    bounds. The Rust panics there (`ub()`), so no reference case reaches it.
  - Likewise an unfilled `ishortest(it)` is 0 here. `atype[itype[0]]` is '' (written as 3
    blanks, not found).
  - `iabondcc[0, ja]` (removal pass with `nbondcc(ja)=0`) reads the unused row 0.
  - `save_coord_energies` being called only after rotations is molfit's business, not
    changed here.
- `en_coordcon_total` divides by `tb` (not `abs(tb)`) and scales `distort` by `sqrt(wtms)`,
  unlike `en_coordcon`. Reproduced.
- plot_routines quirks reproduced:
  - **The `plot` command is effectively broken.** `setup_plot_filenames` and `gnuplot`
    address symbols by fixed offsets from `ext` (`iref+1` = bspace, `nsymboluptodata+5` =
    xmin, ...). These offsets ignore the `f` and `b` symbols that `setup_symbols_plot`
    inserts after `ext`, so `bspace` is read from `f`, `xmin` from `bspace`, and so on.
    `xcolumn` is read from `logx` (`'no'`), which aborts the list-directed read in both the
    Fortran and the Python.
  - `check_plot_symbol`'s setone test `icheck.ne.0.or.icheck.ne.1` is always true.
  - Format 127 has one `a` for two items, so format reversion prints
    "setup_plot_filenames> setup_plot_filenames> ", then a blank line, then the name.
  - `find_plot_par` calls `get_nterms` with 4 of its 5 arguments (no `ierr`). The Python
    passes a dummy 0 and discards the returned status.
  - gfortran `-fno-backslash`: `'plot \\'` and `', \\'` write two backslashes.
- gnuplot is launched in the **foreground** with `runtime.system(plot_command)`. The command
  is `system_gnu "<file>.gnu" -` with no `&`, so nothing needs re-routing to
  `shell/processes.spawn` (it is interactive: it reads stdin until gnuplot exits).

## Uninitialised locals (value chosen)

All of these are 0 or '' and are always set before use, except where noted:

- coordconstraint: the ishortest/dshortest contents (see above).
- plot_routines:
  - `ierrout` in plot_new;
  - `iplottype` and `yorn` in process_plots;
  - `iaccept` in setup_plot_filenames and check_plot_symbol;
  - `ifound` in find_ncolumn. Any start value gives the same `search_for_symbol` result.
  - `nfirst`/`nlast` in setup_blocknumbers. They are used after a failed read only on
    malformed input.
  - `ierr` and `ndata` in process_plot_commands;
  - `icheck`, `rcheck`, `nfilet`, `nvary` and `nvaryt` in check_plot_symbol;
  - `nvaryt` in read_plot2d/3d;
  - `ncoeffs` in find_plot_par.
- On end of file, a `read(u,'(a)',iostat=ierr) sometext` leaves `sometext` unchanged, as
  gfortran does. In read_plot_input_file this means the last line is processed once more.

## Structural liberties (each statement-equivalent)

- `_check_same_molecules`: the "Neighbouring atoms in same molecule" block that
  `find_constraints` repeats 4 times, with only the pass number changed.
- The two donor/acceptor passes of the `anglesintra` section run as one body executed twice.
  A diff of coordconstraint_mod.f90:759-861 against 862-965 shows they are identical apart
  from the form of the entry test.
- `_set_symbol`: the five assignments per symbol in `setup_symbols_plot`.
- Private I/O helpers in plot_routines:
  - `_parse(nspc)`: the COMMON `parse` call;
  - `_read_list`: internal list-directed read with iostat;
  - `_read_line`: `read(u,'(a)',iostat=)`;
  - `_open`: open with iostat.
  The compat layer could usefully provide iostat-style wrappers like these.

## Verification

Each comparison was against a gfortran build (`-O2 -fno-backslash`) of the same Fortran,
using the real Python `file_and_directory_commands`, `parse`, `minimagec` and
`reallocation_routines`. All matched byte for byte:

- coordconstraint (Fortran patched only at line 623 for FIX 1.2): 10 three-atom molecules,
  constraint `Si 3 ... angles ... 109.5 100.0 anglesintra ... O 100.0`, with
  `find_constraints` called twice to exercise the stale slots.
  - stdout and `coordcon.txt` are identical.
  - The unpatched Fortran differs in 52 lines of `coordcon.txt`, confirming that the fixed
    angle pass is exercised.
- plot_routines: two sessions.
  - `gnuplot` + `process_plots('l')` + `save_files_plots`: stdout and `plot012.gnu` are
    identical.
  - An interactive `process_plot_symbols` session (help, value change, d/u, npt 2, pt 2,
    unknown command, bad real, l, default, e, file name), followed by
    `save_files_plot2d`/`save_files_plot3d`: stdout, the defaults file and both output
    files are identical.
- The energy routines (`en_coordcon*`, `en_mole_coord`, ...) were not run, because
  `molfit.py` does not exist yet.
