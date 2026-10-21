#!/usr/bin/env python3
"""Generate the keyword-table routines of src/epsr26/library/setup_symbols_generic.py.

Translates every ``subroutine setup_symbols_<program>()`` of
fortran_src/setup_symbols_generic.f (all routines except the hand-written
``setup_symbols_generic``) statement by statement into Python.  Only the statement
forms that occur in these tables are recognised; anything else stops the script, so a
new kind of statement can never be dropped silently.  Comments are kept.

The output file keeps everything above the line ``# --- GENERATED BELOW ...`` (module
docstring, imports, the hand-written setup_symbols_generic) and replaces the rest.

    python tools/gen/gen_setup_symbols_generic.py
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "../../fortran_src/setup_symbols_generic.f")
OUT = os.path.join(HERE, "../../src/epsr26/library/setup_symbols_generic.py")
MARKER = "# --- GENERATED BELOW by tools/gen/gen_setup_symbols_generic.py; do not edit by hand ---"

# declared lengths of the symbols.inc character variables
CHLEN = {"symbol": 11, "symbol_type": 1, "symbol_text": 256, "symbol_default": 256,
         "symbol_value": 256, "prompt": 80, "dataprompt": 80, "dataprompttext": 80,
         "datatext": 80, "extinp": 20}
# include file -> common module (only those whose names are used in the tables)
INCMOD = {"chains.inc": "chains", "lsi.inc": "lsi"}
PYNAME = {"is": "is_"}


def split_comment(s):
    """Split a statement into (code, trailing ! comment) outside character constants."""
    q = None
    for k, c in enumerate(s):
        if q:
            if c == q:
                q = None
        elif c in "'\"":
            q = c
        elif c == "!":
            return s[:k].rstrip(), s[k + 1:].strip()
    return s.rstrip(), None


def read_units():
    """Return [(name, [items])]; item = ('c', text) | ('b',) | ('s', label, stmt, comment, line)."""
    raw = open(SRC, encoding="latin-1").read().split("\n")
    units, cur = [], None
    for ln, line in enumerate(raw, 1):
        if not line.strip():
            if cur is not None:
                cur[1].append(("b",))
            continue
        if line[:1] in "cC*!" or line.strip().startswith("!"):
            if cur is not None:
                cur[1].append(("c", line, ln))
            continue
        if len(line) > 5 and line[5] not in " 0" and line[:5].strip() == "":
            # continuation: append to the last statement
            for k in range(len(cur[1]) - 1, -1, -1):
                if cur[1][k][0] == "s":
                    it = cur[1][k]
                    code, com = split_comment(it[2] + line[6:])
                    cur[1][k] = ("s", it[1], code, com or it[3], it[4])
                    break
            continue
        label = line[:5].strip()
        code, com = split_comment(line[6:].strip() if len(line) > 6 else "")
        m = re.match(r"\s*subroutine\s+(\w+)\s*\(", code)
        if m:
            cur = (m.group(1), [])
            units.append(cur)
            continue
        cur[1].append(("s", label, code.strip(), com, ln))
    return units


def fstring_expr(e):
    """Fortran character expression (constants joined by //) -> Python expression."""
    parts = []
    rest = e.strip()
    while rest:
        m = re.match(r"'((?:[^']|'')*)'\s*", rest) or re.match(r'"((?:[^"]|"")*)"\s*', rest)
        if not m:
            raise SystemExit(f"cannot translate character expression: {e!r}")
        q = rest[0]
        parts.append(repr(m.group(1).replace(q + q, q)))
        rest = rest[m.end():]
        if rest.startswith("//"):
            rest = rest[2:].strip()
        elif rest:
            raise SystemExit(f"cannot translate character expression: {e!r}")
    return " + ".join(parts)


def comment(text):
    t = text
    if re.fullmatch(r"!+\s*", t):
        return None
    for p in ("c$$$", "C$$$", "c", "C", "*", "!"):
        if t.startswith(p):
            t = t[len(p):]
            break
    t = t.rstrip()
    return "#" + (" " + t.strip() if t.strip() else "")


def translate(name, items, formats, mods):
    out = []
    loopvar = None
    ind = "    "
    pending_blank = False

    def emit(s, com=None):
        nonlocal pending_blank
        if pending_blank and out and out[-1] != "":
            out.append("")
        pending_blank = False
        out.append(ind + s + (f"  # {com}" if com else ""))

    for it in items:
        if it[0] == "b":
            pending_blank = True
            continue
        if it[0] == "c":
            c = comment(it[1])
            if c is not None:
                emit(c)
            continue
        _, label, s, com, ln = it
        sl = s.replace(" ", "")
        m = re.match(r"(\d+)?\s*format\s*(\(.*\))$", (label + " " + s).strip())
        if m:
            formats[(name, label)] = m.group(2)
            emit(f"# {label} format{m.group(2)}")
            continue
        if label:
            raise SystemExit(f"{name}:{ln}: unexpected label {label}")
        m = re.match(r"include'(\w+)\.inc'$", sl)
        if m:
            if m.group(1) + ".inc" in INCMOD:
                mods.add(INCMOD[m.group(1) + ".inc"])
            continue
        if sl in ("iref=0", "iref=iref+1"):
            emit(s.replace("=", " = ").replace("+", " + "), com)
            continue
        m = None
        if "=" in s:
            k = s.index("=")
            m = re.match(r"(symbol|symbol_type|symbol_text|symbol_default)\(iref\)=(.*)$",
                         s[:k].replace(" ", "") + "=" + s[k + 1:].strip())
        if m:
            emit(f"symbols.{m.group(1)}[iref] = fstr({fstring_expr(m.group(2))}, {CHLEN[m.group(1)]})", com)
            continue
        m = re.match(r"iflag\(iref\)=(-?\d+)$", sl)
        if m:
            emit(f"symbols.iflag[iref] = {m.group(1)}", com)
            continue
        m = re.match(r"(nsymbol|nsymboldata|nsymboluptodata|irefrho)=iref(\+1)?$", sl)
        if m:
            emit(f"symbols.{m.group(1)} = iref" + (" + 1" if m.group(2) else ""), com)
            continue
        if sl == "nsymbolplus=nsymboldata-nsymboluptodata":
            emit("symbols.nsymbolplus = symbols.nsymboldata - symbols.nsymboluptodata", com)
            continue
        m = re.match(r"nfiles=(\d+)$", sl)
        if m:
            emit(f"symbols.nfiles = {m.group(1)}", com)
            continue
        m = re.match(r"do(\w+)=1,nsymbol$", sl)
        if m:
            loopvar = PYNAME.get(m.group(1), m.group(1))
            emit(f"for {loopvar} in range(1, symbols.nsymbol + 1):", com)
            ind = "        "
            continue
        if loopvar:
            v = re.escape(loopvar.rstrip("_"))
            if re.fullmatch(rf"symbol_value\({v}\)=symbol_default\({v}\)", sl):
                emit(f"symbols.symbol_value[{loopvar}] = symbols.symbol_default[{loopvar}]", com)
                continue
            m = re.fullmatch(rf"write\(6,(\d+)\){v},symbol\({v}\),symbol_value\({v}\)", sl)
            if m:
                emit(f"units.write(6, fmt(_F{m.group(1)}, {loopvar}, Ch(symbols.symbol[{loopvar}], 11), "
                     f"Ch(symbols.symbol_value[{loopvar}], 256)))", com)
                formats.setdefault(("use", m.group(1)), name)
                continue
            if sl in ("enddo",):
                loopvar = None
                ind = "    "
                continue
        if sl == "read(symbol_value(nsymboluptodata),*,iostat=ierr)nfiles":
            emit("ierr = 0", com)
            emit("try:")
            emit("    symbols.nfiles = read_internal(symbols.symbol_value[symbols.nsymboluptodata]).int(symbols.nfiles)")
            emit("except (FortranEOF, FortranIOError) as e:")
            emit("    ierr = iostat_of(e)")
            continue
        if sl == "if(nsymboluptodata.ge.nsymboldata.or.ierr.ne.0)nfiles=1":
            emit("if symbols.nsymboluptodata >= symbols.nsymboldata or ierr != 0:", com)
            emit("    symbols.nfiles = 1")
            continue
        m = re.match(r"(prompt|dataprompt|dataprompttext|extinp)\s*=\s*(['\"].*)$", s)
        if m:
            emit(f"symbols.{m.group(1)} = fstr({fstring_expr(m.group(2))}, {CHLEN[m.group(1)]})", com)
            continue
        m = re.match(r"(ilenprompt|ilendataprompt|ilenextinp)=index\((prompt|dataprompt|extinp),'(.)'\)(-1)?$", sl)
        if m:
            ch = m.group(3)
            if ch == ">":
                expr = f"findex(symbols.{m.group(2)}, '>')"
            else:
                expr = f"findex_padded(symbols.{m.group(2)}, {ch!r}, {CHLEN[m.group(2)]})"
            emit(f"symbols.{m.group(1)} = {expr}" + (" - 1" if m.group(4) else ""), com)
            continue
        # index(...,'>') with a blank in the separator was removed above by the
        # .replace(" ", ""): redo index(x,' ') on the original statement
        m = re.match(r"(ilenprompt|ilendataprompt|ilenextinp)\s*=\s*index\((prompt|dataprompt|extinp),' '\)-1$", s)
        if m:
            emit(f"symbols.{m.group(1)} = findex_padded(symbols.{m.group(2)}, ' ', {CHLEN[m.group(2)]}) - 1", com)
            continue
        m = re.match(r"(ilenprompt|ilendataprompt)=len_trim\((prompt|dataprompt)\)$", sl)
        if m:
            emit(f"symbols.{m.group(1)} = len(symbols.{m.group(2)})", com)
            continue
        m = re.match(r"(ilenprompt|ilendataprompt)=(\d+)$", sl)
        if m:
            emit(f"symbols.{m.group(1)} = {m.group(2)}", com)
            continue
        m = re.match(r"write\((symbol_text|symbol_default)\(iref\),(\d+|\*)\)(\w+)$", sl)
        if m:
            par = m.group(3)
            mod = [mm for mm in mods if par in ("mring", "mchain")]
            if len(mod) != 1:
                raise SystemExit(f"{name}:{ln}: cannot resolve {par}")
            arg = f"{mod[0]}.{par}"
            if m.group(2) == "*":
                emit(f"symbols.{m.group(1)}[iref] = fstr(lw({arg}), 256)", com)
            else:
                formats.setdefault(("use", m.group(2)), name)
                emit(f"symbols.{m.group(1)}[iref] = fstr(fmt(_F{m.group(2)}_{name[len('setup_symbols_'):]}, {arg}), 256)", com)
            continue
        if sl == "return":
            emit("return", com)
            continue
        if sl == "end":
            continue
        raise SystemExit(f"{name}:{ln}: statement not recognised: {s!r}")
    while out and out[-1] == "":
        out.pop()
    return out


def docstring(name, items):
    """First header comment lines of the routine (after the includes)."""
    lines = []
    started = False
    for it in items:
        if it[0] == "s":
            if it[2].replace(" ", "").startswith("include"):
                continue
            break
        if it[0] == "c":
            t = comment(it[1])
            if t is None:
                continue
            t = t[1:].strip()
            if t:
                lines.append(t)
                started = True
            elif started and lines and lines[-1] != "":
                lines.append("")
    while lines and lines[-1] == "":
        lines.pop()
    return lines


def main():
    units = read_units()
    names = [u[0] for u in units]
    assert names[0] == "setup_symbols_generic", names[0]
    formats = {}
    body = []
    allmods = set()
    for name, items in units[1:]:
        mods = set()
        code = translate(name, items, formats, mods)
        allmods |= mods
        doc = docstring(name, items)
        body.append("")
        body.append("")
        body.append("@routine")
        body.append(f"def {name}():")
        body.append(f'    """setup_symbols_generic.f:{name} -- keyword table and defaults.')
        body.append("")
        if doc:
            body.append("    " + doc[0])
        body.append("    Returns nothing.")
        body.append('    """')
        body += code
    # module-level FORMAT constants
    fdefs = []
    f100 = {v for (n, lab), v in formats.items() if n != "use" and lab == "100"}
    assert len(f100) == 1, f100
    fdefs.append(f"_F100 = {f100.pop()!r}")
    for (n, lab), v in formats.items():
        if n == "use" or lab == "100":
            continue
        fdefs.append(f"_F{lab}_{n[len('setup_symbols_'):]} = {v!r}  # {n}")
    head = open(OUT).read().split(MARKER)[0] if os.path.exists(OUT) else ""
    if not head:
        raise SystemExit(f"{OUT} has no head (write the module header and marker first)")
    text = head + MARKER + "\n" + "\n".join(fdefs) + "\n" + "\n".join(body) + "\n"
    open(OUT, "w").write(text)
    print(f"wrote {len(units) - 1} table routines to {os.path.relpath(OUT)}; modules used: {sorted(allmods)}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
