#!/usr/bin/env python3
"""Generate src/epsr26/common/<name>.py from the Fortran include files.

usage: tools/gen/gen_common.py            (regenerates every include used by ACTIVE_SOURCES)

Each include file becomes one module:

* ``parameter`` constants -> module constants (evaluated as gfortran does: integer
  arithmetic truncates, real constants are folded in single precision);
* variables in a COMMON block -> module globals created by ``reset()`` (zero / blank /
  ``.false.``, as in the Fortran's zero-initialised .bss);
* variables declared but in no COMMON block -> these are *locals* of every routine that
  includes the file.  They are listed in ``LOCALS`` (name -> description) and the
  translator must create them in the routine.  If they have a DATA statement they are
  module constants instead (DATA implies SAVE and the values are never changed).

Translated code accesses common variables through the module, never with ``from``:
``molecom.xa[i]``, ``clustercom.nmol = 5``.  ``common.reset()`` re-runs every reset().
"""

import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "fortran_src")
OUT = os.path.join(ROOT, "src", "epsr26", "common")

# epsr_dimension.inc must come first: everything else uses its parameters.
ORDER_FIRST = ["epsr_dimension.inc"]


def modname(inc):
    return re.sub(r"[^a-z0-9_]", "_", inc[: -len(".inc")].lower())


def strip_comment(line):
    out = []
    q = None
    for c in line:
        if q:
            out.append(c)
            if c == q:
                q = None
        elif c in "'\"":
            q = c
            out.append(c)
        elif c == "!":
            break
        else:
            out.append(c)
    return "".join(out)


def statements(path):
    """Fixed-form (with some free-form '&') statements, lower-cased outside strings."""
    stmts = []
    with open(path, encoding="latin-1") as fh:
        raw = fh.read().split("\n")
    for line in raw:
        line = line.rstrip("\r")
        if not line.strip():
            continue
        if line[0] in "cC*!":
            continue
        if line.lstrip().startswith("!"):
            continue
        body = strip_comment(line)
        if not body.strip():
            continue
        cont = len(body) > 5 and body[:5].strip() == "" and body[5] not in " 0"
        if cont and stmts:
            stmts[-1] += body[6:]
        elif stmts and stmts[-1].rstrip().endswith("&"):
            stmts[-1] = stmts[-1].rstrip()[:-1] + body.lstrip().lstrip("&")
        else:
            stmts.append(body)
    return [lower_outside_strings(s.strip()) for s in stmts]


def lower_outside_strings(s):
    out = []
    q = None
    for c in s:
        if q:
            out.append(c)
            if c == q:
                q = None
        elif c in "'\"":
            q = c
            out.append(c)
        else:
            out.append(c.lower())
    return "".join(out)


def split_top(s, sep=","):
    parts, depth, cur, q = [], 0, [], None
    for c in s:
        if q:
            cur.append(c)
            if c == q:
                q = None
            continue
        if c in "'\"":
            q = c
        elif c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        if c == sep and depth == 0:
            parts.append("".join(cur).strip())
            cur = []
        else:
            cur.append(c)
    if "".join(cur).strip():
        parts.append("".join(cur).strip())
    return parts


TYPE_RE = re.compile(
    r"^(real\s*\*\s*8|real\s*\*\s*4|real|double\s+precision|integer\s*\*\s*8|integer\s*\*\s*4|integer|"
    r"logical|character\s*\*\s*\(\s*[^)]*\)|character\s*\*\s*\d+|character\s*\(\s*len\s*=\s*[^)]*\)|character)"
    r"(\s*,\s*[a-z]+)*\s*(::)?\s*(.*)$"
)


def norm_type(t):
    t = re.sub(r"\s+", "", t)
    if t in ("real", "real*4"):
        return ("real", None)
    if t in ("real*8", "doubleprecision"):
        return ("dble", None)
    if t in ("integer", "integer*4"):
        return ("int", None)
    if t == "integer*8":
        return ("int8", None)
    if t == "logical":
        return ("logical", None)
    m = re.match(r"character\*\(?(?:len=)?([^)]*)\)?$", t) or re.match(r"character\(len=([^)]*)\)$", t)
    if m:
        return ("char", m.group(1))
    if t == "character":
        return ("char", "1")
    raise ValueError(t)


_NUM = re.compile(r"\d+\.\d*(?:[ed][-+]?\d+)?|\.\d+(?:[ed][-+]?\d+)?|\d+[ed][-+]?\d+|\d+")


def names_in(expr):
    """Names referenced by an expression (numeric literals such as 1.2e-7 removed)."""
    return re.findall(r"[a-z_][a-z0-9_]*", _NUM.sub(" ", re.sub(r"'[^']*'", " ", expr)))


def implicit_type(name):
    return ("int", None) if name[0] in "ijklmn" else ("real", None)


class Inc:
    def __init__(self, fname):
        self.fname = fname
        self.mod = modname(fname)
        self.params = []  # (name, expr, type)
        self.decl = {}  # name -> [type, dims or None]
        self.order = []
        self.commons = []  # (block, [names])
        self.data = {}  # name -> list of python value strings
        self.notes = []


def data_groups(s):
    """Split ``a/v1,v2/, b/v3/`` into [(name, 'v1,v2'), ...], honouring quotes."""
    out = []
    i = 0
    n = len(s)
    while i < n:
        while i < n and s[i] in " ,":
            i += 1
        st = i
        while i < n and s[i] != "/":
            i += 1
        name = s[st:i].strip()
        i += 1
        st = i
        q = None
        while i < n:
            c = s[i]
            if q:
                if c == q:
                    q = None
            elif c in "'\"":
                q = c
            elif c == "/":
                break
            i += 1
        if name:
            out.append((name, s[st:i]))
        i += 1
    return out


def parse_entity(ent):
    m = re.match(r"^([a-z_][a-z0-9_]*)\s*(?:\((.*)\))?\s*(?:\*\s*(\d+))?$", ent.strip())
    if not m:
        raise ValueError(f"cannot parse entity {ent!r}")
    name, dims, clen = m.group(1), m.group(2), m.group(3)
    return name, (split_top(dims) if dims else None), clen


def parse_inc(fname):
    inc = Inc(fname)
    for st in statements(os.path.join(SRC, fname)):
        if st.startswith("parameter"):
            body = st[st.index("(") + 1 : st.rindex(")")]
            for p in split_top(body):
                n, e = p.split("=", 1)
                inc.params.append((n.strip(), e.strip()))
            continue
        if st.startswith("common"):
            # common /blk/ a,b /blk2/ c ...
            rest = st[len("common") :].strip()
            for m in re.finditer(r"/\s*([a-z0-9_]*)\s*/([^/]*)", rest):
                names = []
                for ent in split_top(m.group(2)):
                    n, dims, _ = parse_entity(ent)
                    names.append(n)
                    if dims:
                        inc.decl.setdefault(n, [None, None])[1] = dims
                inc.commons.append((m.group(1), names))
            continue
        if st.startswith("data"):
            for name, vals in data_groups(st[4:]):
                inc.data[name] = split_top(vals)
            continue
        if st.startswith("dimension"):
            for ent in split_top(st[len("dimension") :]):
                n, dims, _ = parse_entity(ent)
                inc.decl.setdefault(n, [None, None])[1] = dims
            continue
        if st.startswith(("save", "external", "implicit")):
            inc.notes.append(st)
            continue
        m = TYPE_RE.match(st)
        if m:
            base = norm_type(m.group(1))
            ents = split_top(m.group(4))
            for ent in ents:
                n, dims, clen = parse_entity(ent)
                t = ("char", clen) if (clen and base[0] == "char") else base
                d = inc.decl.setdefault(n, [None, None])
                d[0] = t
                if dims:
                    d[1] = dims
                if n not in inc.order:
                    inc.order.append(n)
            continue
        raise ValueError(f"{fname}: unrecognised statement {st!r}")
    return inc


# --------------------------------------------------------------------------- evaluation
class Env:
    def __init__(self):
        self.vals = {}  # name -> python value
        self.owner = {}  # name -> module


def eval_param(expr, typ, env):
    """Evaluate a parameter expression with Fortran semantics."""
    toks = re.findall(r"\d+\.\d*(?:[ed][-+]?\d+)?|\.\d+(?:[ed][-+]?\d+)?|\d+[ed][-+]?\d+|\d+|[a-z_][a-z0-9_]*|'[^']*'|\*\*|[-+*/()]", expr)
    py = []
    for t in toks:
        if t[0].isdigit() or (t[0] == "." and len(t) > 1):
            if re.match(r"^\d+$", t):
                py.append(f"_I({t})")
            else:
                py.append(f"_R('{t.replace('d', 'e')}')")
        elif t[0] == "'":
            py.append(repr(t[1:-1]))
        elif re.match(r"[a-z_]", t):
            env.vals[t]  # KeyError if the parameter is not defined yet
            py.append(f"_V({t!r})")
        else:
            py.append(t)
    code = " ".join(py)

    class _I(int):
        def _bin(self, o, op):
            if isinstance(o, _I):
                a, b = int(self), int(o)
                if op == "/":
                    q = abs(a) // abs(b)
                    return _I(q if (a >= 0) == (b >= 0) else -q)
                return _I({"+": a + b, "-": a - b, "*": a * b}[op])
            return _R(np.float32(int(self)))._bin(o, op)

        def __add__(s, o): return s._bin(o, "+")
        def __sub__(s, o): return s._bin(o, "-")
        def __mul__(s, o): return s._bin(o, "*")
        def __truediv__(s, o): return s._bin(o, "/")
        def __radd__(s, o): return _wrap(o)._bin(s, "+")
        def __rsub__(s, o): return _wrap(o)._bin(s, "-")
        def __rmul__(s, o): return _wrap(o)._bin(s, "*")
        def __rtruediv__(s, o): return _wrap(o)._bin(s, "/")
        def __neg__(s): return _I(-int(s))
        def __pos__(s): return s

    class _R:
        def __init__(self, v):
            self.v = np.float32(v)

        def _bin(self, o, op):
            b = np.float32(int(o)) if isinstance(o, _I) else o.v
            a = self.v
            return _R({"+": a + b, "-": a - b, "*": a * b, "/": a / b}[op])

        def __add__(s, o): return s._bin(o, "+")
        def __sub__(s, o): return s._bin(o, "-")
        def __mul__(s, o): return s._bin(o, "*")
        def __truediv__(s, o): return s._bin(o, "/")
        def __radd__(s, o): return _wrap(o)._bin(s, "+")
        def __rsub__(s, o): return _wrap(o)._bin(s, "-")
        def __rmul__(s, o): return _wrap(o)._bin(s, "*")
        def __rtruediv__(s, o): return _wrap(o)._bin(s, "/")
        def __neg__(s): return _R(-s.v)
        def __pos__(s): return s

    def _wrap(o):
        return o if isinstance(o, (_I, _R)) else (_I(o) if isinstance(o, int) else _R(o))

    def _V(n):
        v = env.vals[n]
        if isinstance(v, (int, np.integer)):
            return _I(int(v))
        if isinstance(v, str):
            return v
        return _R(v)

    r = eval(code, {"_I": _I, "_R": _R, "_V": _V})
    if typ[0] == "int":
        return int(r) if isinstance(r, _I) else int(r.v)
    if typ[0] == "real":
        return np.float32(int(r)) if isinstance(r, _I) else r.v
    if typ[0] == "dble":
        return float(int(r)) if isinstance(r, _I) else float(r.v)
    if typ[0] == "char":
        return r.rstrip(" ")
    return r


def dims_code(dims, env, inc, imports):
    out = []
    for d in dims:
        if ":" in d:
            lo, hi = d.split(":")
            if lo.strip() != "0":
                raise ValueError(f"lower bound {lo} not supported ({inc.fname})")
            out.append(f"(0, {expr_code(hi, env, inc, imports)})")
        else:
            out.append(expr_code(d, env, inc, imports))
    return ", ".join(out)


def expr_code(e, env, inc, imports):
    e = e.strip()
    names = names_in(e)
    for n in names:
        if n in env.owner and env.owner[n] != inc.mod:
            imports.setdefault(env.owner[n], set()).add(n)
    if re.search(r"[a-z_]", e) is None and re.fullmatch(r"[\d\s+*()-]+", e):
        return e
    # integer expression of parameters: keep it readable, evaluate with idiv semantics
    v = eval_param(e, ("int", None), env)
    if re.fullmatch(r"[a-z_][a-z0-9_]*", e):
        return e
    return f"{v}  # {e}" if False else str(v)


PYTYPE = {"real": "real", "dble": "dble", "int": "int", "int8": "int8", "logical": "logical", "char": "char"}


def scalar_init(t):
    return {"real": "f32(0.0)", "dble": "0.0", "int": "0", "int8": "0", "logical": "False", "char": "''"}[t[0]]


def array_init(t, dimc):
    fn = {"real": "real_array", "dble": "dble_array", "int": "int_array", "int8": "int8_array",
          "logical": "logical_array", "char": "char_array"}[t[0]]
    return f"{fn}({dimc})"


def data_value(tok, t):
    tok = tok.strip()
    if tok.startswith("'"):
        return repr(tok[1:-1].rstrip(" "))
    if t[0] == "real":
        return f"f32({tok!r})" if False else f"f32('{tok}')"
    if t[0] == "dble":
        return f"float('{tok.replace('d', 'e')}')"
    if t[0] == "logical":
        return "True" if "t" in tok else "False"
    return str(int(tok))


def generate(inc, env):
    lines = []
    imports = {}
    pconst = []
    for n, e in inc.params:
        t = inc.decl.get(n, [None])[0] or implicit_type(n)
        v = eval_param(e, t, env)
        env.vals[n] = v
        env.owner[n] = inc.mod
        for nn in names_in(e):
            if nn in env.owner and env.owner[nn] != inc.mod:
                imports.setdefault(env.owner[nn], set()).add(nn)
        if t[0] == "real":
            pconst.append(f"{n} = f32({float(v)!r})  # {e}")
        elif t[0] == "char":
            pconst.append(f"{n} = {v!r}")
        elif t[0] == "dble":
            pconst.append(f"{n} = {v!r}  # {e}")
        else:
            pconst.append(f"{n} = {v}" + (f"  # {e}" if not re.fullmatch(r"\d+", e) else ""))
    in_common = {}
    for blk, names in inc.commons:
        for n in names:
            in_common[n] = blk
    params = {n for n, _ in inc.params}
    body = []
    for blk, names in inc.commons:
        body.append(f"    # common /{blk}/")
        for n in names:
            t, dims = inc.decl.get(n, [None, None])
            t = t or implicit_type(n)
            if dims:
                body.append(f"    {n} = {array_init(t, dims_code(dims, env, inc, imports))}")
            else:
                body.append(f"    {n} = {scalar_init(t)}")
            if n in inc.data:
                body.append(f"    # NB: DATA for common variable {n}")
                body.extend(data_assign(n, inc.data[n], t, dims))
    locals_ = []
    consts = []
    for n in inc.order + [k for k in inc.decl if k not in inc.order]:
        if n in in_common or n in params:
            continue
        t, dims = inc.decl[n]
        t = t or implicit_type(n)
        desc = f"{t[0]}" + (f"*{t[1]}" if t[0] == "char" else "")
        if dims:
            desc += f"({', '.join(dims)})"
        if n in inc.data:
            consts.append((n, t, dims))
        else:
            locals_.append((n, desc))
    hdr = [
        f'"""Generated from fortran_src/{inc.fname} by tools/gen/gen_common.py -- do not edit.',
        "",
        "Common variables are module globals created by ``reset()``; access them as",
        f"``{inc.mod}.<name>``.",
    ]
    if locals_:
        hdr += ["", "Declared here but in no COMMON block (locals of each including routine, see",
                "``LOCALS``): " + ", ".join(n for n, _ in locals_) + "."]
    hdr.append('"""')
    lines += hdr
    lines.append("")
    lines.append("import numpy as np  # noqa: F401")
    lines.append("")
    lines.append("from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401")
    lines.append("from ..fortran.fmath import f32  # noqa: F401")
    for m in sorted(imports):
        lines.append(f"from .{m} import {', '.join(sorted(imports[m]))}  # noqa: F401")
    lines.append("")
    if pconst:
        lines += pconst
        lines.append("")
    for n, t, dims in consts:
        if dims:
            lines.append(f"{n} = {array_init(t, dims_code(dims, env, inc, imports))}")
            lines += [ln[4:] for ln in data_assign(n, inc.data[n], t, dims)]
        else:
            lines.append(f"{n} = {data_value(inc.data[n][0], t)}")
        lines.append("")
    lines.append(f"LOCALS = {{{', '.join(repr(n) + ': ' + repr(d) for n, d in locals_)}}}")
    lines.append("")
    names = [n for _, ns in inc.commons for n in ns]
    lines.append("")
    lines.append("def reset():")
    lines.append(f'    """(Re)create the COMMON storage of {inc.fname}."""')
    if names:
        for k in range(0, len(names), 8):
            lines.append("    global " + ", ".join(names[k : k + 8]))
        lines += body
    else:
        lines.append("    pass")
    lines.append("")
    lines.append("")
    lines.append("reset()")
    lines.append("")
    # imports may have grown while generating consts: regenerate import block
    text = "\n".join(lines)
    imp = "".join(f"from .{m} import {', '.join(sorted(imports[m]))}  # noqa: F401\n" for m in sorted(imports))
    text = re.sub(r"(from \.\.fortran\.fmath import f32  # noqa: F401\n)(?:from \.[a-z0-9_]+ import [^\n]*\n)*", lambda m: m.group(1) + imp, text, count=1)
    return text


def data_assign(n, toks, t, dims):
    vals = []
    for tok in toks:
        m = re.match(r"^(\d+)\s*\*\s*(.*)$", tok)
        if m and not tok.startswith("'"):
            vals += [m.group(2)] * int(m.group(1))
        else:
            vals.append(tok)
    if not dims:
        return [f"    {n} = {data_value(vals[0], t)}"]
    if len(dims) != 1:
        raise ValueError(f"DATA for multi-dimensional {n}")
    return [f"    {n}[1:] = [{', '.join(data_value(v, t) for v in vals)}]"]


def used_includes():
    files = [ln.strip() for ln in open(os.path.join(SRC, "ACTIVE_SOURCES.txt")) if ln.strip() and not ln.startswith("#")]
    incs = set()
    for f in files:
        p = os.path.join(SRC, f)
        if not os.path.exists(p):
            continue
        for m in re.finditer(r"^\s+include\s+'([^']+)'", open(p, encoding="latin-1").read(), re.M | re.I):
            incs.add(m.group(1))
    return sorted(incs, key=lambda x: (x not in ORDER_FIRST, x))


def main():
    env = Env()
    incs = [parse_inc(f) for f in used_includes()]
    # resolve parameter dependencies: process in an order where every referenced
    # parameter is already known
    done = []
    pending = incs[:]
    while pending:
        progress = False
        for inc in pending[:]:
            needed = set()
            for _, e in inc.params:
                needed |= set(names_in(e))
            for n, (t, dims) in inc.decl.items():
                for d in dims or []:
                    needed |= set(names_in(d))
            own = {n for n, _ in inc.params}
            if all(n in env.vals or n in own or n == "len" for n in needed):
                text = generate(inc, env)
                with open(os.path.join(OUT, inc.mod + ".py"), "w") as fh:
                    fh.write(text)
                done.append(inc)
                pending.remove(inc)
                progress = True
        if not progress:
            sys.exit("unresolved parameters in: " + ", ".join(i.fname for i in pending))
    # common/__init__.py
    mods = [i.mod for i in done]
    init = ['"""COMMON blocks, include files and Fortran module data (global state, as in the Fortran).',
            "",
            "``reset()`` recreates all of it; ``programs.run_program`` calls it before every program,",
            "because every Fortran program started in a fresh process.  The include-file modules",
            "are generated by tools/gen/gen_common.py; the Fortran-module data modules (``*_mod``)",
            "are hand-written and register themselves with ``register()``.",
            '"""',
            "",
            "import importlib",
            "",
            "INCLUDE_MODULES = [",
            *[f"    {m!r}," for m in mods],
            "]",
            "",
            "_extra = []",
            "",
            "",
            "def register(module):",
            '    """Add a hand-written data module (anything with a reset() method) to reset()."""',
            "    if module not in _extra:",
            "        _extra.append(module)",
            "",
            "",
            "class SaveBlock:",
            '    """SAVEd local variables of a routine (``save``, DATA, or initialised in the',
            "    declaration).  They keep their values between calls within one program run and",
            "    are restored to the initial values by ``reset()``, i.e. at the start of each program.",
            "",
            "        _save = SaveBlock('molfit', icall=0, first=True)",
            "        ...",
            "        _save.icall += 1",
            '    """',
            "",
            "    def __init__(self, routine, **initial):",
            "        object.__setattr__(self, '_routine', routine)",
            "        object.__setattr__(self, '_initial', initial)",
            "        self.reset()",
            "        register(self)",
            "",
            "    def reset(self):",
            "        import copy",
            "",
            "        for k, v in self._initial.items():",
            "            object.__setattr__(self, k, copy.deepcopy(v))",
            "",
            "",
            "def reset():",
            "    for m in INCLUDE_MODULES:",
            "        importlib.import_module(f'{__name__}.{m}').reset()",
            "    for m in _extra:",
            "        m.reset()",
            "",
            "",
            "def _state_holders():",
            "    import types",
            "",
            "    for m in INCLUDE_MODULES:",
            "        yield importlib.import_module(f'{__name__}.{m}')",
            "    for m in _extra:",
            "        yield m",
            "",
            "",
            "def snapshot():",
            '    """All global state, by reference (reset() rebinds names, so the old objects survive).',
            "",
            "    The shell runs programs in its own process; it takes a snapshot before a program",
            "    and restores it afterwards, so its own COMMON variables are untouched, as they were",
            "    when every Fortran program ran in a separate process.",
            '    """',
            "    import types",
            "",
            "    snap = []",
            "    for m in _state_holders():",
            "        d = {k: v for k, v in vars(m).items() if not k.startswith('_') and not callable(v)",
            "             and not isinstance(v, types.ModuleType)}",
            "        snap.append((m, d))",
            "    return snap",
            "",
            "",
            "def restore(snap):",
            "    for m, d in snap:",
            "        for k, v in d.items():",
            "            if isinstance(m, SaveBlock):",
            "                object.__setattr__(m, k, v)",
            "            else:",
            "                setattr(m, k, v)",
            ""]
    with open(os.path.join(OUT, "__init__.py"), "w") as fh:
        fh.write("\n".join(init))
    print(f"generated {len(done)} modules in {OUT}")


if __name__ == "__main__":
    main()
