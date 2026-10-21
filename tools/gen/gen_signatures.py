#!/usr/bin/env python3
"""Static index of every routine in the active Fortran sources -> docs/SIGNATURES.md
(Python calling convention) and docs/call_graph.md (static call graph).

Calling convention of the translation (docs/TRANSLATION_GUIDE.md):

* the Python function takes the same arguments as the Fortran routine, in order;
* arrays are passed as NumPy arrays and modified in place;
* the routine *returns* every scalar dummy argument it may assign, in argument order
  (a single value, a tuple, or nothing).  A FUNCTION returns its result first, then
  any assigned scalar dummies.

"May assign" is found statically: the dummy appears on the left of an assignment
(including a substring assignment of a character scalar), in a READ list, as a DO
variable, as ``iostat=``/``exist=``/... in an I/O statement, has ``intent(out|inout)``,
or is passed as an actual argument in a position the callee may assign (iterated to a
fixed point).
"""

import collections
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "fortran_src")


def strip_comment(line):
    out, q = [], None
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


def lower_outside_strings(s):
    out, q = [], None
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


def blank_strings(s):
    """Replace string contents by x so that regexes do not look inside them."""
    out, q = [], None
    for c in s:
        if q:
            if c == q:
                q = None
                out.append(c)
            else:
                out.append("x")
        elif c in "'\"":
            q = c
            out.append(c)
        else:
            out.append(c)
    return "".join(out)


def statements(path):
    """[(line_number, statement)] for fixed (.f) or free (.f90) form."""
    free = path.endswith(".f90")
    with open(path, encoding="latin-1") as fh:
        raw = fh.read().split("\n")
    out = []
    pending_amp = False
    for k, line in enumerate(raw, 1):
        line = line.rstrip("\r").expandtabs(6) if not free else line.rstrip("\r")
        if not free:
            if not line.strip() or line[0] in "cC*!":
                continue
        if line.lstrip().startswith("!") or not line.strip():
            continue
        body = strip_comment(line)
        if not body.strip():
            continue
        if free:
            b = body.strip()
            if pending_amp and out:
                out[-1] = (out[-1][0], out[-1][1] + " " + b.lstrip("&"))
            else:
                out.append((k, b))
            pending_amp = out[-1][1].rstrip().endswith("&")
            if pending_amp:
                out[-1] = (out[-1][0], out[-1][1].rstrip()[:-1])
        else:
            cont = len(body) > 5 and body[:5].strip() == "" and body[5] not in " 0"
            if cont and out:
                out[-1] = (out[-1][0], out[-1][1] + body[6:])
            elif out and out[-1][1].rstrip().endswith("&"):
                out[-1] = (out[-1][0], out[-1][1].rstrip()[:-1] + body.strip())
            else:
                out.append((k, body[6:] if len(body) > 6 and body[:6].strip().isdigit() else body))
    res = []
    for k, s in out:
        s = lower_outside_strings(s.strip())
        # split on ';' outside strings
        for part in split_top(s, ";"):
            if part:
                res.append((k, part))
    return res


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


def match_paren(s, i):
    """Index after the parenthesis group starting at s[i] == '('."""
    depth = 0
    for j in range(i, len(s)):
        if s[j] == "(":
            depth += 1
        elif s[j] == ")":
            depth -= 1
            if depth == 0:
                return j + 1
    return len(s)


UNIT_RE = re.compile(
    r"^(?:(?:recursive|pure|elemental)\s+)*(?:(?:real\s*\*\s*\d+|real\s*\(\s*kind\s*=\s*\d+\s*\)|real|integer\s*\*\s*\d+|integer|logical|double\s+precision|character\s*\*\s*\(?[^ (]*\)?)\s+)?"
    r"(subroutine|function|program)\s+([a-z_][a-z0-9_]*)\s*(\(([^)]*)\))?\s*(result\s*\(\s*([a-z0-9_]+)\s*\))?"
)
DECL_RE = re.compile(
    r"^(real|integer|logical|character|double\s+precision|complex|type\s*\()"
)


class Unit:
    def __init__(self, kind, name, args, file, line, module, result):
        self.kind = kind
        self.name = name
        self.args = args
        self.file = file
        self.line = line
        self.module = module
        self.result = result or name
        self.arrays = set()
        self.intent = {}
        self.assigned = set()
        self.calls = []  # (callee, [actual arg strings], line)
        self.fcalls = []  # possible function references (name, [actual args], line)
        self.body = []


def parse_units(path):
    units = []
    module = None
    cur = None
    depth_stack = []
    for k, s in statements(path):
        sb = blank_strings(s)
        m = re.match(r"^module\s+([a-z_][a-z0-9_]*)\s*$", sb)
        if m and not sb.startswith("module procedure"):
            module = m.group(1)
            continue
        if re.match(r"^end\s*module", sb):
            module = None
            continue
        m = UNIT_RE.match(sb)
        if m and not sb.startswith("end"):
            if cur is not None:
                depth_stack.append(cur)
            args = [a.strip() for a in (m.group(4) or "").split(",") if a.strip()]
            cur = Unit(m.group(1), m.group(2), args, os.path.basename(path), k, module, m.group(6))
            units.append(cur)
            continue
        if cur is None:
            continue
        if re.match(r"^end\s*(subroutine|function|program)?\b", sb) and not re.match(r"^end\s*(if|do|select|where|type|interface|module|forall|associate|block)", sb):
            cur = depth_stack.pop() if depth_stack else None
            continue
        cur.body.append((k, s, sb))
    return units


def analyse_decls(u):
    argset = set(u.args)
    for k, s, sb in u.body:
        if DECL_RE.match(sb) or sb.startswith("dimension"):
            attrs = ""
            rest = sb
            if "::" in sb:
                attrs, rest = sb.split("::", 1)
            else:
                m = re.match(r"^(real\s*\*\s*\d+|integer\s*\*\s*\d+|character\s*\*\s*\(\s*\*\s*\)|character\s*\*\s*\d+|character\s*\([^)]*\)|real\s*\([^)]*\)|integer\s*\([^)]*\)|double\s+precision|real|integer|logical|character|dimension)\s*(.*)$", sb)
                if not m:
                    continue
                rest = m.group(2)
            isdim = "dimension" in attrs
            m = re.search(r"intent\s*\(\s*(in|out|inout|in\s+out)\s*\)", attrs)
            intent = m.group(1).replace(" ", "") if m else None
            for ent in split_top(rest):
                mm = re.match(r"^([a-z_][a-z0-9_]*)\s*(\()?", ent)
                if not mm:
                    continue
                n = mm.group(1)
                if n in argset:
                    if mm.group(2) or isdim:
                        u.arrays.add(n)
                    if intent:
                        u.intent[n] = intent


def actual_args(s, start):
    end = match_paren(s, start)
    return split_top(s[start + 1 : end - 1]), end


IO_KW = ("iostat", "exist", "opened", "number", "size", "nextrec", "named", "name")


def analyse_body(u):
    args = set(u.args)
    for k, s, sb in u.body:
        st = sb
        # strip 'if (...)' prefixes of logical IFs (and labels)
        st = re.sub(r"^\d+\s+", "", st)
        while True:
            m = re.match(r"^(else\s*)?if\s*\(", st)
            if not m:
                break
            e = match_paren(st, st.index("("))
            st = st[e:].strip()
            if st.startswith("then") or not st:
                st = ""
                break
        if not st:
            continue
        # function references f(args) anywhere in the statement (resolved against the
        # known functions after parsing; arguments of functions may be assigned too)
        for m in re.finditer(r"(?<![a-z0-9_%])([a-z_][a-z0-9_]*)\s*\(", st):
            if m.group(1) in ("if", "call", "while", "write", "read", "open", "close", "print", "format"):
                continue
            acts, _ = actual_args(st, m.end() - 1)
            u.fcalls.append((m.group(1), acts, k))
        # calls
        for m in re.finditer(r"\bcall\s+([a-z_][a-z0-9_]*)\s*(\()?", st):
            if m.group(2):
                acts, _ = actual_args(st, m.end() - 1)
            else:
                acts = []
            u.calls.append((m.group(1), acts, k))
        # DO variable
        m = re.match(r"^do\s+(\d+\s*,?\s*)?([a-z_][a-z0-9_]*)\s*=", st)
        if m:
            if m.group(2) in args:
                u.assigned.add(m.group(2))
            continue
        # assignment
        m = re.match(r"^([a-z_][a-z0-9_]*)\s*(\(.*?\))?\s*(%[a-z0-9_%]+)?\s*=(?!=)", st)
        if m and not re.match(r"^(if|do|while|call|print|write|read|open|close|return|goto|go)\b", st):
            if m.group(1) in args:
                u.assigned.add(m.group(1))
            continue
        # I/O statements
        m = re.match(r"^(read|inquire|open|close|write|backspace|rewind|allocate|deallocate)\s*\(", st)
        if m:
            e = match_paren(st, st.index("("))
            ctl = st[st.index("(") + 1 : e - 1]
            for item in split_top(ctl):
                if "=" in item:
                    key, val = item.split("=", 1)
                    if key.strip() in IO_KW or key.strip() == "stat":
                        v = val.strip()
                        if v in args:
                            u.assigned.add(v)
            if m.group(1) == "read":
                lst = st[e:]
                for n in re.findall(r"[a-z_][a-z0-9_]*", lst):
                    if n in args:
                        u.assigned.add(n)
            if m.group(1) == "write" and ctl.split(",")[0].strip() in args:
                u.assigned.add(ctl.split(",")[0].strip())  # internal write to a character dummy
        m = re.match(r"^read\s*\*\s*,(.*)$", st) or re.match(r"^read\s+[^(]*,(.*)$", st)
        if m:
            for n in re.findall(r"[a-z_][a-z0-9_]*", m.group(1)):
                if n in args:
                    u.assigned.add(n)
    for n, it in u.intent.items():
        if it in ("out", "inout"):
            u.assigned.add(n)


INTRINSIC_ASSIGNING = {"getarg": [1], "get_command_argument": [1], "date_and_time": [0, 1, 2, 3],
                       "cpu_time": [0], "system_clock": [0, 1, 2], "random_number": [0], "get_environment_variable": [1],
                       "random_seed": [], "getenv": [1], "move_alloc": [0, 1]}


def main():
    files = [ln.strip() for ln in open(os.path.join(SRC, "ACTIVE_SOURCES.txt")) if ln.strip() and not ln.startswith("#")]
    units = []
    for f in files:
        p = os.path.join(SRC, f)
        if os.path.exists(p):
            units += parse_units(p)
    for u in units:
        analyse_decls(u)
        analyse_body(u)
    byname = collections.defaultdict(list)
    for u in units:
        byname[u.name].append(u)
    # fixed point over calls
    changed = True
    while changed:
        changed = False
        for u in units:
            argset = set(u.args)
            fref = [(c, a, k) for c, a, k in u.fcalls if any(t.kind == "function" for t in byname.get(c, []))
                    and c not in u.arrays and c not in u.args]
            for callee, acts, _ in u.calls + fref:
                targets = byname.get(callee, [])
                for t in targets[:1]:
                    for i, a in enumerate(acts):
                        a = a.strip()
                        if i < len(t.args) and a in argset and t.args[i] in t.assigned and a not in u.assigned:
                            u.assigned.add(a)
                            changed = True
                if callee in INTRINSIC_ASSIGNING:
                    for i in INTRINSIC_ASSIGNING[callee]:
                        if i < len(acts) and acts[i].strip() in argset and acts[i].strip() not in u.assigned:
                            u.assigned.add(acts[i].strip())
                            changed = True
    dups = {n for n, us in byname.items() if len({(x.file) for x in us}) > 1}
    lines = ["# Routine signatures (generated by tools/gen/gen_signatures.py)", "",
             "Python calling convention: same arguments as the Fortran, in order; arrays are",
             "modified in place; the function **returns the scalar dummies it may assign**, in",
             "argument order (a FUNCTION returns its result first).  `[]` marks array dummies.",
             "", "Routines defined in more than one active file are marked **(duplicate)**: check",
             "which file the program links (tools/fortran_link_lists.txt).", ""]
    byfile = collections.defaultdict(list)
    for u in units:
        byfile[u.file].append(u)
    for f in files:
        if f not in byfile:
            continue
        lines.append(f"## {f}")
        lines.append("")
        lines.append("| line | routine | arguments | returns |")
        lines.append("|---:|---|---|---|")
        for u in byfile[f]:
            argtxt = ", ".join(a + ("[]" if a in u.arrays else "") for a in u.args)
            rets = [a for a in u.args if a in u.assigned and a not in u.arrays]
            if u.kind == "function":
                rets = [f"*{u.result}*"] + rets
            name = u.name + (f" (module {u.module})" if u.module else "") + (" **(duplicate)**" if u.name in dups else "")
            kind = {"subroutine": "", "function": "function ", "program": "program "}[u.kind]
            lines.append(f"| {u.line} | {kind}`{name}` | {argtxt} | {', '.join(rets)} |")
        lines.append("")
    with open(os.path.join(ROOT, "docs", "SIGNATURES.md"), "w") as fh:
        fh.write("\n".join(lines))
    # call graph
    cg = ["# Static call graph (generated by tools/gen/gen_signatures.py)", "",
          "`caller (file:line) -> callees` from the Fortran `call` statements (function calls",
          "are not listed).  Compare with a runtime trace: `EPSR_TRACE=calls`.", ""]
    for u in units:
        callees = []
        for c, _, _ in u.calls:
            if c not in callees:
                callees.append(c)
        if callees:
            cg.append(f"- `{u.name}` ({u.file}:{u.line}) -> " + ", ".join(f"`{c}`" for c in callees))
    with open(os.path.join(ROOT, "docs", "call_graph.md"), "w") as fh:
        fh.write("\n".join(cg) + "\n")
    print(f"{len(units)} routines in {len(byfile)} files; {len(dups)} names defined in several files")


if __name__ == "__main__":
    main()
