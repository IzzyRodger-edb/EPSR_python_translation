#!/usr/bin/env python3
"""Merge translation part files (written in parallel) into one module.

usage: tools/gen/merge_parts.py OUT PART1 [PART2 ...]

The module docstrings are combined, every top-level import is collected into one block
(a name bound to two different things in different parts is an error), and the rest of
each part follows in order, verbatim.  The part files are deleted afterwards.
"""

import ast
import os
import sys


def binding(node):
    out = []
    if isinstance(node, ast.Import):
        for a in node.names:
            out.append((a.asname or a.name.split(".")[0], ("import", a.name, a.asname)))
    else:
        for a in node.names:
            out.append((a.asname or a.name, ("from", node.module, node.level, a.name)))
    return out


def main():
    out, parts = sys.argv[1], sys.argv[2:]
    docs, imports, bodies = [], {}, []
    order = []
    for p in parts:
        src = open(p).read()
        tree = ast.parse(src)
        lines = src.split("\n")
        body = tree.body
        k = 0
        if body and isinstance(body[0], ast.Expr) and isinstance(getattr(body[0], "value", None), ast.Constant):
            docs.append(body[0].value.value.strip())
            k = 1
        rest_start = None
        for node in body[k:]:
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                if getattr(node, "module", None) == "__future__":
                    continue
                for name, b in binding(node):
                    if name == "*":
                        key = ("*", b[1], b[2])
                        if key not in imports:
                            imports[key] = b
                            order.append(key)
                        continue
                    if name in imports and imports[name] != b:
                        sys.exit(f"conflict for {name!r}: {imports[name]} vs {b} ({p})")
                    if name not in imports:
                        imports[name] = b
                        order.append(name)
            else:
                rest_start = node.lineno - 1
                # include decorators / comments directly above
                if hasattr(node, "decorator_list") and node.decorator_list:
                    rest_start = node.decorator_list[0].lineno - 1
                while rest_start > 0 and lines[rest_start - 1].startswith("#"):
                    rest_start -= 1
                break
        if rest_start is not None:
            bodies.append(f"# {'=' * 77}\n# from {os.path.basename(p)}\n# {'=' * 77}\n\n" + "\n".join(lines[rest_start:]).rstrip() + "\n")
    # emit imports
    plain, froms = [], {}
    for key in order:
        b = imports[key]
        if b[0] == "import":
            plain.append(f"import {b[1]}" + (f" as {b[2]}" if b[2] else ""))
        else:
            mod = "." * b[2] + (b[1] or "")
            froms.setdefault(mod, [])
            nm = b[3] if key == b[3] or key[0] == "*" else f"{b[3]} as {key}"
            if key[0] == "*" if isinstance(key, tuple) else False:
                nm = "*"
            if nm not in froms[mod]:
                froms[mod].append(nm)
    imp = sorted(set(plain))
    for mod in sorted(froms, key=lambda m: (not m.startswith(".."), m)):
        names = froms[mod]
        if "*" in names:
            imp.append(f"from {mod} import *  # noqa: F403")
            names = [n for n in names if n != "*"]
        for n in names:
            imp.append(f"from {mod} import {n}")
    doc = '"""' + "\n\n".join(docs) + '\n"""'
    text = doc + "\n\n" + "\n".join(imp) + "\n\n\n" + "\n\n".join(bodies)
    with open(out, "w") as fh:
        fh.write(text)
    for p in parts:
        if os.path.abspath(p) != os.path.abspath(out):
            os.remove(p)
    print(f"merged {len(parts)} parts into {out}")


if __name__ == "__main__":
    main()
