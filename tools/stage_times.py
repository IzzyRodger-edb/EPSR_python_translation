#!/usr/bin/env python3
"""Time each stage of an EPSR run in the Fortran and the Python from their stdout.

usage: tools/stage_times.py FORTRAN_LOG PYTHON_LOG [--markdown]

EPSR prints time-stamped progress lines ("Start molfit iterations 13:51:58.152 ...",
"fmolec1> Finish energy calculation ..."); the translation prints the same lines at the
same points.  The interval from each time-stamped line to the next is a stage; the
stages of the two logs are paired in order and their durations compared.
"""

import re
import sys

TS = re.compile(r"^(.*?)\s+(\d\d):(\d\d):(\d\d)\.(\d\d\d)\s+\d\d-[A-Z]{3}-\d{4}\s*$")


def stamps(path):
    out = []
    for line in open(path, encoding="latin-1"):
        m = TS.match(line.rstrip("\n"))
        if m:
            h, mi, s, ms = (int(x) for x in m.groups()[1:])
            out.append((m.group(1).strip(), h * 3600 + mi * 60 + s + ms / 1000.0))
    return out


def stages(st):
    res = []
    for (a, ta), (b, tb) in zip(st, st[1:]):
        d = tb - ta
        if d < 0:
            d += 86400
        res.append((f"{a}  ->  {b}", d))
    return res


def main():
    f, p = stamps(sys.argv[1]), stamps(sys.argv[2])
    if [a for a, _ in f] != [a for a, _ in p]:
        print("warning: the logs do not have the same time-stamped lines", file=sys.stderr)
    sf, sp = stages(f), stages(p)
    md = "--markdown" in sys.argv
    if md:
        print("| stage (from line -> to line) | Fortran [s] | Python [s] | ratio |")
        print("|---|---:|---:|---:|")
    agg = {}
    for (name, tf), (_, tp) in zip(sf, sp):
        key = name
        agg.setdefault(key, [0.0, 0.0, 0])
        agg[key][0] += tf
        agg[key][1] += tp
        agg[key][2] += 1
    tot_f = sum(v[0] for v in agg.values())
    tot_p = sum(v[1] for v in agg.values())
    for name, (tf, tp, n) in sorted(agg.items(), key=lambda kv: -kv[1][1]):
        r = f"{tp / tf:.0f}x" if tf > 0.0005 else "-"
        cnt = f" (x{n})" if n > 1 else ""
        if md:
            print(f"| {name}{cnt} | {tf:.3f} | {tp:.1f} | {r} |")
        else:
            print(f"{tf:9.3f} {tp:10.1f} {r:>6}  {name}{cnt}")
    if md:
        print(f"| **total between first and last stamp** | {tot_f:.2f} | {tot_p:.1f} | {tot_p / tot_f:.0f}x |")
    else:
        print(f"{tot_f:9.3f} {tot_p:10.1f}         total")


if __name__ == "__main__":
    main()
