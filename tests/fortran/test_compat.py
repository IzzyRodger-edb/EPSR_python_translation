"""The compatibility layer reproduces gfortran output exactly.

``data/compat_expected.txt`` was produced by ``data/gen_compat.f90`` (gfortran -O2
-fno-tree-vectorize, glibc 2.43).  This test regenerates the same text in Python.
"""

import os

import numpy as np
import pytest

from epsr26.fortran import (
    Ch,
    atan2f,
    cosf,
    expf,
    f32,
    fmt,
    idiv,
    ifix,
    imod,
    logf,
    lw,
    nint,
    powf,
    powi,
    sinf,
    sqrtf,
)
from epsr26.fortran.listio import ListRead, parse_int, parse_real, read_formatted

HERE = os.path.dirname(__file__)

V = [1.0, -0.123456789, 0.0, 1.82, -2.44, 0.1, 0.0999999, 1.0e8, 1.0e9, -0.000123,
     -40.6027222, 99999999.0, -0.0005, 12.3456, 3.3333333e-5, 1.0e-30, 7.5e35, -1.5,
     0.5, 2.5, -2.5, 0.49999997, 123456.789, -9.87654e-12]
D = [1.0, 0.1, -3.3333333333333e-5, 1.0e20, -2.5, 123456789.123, 1.0e-300, 0.0]


def bits(x):
    return int(np.array([x], dtype=np.float32).view(np.int32)[0])


def generate():
    out = []
    v = [np.float32(x) for x in V]
    for i, x in enumerate(v, 1):
        out.append(fmt("(a,i0)", "#v ", i))
        out.append(fmt("(e14.7,1x,e13.6,1x,e12.5,1x,e9.3,1x,e10.3e3)", x, x, x, x, x))
        out.append(fmt("(f10.5,1x,f8.3,1x,f6.3,1x,f14.8,1x,f10.2,1x,f7.0)", x, x, x, x, x, x))
        out.append(lw(x))
        out.append(lw(i, x, "ab", "cd", x))
    for i, d in enumerate(D, 1):
        out.append(fmt("(a,i0)", "#d ", i))
        out.append(lw(d))
        out.append(fmt("(e25.16,1x,f30.10)", d, d))
    out.append("#i")
    iv = [0, 7, -7, 12345, -2147483647, 42]
    out.append(fmt("(i0,1x,i5,1x,i2.2,1x,i3.3,1x,i4,1x,i11)", iv[0], iv[1], iv[2], iv[3], iv[3], iv[4]))
    out.append(lw(*iv))
    vals = [7]
    for i in range(1, 8):
        vals += [i, f32(i) + f32(0.5)]
    out.append(fmt("(1x,i4,(5(1x,i4,1x,e9.3)))", *vals))
    out.append(fmt("(a11,1x,a,a5,a2)", Ch("abc", 11), "x", "yz", "pqr"))
    out.append(fmt("(3x,a,1h-,a,09x,a,1h-,a,05x)", "S", "C", "C", "S"))
    out.append(fmt("(l1,1x,l3)", True, False))
    out.append(lw(True, 3, "text"))
    out.append("#math")
    for i, x in enumerate(v, 1):
        y = abs(x) + f32(0.25)
        out.append(fmt("(10(1x,i0))", bits(expf(x / f32(100.0))), bits(logf(y)), bits(sinf(x)),
                       bits(cosf(x)), bits(atan2f(x, y)), bits(sqrtf(y)), bits(powf(y, f32(1.37))),
                       bits(powi(y, i - 12)), nint(x), ifix(x)))
    out.append("#imod")
    out.append(fmt("(8(1x,i0))", imod(7, 3), imod(-7, 3), imod(7, -3), imod(-7, -3),
                   idiv(7, 3), idiv(-7, 3), idiv(7, -3), idiv(-7, -3)))
    return "\n".join(out) + "\n"


def test_matches_gfortran():
    with open(os.path.join(HERE, "data", "compat_expected.txt")) as fh:
        expected = fh.read()
    got = generate()
    exp_lines = expected.split("\n")
    got_lines = got.split("\n")
    for k, (e, g) in enumerate(zip(exp_lines, got_lines), 1):
        assert g == e, f"line {k}:\n expected {e!r}\n got      {g!r}"
    assert len(exp_lines) == len(got_lines)


def test_list_read_multi_record():
    lines = iter(["3 1 0.5 2", " 0.25 3 0.125 extra", "next"])
    r = ListRead(source=lambda: next(lines))
    n = r.int()
    v = [(r.int(), r.real()) for _ in range(n)]
    assert v == [(1, 0.5), (2, 0.25), (3, 0.125)]
    assert next(lines) == "next"


def test_parse_numbers():
    assert parse_real("1.0D+01") == 10.0
    assert parse_real("1.0+01") == 10.0
    assert parse_real("-.244E+01") == np.float32(-2.44)
    assert parse_real("5.") == 5.0
    with pytest.raises(Exception):
        parse_int("1.0")
    r = ListRead(text="1, ,3*2.5 /")
    assert r.int() == 1
    assert r.real(np.float32(9)) == 9  # null value keeps the old value
    assert [r.real() for _ in range(3)] == [2.5, 2.5, 2.5]
    assert r.real(np.float32(7)) == 7  # after '/'


def test_formatted_read():
    assert read_formatted("(i7,3g14.7)", "      1" + " 0.1000000E+01" + "      2.500000" + "         -3.25") == [
        1, np.float32(1.0), np.float32(2.5), np.float32(-3.25)]
    assert read_formatted("(1x,a3,1x,a3)", " OW  HW") == ["OW", "HW"]
