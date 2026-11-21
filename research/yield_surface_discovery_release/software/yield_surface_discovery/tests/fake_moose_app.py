#!/usr/bin/env python3
"""FAKE MOOSE application for unit tests (NOT MOOSE).  Reads the filled input deck (-i file), extracts the traction directions, J2 material
parameters and time levels, and writes <file_base>.csv with the columns MOOSE's CSV output would have, using the exact J2/linear-hardening
response of a homogeneous element under stress control.  Exit code 3 if FAKE_MOOSE_FAIL is set (to test failure handling)."""
import os, re, sys
if os.environ.get("FAKE_MOOSE_FAIL"):
    print("fake solver failure", file=sys.stderr); sys.exit(3)
text = open(sys.argv[sys.argv.index("-i") + 1]).read()
num = lambda name: float(re.search(rf"{name}\s*=\s*([-+0-9.eE]+)", text).group(1))
func = lambda name: float(re.search(rf"\[{name}\][^\[]*expression = '([-+0-9.eE]+)\*t'", text).group(1))
sxx, syy, sxy = func("txx"), func("tyy"), func("txy"); sy, H = num("yield_stress"), num("hardening_constant")
levels = [float(v) for v in re.search(r"time_sequence = '([^']*)'", text).group(1).split()]
base = re.search(r"file_base = (\S+)", text).group(1)
vm = (sxx * sxx + syy * syy - sxx * syy + 3 * sxy * sxy) ** 0.5
with open(base + ".csv", "w") as fh:
    fh.write("time,eqps,sxx,syy,sxy,szz,sxz,syz\n0,0,0,0,0,0,0,0\n")
    for t in levels:
        fh.write(f"{t!r},{max(0.0, (t * vm - sy) / H)!r},{t * sxx!r},{t * syy!r},{t * sxy!r},0,0,0\n")
print(f"fake MOOSE: {len(levels)} steps written to {base}.csv")
