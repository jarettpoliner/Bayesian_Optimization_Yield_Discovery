"""ANALYTIC STAND-IN for a simulation driver (NOT MOOSE): writes result.json for one loading direction (plane-stress von Mises, 250 MPa).
A real driver would generate a MOOSE input for the direction, run it, extract the yield radius and write the same JSON."""
import json, sys, time
from pathlib import Path
d = json.loads(sys.argv[1]); wd = Path(sys.argv[2]); t0 = time.process_time()
vm = (d[0] ** 2 + d[1] ** 2 - d[0] * d[1] + 1.5 * d[2] ** 2) ** 0.5; r = 250.0 / vm
(wd / "result.json").write_text(json.dumps(dict(yield_radius=r, accepted=True, stress_mandel=[r * x for x in d], cpu_s=time.process_time() - t0,
                                                 note="analytic stand-in, not a MOOSE result")))
