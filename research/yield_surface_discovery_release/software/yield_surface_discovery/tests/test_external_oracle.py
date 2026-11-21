import sys
import numpy as np
from _helpers import ROOT
from yield_surface_discovery.solvers.external import ExternalCommandOracle


def test_stub_solver_contract(tmp_path):
    o = ExternalCommandOracle([sys.executable, str(ROOT / "examples/external_solver_stub/stub_solver.py"), "--z1", "{z1}", "--z2", "{z2}", "--z3", "{z3}", "--out", "{run_dir}/result.json"], tmp_path)
    d = np.array([0.6, 0.8, 0.0]); ob = o.evaluate(d)
    assert ob.success and ob.yield_radius > 0 and np.allclose(ob.achieved_direction, d) and (tmp_path / "test_0001" / "stdout.log").exists()
    bad = ExternalCommandOracle([sys.executable, "-c", "import sys; sys.exit(2)"], tmp_path / "bad"); assert bad.evaluate(d).status == "solver_failed"
