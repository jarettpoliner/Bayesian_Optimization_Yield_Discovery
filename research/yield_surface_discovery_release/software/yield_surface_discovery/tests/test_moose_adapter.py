"""MOOSE adapter tested against a FAKE executable (tests/fake_moose_app.py); real MOOSE is not required and was not available."""
import os, stat, sys
from pathlib import Path
import numpy as np, pytest
from _helpers import ROOT
from yield_surface_discovery.solvers.moose import MooseYieldOracle
from yield_surface_discovery.analytical import von_mises_radius

TEMPLATE = ROOT / "examples" / "moose_adapter" / "single_element_yield.i"; MAT = dict(E=70000.0, nu=0.3, sigma_y=250.0, H=700.0)


def fake_exe(d):
    p = Path(d) / "fake-moose-opt"; p.write_text(f"#!/bin/sh\nexec \"{sys.executable}\" \"{ROOT / 'tests' / 'fake_moose_app.py'}\" \"$@\"\n")
    p.chmod(p.stat().st_mode | stat.S_IEXEC); return p


def run_fake_check(d):
    o = MooseYieldOracle(fake_exe(d), TEMPLATE, Path(d) / "runs", MAT); worst = 0.0
    for z in ([1, 0, 0], [-1, 0, 0], [0.6, 0.8, 0], [0, 0, 1], [0.5, -0.5, 0.7071067811865476]):
        z = np.array(z, float) / np.linalg.norm(z); ob = o.evaluate(z); assert ob.success, ob.status
        worst = max(worst, abs(ob.yield_radius / von_mises_radius(z[None], 250.0)[0] - 1))
        assert ob.metadata["achieved_vs_requested_rad"] < 1e-9 and Path(ob.metadata["runs"][0]["run_dir"], "stdout.log").exists()
    assert worst < 2e-4; return worst


def test_fake_moose_first_yield(tmp_path):
    run_fake_check(tmp_path)


def test_fake_moose_failure_and_capacity(tmp_path, monkeypatch):
    o = MooseYieldOracle(fake_exe(tmp_path), TEMPLATE, tmp_path / "runs", dict(MAT, sigma_y=1e9), capacity_escalations=1)
    assert o.evaluate(np.array([1.0, 0, 0])).status == "no_yield_within_capacity"
    monkeypatch.setenv("FAKE_MOOSE_FAIL", "1"); o2 = MooseYieldOracle(fake_exe(tmp_path), TEMPLATE, tmp_path / "runs2", MAT)
    ob = o2.evaluate(np.array([1.0, 0, 0])); assert not ob.success and ob.status == "solver_failed"
