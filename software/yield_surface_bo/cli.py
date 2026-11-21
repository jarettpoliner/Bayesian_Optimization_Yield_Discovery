"""Command-line interface.   python -m yield_surface_bo <command> ...   (or software/ysbo <command>)
  methods                         list the shipped methods, their maturity and default
  libraries                       list built-in candidate libraries
  check   CONFIG                  validate a campaign configuration (symmetry, signed directions, library, method constraints); no oracle call
  run     CONFIG [--outdir DIR]   run or resume a campaign (Mode A analytic, Mode B replay, Mode C live MOOSE / external command)
  verify-frozen                   compare frozen research sources with FROZEN_SOURCES.json
  moose-check [--exe PATH] [--check-input]   detect and validate an external MOOSE executable; optionally parse a generated TL1 input (no solve)
  truss-workspace LATTICE.json DIR          validate a truss specimen and build its workspace (no MOOSE run)"""
import argparse, json, sys
from pathlib import Path


def main(argv=None):
    ap = argparse.ArgumentParser(prog="ysbo", description="Yield Surface BO: MOOSE-connected Bayesian experimental design for yield surfaces")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("methods"); sub.add_parser("libraries"); sub.add_parser("verify-frozen")
    p = sub.add_parser("check"); p.add_argument("config")
    p = sub.add_parser("run"); p.add_argument("config"); p.add_argument("--outdir")
    p = sub.add_parser("moose-check"); p.add_argument("--exe"); p.add_argument("--check-input", action="store_true"); p.add_argument("--workspace", default="moose_check_workspace")
    p = sub.add_parser("truss-workspace"); p.add_argument("lattice"); p.add_argument("dir")
    a = ap.parse_args(argv)
    from . import core, campaign, _frozen as F
    if a.cmd == "methods":
        for k, v in core.METHODS.items():
            print(f"{k:17s} [{v['maturity']}{', DEFAULT for live MOOSE' if k == core.DEFAULT_LIVE_METHOD else ''}]  {v['label']}\n{'':19s}{v['note']}")
    elif a.cmd == "libraries":
        for k, v in core.BUILTIN_LIBRARIES.items():
            print(f"{k:20s} {v}")
    elif a.cmd == "verify-frozen":
        r = F.verify(strict=False); print(json.dumps(r, indent=1)); return 0 if r["ok"] else 1
    elif a.cmd == "check":
        cfg = campaign.load_config(a.config); lib = core.load_library(cfg.get("library", {"builtin": "tl1_signed_8192"})); c = core.check_compat(cfg, lib)
        print(json.dumps(dict(ok=c.ok, errors=c.errors, warnings=c.warnings, info=c.info), indent=1, default=str)); return 0 if c.ok else 2
    elif a.cmd == "run":
        cfg = campaign.load_config(a.config)
        if a.outdir:
            cfg["outdir"] = a.outdir
        r = campaign.run_campaign(cfg); return 0 if r["info"]["status"] == "complete" else 3
    elif a.cmd == "moose-check":
        from . import moose_truss as MT
        f = MT.find_moose(a.exe); out = dict(detection=f)
        if f["found"]:
            out["version"] = MT.moose_version(f["exe"], MT._env())
            if a.check_input:
                orc = MT.TrussMooseOracle(a.workspace, None, f["exe"]); inp = Path(a.workspace) / "check" / "run.i"; inp.parent.mkdir(parents=True, exist_ok=True)
                out["input"] = orc.write_input([1.0, 0.0, 0.0], inp)
                import subprocess
                p = subprocess.run([f["exe"], "-i", "run.i", "--check-input", "--color", "off"], cwd=inp.parent, capture_output=True, text=True, timeout=600, env=MT._env())
                out["check_input"] = dict(returncode=p.returncode, tail=(p.stdout + p.stderr)[-800:])
        print(json.dumps(out, indent=1, default=str)); return 0 if f["found"] and out.get("check_input", {}).get("returncode", 0) == 0 else 4
    elif a.cmd == "truss-workspace":
        from . import moose_truss as MT
        v = MT.prepare_workspace(a.dir, a.lattice); print(json.dumps(v, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
