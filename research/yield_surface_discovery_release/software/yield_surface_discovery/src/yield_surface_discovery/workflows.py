"""High-level workflows used by the command-line scripts (benchmark, closed loop, minimal example)."""
import json, time
from pathlib import Path
import numpy as np, pandas as pd
from .io import load_config, build_controller, load_eval_directions, environment_info, write_json
from .metrics import evaluate, region_masks, coverage_gaps
from .reconstruction import Reconstruction
from .analytical import MODELS


def reference_for(cfg):
    """reference radii on the evaluation directions -- only for analytical oracles; used for scoring, never by the controller."""
    ev = cfg.get("evaluation"); o = cfg["oracle"]
    if not ev or o["type"] != "analytical":
        return None, None, None
    D = load_eval_directions(ev["directions_path"]); r = MODELS[o["model"]](D, **o.get("parameters", {})); return D, r, region_masks(r, D)


def _summarise(per_run):
    m = per_run.groupby("n_real_tests")
    return pd.DataFrame(dict(n_runs=m.global_p95.size(), global_p95_median=m.global_p95.median(), global_p95_min=m.global_p95.min(), global_p95_max=m.global_p95.max(),
                             weak10_max_over_median=m.weak10_max_overprediction.median(), sharp10_max_over_median=m.sharp10_max_overprediction.median(),
                             high_shear_p95_median=m.high_shear_p95.median())).reset_index()


def run_config(config_path, output=None, starts=None, export_mesh=None, quiet=False):
    cfg = load_config(config_path); out = Path(output or cfg["output"]["directory"]); out.mkdir(parents=True, exist_ok=True)
    ctrl = build_controller(cfg); EV, rref, M = reference_for(cfg); starts = cfg.get("starts", [0]) if starts is None else starts
    report = cfg.get("budget", {}).get("report_at", [ctrl.settings.budget]); tests, steps, rows, t0 = [], [], [], time.time()
    for s in starts:
        R = ctrl.run(s); tests.append(R.tests_frame()); steps.append(pd.DataFrame(R.steps).assign(start=s))
        for n in report:
            if n > len(R.observations):
                continue
            gp = R.fit(ctrl, n); hp = gp.hyperparameters(); row = dict(start=s, n_real_tests=n, symmetry=ctrl.symmetry.name, **{f"hyp_{k}": v for k, v in hp.items() if k != "kernel"})
            if EV is not None:
                row.update(evaluate(gp.predict_radius(EV), rref, EV, M)); row.update(coverage_gaps(EV, ctrl.library.directions[R.sequence[:n]], ctrl.symmetry))
            rows.append(row)
        if (export_mesh if export_mesh is not None else cfg.get("output", {}).get("export_mesh", False)):
            rec = Reconstruction(R.fit(ctrl, len(R.observations))); rec.export(out / f"surface_start{s}.obj", level=4)
            th, Ds, rs, ss = rec.plane_section(); pd.DataFrame(dict(theta_deg=th, radius=rs, log_sd=ss)).to_csv(out / f"section_sxy0_start{s}.csv", index=False)
        write_json([o.to_dict() for o in R.observations], out / f"observations_start{s}.json")
    pd.concat(tests).to_csv(out / "tests.csv", index=False); pd.concat(steps).to_csv(out / "steps.csv", index=False)
    per_run = pd.DataFrame(rows); per_run.to_csv(out / "per_run.csv", index=False)
    summ = _summarise(per_run) if EV is not None and len(per_run) else None
    if summ is not None:
        summ.to_csv(out / "summary.csv", index=False)
    write_json(dict(config=cfg, environment=environment_info(), starts=starts, wall_s=time.time() - t0, library_size=len(ctrl.library),
                    removed_duplicates=list(ctrl.library.removed_duplicates), oracle_calls=getattr(ctrl.oracle, "calls", None)), out / "run_info.json")
    if not quiet and summ is not None:
        print(summ.to_string(index=False))
    return dict(per_run=per_run, summary=summ, tests=pd.concat(tests), output=out, controller=ctrl)


def plot_run(output_dir, start=0, title=""):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from .geometry import lonlat_deg
    out = Path(output_dir); sec = pd.read_csv(out / f"section_sxy0_start{start}.csv"); t = pd.read_csv(out / "tests.csv"); t = t[t.start == start]
    fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.4), constrained_layout=True)
    ax = axs[0]; th = sec.theta_deg.values; r = sec.radius.values; sd = sec.log_sd.values
    ax.fill_between(th, r * np.exp(-1.96 * sd), r * np.exp(1.96 * sd), color="#eb6834", alpha=0.18, lw=0, label="±1.96 sd (log radius)")
    ax.plot(th, r, color="#eb6834", lw=1.6, label="reconstruction (s_xy = 0 section)")
    near = np.abs(t.req_z3.values) < np.sin(np.radians(10)); ang = np.degrees(np.arctan2(t.req_z2.values, t.req_z1.values)) % 360
    ax.scatter(ang[near], t.yield_radius.values[near], s=22, color="#0b0b0b", zorder=4, label="tests within 10° of the section")
    ax.set_xlabel("angle in the (s_xx, s_yy) plane (deg)"); ax.set_ylabel("yield radius |z|"); ax.legend(fontsize=7, frameon=False); ax.set_xlim(0, 360)
    ax = axs[1]; lo, la = lonlat_deg(t[["req_z1", "req_z2", "req_z3"]].values)
    sc = ax.scatter(lo, la, c=t.test_number, cmap="Blues", s=30, edgecolors="#0b0b0b", lw=0.5, vmin=1, vmax=t.test_number.max())
    fig.colorbar(sc, ax=ax, label="test number"); ax.set_xlim(-180, 180); ax.set_ylim(-90, 90); ax.set_xlabel("in-plane angle atan2(z2, z1) (deg)"); ax.set_ylabel("shear latitude asin(z3) (deg)")
    for a in axs:
        a.spines[["top", "right"]].set_visible(False)
    fig.suptitle(title, x=0.01, ha="left", fontsize=9); p = out / f"run_start{start}.png"; fig.savefig(p, dpi=150); fig.savefig(out / f"run_start{start}.pdf"); plt.close(fig); return p
