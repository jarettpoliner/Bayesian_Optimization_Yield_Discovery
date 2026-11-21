"""Release tables (CSV + LaTeX) rebuilt from the bundled frozen results.  'Tests' always = real physical experiments / oracle calls."""
from pathlib import Path
import numpy as np, pandas as pd
from .figures import Bundle, effects

BUD = [8, 12, 16, 20, 24, 28]


def _distinct(g):
    """number of distinct designs among the starts (identical starts give identical metric sequences)."""
    return int(len({tuple(np.round(x.sort_values("n_tests").abs_p95.values, 12)) for _, x in g.groupby("start")}))


def _row(g, **meta):
    r = dict(**meta, starts=int(g.start.nunique()), distinct_result_sequences=_distinct(g), statistic="median [min-max] over starts", score="pinned mesh",
             note="starts 2 and 3 share the same initial 8-test set" if str(meta.get("source_arm", "")).startswith(("same_surrogate", "blind_study:blind")) and meta.get("loading_library") == "expanded physical" else "")
    for n in BUD:
        x = g[g.n_tests == n].abs_p95
        r[f"p95_{n}_median"], r[f"p95_{n}_min"], r[f"p95_{n}_max"] = (float(x.median()), float(x.min()), float(x.max())) if len(x) else (np.nan,) * 3
    x = g[g.n_tests == 28]
    for c, k in (("weak10_max_over", "weak10_max_overprediction_28_median"), ("sharp10_max_over", "sharp10_max_overprediction_28_median"),
                 ("high_shear_0.85_1_abs_p95", "high_shear_p95_28_median")):
        r[k] = float(x[c].median()) if c in x and len(x) else np.nan
    return r


def main_results(B):
    rows = []
    M = [("same_surrogate:C", "expanded physical", "FROZEN PRACTICAL METHOD", "physics-informed", "none", "shear mirror", "adaptive", "8 angular maximin", "predictive uncertainty (same GP)", "registered"),
         ("same_surrogate:A", "expanded physical", "generic blind baseline", "generic blind", "none", "none", "adaptive", "8 angular maximin", "predictive uncertainty (same GP)", "registered")]
    for lib in ("conventional", "expanded physical", "arbitrary-direction reference"):
        M.append(("blind_study:directGP_mirror_on_frozen_sequences", lib, "direct GP + mirror on frozen symmetric test sequences (library comparison)", "physics-informed (reconstruction)",
                  "none", "shear mirror", "adaptive (frozen sequences)", "4 Hill-calibration + 4 symmetric maximin", "predictive uncertainty (separate symmetric RBF GP)", "registered decomposition"))
    for lib in ("conventional", "expanded physical", "arbitrary-direction reference"):
        M.append(("pool_experiment:template_informed", lib, "template-informed residual GP (earlier method)", "template-informed", "Hill-1948 (4 oracle tests)", "shear mirror", "adaptive",
                  "4 Hill-calibration + 4 symmetric maximin", "predictive uncertainty (separate symmetric RBF GP)", "frozen existing"))
    for lib in ("conventional", "expanded physical", "arbitrary-direction reference"):
        M.append(("blind_study:blind", lib, "generic blind with separate RBF selector (blind study)", "generic blind", "none", "none", "adaptive", "8 angular maximin",
                  "predictive uncertainty (separate RBF GP)", "registered"))
    for arm, lib, name, mode, tpl, sym, des, ini, acq, reg in M:
        rows.append(_row(B.sel(arm, lib), method=name, mode=mode, analytical_template=tpl, symmetry=sym, design=des, initialization=ini, acquisition=acq,
                         loading_library=lib, registration=reg, source_arm=arm))
    F = B.csv("fuhg_same_points.csv"); F = F[~F.model.str.contains("continuous")]
    for model, mode, tpl, sym in (("template-informed GP + shear mirror", "template-informed", "Hill-1948 (4 oracle tests)", "shear mirror"),
                                  ("direct Matern-3/2 GP (generic blind)", "generic blind", "none", "none")):
        g = F[F.model == model]; r = dict(method=f"{model}, published fixed design", mode=mode, analytical_template=tpl, symmetry=sym, design="fixed (published D_N)",
                                          initialization="-", acquisition="-", loading_library="conventional", registration="frozen existing" if "template" in model else "registered",
                                          source_arm="fuhg_same_points", starts=1, distinct_result_sequences=1, statistic="single run", score="pinned mesh", note="")
        for n in BUD:
            x = g[g.N == n].abs_p95; r[f"p95_{n}_median"] = float(x.iloc[0]) if len(x) else np.nan; r[f"p95_{n}_min"] = r[f"p95_{n}_max"] = r[f"p95_{n}_median"]
        x = g[g.N == 28]; r["weak10_max_overprediction_28_median"] = float(x.weak10_max_over.iloc[0]); r["sharp10_max_overprediction_28_median"] = float(x.sharp10_max_over.iloc[0])
        r["high_shear_p95_28_median"] = np.nan; rows.append(r)
    return pd.DataFrame(rows)


def method_ablation(B):
    rows = []
    def med(arm, lib):
        return B.sel(arm, lib).groupby("n_tests").abs_p95.median()
    def paired(a, b, lib="expanded physical"):
        x, y = B.sel(a, lib), B.sel(b, lib); out = {}
        for n in BUD:
            xa = x[x.n_tests == n].set_index("start").abs_p95.sort_index(); yb = y[y.n_tests == n].set_index("start").abs_p95.sort_index(); out[n] = float((xa - yb).median())
        return out
    pairs = {"no shear-mirror symmetry": paired("same_surrogate:A", "same_surrogate:C"), "IVR instead of predictive uncertainty": paired("same_surrogate:D", "same_surrogate:C"),
             "Hill-1948 template added (residual GP)": paired("pool_experiment:template_informed", "blind_study:directGP_mirror_on_frozen_sequences")}
    for lab, s, src in effects(B):
        r = dict(change_from_frozen_method=lab, comparison=src, delta_definition="alternative - frozen practical method; median global p95 (pp)", registration="registered")
        for n in BUD:
            r[f"delta_{n}_pp_difference_of_medians"] = 100 * float(s.get(n, np.nan)); r[f"delta_{n}_pp_paired_median"] = 100 * pairs[lab][n] if lab in pairs else np.nan
        rows.append(r)
    extra = [("separate RBF selector instead of the same-GP predictive uncertainty (no symmetry)", med("blind_study:blind", "expanded physical") - med("same_surrogate:A", "expanded physical"),
              "blind study (registered) vs same-surrogate Arm A (registered); identical first 8 tests", "registered"),
             ("symmetry-aware initialisation + deduplication + deterministic ties (cleanup)", B.runs[B.runs.arm_id == "cleanup:cleaned"].groupby("n_tests").abs_p95.median() - med("same_surrogate:C", "expanded physical"),
              "symmetry-aware cleanup (registered; 10 vs 5 starts; NOT adopted: regional rule failed at 16 tests)", "registered (not adopted)")]
    for lab, s, src, reg in extra:
        r = dict(change_from_frozen_method=lab, comparison=src, delta_definition="alternative - reference; median global p95 (pp)", registration=reg)
        for n in BUD:
            r[f"delta_{n}_pp_difference_of_medians"] = 100 * float(s.get(n, np.nan)); r[f"delta_{n}_pp_paired_median"] = np.nan
        rows.append(r)
    return pd.DataFrame(rows)


def fuhg_same_points(B):
    F = B.csv("fuhg_same_points.csv"); F = F[~F.model.str.contains("continuous")].copy()
    F["mode"] = np.where(F.model.str.startswith("Fuhg"), "Fuhg et al. (2023) released model", np.where(F.model.str.contains("template"), "template-informed (ours, earlier)", "generic blind (ours)"))
    F["comparison_status"] = np.where(F.information_matched_no_symmetry, "information-matched (no symmetry)", "symmetry-assisted (not information-matched)")
    F["statistic"] = "single run on the published design"
    cols = ["N", "model", "mode", "symmetry", "comparison_status", "design", "statistic", "abs_p95", "weak10_signed_p95", "weak10_max_over", "sharp10_signed_p95", "sharp10_max_over",
            "n_fail_eval", "chord_viol"]
    return F[cols].rename(columns={"N": "n_real_tests", "abs_p95": "global_p95", "weak10_max_over": "weak10_max_overprediction", "sharp10_max_over": "sharp10_max_overprediction",
                                   "n_fail_eval": "ray_failures", "chord_viol": "chord_convexity_violations"}).sort_values(["comparison_status", "n_real_tests", "mode"])


def _pct(v):
    return "--" if not np.isfinite(v) else f"{100 * v:.2f}"


def _tex(df, cols, header, caption, label):
    lines = ["\\begin{table}[t]", "\\centering", "\\scriptsize", f"\\caption{{{caption}}}", f"\\label{{{label}}}", "\\begin{tabular}{" + "l" * len(cols) + "}", "\\hline",
             " & ".join(header) + " \\\\", "\\hline"]
    for _, r in df.iterrows():
        lines.append(" & ".join(str(c(r)) if callable(c) else str(r[c]) for c in cols).replace("%", "\\%").replace("_", "\\_") + " \\\\")
    return "\n".join(lines + ["\\hline", "\\end{tabular}", "\\end{table}"]) + "\n"


def make_all(results_dir, out_dir):
    B = Bundle(results_dir); out = Path(out_dir); out.mkdir(parents=True, exist_ok=True); paths = []
    MR = main_results(B); MR.to_csv(out / "main_results.csv", index=False)
    cell = lambda n: (lambda r: _pct(r[f"p95_{n}_median"]) if r["statistic"] == "single run" else f"{_pct(r[f'p95_{n}_median'])} [{_pct(r[f'p95_{n}_min'])}--{_pct(r[f'p95_{n}_max'])}]")
    (out / "main_results.tex").write_text(_tex(MR, ["method", "mode", "symmetry", "loading_library", lambda r: f"{r.starts} ({r.distinct_result_sequences})"] + [cell(n) for n in (8, 16, 20, 28)],
        ["Method", "Mode", "Symmetry", "Library", "Starts (distinct runs)", "8 tests", "16 tests", "20 tests", "28 tests"],
        "Global p95 relative radial error (\\%) on 3183 held-out directions, Mg--4\\%Li; median [min--max] over starts. Test counts are real experiments; shear-mirror pseudo-observations are not tests.",
        "tab:main_results")); paths += [out / "main_results.csv", out / "main_results.tex"]
    MA = method_ablation(B); MA.to_csv(out / "method_ablation.csv", index=False)
    (out / "method_ablation.tex").write_text(_tex(MA, ["change_from_frozen_method"] + [(lambda n: (lambda r: f"{r[f'delta_{n}_pp_difference_of_medians']:+.2f}"))(n) for n in (12, 16, 20, 28)] + ["registration"],
        ["Change", "12 tests", "16 tests", "20 tests", "28 tests", "Status"],
        "Change in median global p95 (percentage points) caused by one change to the frozen practical method; each row from one registered comparison (sources in the CSV).",
        "tab:method_ablation")); paths += [out / "method_ablation.csv", out / "method_ablation.tex"]
    FS = fuhg_same_points(B); FS.to_csv(out / "fuhg_same_points.csv", index=False)
    (out / "fuhg_same_points.tex").write_text(_tex(FS, ["n_real_tests", "model", "comparison_status", lambda r: _pct(r.global_p95), lambda r: _pct(r.weak10_max_overprediction),
                                                         lambda r: _pct(r.sharp10_max_overprediction), "ray_failures"],
        ["Tests", "Model", "Status", "Global p95", "Weak max over", "Sharp max over", "Ray failures"],
        "Reconstruction on the identical published test points (single run per design); only the information-matched rows compare like with like.", "tab:fuhg_same_points"))
    paths += [out / "fuhg_same_points.csv", out / "fuhg_same_points.tex"]; return paths
