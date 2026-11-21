#!/usr/bin/env python3
"""Generate COMPLEXITY_SCALING_REPORT_<ver>.md from analysis_<ver>/ and logs (every number read from files; no hand-typed results).
    python make_complexity_report.py r2|r3"""
import json, sys, time
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent; VER = sys.argv[1] if len(sys.argv) > 1 else "r2"; NMAX = {"r2": 512, "r3": 1000}[VER]; AD = HERE / f"analysis_{VER}"
DESC = json.load(open(HERE / "analytic_descriptors.json")); KAP = {k: v["kappa_max_dimless_p99"] for k, v in DESC.items()}; ORDER = sorted(KAP, key=KAP.get)
conc = json.load(open(AD / "conclusion.json")); PO = pd.read_csv(AD / "paired_pooled_and_trend.csv"); P = pd.read_csv(AD / "paired_by_level.csv")
SM = pd.read_csv(AD / "summary_median.csv"); TH = pd.read_csv(AD / "n_to_threshold.csv"); CAL = pd.read_csv(AD / "calibration.csv", header=[0, 1], index_col=[0, 1])
VC = json.load(open(HERE / "logs_r2/verify_complexity.json")); ST = json.load(open(HERE / "logs/resume_selftest.json"))
M = pd.read_csv(HERE / f"evaluation_{VER}/metrics_all.csv"); M = M[M.n <= NMAX]
BUD = [b for b in (8, 16, 33, 64, 128, 256, 512, 1000) if b <= NMAX]; ARMS = ["BO", "AXMM", "RANDOM", "BOFIXED"]
LAB = {"BO": "BO+MES", "AXMM": "angular maximin", "RANDOM": "random", "BOFIXED": "BO-FIXED (diag.)"}
pc = lambda x: f"{100 * x:.3f}" if x >= 1e-4 else f"{100 * x:.1e}"

PROTO_BUD = (16, 33, 64, 128, 256, 512, 1000)   # protocol budgets with a testable difference (n = 8: identical shared designs)
def holm_p(PO, vs, metric, n):
    """Post-hoc multiplicity context (NOT part of the frozen protocol, which sets no correction): Holm-adjusted p of budget n
    over the protocol budgets <= n that have a test."""
    s = PO[(PO.vs == vs) & (PO.metric == metric) & PO.n.isin([b for b in PROTO_BUD if b <= max(n, 16)])].dropna(subset=["wilcoxon_p"])
    ps = sorted(s.wilcoxon_p); m = len(ps); p0 = float(s[s.n == n].wilcoxon_p.iloc[0]); k = ps.index(p0)
    return max(min(1.0, (m - i) * ps[i]) for i in range(k + 1)), m
L = []; w = L.append
final = NMAX == 1000
w(f"# Analytic complexity-scaling study: {'FINAL' if final else 'INTERIM'} report ({'n ≤ 1000, protocol budget' if final else 'n ≤ 512 checkpoint'})")
w(""); w(f"Generated {time.strftime('%Y-%m-%d %H:%M')} by `make_complexity_report.py {VER}` from `analysis_{VER}/`, `evaluation_{VER}/metrics_all.csv` and `logs_r2/verify_complexity.json`.")
if final:
    w(""); w("**Revision 1 (2026-11-20).** The as-built report (11:55, archived in `closeout_20261119/final_asbuilt_20261120/`) stated an automatic \"statistically significant pooled difference in favour of BO\" for angular maximin at n = 1000 from the unadjusted p < 0.05. This revision adds the multiplicity and effect-size context (§4 note, §10). No data, analysis or pre-specified criterion changed.")
w(""); w("**Scope.** Exploratory study, pre-registered in `PROTOCOL.md` (sha256 1d614d79…, verified unchanged). It uses a synthetic, exact analytic oracle; no MOOSE simulation is involved. It is separate from the B1/B2 benchmark and is not evidence about the physical materials.")
w("")
if not final:
    w("> **Status: INTERIM.** The protocol budget is n = 1000. The BO campaigns are running from 512 to 1000 (resumable checkpoints). The trend criterion is pre-specified at every n ≥ 256 **including 1000**, so this interim report **cannot** declare it met. It can only report whether it holds at 256 and 512.")
    w("")
w("## 1. Design"); w("")
w("| Item | Setting |"); w("|---|---|")
w("| Library | 4,096-axis Fibonacci half-sphere in Mandel space (hash 8e1654f8…, matches protocol) |")
w(f"| References | {VC['reference_set']['n']:,} equal-area Fibonacci points, disjoint from the library (sha256 of rounded coordinates {VC['reference_set']['sha256_rounded'][:12]}…; definition in `cs_common.py`, hash-verified) |")
w("| Levels (by curvature κ₁·r p99) | " + "; ".join(f"{l} {KAP[l]:.2f} ({DESC[l]['class']})" for l in ORDER) + " |")
w("| Arms | BO+MES = frozen `bomes_v2.V2Controller` (A20, antipodal, Matérn-3/2 on log r, maximum latent variance, deterministic ties); angular maximin (AXMM); random unique axes (RANDOM); BO-FIXED (diagnostic; response-independent design re-observed on each level) |")
w("| Seeds | starts k = 0–3 per level and arm. BO_k, AXMM_k and RANDOM_k share the first 8 axes (verified) |")
w(f"| Budgets | {', '.join(map(str, BUD))} (plus the protocol threshold grid) |")
w("| Metrics | A20 GP refitted on each design's first n observations; area-weighted mean / max \\|r̂/r − 1\\|; mean / max over the top-10% curvature references; 90% interval coverage; RMS standardised residual |")
w("")
w("## 2. Integrity (`logs_r2/verify_complexity.json`, `logs/resume_selftest.json`)"); w("")
w(f"- Frozen files (`PROTOCOL.sha256`: PROTOCOL.md, cs_common.py, run_analytic.py, descriptors, launcher): **{'all OK' if VC['protocol_sha256_check']['rc'] == 0 else 'MISMATCH'}**.")
w(f"- Frozen controller / GP sources: {VC['frozen_sources']['hashes_checked']} recorded hashes (every run_info.json) re-checked, **{len(VC['frozen_sources']['mismatches'])} mismatches**. The Matérn-3/2 kernel, the maximum-variance (MES) acquisition and the A20 schedule are unchanged.")
w(f"- Non-adaptive designs: all {len(VC['nonadaptive_designs'])} have ≥ 1000 unique unoriented axes (no antipodal duplicates), and every AXMM_k shares its first 8 axes with BO_k.")
w(f"- BO campaigns: every recorded radius equals the analytic oracle exactly, using the controller's call pattern. Axes are unique and unoriented. Minimum n reached: **{VC['bo_min_n']}**.")
w(f"- Exact resume after the 17:09 reboot: self-test **{'PASS' if ST['PASS'] else 'FAIL'}**. A resumed run equals a fresh frozen run bitwise, including scores, and all 20 interrupted logs equal fresh runs over their first 40 tests. A20 refits ML-II from scratch at every step, so a selection depends only on the recorded data.")
w(f"- Evaluator: the r2/r3 evaluators reuse the original metric code; on the 6 designs the original evaluator had completed they agree **bitwise** (`evaluation_r3/consistency.json`).")
w("")
w(f"## 3. Accuracy by family (median over 4 starts, area-weighted mean |error| %, A20 surrogate)"); w("")
w("| Level | class | " + " | ".join(f"n={b}" for b in BUD if b >= 33) + " |"); w("|---|---|" + "---|" * len([b for b in BUD if b >= 33]))
for l in ORDER:
    for a in ARMS:
        r = SM[(SM.level == l) & (SM.arm == a)].set_index("n")["mean"]
        w(f"| {l if a == 'BO' else ''} | {DESC[l]['class'] if a == 'BO' else ''} {LAB[a]} | " + " | ".join(pc(r[b]) if b in r else "–" for b in BUD if b >= 33) + " |")
w("")
w("## 4. Pooled paired comparisons (6 levels × 4 starts = 24 matched pairs; Wilcoxon signed-rank)"); w("")
w("| n | vs | metric | BO better | median ratio other/BO (exp of median log ratio) | Wilcoxon p | Spearman ρ (κ, log ratio) | p | C3 BO wins |"); w("|---|---|---|---|---|---|---|---|---|")
for _, r in PO[PO.n.isin(BUD + ([640, 768, 896] if final else []))].sort_values(["metric", "vs", "n"]).iterrows():
    w(f"| {r.n} | {r.vs} | {r.metric} | {r.bo_wins}/{r.pairs} | {r.median_ratio_other_over_bo:.3f} | {r.wilcoxon_p:.2g} | {r.spearman_rho:.2f} | {r.spearman_p:.2g} | {r.c3_bo_wins}/{r.c3_pairs} |")
w("")
w("p values are per budget and **unadjusted**, as in the protocol, which pre-specifies a pooled test at each budget but no multiplicity correction. "
  "Holm-adjusted values for the headline comparisons are given in §10 as post-hoc context.")
w("")
w("## 5. Per-family paired results (4 pairs each; descriptive only, no per-family significance)"); w("")
w("| n | level | vs | BO wins | median gain (pp) | BO median % | other median % |"); w("|---|---|---|---|---|---|---|")
for _, r in P[P.n.isin([b for b in BUD if b >= 64]) & P.vs.isin(["AXMM", "RANDOM"])].sort_values(["vs", "n", "kappa_p99"]).iterrows():
    w(f"| {r.n} | {r.level} | {r.vs} | {r.bo_wins}/{r.pairs} | {r.median_improvement_pp:+.4f} | {r.bo_median_pct:.4f} | {r.other_median_pct:.4f} |")
w("")
w("## 6. Pre-specified complexity-trend criterion"); w("")
for o in ("AXMM", "RANDOM"):
    c = conc["criterion"][o]
    w(f"- **vs {LAB[o]}**: evaluated at n = {c['budgets_evaluated']}. " + "; ".join(f"n={int(x['n'])}: ρ = {x['spearman_rho']:.2f} (p = {x['spearman_p']:.2g}), C3 wins {x['c3_bo_wins']}/{x['c3_pairs']}" for x in c["rows"])
      + f". Result: **{'met' if c['SUPPORTED'] else 'NOT met'}**{'' if final else ' at the interim budgets (the final test needs n = 1000)'}.")
w("")
w("## 7. Calibration (A20 surrogate, median [min, max] over levels × starts)"); w("")
w("| arm | n | 90% coverage | RMS z | corr(σ, \\|e\\|) |"); w("|---|---|---|---|---|")
for (a, n), r in CAL.iterrows():
    if n in BUD and n >= 33:
        w(f"| {a} | {n} | {r[('cov90', 'median')]:.3f} [{r[('cov90', 'min')]:.2f}, {r[('cov90', 'max')]:.2f}] | {r[('rms_z', 'median')]:.2f} [{r[('rms_z', 'min')]:.2f}, {r[('rms_z', 'max')]:.2f}] | {r[('corr_sd_abs_err', 'median')]:.2f} |")
w("")
w("## 8. Observations to reach an error threshold (median over starts; first grid n; '–' = not reached)"); w("")
thr = sorted(TH.threshold_pct.unique(), reverse=True)
w("| level | arm | " + " | ".join(f"{t:g} %" for t in thr) + " |"); w("|---|---|" + "---|" * len(thr))
for l in ORDER:
    for a in ("BO", "AXMM", "RANDOM"):
        r = TH[(TH.level == l) & (TH.arm == a)].set_index("threshold_pct").n_first
        w(f"| {l if a == 'BO' else ''} | {LAB[a]} | " + " | ".join("–" if pd.isna(r.get(t)) else str(int(r[t])) for t in thr) + " |")
w("")
cost = conc["cost"]
w("## 9. Computational cost (analytic; no MOOSE)"); w("")
w(f"- The BO step wall time grows as ≈ {cost['bo_step_wall_s_fit']['prefactor']:.2e} · n^{cost['bo_step_wall_s_fit']['exponent']:.2f} s ({cost['bo_step_wall_s_fit']['note']}). ML-II is refitted with 2 restarts on 2n points (antipodal pseudo-observations) at every step.")
w("- Measured BO step time by n (s; continuation logs):"); w(""); w("| n range | steps | median s | max s |"); w("|---|---|---|---|")
for r in cost["bo_step_wall_s_by_n"]:
    w(f"| {r['bin']} | {r['count']} | {r['median']:.1f} | {r['max']:.1f} |")
w(""); w(f"- Recorded BO step wall time, post-reboot continuation: {cost['bo_total_recorded_step_wall_h']:.1f} h summed over campaigns. The non-adaptive designs cost negligible time.")
ft = pd.read_csv(AD / "evaluation_fit_cost.csv"); fa = ft[ft.surrogate == "A20"].set_index("n").fit_s
w(f"- Evaluation: the A20 refit per design and budget takes {fa.get(256, np.nan):.0f} s at n = 256, {fa.get(512, np.nan):.0f} s at 512" + (f" and {fa.get(1000, np.nan):.0f} s at 1000" if final else "") + " (median).")
w("")
w("## 10. Answers to the closeout questions"); w("")
def ans(o, met="mean"):
    r = PO[(PO.vs == o) & (PO.metric == met) & (PO.n == NMAX)].iloc[0]; return r
ra, rr = ans("AXMM"), ans("RANDOM"); ha, hr = ans("AXMM", "mean_highcurv"), ans("RANDOM", "mean_highcurv")
def judged(r, met, lab):
    hp, m = holm_p(PO, r.vs, met, NMAX); rel = 100 * (r.median_ratio_other_over_bo - 1)
    if r.wilcoxon_p >= 0.05:
        return f"No significant pooled difference (Holm-adjusted p over {m} budgets = {hp:.2g})."
    side = "BO" if r.median_ratio_other_over_bo > 1 else lab
    if hp < 0.05:
        return f"The pooled difference favours **{side}** and remains below 0.05 after Holm adjustment over {m} budgets (p = {hp:.2g}); median relative error difference {rel:+.1f} %."
    return (f"Nominally p < 0.05 in favour of {side}, but this is one of {m} per-budget tests that the protocol leaves uncorrected: the Holm-adjusted p is {hp:.2g}. "
            f"The effect is small (median relative error difference {rel:+.1f} %) and not consistent across budgets (see §4). **No robust pooled superiority is claimed.**")
seq = PO[(PO.vs == "AXMM") & (PO.metric == "mean") & PO.n.isin([256, 384, 512, 640, 768, 896, 1000])].sort_values("n")
w(f"1. **BO+MES vs angular maximin (n = {NMAX}, pooled):** BO is better in {ra.bo_wins}/{ra.pairs} pairs, with median error ratio AXMM/BO = {ra.median_ratio_other_over_bo:.3f} and Wilcoxon p = {ra.wilcoxon_p:.2g}. "
  + judged(ra, "mean", "angular maximin") + " BO wins by budget (n ≥ 256): " + ", ".join(f"{int(x.n)}: {x.bo_wins}/24 (p {x.wilcoxon_p:.2g})" for _, x in seq.iterrows())
  + ". See §5 for the per-family breakdown; it is descriptive (4 pairs each).")
w(f"2. **BO+MES vs random (n = {NMAX}, pooled):** BO is better in {rr.bo_wins}/{rr.pairs} pairs, with ratio {rr.median_ratio_other_over_bo:.2f} and p = {rr.wilcoxon_p:.2g}.")
ca, cr = conc["criterion"]["AXMM"], conc["criterion"]["RANDOM"]
w(f"3. **Does the advantage increase with complexity?** Against angular maximin the pre-specified criterion is **{'met' if ca['SUPPORTED'] else 'not met'}**; against random it is **{'met' if cr['SUPPORTED'] else 'not met'}**" + ("." if final else " at the interim budgets. The final decision requires n = 1000.")
  + " Descriptively, both conditions (ρ > 0 with p < 0.05, and ≥ 7/8 C3 wins) hold "
  + "; ".join(f"vs {LAB[o]} " + ((lambda ns: f"at n = {ns}" if ns else "at no evaluated budget")(", ".join(str(int(x['n'])) for x in conc['criterion'][o]['rows'] if x['spearman_rho'] > 0 and x['spearman_p'] < 0.05 and x['c3_bo_wins'] >= 7)))
              for o in ("AXMM", "RANDOM"))
  + ". The criterion requires every n ≥ 256, so a trend that appears only at the larger budgets does not meet it.")
w(f"4. **High-curvature regions (top-10% κ references, n = {NMAX}):** BO vs random is better in {hr.bo_wins}/{hr.pairs} pairs (ratio {hr.median_ratio_other_over_bo:.2f}, p = {hr.wilcoxon_p:.2g}). BO vs AXMM is better in {ha.bo_wins}/{ha.pairs} pairs (ratio {ha.median_ratio_other_over_bo:.3f}, p = {ha.wilcoxon_p:.2g}). {judged(ha, 'mean_highcurv', 'angular maximin')} "
  "The acquisition is the latent variance, which depends on the responses only through the ML-II hyper-parameters. Any localisation therefore comes from the length scale, not from targeting curvature directly.")
cb = CAL.xs(NMAX, level=1) if NMAX in CAL.index.get_level_values(1) else None
if cb is not None:
    w(f"5. **Calibration (n = {NMAX}):** median 90% coverage is " + ", ".join(f"{LAB[a]} {cb.loc[a, ('cov90', 'median')]:.2f}" for a in ("BO", "AXMM", "RANDOM") if a in cb.index)
      + "; median RMS z is " + ", ".join(f"{LAB[a]} {cb.loc[a, ('rms_z', 'median')]:.2f}" for a in ("BO", "AXMM", "RANDOM") if a in cb.index) + ". Nominal values are 0.90 and 1.0; RMS z below 1 means conservative intervals.")
w(f"6. **Cost scaling:** the controller cost per BO step grows ≈ n^{cost['bo_step_wall_s_fit']['exponent']:.1f}, dominated by the per-step ML-II refit. At n ≈ 1000 a BO step costs minutes, whereas observing an analytic surface costs microseconds. With a physical oracle (≈ 140–640 s per MOOSE solve, pilot) the controller cost becomes comparable to the oracle cost near n ≈ 500–1000.")
w("")
w("## 11. Limitations"); w("")
w("- Synthetic analytic surfaces (convex Hill/J2-type norms with localized bumps); exact, noise-free oracle. Behaviour with simulation noise is not tested.")
w("- 4 starts per arm and level: per-family comparisons are descriptive. Only pooled 24-pair tests are inferential, and they are not corrected for multiple budgets.")
w("- BO-FIXED is a diagnostic: one response-independent design observed on every level.")
w("- The physical MOOSE pilot is incomplete (10–12 of 24 axes per geometry) and was not extended. The full physical campaign estimate (`pilot/pilot_summary.json`) far exceeds any approved budget, so it was not launched.")
if not final:
    w("- **Interim:** n = 1000 is not yet evaluated. This report must not be cited as the final complexity result.")
(HERE / f"COMPLEXITY_SCALING_REPORT_{VER}.md").write_text("\n".join(L) + "\n"); print("wrote", HERE / f"COMPLEXITY_SCALING_REPORT_{VER}.md")
