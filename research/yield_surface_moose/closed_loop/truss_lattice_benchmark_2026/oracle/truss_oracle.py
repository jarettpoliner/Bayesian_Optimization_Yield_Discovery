"""MOOSE periodic-truss yield-radius oracle for lattice TL1 (stress-controlled, exact periodic multi-point constraints).

Kinematics (all exact linear multi-point constraints, libMesh constraint rows via FileMeshGenerator/constraint_matrix; no penalty):
  * real lattice nodes N0..N3 (unconstrained DOFs);
  * control nodes C1, C2 (unconstrained DOFs) carry the macroscopic jumps u_Ck = H L_k (H = macroscopic displacement gradient);
  * every bar that crosses the cell boundary ends on a GHOST node g = image of N_j shifted by n1 L1 + n2 L2, constrained by
        u_g = u_j + n1 u_C1 + n2 u_C2                     (periodic fluctuation + affine jump),
    so the periodic fluctuation is free and the affine part is set by (H11, H12, H22).
  * rigid translation: u(N0) = 0;  rigid rotation: u_y(C1) = 0  (i.e. H21 = 0; H12 carries the full shear, small-strain E12 = H12/2).
Stress control (Hill-Mandel conjugate loads): nodal forces F_k = A_cell * Sigma . g^k on C_k (g^k reciprocal basis, g^k.L_j = delta_kj);
  with the gauge H21 = 0 the y-component on C1 is carried by the support (moment balance).  Sigma(t) = t * Sigma_hat, Sigma_hat from the unit
  Mandel direction d = (s_xx, s_yy, sqrt2 s_xy), so t is the Mandel stress radius in MPa.
Kernels/materials: MOOSE StressDivergenceTensorsTruss + PlasticTruss (linear hardening), one material per strut group, area per family.
Observation (primary): offset yield radius r = t at which the macroscopic plastic strain  E_p = E - C_el Sigma  reaches ||E_p||_Mandel = 1e-4,
  C_el = macroscopic elastic compliance measured by MOOSE itself (3 elastic probes, `calibrate`).
"""
import json, os, re, subprocess, sys, time, hashlib, resource, shutil
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lattice"))
import lattice_def as LD                                                          # noqa: E402

MOOSE_EXE = os.environ.get("TRUSS_MOOSE_EXE", str(Path.home() / "projects/moose/modules/combined/combined-opt"))
MESH_DIR = ROOT / "oracle" / "mesh_TL1"
S2 = np.sqrt(2.0)
DEFAULTS = dict(dt=0.02, threshold=1e-4, stop_factor=1.5, t_max=200.0, nl_rel_tol=1e-10, nl_abs_tol=1e-9, exodus=False,
                n_coarse=10, switch_frac=0.95, fine_per_r1=500)   # time grid: 10 coarse elastic steps to 0.95 r1, then dt = r1/500 (see observe)


# ------------------------------------------------------------------ mesh layout
def layout():
    bars = LD.bars(); nodes = [list(x) for x in LD.NODES]; ghost = {}; elems = []
    for b in bars:
        if tuple(b["shift"]) == (0, 0):
            elems.append((b["i"], b["j"]))
        else:
            key = (b["j"],) + tuple(b["shift"])
            if key not in ghost:
                n1, n2 = b["shift"]; ghost[key] = len(nodes); nodes.append(list(LD.NODES[b["j"]] + n1 * LD.L1 + n2 * LD.L2))
            elems.append((b["i"], ghost[key]))
    c1 = len(nodes); nodes.append([-1.0, -1.0]); c2 = len(nodes); nodes.append([-2.0, -1.0])
    elems.append((c1, c2))                                                        # dummy element (block 99, no kernels) to host the C DOFs
    blocks = [LD.BLOCK_ID[b["family"]] for b in bars] + [99]
    return dict(nodes=nodes, elems=elems, blocks=blocks, ghost={f"{k[0]},{k[1]},{k[2]}": v for k, v in ghost.items()}, c1=c1, c2=c2, pin=0,
                bars=[dict(b, shift=list(b["shift"])) for b in bars], n_real=len(LD.NODES))


def constraint_matrix_text(L):
    """PETSc-ASCII (Matlab) file of the TRANSPOSED constraint operator (FileMeshGenerator transposes it back): rows = free DOF nodes,
    columns = mesh nodes (libMesh ids).  Free nodes have unit rows; ghost g: u_g = u_j + n1 u_C1 + n2 u_C2."""
    fid = L["file_id"]; free = list(range(L["n_real"])) + [L["c1"], L["c2"]]; col = {n: k for k, n in enumerate(free)}; nn = len(L["nodes"]); ent = []
    for n in free:
        ent.append((col[n], fid[n], 1.0))
    for key, g in L["ghost"].items():
        j, n1, n2 = map(int, key.split(","))
        ent.append((col[j], fid[g], 1.0))
        if n1:
            ent.append((col[L["c1"]], fid[g], float(n1)))
        if n2:
            ent.append((col[L["c2"]], fid[g], float(n2)))
    ent.sort(key=lambda e: (e[0], e[1]))
    s = ["%Mat Object: () 1 MPI process", "%  type: seqaij", f"% Size = {len(free)} {nn}", f"% Nonzeros = {len(ent)} ",
         f"zzz = zeros({len(ent)},3);", "zzz = ["]
    s += [f"{r + 1} {c + 1}  {v:.16e}" for r, c, v in ent]
    s += ["];", "Mat_0x84000001_0 = spconvert(zzz);", ""]
    return "\n".join(s)


def _moose_env():
    """runtime environment of the project's `moose` conda env (the env MOOSE trials run in, per closed_loop/README.md)."""
    pre = str(Path.home() / "micromamba/envs/moose")
    return dict(os.environ, CONDA_PREFIX=pre, PATH=pre + "/bin:" + os.environ.get("PATH", ""), OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                MKL_NUM_THREADS="1", FI_PROVIDER="tcp", MPICH_CH4_NETMOD="ofi")


def _ncdump_var(path, var):
    txt = subprocess.run([str(Path.home() / "micromamba/envs/moose/bin/ncdump"), "-v", var, str(path)], capture_output=True, text=True, check=True, env=_moose_env()).stdout
    body = txt.split(f" {var} =", 1)[1].split(";")[0]; return np.array([float(x) for x in re.findall(r"[-+0-9.eE]+", body)])


def build_mesh(force=False):
    """MOOSE builds the mesh (one ElementGenerator per bar, overlapping nodes merged, nodesets by coordinate); the node numbering of the written
    Exodus file is then READ BACK (ncdump) and mapped to the logical layout by coordinates, so the constraint matrix matches the real numbering."""
    MESH_DIR.mkdir(parents=True, exist_ok=True); L = layout()
    if (MESH_DIR / "mesh.e").exists() and (MESH_DIR / "layout.json").exists() and not force:
        return json.loads((MESH_DIR / "layout.json").read_text())
    g = ["[Mesh]", "  allow_renumbering = true"]; prev = None
    for e, ((a, b), blk) in enumerate(zip(L["elems"], L["blocks"])):
        pa, pb = L["nodes"][a], L["nodes"][b]
        g += [f"  [e{e}]", "    type = ElementGenerator"] + ([f"    input = {prev}"] if prev else []) + [
              f"    nodal_positions = '{pa[0]!r} {pa[1]!r} 0  {pb[0]!r} {pb[1]!r} 0'", "    element_connectivity = '0 1'", "    elem_type = EDGE2",
              f"    subdomain_id = {blk}", "  []"]; prev = f"e{e}"
    g += ["  [merge]", "    type = MeshRepairGenerator", f"    input = {prev}", "    fix_node_overlap = true", "  []"]; prev = "merge"
    for nm, n in (("pin", L["pin"]), ("c1", L["c1"]), ("c2", L["c2"])):
        x, y = L["nodes"][n]
        g += [f"  [ns_{nm}]", "    type = ExtraNodesetGenerator", f"    input = {prev}", f"    new_boundary = '{nm}'", f"    coord = '{x!r} {y!r} 0'", "  []"]; prev = f"ns_{nm}"
    g += ["[]", ""]
    (MESH_DIR / "mesh_gen.i").write_text("\n".join(g))
    p = subprocess.run([MOOSE_EXE, "-i", "mesh_gen.i", "--mesh-only", "mesh.e", "--color", "off"], cwd=MESH_DIR, capture_output=True, text=True, env=_moose_env())
    (MESH_DIR / "mesh_gen.log").write_text(p.stdout + p.stderr)
    if p.returncode:
        raise RuntimeError("mesh generation failed; see mesh_gen.log")
    X = np.column_stack([_ncdump_var(MESH_DIR / "mesh.e", "coordx"), _ncdump_var(MESH_DIR / "mesh.e", "coordy")])
    if len(X) != len(L["nodes"]):
        raise RuntimeError(f"node count {len(X)} != layout {len(L['nodes'])}")
    nmap = _ncdump_var(MESH_DIR / "mesh.e", "node_num_map").astype(int) - 1       # libMesh ids come from the Exodus number maps
    perm_pos = [int(np.argmin(np.linalg.norm(X - np.array(p_), axis=1))) for p_ in L["nodes"]]   # logical node -> position in file
    perm = [int(nmap[k]) for k in perm_pos]                                                       # logical node -> libMesh id
    assert len(set(perm)) == len(perm) and max(np.linalg.norm(X[perm_pos] - np.array(L["nodes"]), axis=1)) < 1e-12
    L["file_id"] = perm; L["file_coords"] = X.tolist()
    # element ids: libMesh numbers Exodus elements sequentially over the blocks in file order (connect1, connect2, ...)
    hdr = subprocess.run([str(Path.home() / "micromamba/envs/moose/bin/ncdump"), "-h", str(MESH_DIR / "mesh.e")], capture_output=True, text=True, env=_moose_env()).stdout
    nblk = int(re.search(r"num_el_blk = (\d+)", hdr).group(1)); conn = []
    for k in range(1, nblk + 1):
        c = _ncdump_var(MESH_DIR / "mesh.e", f"connect{k}").astype(int).reshape(-1, 2) - 1; conn += [tuple(sorted(x)) for x in c]
    emap = _ncdump_var(MESH_DIR / "mesh.e", "elem_num_map").astype(int) - 1
    inv = {tuple(sorted((perm_pos[a], perm_pos[b]))): k for k, (a, b) in enumerate(L["elems"])}      # connectivity is in file positions
    eid = [None] * len(L["elems"])
    for file_e, pair in enumerate(conn):
        eid[inv[pair]] = int(emap[file_e])
    assert None not in eid
    L["file_elem_id"] = eid
    (MESH_DIR / "layout.json").write_text(json.dumps(L, indent=1)); (MESH_DIR / "constraint_matrix.m").write_text(constraint_matrix_text(L))
    return L


# ------------------------------------------------------------------ loads
def reciprocal():
    M = np.column_stack([LD.L1, LD.L2]); G = np.linalg.inv(M)                    # rows of G are g^1, g^2  (g^k . L_j = delta)
    return G[0], G[1]


def unit_stress_tensor(d):
    d = np.asarray(d, float); d = d / np.linalg.norm(d); return np.array([[d[0], d[2] / S2], [d[2] / S2, d[1]]]), d


def control_forces(d):
    """forces per unit t on C1 (x only) and C2 (x, y)."""
    S, _ = unit_stress_tensor(d); g1, g2 = reciprocal(); A = LD.cell_area(); F1 = A * S @ g1; F2 = A * S @ g2
    return dict(F1x=F1[0], F1y_reaction=F1[1], F2x=F2[0], F2y=F2[1])


# ------------------------------------------------------------------ input file
def write_input(d, path, C=None, Sbar=None, **opt):
    o = dict(DEFAULTS, **opt); L = build_mesh(); F = control_forces(d); C = np.zeros((3, 3)) if C is None else np.asarray(C)
    _, du = unit_stress_tensor(d); nb = len(L["bars"]); blocks = sorted(set(L["blocks"]) - {99}); bl = " ".join(map(str, blocks))
    fam_by_block = {v: k for k, v in LD.BLOCK_ID.items()}
    grp = {g: " ".join(str(LD.BLOCK_ID[f]) for f, (gg, _, _) in LD.FAMILY.items() if gg == g) for g in LD.MATERIALS}
    md = MESH_DIR
    s = []
    s.append(f"""# TL1 periodic truss RVE, stress-controlled proportional loading. direction (Mandel unit) = {du[0]:.17g} {du[1]:.17g} {du[2]:.17g}
[Mesh]
  displacements = 'disp_x disp_y'
  allow_renumbering = false
  [fmg]
    type = FileMeshGenerator
    file = '{md / "mesh.e"}'
    constraint_matrix = '{md / "constraint_matrix.m"}'
    allow_renumbering = false
  []
[]
[GlobalParams]
  displacements = 'disp_x disp_y'
[]
[Problem]
  # block 99 is the kernel-free dummy element that only hosts the control-node DOFs C1, C2
  material_coverage_check = false
[]
[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]
[AuxVariables]
  [area]
    order = CONSTANT
    family = MONOMIAL
    block = '{bl}'
  []
""")
    for q in ("axial_stress", "plastic_stretch", "total_stretch", "hardening_variable"):
        s.append(f"  [{q}]\n    order = CONSTANT\n    family = MONOMIAL\n    block = '{bl}'\n  []\n")
    s.append("[]\n[AuxKernels]\n")
    for b in blocks:
        s.append(f"  [area_{b}]\n    type = ConstantAux\n    variable = area\n    block = {b}\n    value = {LD.FAMILY[fam_by_block[b]][1]!r}\n    execute_on = 'initial timestep_begin'\n  []\n")
    for q in ("axial_stress", "plastic_stretch", "total_stretch", "hardening_variable"):
        s.append(f"  [{q}]\n    type = MaterialRealAux\n    property = {q}\n    variable = {q}\n    block = '{bl}'\n    execute_on = 'initial timestep_end'\n  []\n")
    s.append("[]\n[Kernels]\n")
    for c, v in ((0, "disp_x"), (1, "disp_y")):
        s.append(f"  [solid_{v}]\n    type = StressDivergenceTensorsTruss\n    block = '{bl}'\n    component = {c}\n    variable = {v}\n    area = area\n  []\n")
    s.append("[]\n[Materials]\n")
    for g, m in LD.MATERIALS.items():
        s.append(f"  [mat_{g}]\n    type = PlasticTruss\n    block = '{grp[g]}'\n    youngs_modulus = {m['E']!r}\n    yield_stress = {m['sigma_y']!r}\n"
                 f"    hardening_constant = {m['H']!r}\n    absolute_tolerance = 1e-12\n    relative_tolerance = 1e-12\n  []\n")
    s.append(f"""[]
[BCs]
  [pin_x]
    type = DirichletBC
    variable = disp_x
    boundary = pin
    value = 0
  []
  [pin_y]
    type = DirichletBC
    variable = disp_y
    boundary = pin
    value = 0
  []
  [rot_gauge]
    type = DirichletBC
    variable = disp_y
    boundary = c1
    value = 0
  []
[]
[Functions]
  [f1x]
    type = ParsedFunction
    expression = '{F["F1x"]!r}*t'
  []
  [f2x]
    type = ParsedFunction
    expression = '{F["F2x"]!r}*t'
  []
  [f2y]
    type = ParsedFunction
    expression = '{F["F2y"]!r}*t'
  []
  [dtfun]
    type = ParsedFunction
    expression = 'if(t < {o.get('t_switch', 0.0) - 1e-9 * max(1.0, o.get('t_switch', 0.0))!r}, {o.get('dt_coarse', o['dt'])!r}, {o['dt']!r})'
  []
[]
[NodalKernels]
  [load_c1x]
    type = UserForcingFunctorNodalKernel
    variable = disp_x
    boundary = c1
    functor = f1x
  []
  [load_c2x]
    type = UserForcingFunctorNodalKernel
    variable = disp_x
    boundary = c2
    functor = f2x
  []
  [load_c2y]
    type = UserForcingFunctorNodalKernel
    variable = disp_y
    boundary = c2
    functor = f2y
  []
[]
[Postprocessors]
""")
    for n in range(len(L["nodes"])):
        for v in ("disp_x", "disp_y"):
            s.append(f"  [u{n}_{v[-1]}]\n    type = NodalVariableValue\n    variable = {v}\n    nodeid = {L["file_id"][n]}\n  []\n")
    for e in range(nb):
        for q, tag in (("axial_stress", "s"), ("plastic_stretch", "ep"), ("total_stretch", "eps"), ("hardening_variable", "hv")):
            s.append(f"  [{tag}{e}]\n    type = ElementalVariableValue\n    variable = {q}\n    elementid = {L["file_elem_id"][e]}\n  []\n")
    l1x, l2x, l2y = LD.L1[0], LD.L2[0], LD.L2[1]
    c1, c2 = L["c1"], L["c2"]
    # macro strain (Mandel) and plastic strain norm, with the MOOSE-calibrated compliance C
    e11 = f"(u{c1}_x/{l1x!r})"; e12 = f"((u{c2}_x-{l2x!r}*u{c1}_x/{l1x!r})/{l2y!r})"; e22 = f"(u{c2}_y/{l2y!r})"
    E = [e11, e22, f"({S2!r}*0.5*{e12})"]
    Ep = [f"({E[i]}-t*({C[i,0]!r}*{du[0]!r}+{C[i,1]!r}*{du[1]!r}+{C[i,2]!r}*{du[2]!r}))" for i in range(3)]
    pps = " ".join(f"u{c1}_x u{c2}_x u{c2}_y".split())
    for nm, ex in (("E11", E[0]), ("E22", E[1]), ("E12m", E[2]), ("ep_norm", f"sqrt({Ep[0]}^2+{Ep[1]}^2+{Ep[2]}^2)")):
        s.append(f"  [{nm}]\n    type = ParsedPostprocessor\n    pp_names = '{pps}'\n    expression = '{ex}'\n    use_t = true\n  []\n")
    if Sbar is not None:
        # Levin macro plastic strain  Ep_k = (1/A) sum_b V_b Sbar[b,k] eps_p,b  (Mandel); Sbar = MOOSE-measured bar-stress influence matrix
        Sb = np.asarray(Sbar); V = AREAS * LD.A_BAR / LD.cell_area(); comp = []
        for k in range(3):
            comp.append("(" + "+".join(f"{V[b] * Sb[b, k]!r}*ep{b}" for b in range(nb)) + ")")
        s.append(f"  [epL_norm]\n    type = ParsedPostprocessor\n    pp_names = '{' '.join(f'ep{b}' for b in range(nb))}'\n"
                 f"    expression = 'sqrt({comp[0]}^2+{comp[1]}^2+{comp[2]}^2)'\n  []\n")
    s.append("[]\n")
    stop_pp = "epL_norm" if Sbar is not None else "ep_norm"
    if o["stop_factor"]:
        s.append(f"[UserObjects]\n  [stop]\n    type = Terminator\n    expression = '{stop_pp} > {o['stop_factor'] * o['threshold']!r}'\n    fail_mode = HARD\n"
                 f"    error_level = NONE\n    execute_on = 'timestep_end'\n  []\n[]\n")
    s.append(f"""[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]
[Executioner]
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu petsc'
  nl_rel_tol = {o['nl_rel_tol']!r}
  nl_abs_tol = {o['nl_abs_tol']!r}
  nl_max_its = 30
  l_tol = 1e-12
  l_max_its = 50
  end_time = {o['t_max']!r}
  [TimeStepper]
    type = FunctionDT
    function = dtfun
    cutback_factor_at_failure = 0.5
  []
[]
[Outputs]
  csv = true
  exodus = {str(bool(o['exodus'])).lower()}
  [pc]
    type = PerfGraphOutput
    level = 0
    execute_on = final
  []
[]
""")
    Path(path).write_text("".join(s)); return L


# ------------------------------------------------------------------ run + parse
def _run_moose(workdir, inp="run.i", timeout=600):
    env = _moose_env()
    r0 = resource.getrusage(resource.RUSAGE_CHILDREN); t0 = time.time()
    p = subprocess.run([MOOSE_EXE, "-i", inp, "--color", "off", "--no-gdb-backtrace"], cwd=workdir, capture_output=True, text=True, timeout=timeout, env=env)
    r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    (Path(workdir) / "moose.log").write_text(p.stdout + "\n--- stderr ---\n" + p.stderr)
    return p.returncode, time.time() - t0, (r1.ru_utime - r0.ru_utime) + (r1.ru_stime - r0.ru_stime), p.stdout


def read_csv(workdir):
    import pandas as pd
    f = Path(workdir) / "run_out.csv"; return pd.read_csv(f)


def bar_geometry():
    out = []
    for b in LD.bars():
        n1, n2 = b["shift"]; X = LD.NODES[b["j"]] + n1 * LD.L1 + n2 * LD.L2 - LD.NODES[b["i"]]; out.append(X / np.linalg.norm(X))
    return np.array(out)


AREAS = np.array([LD.FAMILY[b["family"]][1] for b in LD.bars()])
SY = np.array([LD.MATERIALS[LD.FAMILY[b["family"]][0]]["sigma_y"] for b in LD.bars()])
FAMS = [b["family"] for b in LD.bars()]


def homogenised_stress(sig):
    """Mandel Sigma from bar stresses (reference geometry): (1/A) sum N_b a t (x) t."""
    T = bar_geometry(); S = np.einsum("...b,bi,bj->...ij", np.asarray(sig) * AREAS * LD.A_BAR, T, T) / LD.cell_area()
    return np.stack([S[..., 0, 0], S[..., 1, 1], S2 * S[..., 0, 1]], -1)


def nominal_stress_from_bars(df, L):
    """controlled components of the nominal (1st Piola-Kirchhoff) macro stress from MOOSE bar forces in the CURRENT bar orientation:
    P = (1/A0) sum_b N_b t_cur (x) X_ref  ->  (P11, P22, sqrt2 P12)  [P12 = row x, column y; P21 is the rotation-gauge reaction].
    Also returns the small-strain symmetric estimate with reference orientation (difference = geometric nonlinearity)."""
    nb = len(FAMS); bars = L["bars"]; c1, c2 = L["c1"], L["c2"]; P = np.zeros((len(df), 2, 2))
    U = lambda n: np.column_stack([df[f"u{n}_x"].to_numpy(), df[f"u{n}_y"].to_numpy()])
    for e, b in enumerate(bars):
        n1, n2 = b["shift"]; X = LD.NODES[b["j"]] + n1 * LD.L1 + n2 * LD.L2 - LD.NODES[b["i"]]
        du = U(b["j"]) + n1 * U(c1) + n2 * U(c2) - U(b["i"]); x = X[None, :] + du; t = x / np.linalg.norm(x, axis=1, keepdims=True)
        N = df[f"s{e}"].to_numpy() * AREAS[e]; P += N[:, None, None] * t[:, :, None] * X[None, None, :]
    P /= LD.cell_area(); return np.stack([P[:, 0, 0], P[:, 1, 1], S2 * P[:, 0, 1]], -1)


def levin_plastic_strain(ep, Sbar):
    """macroscopic plastic strain (Mandel) from bar plastic strains by Levin/reciprocity: Ep = (1/A) Sbar^T (V * eps_p).
    Exact for the small-strain truss (identical to E - C Sigma); unlike E - C Sigma it contains no nonlinear-elastic (geometric) strain."""
    V = AREAS * LD.A_BAR / LD.cell_area(); return (np.asarray(ep) * V) @ np.asarray(Sbar)


def analyse(df, d, C, threshold=1e-4, L=None, Sbar=None):
    """offset-yield extraction + all per-run diagnostics from the MOOSE CSV."""
    _, du = unit_stress_tensor(d); nb = len(FAMS)
    df = df[df.time > 0].reset_index(drop=True); t = df.time.to_numpy()
    sig = df[[f"s{e}" for e in range(nb)]].to_numpy(); ep = df[[f"ep{e}" for e in range(nb)]].to_numpy()
    E = df[["E11", "E22", "E12m"]].to_numpy(); Ep_tot = E - t[:, None] * (C @ du)[None, :]; n_tot = np.linalg.norm(Ep_tot, axis=1)
    Ep = levin_plastic_strain(ep, Sbar); n = np.linalg.norm(Ep, axis=1)
    rec = dict(n_steps=int(len(t)), t_last=float(t[-1]), criterion="levin_macro_plastic_strain_mandel_norm")
    kk = int(np.argmax(n_tot >= threshold))                       # diagnostic only: superseded total-strain offset definition (E - C Sigma)
    if n_tot[kk] >= threshold and kk > 0:
        rec["r_diag_total_strain_offset"] = float(t[kk - 1] + (threshold - n_tot[kk - 1]) / (n_tot[kk] - n_tot[kk - 1]) * (t[kk] - t[kk - 1]))
    # first local yield (elastic linearity: first step is elastic; r1 = t1 * min(sigma_y/|s|))
    k0 = 0; r1_lin = float(t[k0] * np.min(SY / np.abs(sig[k0]))); fb = int(np.argmin(SY / np.abs(sig[k0])))
    first_plastic_step = int(np.argmax(np.any(np.abs(ep) > 0, axis=1))) if np.any(np.abs(ep) > 0) else None
    rec.update(r1_first_yield=r1_lin, r1_bar=fb, r1_family=FAMS[fb], first_plastic_step_t=None if first_plastic_step is None else float(t[first_plastic_step]))
    k = int(np.argmax(n >= threshold))
    if n[k] < threshold or k == 0:
        rec.update(yield_radius=float("nan"), crossed=False); return rec
    w = (threshold - n[k - 1]) / (n[k] - n[k - 1]); r = float(t[k - 1] + w * (t[k] - t[k - 1]))
    at = lambda X: X[k - 1] + w * (X[k] - X[k - 1])
    Sh = nominal_stress_from_bars(df, L); Sa = t[:, None] * du[None, :]; Ss = homogenised_stress(sig)
    eq_err = float(np.max(np.linalg.norm(Sh - Sa, axis=1) / np.linalg.norm(Sa, axis=1)))
    geo_nl = float(np.max(np.linalg.norm(Ss - Sa, axis=1) / np.linalg.norm(Sa, axis=1)))
    sh_y = at(Sh); dir_err = float(np.degrees(np.arccos(np.clip(sh_y @ du / np.linalg.norm(sh_y), -1, 1))))
    # plastic work density (per cell area) and solid-averaged |eps_p|
    vol = AREAS * LD.A_BAR; dW = 0.5 * (sig[1:] + sig[:-1]) * np.diff(ep, axis=0); Wp = np.concatenate([[0.0], np.cumsum(dW @ vol)]) / LD.cell_area()
    ep_solid = (np.abs(ep) @ vol) / vol.sum()
    epy = at(ep); sy = at(sig)
    rec.update(yield_radius=r, crossed=True, step_at_yield=k, Ep_at_yield=at(Ep).tolist(), Ep_total_strain_def_at_yield=at(Ep_tot).tolist(), E_at_yield=at(E).tolist(),
               stress_mandel=(r * du).tolist(), stress_homog_from_bars=sh_y.tolist(), achieved_dir_err_deg=dir_err, equilibrium_rel_err_max=eq_err,
               small_strain_homog_vs_applied_max=geo_nl,
               Wp_at_yield=float(at(Wp)), ep_solid_avg_at_yield=float(at(ep_solid)), ep_cell_avg_at_yield=float(at(ep_solid) * vol.sum() / LD.cell_area()),
               bar_sigma_at_yield=sy.tolist(), bar_ep_at_yield=epy.tolist(), yielded_bars=[int(b) for b in np.where(np.abs(epy) > 0)[0]],
               yielded_families=sorted({FAMS[b] for b in np.where(np.abs(epy) > 0)[0]}),
               plastic_share_by_family={f: float(sum(abs(epy[b]) * vol[b] for b in range(nb) if FAMS[b] == f) / max(np.abs(epy) @ vol, 1e-300)) for f in LD.FAMILY},
               offset_to_first_yield_ratio=float(r / r1_lin))
    return rec


def periodic_kinematics_check(df, L):
    """max |u_g - (u_j + n1 u_C1 + n2 u_C2)| over ghosts and steps, relative to max |u|."""
    c1, c2 = L["c1"], L["c2"]; err, umax = 0.0, 1e-300
    for key, g in L["ghost"].items():
        j, n1, n2 = map(int, key.split(","))
        for c in "xy":
            pred = df[f"u{j}_{c}"] + n1 * df[f"u{c1}_{c}"] + n2 * df[f"u{c2}_{c}"]; err = max(err, float(np.max(np.abs(df[f"u{g}_{c}"] - pred))))
    for n in range(len(L["nodes"])):
        for c in "xy":
            umax = max(umax, float(np.max(np.abs(df[f"u{n}_{c}"]))))
    return dict(ghost_constraint_abs_err=err, rel=err / umax, pin_abs=float(max(np.abs(df["u0_x"]).max(), np.abs(df["u0_y"]).max())),
                c1_y_abs=float(np.abs(df[f"u{c1}_y"]).max()))


def schedule(d, cal, o):
    """time grid from the MOOSE elastic calibration: r1_pred = exact first-yield radius by elastic linearity (bar influence matrix measured by
    MOOSE); n_coarse equal steps up to switch_frac*r1_pred, then dt = r1_pred/fine_per_r1."""
    _, du = unit_stress_tensor(d); sb = np.asarray(cal["Sbar"]) @ du; r1 = float(np.min(SY / np.abs(sb)))
    ts = o["switch_frac"] * r1; return dict(t_switch=ts, dt_coarse=ts / o["n_coarse"], dt=r1 / o["fine_per_r1"], t_max=50.0 * r1), r1


def observe(d, workdir, cal, keep=False, **opt):
    """one MOOSE yield-radius experiment.  cal = calibrate() output (C, Sbar).  returns record (json-able)."""
    workdir = Path(workdir); workdir.mkdir(parents=True, exist_ok=True); o = dict(DEFAULTS, **opt); C = np.asarray(cal["C"])
    sch, r1_pred = schedule(d, cal, o); o.update(sch); opt = dict(opt, **sch)
    L = write_input(d, workdir / "run.i", C=C, Sbar=cal["Sbar"], **opt)
    rc, wall, cpu, out = _run_moose(workdir)
    _, du = unit_stress_tensor(d)
    rec = dict(direction=du.tolist(), r1_predicted_elastic=r1_pred, returncode=rc, wall_s=wall, cpu_s=cpu, options=o, input_sha256=hashlib.sha256((workdir / "run.i").read_bytes()).hexdigest())
    nl = [int(x) for x in re.findall(r"(\d+) Nonlinear \|R\|", out)]
    rec["nonlinear_residual_evals"] = len(nl); rec["failed_steps"] = len(re.findall(r"Solve Did NOT Converge|Solve Failed", out))
    if rc != 0 or not (workdir / "run_out.csv").exists():
        rec.update(accepted=False, reason="moose_error"); return rec
    df = read_csv(workdir); rec.update(analyse(df, d, np.asarray(C), o["threshold"], L, cal["Sbar"])); rec["periodic_check"] = periodic_kinematics_check(df, L)
    flags = []
    if not rec.get("crossed"):
        flags.append("threshold_not_reached")
    if rec["failed_steps"]:
        flags.append("failed_steps")
    if rec.get("equilibrium_rel_err_max", 1) > 1e-6:          # nominal-stress equilibrium is exact up to the nonlinear tolerance
        flags.append("equilibrium")
    if rec.get("achieved_dir_err_deg", 90) > 1e-3:
        flags.append("direction")
    if rec["periodic_check"]["rel"] > 1e-8:
        flags.append("periodicity")
    r = rec.get("yield_radius", np.nan)
    rec["flags"] = flags; rec["accepted"] = (not flags) and np.isfinite(r) and r > 0
    if not keep:
        for f in workdir.glob("*.e"):
            f.unlink()
    return rec


def calibrate(workdir, dt=1.0):
    """MOOSE elastic compliance C (Mandel) and bar-stress influence matrix Sbar (12x3, bar stress per unit Mandel stress) from three
    single-step elastic probes at Mandel radius dt (must be below first yield; checked: max |eps_p| = 0 recorded)."""
    C = np.zeros((3, 3)); Sb = np.zeros((12, 3)); recs = []
    for k in range(3):
        d = np.zeros(3); d[k] = 1.0; w = Path(workdir) / f"cal_{k}"; w.mkdir(parents=True, exist_ok=True)
        write_input(d, w / "run.i", C=None, dt=dt, t_max=dt, stop_factor=0); rc, wall, cpu, _ = _run_moose(w)
        df = read_csv(w); row = df[df.time > 0].iloc[-1]; C[:, k] = np.array([row.E11, row.E22, row.E12m]) / row.time
        Sb[:, k] = np.array([row[f"s{e}"] for e in range(12)]) / row.time
        recs.append(dict(k=k, rc=rc, wall_s=wall, cpu_s=cpu, max_abs_ep=float(np.max(np.abs(row[[f"ep{e}" for e in range(12)]].to_numpy())))))
    return dict(C=C.tolist(), Sbar=Sb.tolist(), probes=recs, probe_radius=dt)
