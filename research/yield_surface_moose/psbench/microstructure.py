"""2D periodic two-phase RVE generator (gmsh 4.4 .geo scripting, OpenCASCADE kernel).

Inclusions are kept a fixed margin away from the RVE boundary, so every boundary
edge is a single straight matrix curve; this makes the periodic node matching
required by MOOSE periodic BCs trivial and robust (gmsh `Periodic Curve`). In the
periodic tiling this leaves a matrix ligament of width 2*margin between cells.
"""

import json
import math
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

GMSH = "gmsh"


@dataclass
class Inclusion:
    """Elliptical inclusion; a circle has a == b. theta in radians, a >= b."""

    cx: float
    cy: float
    a: float
    b: float
    theta: float = 0.0

    @property
    def area(self):
        return math.pi * self.a * self.b

    def boundary(self, inflate=0.0, n=96):
        t = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
        x = (self.a + inflate) * np.cos(t)
        y = (self.b + inflate) * np.sin(t)
        c, s = math.cos(self.theta), math.sin(self.theta)
        return np.c_[self.cx + c * x - s * y, self.cy + s * x + c * y]

    def contains(self, P, inflate=0.0):
        c, s = math.cos(self.theta), math.sin(self.theta)
        dx, dy = P[:, 0] - self.cx, P[:, 1] - self.cy
        u, v = c * dx + s * dy, -s * dx + c * dy
        return (u / (self.a + inflate)) ** 2 + (v / (self.b + inflate)) ** 2 <= 1.0

    def bbox_half(self):
        c, s = math.cos(self.theta), math.sin(self.theta)
        return (math.hypot(self.a * c, self.b * s), math.hypot(self.a * s, self.b * c))


def _fits(inc, placed, L, gap, margin):
    hx, hy = inc.bbox_half()
    if not (margin + hx <= inc.cx <= L - margin - hx and margin + hy <= inc.cy <= L - margin - hy):
        return False
    # Two inclusions are separated by at least `gap` if neither one inflated by gap/2
    # contains a boundary point of the other inflated by gap/2 (sampled test; 96
    # points resolve the boundary well below the minimum gap used here).
    for o in placed:
        if math.hypot(inc.cx - o.cx, inc.cy - o.cy) > max(inc.a, o.a) + min(inc.a, o.a) + gap:
            continue
        if o.contains(inc.boundary(gap / 2), gap / 2).any() or inc.contains(
            o.boundary(gap / 2), gap / 2
        ).any():
            return False
    return True


def place_inclusions(family, params, seed, L=1.0):
    """Random sequential adsorption for the three microstructure families.

    family: 'particulate' (equal circles, uniform positions),
            'aligned' (equal ellipses, common orientation, uniform positions),
            'clustered' (lognormal circle radii, positions drawn around cluster centers).
    Inclusions are added until the target area fraction would be exceeded.
    """
    rng = np.random.default_rng(seed)
    vf, gap, margin = params["vf_target"], params["gap"], params["margin"]
    placed, area, tries = [], 0.0, 0
    if family == "clustered":
        centers = rng.uniform(0.25 * L, 0.75 * L, size=(params["n_clusters"], 2))
    while tries < 200000:
        tries += 1
        if family == "particulate":
            a = b = params["radius"]
            th = 0.0
            cx, cy = rng.uniform(0, L, 2)
        elif family == "aligned":
            a, b, th = params["a"], params["b"], math.radians(params["theta_deg"])
            cx, cy = rng.uniform(0, L, 2)
        elif family == "clustered":
            r = float(
                np.clip(rng.lognormal(math.log(params["r_median"]), params["r_sigma"]),
                        params["r_min"], params["r_max"])
            )
            a = b = r
            th = 0.0
            k = rng.integers(len(centers))
            cx, cy = centers[k] + rng.normal(0.0, params["cluster_std"] * L, 2)
        else:
            raise ValueError(family)
        inc = Inclusion(float(cx), float(cy), a, b, th)
        if (area + inc.area) / L**2 > vf + 1e-12:
            # for fixed-size families no further inclusion can fit the budget
            if family != "clustered":
                break
            continue
        if _fits(inc, placed, L, gap, margin):
            placed.append(inc)
            area += inc.area
            if area / L**2 > vf - 0.005:
                break
    return placed


def write_geo(path, inclusions, L, h):
    lines = [
        'SetFactory("OpenCASCADE");',
        "Geometry.OCCBooleanPreserveNumbering = 1;",
        f"Rectangle(1) = {{0, 0, 0, {L}, {L}}};",
    ]
    tags = []
    for k, inc in enumerate(inclusions):
        t = k + 2
        tags.append(t)
        lines.append(f"Disk({t}) = {{{inc.cx!r}, {inc.cy!r}, 0, {inc.a!r}, {inc.b!r}}};")
        if inc.theta != 0.0:
            lines.append(f"Rotate {{{{0,0,1}}, {{{inc.cx!r},{inc.cy!r},0}}, {inc.theta!r}}} {{ Surface{{{t}}}; }}")
    e = 1e-6 * L
    if tags:
        tl = ",".join(map(str, tags))
        lines += [
            f"f() = BooleanFragments{{ Surface{{1}}; Delete; }}{{ Surface{{{tl}}}; Delete; }};",
            "m() = f();",
            f"m() -= {{{tl}}};",
            'Physical Surface("matrix", 1) = {m()};',
            f'Physical Surface("inclusion", 2) = {{{tl}}};',
        ]
    else:
        lines.append('Physical Surface("matrix", 1) = {1};')
    lines += [
        f"lc() = Curve In BoundingBox{{{-e},{-e},{-e},{e},{L + e},{e}}};",
        f"rc() = Curve In BoundingBox{{{L - e},{-e},{-e},{L + e},{L + e},{e}}};",
        f"bc() = Curve In BoundingBox{{{-e},{-e},{-e},{L + e},{e},{e}}};",
        f"tc() = Curve In BoundingBox{{{-e},{L - e},{-e},{L + e},{L + e},{e}}};",
        f"Periodic Curve{{rc()}} = {{lc()}} Translate{{{L},0,0}};",
        f"Periodic Curve{{tc()}} = {{bc()}} Translate{{0,{L},0}};",
        'Physical Curve("left", 11) = {lc()};',
        'Physical Curve("right", 12) = {rc()};',
        'Physical Curve("bottom", 13) = {bc()};',
        'Physical Curve("top", 14) = {tc()};',
        f"Mesh.CharacteristicLengthMin = {h / 3};",
        f"Mesh.CharacteristicLengthMax = {h};",
        # resolve each inclusion boundary with at least ~20 segments
        "Mesh.CharacteristicLengthFromCurvature = 1;",
        "Mesh.MinimumCirclePoints = 20;",
        "Mesh.RecombineAll = 1;",
        # Frontal-Delaunay + recombination; the frontal-quad algorithm (8) crashes in gmsh 4.4.1
        "Mesh.Algorithm = 6;",
    ]
    Path(path).write_text("\n".join(lines) + "\n")


def read_msh_stats(msh):
    """Return node coordinates, per-physical-surface areas and element-type counts (msh 2.2)."""
    L = Path(msh).read_text().split("\n")
    i = L.index("$Nodes")
    n = int(L[i + 1])
    X = {}
    for l in L[i + 2 : i + 2 + n]:
        a = l.split()
        X[int(a[0])] = (float(a[1]), float(a[2]))
    i = L.index("$Elements")
    m = int(L[i + 1])
    area, types = {}, {}
    for l in L[i + 2 : i + 2 + m]:
        a = list(map(int, l.split()))
        t, nt = a[1], a[2]
        tag, nodes = a[3], a[3 + nt :]
        if t in (2, 3):  # TRI3, QUAD4
            types[t] = types.get(t, 0) + 1
            P = np.array([X[k] for k in nodes])
            x, y = P[:, 0], P[:, 1]
            area[tag] = area.get(tag, 0.0) + 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    return np.array(list(X.values())), area, types


def build_mesh(outdir, inclusions, L=1.0, h=1.0 / 40):
    """Write the .geo, mesh it, and verify areas and periodic node matching."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    geo, msh = outdir / "rve.geo", outdir / "rve.msh"
    write_geo(geo, inclusions, L, h)
    # gmsh 4.4.1 exits 1 on an unrelated GUI-option warning, so check the output instead
    res = subprocess.run([GMSH, "-2", "-format", "msh22", str(geo), "-o", str(msh)],
                         capture_output=True, text=True)
    (outdir / "gmsh.log").write_text(res.stdout + res.stderr)
    if not msh.exists() or "Done writing" not in res.stdout + res.stderr:
        raise RuntimeError(f"gmsh failed, see {outdir/'gmsh.log'}")
    X, area, types = read_msh_stats(msh)
    for c in (0, 1):
        lo = np.sort(X[np.abs(X[:, c]) < 1e-9 * L][:, 1 - c])
        hi = np.sort(X[np.abs(X[:, c] - L) < 1e-9 * L][:, 1 - c])
        if lo.size != hi.size or np.max(np.abs(lo - hi)) > 1e-9 * L:
            raise RuntimeError("periodic node matching failed")
    a_incl_exact = sum(i.area for i in inclusions)
    a_incl_mesh = area.get(2, 0.0)
    # The OCC fragment tags are checked indirectly: the meshed inclusion area must match
    # the exact ellipse area up to polygonal discretization of the interface (a few %).
    if inclusions and abs(a_incl_mesh - a_incl_exact) > 0.05 * a_incl_exact:
        raise RuntimeError(f"inclusion area mismatch {a_incl_mesh} vs {a_incl_exact}")
    stats = dict(
        n_elem=int(sum(types.values())),
        n_quad=int(types.get(3, 0)),
        n_tri=int(types.get(2, 0)),
        n_nodes=int(len(X)),
        vf_mesh=a_incl_mesh / L**2,
        vf_exact=a_incl_exact / L**2,
        h=h,
        L=L,
    )
    (outdir / "mesh_stats.json").write_text(json.dumps(stats, indent=2))
    (outdir / "inclusions.json").write_text(json.dumps([asdict(i) for i in inclusions], indent=1))
    return msh, stats
