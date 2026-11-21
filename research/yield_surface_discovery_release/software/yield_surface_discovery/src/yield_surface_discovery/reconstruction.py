"""Reconstructed yield surface from a fitted DirectGP: radius and uncertainty along any direction, sections, and a radial surface mesh.

The mesh is the continuous GP prediction sampled at icosphere vertices (not pinned to the observations).  The registered benchmark additionally
pinned every observed point as an exact mesh vertex; on the frozen results the global-p95 difference between the pinned mesh and the continuous
prediction is at most 0.04 percentage points (data/frozen_results/PROVENANCE.json)."""
import numpy as np
from .geometry import icosphere, mandel_to_cauchy


class Reconstruction:
    def __init__(self, gp):
        self.gp = gp

    def radius(self, directions):
        return self.gp.predict_radius(directions)

    def log_std(self, directions):
        return np.sqrt(self.gp.predict_log(directions)[1])

    def plane_section(self, n=720):
        """directions in the s_xy = 0 plane (Mandel z3 = 0), angle theta in the (s_xx, s_yy) plane."""
        th = np.linspace(0, 2 * np.pi, n, endpoint=False); D = np.stack([np.cos(th), np.sin(th), np.zeros(n)], 1)
        return np.degrees(th), D, self.radius(D), self.log_std(D)

    def mesh(self, level=4):
        U, F = icosphere(level); return U * self.radius(U)[:, None], F

    def export(self, path, level=4):
        V, F = self.mesh(level); path = str(path)
        if path.endswith(".obj"):
            with open(path, "w") as fh:
                fh.write("# yield surface in Mandel coordinates z = [s_xx, s_yy, sqrt2 s_xy]\n")
                fh.writelines(f"v {x:.10g} {y:.10g} {z:.10g}\n" for x, y, z in V); fh.writelines(f"f {a + 1} {b + 1} {c + 1}\n" for a, b, c in F)
        else:
            np.savez_compressed(path, vertices_mandel=V, vertices_cauchy=mandel_to_cauchy(V), faces=F)
        return path
