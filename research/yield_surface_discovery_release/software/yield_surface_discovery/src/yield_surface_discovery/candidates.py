"""Candidate loading libraries: the finite set of physically executable loading directions the controller may choose from.

CSV format (one row per candidate): id, label, family, z1, z2, z3 [unit Mandel direction], optional s_xx_dir, s_yy_dir, s_xy_dir, ...
Rows are normalised to unit length on load.  Exact duplicate directions (Euclidean distance < dedup_tol between unit vectors) can be merged; the
lowest original id is kept and all labels stay traceable."""
from dataclasses import dataclass
from pathlib import Path
import numpy as np, pandas as pd
from .geometry import unit, mandel_to_cauchy


@dataclass
class CandidateLibrary:
    directions: np.ndarray
    ids: np.ndarray
    labels: list
    families: list
    name: str = "candidates"
    removed_duplicates: tuple = ()

    def __len__(self):
        return len(self.directions)

    @classmethod
    def from_csv(cls, path, deduplicate=False, dedup_tol=1e-9, name=None):
        df = pd.read_csv(path, float_precision="round_trip")
        for c in ("z1", "z2", "z3"):
            if c not in df:
                raise ValueError(f"{path}: missing column {c} (unit Mandel direction z = [s_xx, s_yy, sqrt2 s_xy])")
        D = df[["z1", "z2", "z3"]].to_numpy(float); n = np.linalg.norm(D, axis=1)
        if not np.allclose(n, 1.0, atol=1e-12):
            D = unit(D)
        lib = cls(D, df["id"].to_numpy() if "id" in df else np.arange(len(df)), list(df.get("label", pd.Series(map(str, range(len(df)))))),
                  list(df.get("family", pd.Series(["unspecified"] * len(df)))), name or Path(path).stem)
        return lib.deduplicated(dedup_tol) if deduplicate else lib

    def deduplicated(self, tol=1e-9):
        keep, removed = [], []
        for i in range(len(self.directions)):
            hit = next((j for j in keep if np.linalg.norm(self.directions[i] - self.directions[j]) < tol), None)
            if hit is None:
                keep.append(i)
            else:
                removed.append((int(self.ids[i]), self.labels[i], int(self.ids[hit]), self.labels[hit]))
        k = np.array(keep)
        return CandidateLibrary(self.directions[k], self.ids[k], [self.labels[i] for i in k], [self.families[i] for i in k], self.name, tuple(removed))

    def cauchy_directions(self):
        return mandel_to_cauchy(self.directions)

    def to_frame(self):
        df = pd.DataFrame(dict(id=self.ids, label=self.labels, family=self.families)); df[["z1", "z2", "z3"]] = self.directions; return df
