#!/usr/bin/env python3
"""Raw MERFISH runs (docs/specs/raw-merfish-profile.md): recognise, read in
place, sanity-check.

A run folder holds `final_spots/` (one .npz per set: XF, genes, pix_size,
fov_size, pos_fovs_pix, fovs — de-duplicated, score-filtered) and
`final_mats[WithGFP]/` (one .npz per field: Mcells, Mgenes, Mcounts, Xum —
not de-duplicated). The check (`check`) and the conversion read it through
this module only. The spot files are tens of GB: their members are read in
place (memory-mapped) and in pieces, never whole.

    python raw_merfish.py detect RUN                 # JSON: the run, or why not
    python raw_merfish.py check --run RUN --work DIR --args A.json --out O.json

`check` writes a summary and PNG plots into DIR and its result (with the
plot paths) to O.json; progress lines go to stdout, a refusal as "ERROR: …".
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# The producer's names, in one place (the profile is their spec).
SPOTS_DIR = "final_spots"
MATS_DIRS = ("final_matsWithGFP", "final_mats")  # first found per field wins
MAT_SUFFIX = "--mat.npz"
SPOT_MEMBERS = ("XF", "genes", "pix_size", "fov_size", "pos_fovs_pix", "fovs")
MAT_MEMBERS = ("Mcells", "Mgenes", "Mcounts", "Xum")
GFP_MEMBER = "ctsGFP"
XF_Z, XF_X, XF_Y, XF_SCORE, XF_FIELD = range(5)
BLANK_PREFIX = "blank"
SET_NAME = re.compile(r"set[_-]?(\d+)", re.I)

# Light-mode chart colours (the dataviz reference palette).
INK, INK2, GRID, CONTEXT = "#191919", "#666666", "#e6e6e6", "#c9c9c9"
SERIES = ("#257de0", "#f47230", "#00b894")  # a scatter shows at most three
SEQ = ("#cce3fd", "#9dc8f7", "#6aabf1", "#338ced", "#2370cd", "#1757ac", "#0e3d7e")


class NotARun(Exception):
    """The folder is not a run in the baseline shape; the message says why."""


# ---------------------------------------------------------------- reading


class Member:
    """An .npz member read in place: rows by offset with os.pread when the
    member is stored (not deflated), so a 10 GB member is never loaded; a
    deflated one is loaded whole once. (Reading through np.memmap instead is
    slow on a network share: page faults are served one at a time.)"""

    def __init__(self, path: Path, name: str):
        self.path, self.name = path, name
        info = zipfile.ZipFile(path).getinfo(f"{name}.npy")
        self._whole: np.ndarray | None = None
        if info.compress_type != zipfile.ZIP_STORED:
            with np.load(path, allow_pickle=False) as z:
                self._whole = z[name]
            self.shape, self.dtype, self.offset = self._whole.shape, self._whole.dtype, -1
            return
        from numpy.lib import format as npf

        with open(path, "rb") as fh:
            fh.seek(info.header_offset)
            local = fh.read(30)
            n, extra = int.from_bytes(local[26:28], "little"), int.from_bytes(local[28:30], "little")
            fh.seek(info.header_offset + 30 + n + extra)
            version = npf.read_magic(fh)
            read = npf.read_array_header_1_0 if version == (1, 0) else npf.read_array_header_2_0
            shape, fortran, dtype = read(fh)
            self.offset = fh.tell()
        if dtype.hasobject:
            raise NotARun(f"{path.name}: {name} holds Python objects (a pickle), not read")
        if fortran and len(shape) > 1:
            raise NotARun(f"{path.name}: {name} is stored column-major; rows cannot be read in place")
        self.shape, self.dtype = tuple(shape), dtype

    @property
    def stored(self) -> bool:
        return self._whole is None

    def rows(self, start: int, count: int) -> np.ndarray:
        n = self.shape[0]
        start, stop = max(0, start), min(n, start + count)
        if self._whole is not None:
            return self._whole[start:stop]
        row = self.dtype.itemsize * int(np.prod(self.shape[1:], dtype=np.int64))
        fd = os.open(self.path, os.O_RDONLY)
        try:
            want = (stop - start) * row
            buf = bytearray()
            while len(buf) < want:
                got = os.pread(fd, want - len(buf), self.offset + start * row + len(buf))
                if not got:
                    raise NotARun(f"{self.path.name}: {self.name} ends early (cut short while copying?)")
                buf += got
        finally:
            os.close(fd)
        return np.frombuffer(bytes(buf), self.dtype).reshape((stop - start, *self.shape[1:]))


@dataclass
class SpotSet:
    id: str
    path: Path
    n: int
    fovs: list[str]
    pos_px: np.ndarray  # fields' x, y corners, stitched pixels
    pix: np.ndarray  # µm per z plane, x pixel, y pixel
    fov_size: np.ndarray  # planes, pixels, pixels
    stored: bool  # members stored, so read in place

    @property
    def xy_um(self) -> float:
        return float(self.pix[1])

    def centres_um(self) -> np.ndarray:
        return (self.pos_px + self.fov_size[1:3] / 2.0) * self.xy_um

    def extent_um(self) -> tuple[float, float, float, float]:
        lo = self.pos_px.min(0) * self.xy_um
        hi = (self.pos_px.max(0) + self.fov_size[1:3]) * self.xy_um
        return float(lo[0]), float(lo[1]), float(hi[0]), float(hi[1])

    def xf(self) -> Member:
        return Member(self.path, "XF")

    def genes(self) -> Member:
        return Member(self.path, "genes")


@dataclass
class Run:
    root: Path
    sets: list[SpotSet]
    mats: dict[str, Path]  # field name -> mat file
    gfp_fields: int  # fields whose mat has ctsGFP
    notes: list[str] = field(default_factory=list)

    def to_json(self) -> dict:
        return {
            "root": str(self.root),
            "sets": [{"id": s.id, "molecules": s.n, "fields": len(s.fovs),
                      "stored": s.stored} for s in self.sets],
            "mats": len(self.mats),
            "gfp_fields": self.gfp_fields,
            "notes": self.notes,
        }


def _set_id(p: Path) -> str:
    m = SET_NAME.search(p.stem)
    return m.group(1) if m else p.stem


def _set_order(s: SpotSet):
    return (0, int(s.id), "") if s.id.isdigit() else (1, 0, s.id)


def detect(folder: str | os.PathLike) -> Run:
    """The run at `folder` (or raise NotARun saying what is missing)."""
    root = Path(folder).expanduser().resolve()
    if not root.is_dir():
        raise NotARun(f"{root} is not a folder")
    spots = root / SPOTS_DIR
    if not spots.is_dir():
        hint = ""
        if any(root.glob("*/saveDecodedFinalV1*.npz")) or any(root.glob("saveDecodedFinalV1*.npz")):
            hint = (" It holds per-field saveDecodedFinalV1 files: that layout is not"
                    " read yet (it is not de-duplicated or stitched per set).")
        raise NotARun(f"No {SPOTS_DIR}/ folder in {root}: not a raw MERFISH run in the"
                      f" baseline shape (docs/specs/raw-merfish-profile.md).{hint}")
    sets: list[SpotSet] = []
    for p in sorted(spots.glob("*.npz")):
        try:
            names = {n[:-4] for n in zipfile.ZipFile(p).namelist() if n.endswith(".npy")}
        except zipfile.BadZipFile:
            raise NotARun(f"{p.name} is not a readable .npz (cut short while copying?)")
        missing = [m for m in SPOT_MEMBERS if m not in names]
        if missing:
            raise NotARun(f"{p.name} lacks {', '.join(missing)} (a spot file needs"
                          f" {', '.join(SPOT_MEMBERS)})")
        xf = Member(p, "XF")
        g = Member(p, "genes")
        if len(xf.shape) != 2 or xf.shape[1] < 5 or xf.dtype.kind != "f":
            raise NotARun(f"{p.name}: XF is {xf.dtype} {xf.shape}; expected floats, N × 5"
                          " (z, x, y, score, field)")
        if g.shape[0] != xf.shape[0] or g.dtype.kind not in "US":
            raise NotARun(f"{p.name}: genes is {g.dtype} {g.shape}; expected one name per XF row")
        with np.load(p, allow_pickle=False) as z:
            fovs = [str(f) for f in z["fovs"]]
            pos = np.asarray(z["pos_fovs_pix"], dtype=np.float64).reshape(-1, 2)
            pix = np.asarray(z["pix_size"], dtype=np.float64).ravel()
            fov_size = np.asarray(z["fov_size"], dtype=np.float64).ravel()
        if len(fovs) != len(pos) or pix.size != 3 or fov_size.size != 3:
            raise NotARun(f"{p.name}: fovs ({len(fovs)}), pos_fovs_pix ({len(pos)}),"
                          f" pix_size ({pix.size}) or fov_size ({fov_size.size}) disagree")
        stored = xf.stored and g.stored
        sets.append(SpotSet(_set_id(p), p, int(xf.shape[0]), fovs, pos, pix, fov_size, stored))
    if not sets:
        raise NotARun(f"{spots} holds no .npz spot files")
    ids = [s.id for s in sets]
    if len(set(ids)) != len(ids):
        raise NotARun(f"Two spot files name the same set: {sorted(ids)}")
    sets.sort(key=_set_order)

    mats: dict[str, Path] = {}
    for d in reversed(MATS_DIRS):  # the GFP folder overrides, field by field
        for p in (root / d).glob(f"*{MAT_SUFFIX}") if (root / d).is_dir() else []:
            mats[p.name[: -len(MAT_SUFFIX)]] = p
    gfp_fields = sum(1 for p in mats.values() if p.parent.name == MATS_DIRS[0])
    run = Run(root, sets, mats, gfp_fields)
    if not mats:
        run.notes.append(f"No {' or '.join(MATS_DIRS)}/: molecules only, no cells.")
    if not all(s.stored for s in sets):
        run.notes.append("Some spot files are compressed: they are read whole (slower, more memory).")
    return run


# ------------------------------------------------------------ gene names


class GeneCodes:
    """Turns arrays of gene names into small integer codes fast: a hash of
    each name's characters, looked up in the names seen so far (np.unique
    over 10^7 strings would sort; this is a few multiply-adds per name)."""

    def __init__(self, names: list[str] | None = None):
        self.names: list[str] = []
        self._hashes = np.zeros(0, np.uint64)
        self._codes = np.zeros(0, np.int32)
        if names:
            self._add(np.asarray(names, dtype="U"))

    @staticmethod
    def _hash(a: np.ndarray) -> np.ndarray:
        if a.dtype.kind == "S":
            a = a.astype("U")
        a = np.ascontiguousarray(a)
        width = a.dtype.itemsize // 4
        u = a.view(np.uint32).reshape(len(a), width)
        h = np.full(len(a), 1469598103934665603, np.uint64)
        with np.errstate(over="ignore"):
            for j in range(width):  # padding (code point 0) is skipped: any width hashes alike
                c = u[:, j].astype(np.uint64)
                h = np.where(c != 0, (h ^ c) * np.uint64(1099511628211), h)
        return h

    def _add(self, new: np.ndarray) -> None:
        for name in np.unique(new):
            name = str(name)
            if name in self.names:
                continue
            h = self._hash(np.asarray([name]))[0]
            if np.any(self._hashes == h):
                raise RuntimeError(f"gene name hash collision on {name!r}")
            self.names.append(name)
            order = np.argsort(np.append(self._hashes, h), kind="stable")
            self._hashes = np.append(self._hashes, h)[order]
            self._codes = np.append(self._codes, len(self.names) - 1).astype(np.int32)[order]

    def codes(self, a: np.ndarray) -> np.ndarray:
        h = self._hash(a)
        i = np.searchsorted(self._hashes, h)
        i[i >= len(self._hashes)] = 0
        hit = self._hashes[i] == h if len(self._hashes) else np.zeros(len(h), bool)
        if not hit.all():
            self._add(np.asarray(a)[~hit])
            return self.codes(a)
        return self._codes[i]

    def is_blank(self) -> np.ndarray:
        return np.array([n.lower().startswith(BLANK_PREFIX) for n in self.names], bool)


# ------------------------------------------------------------------ cells


def nearest_field(xy_um: np.ndarray, centres_um: np.ndarray) -> np.ndarray:
    """Index of the field whose centre is nearest, by Chebyshev distance —
    the pipeline's overlap rule (Step4 / Step6: KDTree.query(p=inf)). An
    exact tie (a field imaged twice at one position) goes to the lowest
    index: the saved float32 positions cannot reproduce how the lab's run
    split such a field (sub-pixel differences in the stage positions)."""
    from scipy.spatial import cKDTree

    k = min(4, len(centres_um))
    d, i = cKDTree(centres_um).query(xy_um, k=k, p=np.inf)
    if k == 1:
        return i
    i = np.where(d == d[:, :1], i, np.iinfo(i.dtype).max)
    return i.min(axis=1)


@dataclass
class FieldCells:
    set: str
    field: int
    name: str
    cells: np.ndarray  # names
    counts: np.ndarray  # cells × panel
    xum: np.ndarray  # cells × 3 (z, x, y)
    gfp: np.ndarray | None
    keep: np.ndarray  # the overlap rule


def read_cells(run: Run, workers: int = 16, progress=print) -> tuple[list[str], list[FieldCells], list[str]]:
    """Every field's cells, with the overlap rule applied (`keep`). Returns
    the panel (Mgenes), the fields and problems found (as sentences)."""
    jobs = [(s, i, f) for s in run.sets for i, f in enumerate(s.fovs) if f in run.mats]

    def load(job):
        s, i, f = job
        with np.load(run.mats[f], allow_pickle=False) as z:
            missing = [m for m in MAT_MEMBERS if m not in z.files]
            if missing:
                return job, None, f"{run.mats[f].name} lacks {', '.join(missing)}"
            return job, {
                "Mcells": z["Mcells"].astype(str), "Mgenes": z["Mgenes"].astype(str),
                "Mcounts": np.asarray(z["Mcounts"]), "Xum": np.asarray(z["Xum"], np.float64),
                "gfp": np.asarray(z[GFP_MEMBER]) if GFP_MEMBER in z.files else None,
            }, None

    with ThreadPoolExecutor(workers) as ex:
        loaded = list(ex.map(load, jobs))
    problems: list[str] = [e for _, _, e in loaded if e]
    panel: list[str] | None = None
    fields: list[FieldCells] = []
    centres = {s.id: s.centres_um() for s in run.sets}
    for (s, i, f), d, _ in loaded:
        if d is None:
            continue
        genes = list(d["Mgenes"])
        if panel is None:
            panel = genes
        elif genes != panel:
            problems.append(f"{run.mats[f].name}: its gene panel differs from the first mat's")
            continue
        xum = d["Xum"].reshape(-1, 3)
        keep = (nearest_field(xum[:, 1:3], centres[s.id]) == i) if len(xum) else np.zeros(0, bool)
        fields.append(FieldCells(s.id, i, f, d["Mcells"], d["Mcounts"].reshape(len(xum), -1),
                                 xum, d["gfp"], keep))
    progress(f"cells: {len(fields)} fields read")
    return panel or [], fields, problems


# ------------------------------------------------------------- the check


def sample_set(s: SpotSet, codes: GeneCodes, blocks: int = 8, rows: int = 50_000) -> dict:
    """A few large blocks of a set's molecules, evenly spaced through the
    file (rows are grouped by field): enough for the blank fraction, the
    score range and the overlap test in seconds. Not a spatial sample —
    within a field, rows run along x — so it is never drawn."""
    xf, g = s.xf(), s.genes()
    rows = min(rows, max(1, s.n // blocks))
    starts = np.unique(np.linspace(0, max(0, s.n - rows), blocks).astype(np.int64))
    x = np.concatenate([xf.rows(int(a), rows)[:, :5] for a in starts])
    names = np.concatenate([g.rows(int(a), rows) for a in starts])
    return {
        "x": x[:, XF_X] * s.xy_um, "y": x[:, XF_Y] * s.xy_um,
        "gene": codes.codes(names), "score": x[:, XF_SCORE], "field": x[:, XF_FIELD].astype(np.int64),
    }


def scan_genes(s: SpotSet, codes: GeneCodes, wanted: list[int], chunk: int = 4_000_000,
               progress=print) -> dict:
    """Every molecule of the wanted genes in one set: one pass through its
    gene names, and the XF rows of the chunks that hold any."""
    xf, g = s.xf(), s.genes()
    want = np.asarray(wanted)
    xs, ys, cs = [], [], []
    for a in range(0, s.n, chunk):
        c = codes.codes(g.rows(a, chunk))
        hit = np.isin(c, want)
        if hit.any():
            block = xf.rows(a, chunk)[hit]
            xs.append(block[:, XF_X] * s.xy_um)
            ys.append(block[:, XF_Y] * s.xy_um)
            cs.append(c[hit])
        if (a // chunk) % 5 == 4 or a + chunk >= s.n:
            progress(f"set {s.id}: {min(a + chunk, s.n):,} / {s.n:,} molecules scanned")
    cat = lambda v, t: np.concatenate(v) if v else np.zeros(0, t)  # noqa: E731
    return {"x": cat(xs, np.float64), "y": cat(ys, np.float64), "gene": cat(cs, np.int32)}


def check(run: Run, work: Path, genes: list[str], detail_set: str | None,
          field_metric: str, progress=print) -> dict:
    t0 = time.time()
    work.mkdir(parents=True, exist_ok=True)
    flags: list[str] = []
    notes: list[str] = list(run.notes)

    # cells first: the panel names the codes
    panel, fields, problems = read_cells(run, progress=progress) if run.mats else ([], [], [])
    flags += problems
    codes = GeneCodes(panel)

    # fields with no mat, and mats with no field
    known = {f for s in run.sets for f in s.fovs}
    for s in run.sets:
        miss = [f for f in s.fovs if f not in run.mats]
        if run.mats and miss:
            flags.append(f"Set {s.id}: {len(miss)} of {len(s.fovs)} fields have no cell mat"
                         f" ({', '.join(miss[:5])}{' …' if len(miss) > 5 else ''})")
    from scipy.spatial import cKDTree

    # a field imaged twice: two fields within 5 % of a field's width. The
    # lab's runs image their first fields again at the end of each set: that
    # is noted once; any other pair is a flag.
    reimaged = []
    for s in run.sets:
        for i, j in sorted(cKDTree(s.pos_px).query_pairs(0.05 * float(s.fov_size[1]), p=np.inf)):
            (reimaged if i < 2 and j >= len(s.fovs) - 2 else flags).append(
                f"Set {s.id}: {s.fovs[i]} and {s.fovs[j]} were imaged at nearly one position")
    if reimaged:
        notes.append(f"{len(reimaged)} fields were imaged again at the end of their set (the first"
                     " fields, repeated last): each cell goes to the nearer copy's centre, an exact"
                     " tie to the first.")
    orphans = sorted(set(run.mats) - known)
    if orphans:
        flags.append(f"{len(orphans)} cell mats belong to no set's fields"
                     f" ({', '.join(orphans[:5])}{' …' if len(orphans) > 5 else ''})")
    if run.gfp_fields and run.gfp_fields < len(run.mats):
        flags.append(f"{len(run.mats) - run.gfp_fields} of {len(run.mats)} fields have no"
                     f" {MATS_DIRS[0]} file (their cells have no GFP count)")

    # molecules: blocks of every set, in parallel (reads wait on the disk)
    progress(f"molecules: sampling {len(run.sets)} sets")
    with ThreadPoolExecutor(min(8, len(run.sets))) as ex:
        samples = dict(zip([s.id for s in run.sets],
                           ex.map(lambda s: sample_set(s, codes), run.sets)))
    progress("molecules: sampled")
    blank = codes.is_blank()

    # per set
    per_set = []
    field_rows = []
    for s in run.sets:
        sm = samples[s.id]
        fc = [f for f in fields if f.set == s.id]
        n_before = int(sum(len(f.cells) for f in fc))
        n_after = int(sum(int(f.keep.sum()) for f in fc))
        assigned = int(sum(int(f.counts[f.keep].sum()) for f in fc))
        x0, y0, x1, y1 = s.extent_um()
        fov_um = s.fov_size[1] * s.xy_um
        nf = len(s.fovs)
        mism = float(np.mean(nearest_field(np.c_[sm["x"], sm["y"]], s.centres_um()) != sm["field"])) \
            if len(sm["x"]) else 0.0
        bf = float(blank[sm["gene"]].mean()) if len(sm["gene"]) else float("nan")
        per_set.append({
            "set": s.id, "fields": nf, "molecules": s.n,
            "cells_in_mats": n_before, "cells": n_after,
            "assigned_fraction": round(assigned / s.n, 3) if s.n and fc else None,
            "blank_fraction": round(bf, 4),
            "score": [round(float(q), 3) for q in np.quantile(sm["score"], [0.01, 0.5, 0.99])]
            if len(sm["score"]) else None,
            "extent_um": [round(x1 - x0, 1), round(y1 - y0, 1)],
            "molecules_per_field": round(s.n / max(nf, 1)),
            "outside_own_field": round(mism, 4),
        })
        if mism > 0.01:
            flags.append(f"Set {s.id}: {mism:.1%} of sampled molecules sit nearer another field's"
                         " centre — its spots look not de-duplicated in the field overlaps")
        for f in fc:
            tot = f.counts.sum(1)
            real = tot - f.counts[:, blank[: f.counts.shape[1]]].sum(1)
            field_rows.append({
                "set": s.id, "field": f.field, "name": f.name, "cells": len(f.cells),
                "molecules_per_cell": float(np.median(real)) if len(real) else 0.0,
                "blank_fraction": float(f.counts[:, blank[: f.counts.shape[1]]].sum() / max(tot.sum(), 1)),
                "x": float(s.pos_px[f.field, 0] * s.xy_um), "y": float(s.pos_px[f.field, 1] * s.xy_um),
                "w": float(fov_um),
            })

    # field outliers, judged within their set: no cells, under half the
    # set's median molecules per cell, or over 3× its blank fraction
    for s in run.sets:
        rows = [r for r in field_rows if r["set"] == s.id]
        if len(rows) < 5:
            continue
        mpc = np.array([r["molecules_per_cell"] for r in rows])
        bfr = np.array([r["blank_fraction"] for r in rows])
        cells = np.array([r["cells"] for r in rows])
        live = cells > 0
        m_mpc, m_bfr = np.median(mpc[live]), np.median(bfr[live])
        for r in rows:
            why = []
            if r["cells"] == 0:
                why.append("no cells")
            else:
                if r["molecules_per_cell"] < 0.5 * m_mpc:
                    why.append(f"{r['molecules_per_cell']:.0f} molecules per cell, set median {m_mpc:.0f}")
                if r["blank_fraction"] > max(3 * m_bfr, m_bfr + 0.005):
                    why.append(f"blank fraction {r['blank_fraction']:.1%}, set median {m_bfr:.1%}")
            if why:
                r["flag"] = "; ".join(why)
        bad = [r for r in rows if "flag" in r]
        if bad:
            flags.append(f"Set {s.id}: {len(bad)} field(s) stand out — "
                         + "; ".join(f"{r['name']} ({r['flag']})" for r in bad[:6])
                         + (" …" if len(bad) > 6 else ""))

    for p in per_set:
        v = [r["molecules_per_cell"] for r in field_rows if r["set"] == p["set"] and r["cells"]]
        p["molecules_per_cell"] = round(float(np.median(v))) if v else None
    mpc = np.array([p["molecules_per_cell"] for p in per_set if p["molecules_per_cell"]], float)
    if len(mpc) > 2:
        med = float(np.median(mpc))
        for p in per_set:
            v = p["molecules_per_cell"]
            if v and (v < med / 2 or v > med * 2):
                flags.append(f"Set {p['set']}: median {v:,} molecules per cell, {v / med:.1f}× the run's ({med:,.0f})")

    dens = np.array([p["molecules_per_field"] for p in per_set], float)
    if len(dens) > 2:
        med = float(np.median(dens))
        for p, d in zip(per_set, dens):
            if d < med / 2 or d > med * 2:
                flags.append(f"Set {p['set']}: {d / med:.1f}× the median molecules per field")

    # genes against the blanks (cells' counts, after the overlap rule). The
    # floor is the median blank: a gene below it is as rare as a typical
    # misidentification. A blank far above it is a codeword decoded wrongly.
    gene_rows, weakest, hot = [], [], []
    floor = 0
    if fields:
        totals = np.zeros(len(panel), np.int64)
        for f in fields:
            totals += f.counts[f.keep].sum(0).astype(np.int64)
        pb = blank[: len(panel)]
        floor = int(np.median(totals[pb])) if pb.any() else 0
        gene_rows = [{"gene": g, "count": int(c), "blank": bool(b)} for g, c, b in zip(panel, totals, pb)]
        weakest = sorted([r for r in gene_rows if not r["blank"] and r["count"] <= floor],
                         key=lambda r: r["count"])
        hot = sorted([r for r in gene_rows if r["blank"] and floor and r["count"] > 10 * floor],
                     key=lambda r: -r["count"])
        if weakest:
            flags.append(f"{len(weakest)} genes have fewer molecules in cells than the median"
                         f" blank ({floor:,}): "
                         + ", ".join(r["gene"] for r in weakest[:12]) + (" …" if len(weakest) > 12 else ""))
        if hot:
            flags.append(f"{len(hot)} blanks have over 10× the median blank's molecules (codewords"
                         " decoded wrongly?): "
                         + ", ".join(f"{r['gene']} ({r['count']:,})" for r in hot[:6]))

    # chosen genes
    lower = {n.lower(): i for i, n in enumerate(codes.names)}
    chosen, unknown = [], []
    for g in genes:
        (chosen if g.lower() in lower else unknown).append(g)
    chosen_codes = [lower[g.lower()] for g in chosen]
    if unknown:
        near = [n for n in codes.names if any(n.lower().startswith(u.lower()[:3]) for u in unknown)
                and not n.lower().startswith(BLANK_PREFIX)]
        flags.append(f"Not in the panel: {', '.join(unknown)}"
                     + (f" (similar: {', '.join(near[:8])})" if near else ""))

    progress("plotting")
    plots = {}
    if fields:
        plots["overview"] = str(plot_overview(run, fields, codes, chosen_codes[:1], blank,
                                              work / "overview.png"))
        plots["fields"] = str(plot_fields(run, field_rows, field_metric, work / "fields.png"))
        plots["genes"] = str(plot_genes(gene_rows, floor, hot, work / "genes.png"))
    if chosen_codes:
        sid = str(detail_set) if detail_set not in (None, "") else run.sets[0].id
        s = next((s for s in run.sets if s.id == sid), None)
        if s is None:
            flags.append(f"No set {sid} (sets: {', '.join(x.id for x in run.sets)})")
        else:
            progress(f"set {s.id}: reading every molecule of {', '.join(codes.names[c] for c in chosen_codes)}")
            full = scan_genes(s, codes, chosen_codes, progress=progress)
            context = np.concatenate([f.xum[f.keep][:, 1:3] for f in fields if f.set == s.id]) \
                if fields else np.zeros((0, 2))
            plots["molecules"] = str(plot_molecules(s, full, codes, chosen_codes, context,
                                                    work / f"set{s.id}_molecules.png"))
            per_gene = {codes.names[c]: int((full["gene"] == c).sum()) for c in chosen_codes}
            flags.append(f"Set {s.id}, molecules of the chosen genes: "
                         + ", ".join(f"{g} {n:,}" for g, n in per_gene.items()))

    tot_mol = sum(p["molecules"] for p in per_set)
    summary = {
        "ok": True,
        "run": str(run.root),
        "sets": len(run.sets),
        "fields": sum(len(s.fovs) for s in run.sets),
        "molecules": tot_mol,
        "cells_in_mats": sum(p["cells_in_mats"] for p in per_set),
        "cells": sum(p["cells"] for p in per_set),
        "panel": {"genes": int((~blank[: len(panel)]).sum()) if panel else None,
                  "blanks": int(blank[: len(panel)].sum()) if panel else None},
        "blank_fraction": round(float(np.average([p["blank_fraction"] for p in per_set],
                                                 weights=[p["molecules"] for p in per_set])), 4),
        "gfp": f"{run.gfp_fields} of {len(run.mats)} fields" if run.gfp_fields else None,
        "per_set": per_set,
        "flags": flags,
        "notes": notes,
        "weakest_genes": [r["gene"] for r in weakest[:30]],
        "median_blank": floor,
        "plots": plots,
        "seconds": round(time.time() - t0, 1),
        "note": ("Blank fraction, score and the overlap test come from 8 blocks of 50,000"
                 " molecules per set; cell numbers from every mat, after the overlap rule"
                 " (a cell within float precision of two field centres may go either way)."),
    }
    (work / "summary.json").write_text(json.dumps(summary, indent=1))
    return summary


# ------------------------------------------------------------------ plots


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.size": 9, "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": GRID,
        "xtick.color": INK2, "ytick.color": INK2, "axes.titlesize": 9, "axes.titlecolor": INK,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    })
    return plt


def _grid(n: int) -> tuple[int, int]:
    cols = min(n, max(1, math.ceil(math.sqrt(n * 1.3))))
    return math.ceil(n / cols), cols


def _short(n: int) -> str:
    return f"{n / 1e9:.2f} B" if n >= 1e9 else f"{n / 1e6:.1f} M" if n >= 1e6 else f"{n:,}"


def _panel(ax, s: SpotSet, title: str) -> None:
    x0, y0, x1, y1 = s.extent_um()
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)  # image convention: y grows downwards, as the stitched pixels do
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(title, loc="left")


def plot_overview(run: Run, fields: list[FieldCells], codes: GeneCodes, chosen: list[int],
                  blank: np.ndarray, out: Path) -> Path:
    """Every set's cells at their µm position (after the overlap rule),
    coloured by one gene's molecules per cell — or by all genes' when none
    is chosen. From the mats: seconds, whatever the run's size."""
    plt = _plt()
    from matplotlib.colors import LinearSegmentedColormap, Normalize

    cmap = LinearSegmentedColormap.from_list("seq", ("#e9f2fd",) + SEQ[1:])
    per = {}
    for s in run.sets:
        fc = [f for f in fields if f.set == s.id]
        if not fc:
            continue
        xy = np.concatenate([f.xum[f.keep][:, 1:3] for f in fc])
        cnt = np.concatenate([f.counts[f.keep] for f in fc])
        v = cnt[:, chosen[0]] if chosen else cnt[:, ~blank[: cnt.shape[1]]].sum(1)
        per[s.id] = (xy, v.astype(float))
    allv = np.concatenate([v for _, v in per.values()]) if per else np.zeros(1)
    norm = Normalize(0, max(float(np.quantile(allv, 0.99)), 1.0))
    rows, cols = _grid(len(run.sets))
    fig, axes = plt.subplots(rows, cols, figsize=(2.6 * cols, 2.6 * rows + 1.0), squeeze=False)
    for ax in axes.ravel()[len(run.sets):]:
        ax.axis("off")
    for ax, s in zip(axes.ravel(), run.sets):
        if s.id in per:
            xy, v = per[s.id]
            o = np.argsort(v, kind="stable")
            ax.scatter(xy[o, 0], xy[o, 1], s=0.4, c=v[o], cmap=cmap, norm=norm, linewidths=0, rasterized=True)
        x0, _, x1, _ = s.extent_um()
        n_cells = len(per[s.id][1]) if s.id in per else 0
        _panel(ax, s, f"set {s.id}\n{n_cells:,} cells · {_short(s.n)} molecules")
    label = f"{codes.names[chosen[0]]} molecules per cell" if chosen else "molecules per cell (blanks excluded)"
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes.ravel().tolist(),
                      orientation="horizontal", fraction=0.03, pad=0.02, aspect=50, extend="max")
    cb.outline.set_visible(False)
    cb.set_label(label)
    fig.suptitle("Cells per set, at their position in the set's stitched frame (µm)",
                 x=0.01, ha="left", color=INK, fontsize=10)
    fig.savefig(out, dpi=85)
    plt.close(fig)
    return out


def plot_fields(run: Run, rows: list[dict], metric: str, out: Path) -> Path:
    plt = _plt()
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.patches import Rectangle

    cmap = LinearSegmentedColormap.from_list("seq", SEQ)
    vals = np.array([r[metric] for r in rows], float)
    lo, hi = np.quantile(vals, [0.02, 0.98]) if len(vals) else (0, 1)
    norm = Normalize(lo, hi if hi > lo else lo + 1)
    nr, nc = _grid(len(run.sets))
    fig, axes = plt.subplots(nr, nc, figsize=(2.6 * nc, 2.6 * nr + 0.9), squeeze=False)
    for ax in axes.ravel()[len(run.sets):]:
        ax.axis("off")
    for ax, s in zip(axes.ravel(), run.sets):
        for r in (r for r in rows if r["set"] == s.id):
            ax.add_patch(Rectangle((r["x"], r["y"]), r["w"], r["w"], facecolor=cmap(norm(r[metric])),
                                   edgecolor="white", linewidth=0.6))
            if "flag" in r:
                ax.add_patch(Rectangle((r["x"], r["y"]), r["w"], r["w"], fill=False,
                                       edgecolor=SERIES[1], linewidth=1.6))
                ax.text(r["x"] + r["w"] / 2, r["y"] + r["w"] / 2, str(r["field"]), ha="center",
                        va="center", fontsize=6, color=INK)
        _panel(ax, s, f"set {s.id} · {sum(1 for r in rows if r['set'] == s.id)} fields")
    label = {"molecules_per_cell": "median molecules per cell (blanks excluded)",
             "blank_fraction": "blank fraction of cells' molecules", "cells": "cells per field"}[metric]
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=axes.ravel().tolist(), orientation="horizontal", fraction=0.03,
                      pad=0.02, aspect=50)
    cb.outline.set_visible(False)
    cb.set_label(label)
    fig.suptitle("Fields per set (from the cell mats; orange outline + number: a field that stands out)",
                 x=0.01, ha="left", color=INK, fontsize=10)
    fig.savefig(out, dpi=110)
    plt.close(fig)
    return out


def plot_genes(rows: list[dict], floor: int, hot: list[dict], out: Path) -> Path:
    plt = _plt()
    order = sorted(rows, key=lambda r: -r["count"])
    x = np.arange(len(order))
    y = np.array([max(r["count"], 0.5) for r in order], float)
    b = np.array([r["blank"] for r in order])
    fig, ax = plt.subplots(figsize=(9, 3.6))
    ax.scatter(x[~b], y[~b], s=9, c=SERIES[0], linewidths=0, label=f"genes ({(~b).sum()})")
    ax.scatter(x[b], y[b], s=9, c=SERIES[1], linewidths=0, label=f"blanks ({b.sum()})")
    if floor:
        ax.axhline(floor, color=INK2, lw=1, ls=(0, (3, 3)))
        ax.text(len(order) * 0.995, floor * 1.25, "median blank", ha="right", color=INK2, fontsize=8)
    names = {r["gene"] for r in hot[:3]}
    for i, r in enumerate(order):
        if r["gene"] in names:
            ax.annotate(r["gene"], (x[i], y[i]), xytext=(6, 6), textcoords="offset points",
                        fontsize=7, color=INK)
    ax.set_yscale("log")
    ax.set_xlabel("codewords, most counted first")
    ax.set_ylabel("molecules in cells")
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1, 0.92))
    ax.set_title("Molecules per codeword against the blanks (every cell, after the overlap rule)", loc="left")
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)
    return out


def plot_molecules(s: SpotSet, full: dict, codes: GeneCodes, chosen: list[int],
                   context: np.ndarray, out: Path) -> Path:
    """The XF + genes plot: every molecule of each chosen gene in one set,
    one panel per gene, over the set's cells in grey."""
    plt = _plt()
    n = len(chosen)
    rows, cols = _grid(n)
    x0, y0, x1, y1 = s.extent_um()
    w = 4.6 if n > 1 else 7.5
    h = w * (y1 - y0) / max(x1 - x0, 1)
    fig, axes = plt.subplots(rows, cols, figsize=(w * cols, h * rows + 0.5), squeeze=False)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    for ax, c in zip(axes.ravel(), chosen):
        if len(context):
            ax.scatter(context[:, 0], context[:, 1], s=0.3, c=CONTEXT, linewidths=0, rasterized=True)
        m = full["gene"] == c
        ax.scatter(full["x"][m], full["y"][m], s=0.25, c=SERIES[0], linewidths=0, rasterized=True)
        _panel(ax, s, f"{codes.names[c]} · {int(m.sum()):,} molecules")
    fig.suptitle(f"Set {s.id}: every molecule of each gene (XF x, y), cells in grey",
                 x=0.01, ha="left", color=INK, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


# -------------------------------------------------------------------- CLI


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("detect", help="say whether a folder is a run, and what it holds")
    d.add_argument("run")
    c = sub.add_parser("check", help="summary, flags and plots of a run")
    c.add_argument("--run", required=True)
    c.add_argument("--work", required=True)
    c.add_argument("--args", help="JSON: genes, set, field_metric")
    c.add_argument("--out", help="write the result here (JSON)")
    a = ap.parse_args(argv)
    try:
        if a.cmd == "detect":
            print(json.dumps(detect(a.run).to_json(), indent=1))
            return 0
        started = time.time()
        opts = json.loads(Path(a.args).read_text()) if a.args else {}
        run = detect(a.run)
        print(f"run: {len(run.sets)} sets, {sum(s.n for s in run.sets):,} molecules,"
              f" {len(run.mats)} cell mats", flush=True)
        metric = opts.get("field_metric", "molecules_per_cell")
        if metric not in ("molecules_per_cell", "blank_fraction", "cells"):
            raise NotARun(f"field_metric must be molecules_per_cell, blank_fraction or cells, not {metric}")
        result = check(run, Path(a.work), list(opts.get("genes") or []),
                       None if opts.get("set") in (None, "") else str(opts["set"]), metric,
                       progress=lambda m: print(f"[{time.time() - started:5.1f} s] {m}", flush=True))
        text = json.dumps(result, indent=1)
        if a.out:
            Path(a.out).write_text(text)
        else:
            print(text)
        return 0
    except NotARun as e:
        print(f"ERROR: {e}", flush=True)
        return 2


if __name__ == "__main__":
    sys.exit(main())
