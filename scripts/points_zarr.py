"""Points in Zarr — reference writer and reader for docs/points-zarr-spec.md.

A points element is written as levels (level 0 every point, each further
level a uniform quarter of the one below), each level in two row orders —
tile-major and feature-major — with an index of row numbers.

The writer works out of core, so it handles billions of points:

1. ``spill`` (pass A, any number of worker processes): give each point a
   random 32-bit ``r`` and append it to the file of its horizontal band.
2. ``finish``, for each level:
   - pass B: load one band at a time, sort by (tile, feature, r), append to
     ``tile_major``; spill each feature's rows to its own file, and the points
     of the next level (``r < 2^(32 - 2(k+1))``) to that level's bands;
   - pass C: load one feature at a time, keep tile-major order, append to
     ``feature_major``;
   - write the index from the (tile × feature) counts.

Every Zarr write covers whole shards, except the last one of each array.

Usage (checks and quick reads)::

    python points_zarr.py validate STORE.zarr/points/transcripts
    python points_zarr.py info STORE.zarr/points/transcripts
"""

from __future__ import annotations

import json
import math
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import zarr
from numcodecs.zarr3 import Delta
from zarr.codecs import BloscCodec

CHUNK = 65_536
SHARD = CHUNK * 64

SPEC_VERSION = "0.1"
MULTISCALES = {
    "uuid": "d35379db-88df-4056-af3a-620245f8e347",
    "name": "multiscales",
    "schema_url": "https://raw.githubusercontent.com/zarr-conventions/multiscales/refs/tags/v0.1/schema.json",
    "spec_url": "https://github.com/zarr-conventions/multiscales/blob/v0.1/README.md",
}
DATAFRAME = {"uuid": "1f1b36d2-3185-6da8-91f6-9a0823af2d33", "name": "df:"}
ENUM = {"uuid": "c906c423-56ed-413c-9943-7b2ff52d18f2", "name": "enum:"}
POINTS = {
    "name": "points:",
    "description": "Multiscale points: levels, tile- and feature-major orders, row index",
    "spec_url": "https://github.com/kresnajenie/merfisheyes/blob/develop/docs/points-zarr-spec.md",
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _blosc(typesize: int) -> BloscCodec:
    return BloscCodec(cname="zstd", clevel=3, shuffle="shuffle", typesize=typesize)


def _array(group, name, shape, dtype, *, delta=False, chunk=None, shard=-1):
    chunk = chunk or CHUNK
    shard = SHARD if shard == -1 else shard
    dtype = np.dtype(dtype)
    n = max(int(shape), 1)
    if n <= chunk:
        chunks, shards = (n,), None
    else:
        chunks = (chunk,)
        shards = (min(shard, -(-n // chunk) * chunk),) if shard else None
    return group.create_array(
        name,
        shape=(int(shape),),
        dtype=dtype,
        chunks=chunks,
        shards=shards,
        filters=[Delta(dtype=dtype.str)] if delta else None,
        compressors=[_blosc(dtype.itemsize)],
        fill_value=0,
    )


def _strings(group, name, values):
    a = group.create_array(name, shape=(len(values),), dtype=str, chunks=(max(len(values), 1),))
    if len(values):
        a[:] = np.array(values, dtype=object)
    return a


class _Writer:
    """Appends rows to aligned 1-D arrays, writing whole shards."""

    def __init__(self, arrays: dict):
        self.arrays = arrays
        self.buf = {k: [] for k in arrays}
        self.buffered = 0
        self.offset = 0

    def append(self, cols: dict) -> None:
        n = len(next(iter(cols.values())))
        for k, v in cols.items():
            self.buf[k].append(v)
        self.buffered += n
        if self.buffered >= SHARD:
            self._flush(self.buffered // SHARD * SHARD)

    def _flush(self, n: int) -> None:
        for k, parts in self.buf.items():
            whole = np.concatenate(parts) if len(parts) > 1 else parts[0]
            self.arrays[k][self.offset : self.offset + n] = whole[:n]
            self.buf[k] = [whole[n:]] if len(whole) > n else []
        self.offset += n
        self.buffered -= n

    def close(self) -> int:
        if self.buffered:
            self._flush(self.buffered)
        return self.offset


class PointsWriter:
    """Writes one points element. Picklable, so workers can ``spill``.

    columns: extra per-point columns, ``[(name, dtype, enum_values or None)]``.
    origin, extent: the level-0 tile grid covers ``[origin, origin + extent)``.
    """

    def __init__(
        self,
        out: Path,
        tmp: Path,
        *,
        features: list[str],
        feature_key: str = "gene",
        dims: int = 3,
        columns: list[tuple[str, str, list[str] | None]] = (),
        origin: tuple[float, float],
        extent: tuple[float, float],
        tile_size: float = 100.0,
        band_um: float = 800.0,
        n_levels: int = 8,
        units: str = "micrometer",
        coordinate_system: str = "global",
    ):
        self.out, self.tmp = Path(out), Path(tmp)
        self.features = list(features)
        self.feature_key = feature_key
        self.dims = dims
        self.columns = [tuple(c) for c in columns]
        self.origin = (float(origin[0]), float(origin[1]))
        self.extent = (float(extent[0]), float(extent[1]))
        self.tile_size = float(tile_size)
        self.band_um = float(band_um)
        self.n_levels = n_levels
        self.units = units
        self.coordinate_system = coordinate_system
        self.fdtype = np.dtype("u2" if len(self.features) < 2**16 else "u4")
        # leftovers of an earlier, interrupted run would be appended to
        shutil.rmtree(self.tmp, ignore_errors=True)

    # ---- layout helpers -------------------------------------------------
    @property
    def axes(self):
        return ["x", "y", "z"][: self.dims]

    def grid(self, k: int):
        size = self.tile_size * 2**k
        nx = max(1, math.ceil(self.extent[0] / size))
        ny = max(1, math.ceil(self.extent[1] / size))
        band_rows = max(1, round(self.band_um / size))
        return size, nx, ny, band_rows

    def tiles(self, x, y, k):
        size, nx, ny, _ = self.grid(k)
        tx = np.clip(((x - self.origin[0]) / size).astype(np.int64), 0, nx - 1)
        ty = np.clip(((y - self.origin[1]) / size).astype(np.int64), 0, ny - 1)
        return ty * nx + tx, ty

    def record(self, with_id: bool):
        f = [(a, "<f4") for a in self.axes] + [("r", "<u4"), ("f", self.fdtype.str)]
        f += [(name, np.dtype(dt).str) for name, dt, _ in self.columns]
        if with_id:
            f.append(("id", "<u8"))
        return np.dtype(f)

    def feature_record(self):
        return np.dtype([(a, "<f4") for a in self.axes] + [("row", "<u8"), ("id", "<u8")])

    # ---- pass A ---------------------------------------------------------
    def spill(self, cols: dict, seed: int, tag: str) -> dict:
        """Append points to level-0 band files. ``cols`` holds x, y(, z),
        feature codes under "f", and the extra columns. ``tag`` names this
        worker's files so workers never share one."""
        n = len(cols["x"])
        if n == 0:
            return {"n": 0}
        rec = np.empty(n, self.record(False))
        for a in self.axes:
            rec[a] = cols[a]
        rec["f"] = cols["f"]
        for name, _, _ in self.columns:
            rec[name] = cols[name]
        rec["r"] = np.random.default_rng(seed).integers(0, 2**32, n, dtype=np.uint32)
        _, ty = self.tiles(rec["x"], rec["y"], 0)
        band = ty // self.grid(0)[3]
        order = np.argsort(band, kind="stable")
        rec, band = rec[order], band[order]
        d = self.tmp / "L0"
        d.mkdir(parents=True, exist_ok=True)
        cuts = np.flatnonzero(np.diff(band)) + 1
        for lo, hi in zip(np.r_[0, cuts], np.r_[cuts, n]):
            with open(d / f"band_{band[lo]:05d}.{tag}.bin", "ab") as fh:
                rec[lo:hi].tofile(fh)
        lo = [float(rec[a].min()) for a in self.axes]
        hi = [float(rec[a].max()) for a in self.axes]
        return {"n": n, "min": lo, "max": hi}

    # ---- passes B, C and metadata ---------------------------------------
    def finish(self, bounds=None) -> None:
        self._levels(zarr.open_group(self.out, mode="w"), 0, bounds)

    def resume(self, start: int, bounds=None, seed: int = 0) -> None:
        """Rebuild levels ``start`` and up from level ``start - 1`` as written,
        after a run stopped part way. The new levels are fresh uniform samples
        of that level, so they nest as the spec asks."""
        root = zarr.open_group(self.out, mode="a")
        for k in range(start, self.n_levels):
            shutil.rmtree(Path(self.out) / str(k), ignore_errors=True)
        tm = root[str(start - 1)]["tile_major"]
        names = list(tm["columns"][:])
        col = lambda c: tm[f"column_{names.index(c)}"]
        read = lambda c: (col(c)["codes"] if isinstance(col(c), zarr.Group) else col(c))[:]
        n = tm.attrs["df:shape"][0]
        rec = np.empty(n, self.record(True))
        for a in self.axes:
            rec[a] = read(a)
        rec["f"] = read(self.feature_key)
        for name, _, _ in self.columns:
            rec[name] = read(name)
        rec["id"] = read("id") if "id" in names else np.arange(n, dtype=np.uint64)
        # Level start - 1 holds r < 2^(32 - 2(start - 1)); draw r afresh in that
        # range, so r < 2^(32 - 2·start) keeps a quarter, as in a full run
        r = np.random.default_rng(seed).integers(0, 2 ** (32 - 2 * (start - 1)), n, dtype=np.uint64).astype(np.uint32)
        rec = rec[r < 2 ** (32 - 2 * start)]
        rec["r"] = r[r < 2 ** (32 - 2 * start)]
        _, ty = self.tiles(rec["x"], rec["y"], start)
        band = ty // self.grid(start)[3]
        d = self.tmp / f"L{start}"
        d.mkdir(parents=True, exist_ok=True)
        for b in np.unique(band):
            rec[band == b].tofile(d / f"band_{b:05d}.resume.bin")
        log(f"resuming at level {start}: {len(rec):,} of level {start - 1}'s {n:,} points")
        self._levels(root, start, bounds)

    def _levels(self, root, start: int, bounds=None) -> None:
        levels = []
        for k in range(self.n_levels):
            if k < start:
                n = root[str(k)]["tile_major"].attrs["df:shape"][0]
            else:
                n = self._level(root, k)
            if n == 0:
                break
            size, nx, ny, _ = self.grid(k)
            levels.append((k, n, size, nx, ny))
        n0 = levels[0][1]
        layout = []
        for k, n, size, nx, ny in levels:
            e = {
                "asset": str(k),
                "points:n_points": int(n),
                "points:grid": {"origin": list(self.origin), "tile_size": size, "shape": [ny, nx]},
            }
            if k:
                e.update(
                    derived_from=str(k - 1),
                    transform={"scale": [1.0] * self.dims},
                    resampling_method="random",
                )
            layout.append(e)
        cs = {"name": self.coordinate_system, "axes": [{"name": a, "type": "space", "unit": self.units} for a in self.axes]}
        attrs = {
            "encoding-type": "ngff:points",
            "axes": self.axes,
            "coordinateTransformations": [
                {"type": "identity", "input": {"name": "xyz"[: self.dims], "axes": cs["axes"]}, "output": cs}
            ],
            "spatialdata_attrs": {"feature_key": self.feature_key, "version": "0.3-zarr"},
            "zarr_conventions": [MULTISCALES, POINTS],
            "multiscales": {"layout": layout},
            "points:version": SPEC_VERSION,
            "points:units": self.units,
            "points:n_points": int(n0),
        }
        if bounds:
            attrs["points:extent"] = [bounds["min"], bounds["max"]]
        root.attrs.update(attrs)
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _enum_column(self, df, i, n, dtype, values):
        g = df.create_group(f"column_{i}")
        g.attrs.update({"zarr_conventions": [ENUM], "enum:ordered": False, "enum:codes": ["codes"]})
        _strings(g, "enum_values", values)
        return _array(g, "codes", n, dtype)

    def _dataframe(self, level, name, n, cols):
        """cols: [(name, dtype, enum_values|None, delta)] -> {name: array}"""
        df = level.create_group(name)
        df.attrs.update({
            "zarr_conventions": [DATAFRAME],
            "df:shape": [int(n), len(cols)],
            "points:order": name,
        })
        _strings(df, "columns", [c[0] for c in cols])
        out = {}
        for i, (cname, dtype, values, delta) in enumerate(cols):
            if values is not None:
                out[cname] = self._enum_column(df, i, n, dtype, values)
            else:
                out[cname] = _array(df, f"column_{i}", n, dtype, delta=delta)
        return df, out

    def _level(self, root, k: int) -> int:
        d = self.tmp / f"L{k}"
        files = sorted(d.glob("band_*.bin")) if d.exists() else []
        rec_t = self.record(k > 0)
        n = sum(f.stat().st_size for f in files) // rec_t.itemsize
        if n == 0:
            return 0
        size, nx, ny, _ = self.grid(k)
        V = len(self.features)
        log(f"level {k}: {n:,} points, grid {nx}×{ny} of {size:g} {self.units}")
        level = root.create_group(str(k))
        tm_cols = [(a, "f4", None, False) for a in self.axes]
        tm_cols.append((self.feature_key, self.fdtype, self.features, False))
        tm_cols += [(name, dt, vals, False) for name, dt, vals in self.columns]
        if k:
            tm_cols.append(("id", "u8", None, True))
        tdf, tarr = self._dataframe(level, "tile_major", n, tm_cols)
        tw = _Writer(tarr)
        counts = np.zeros(nx * ny * V, np.uint64)
        fdir = self.tmp / f"F{k}"
        fdir.mkdir(parents=True, exist_ok=True)
        ndir = self.tmp / f"L{k + 1}"
        ndir.mkdir(parents=True, exist_ok=True)
        keep_next = np.uint32(2 ** max(0, 32 - 2 * (k + 1))) if k + 1 < self.n_levels else np.uint32(0)
        fbits = max(1, (V - 1).bit_length())
        tbits = max(1, (nx * ny - 1).bit_length())
        rshift = max(0, tbits + fbits + 32 - 64)
        nrec = self.record(True)
        frec = self.feature_record()
        _, _, _, band_rows_next = self.grid(k + 1)
        row = 0
        bands = {}
        for f in files:
            bands.setdefault(f.name.split(".")[0], []).append(f)
        t0 = time.time()
        for bi, (bname, parts) in enumerate(sorted(bands.items())):
            rec = np.concatenate([np.fromfile(p, rec_t) for p in parts])
            for p in parts:
                p.unlink()
            tile, _ = self.tiles(rec["x"], rec["y"], k)
            if k == 0:
                # Within a (tile, feature) run by r, so any prefix of a run is
                # a uniform sample; the order here is what defines the ids
                key = (tile.astype(np.uint64) << np.uint64(fbits + 32 - rshift)) | (
                    rec["f"].astype(np.uint64) << np.uint64(32 - rshift)
                ) | (rec["r"].astype(np.uint64) >> np.uint64(rshift))
                order = np.argsort(key)
                del key
            else:
                # Coarser levels are samples already; by id within a run, so
                # the ids rise and delta-code to almost nothing
                order = np.lexsort((rec["id"], rec["f"], tile))
            rec, tile = rec[order], tile[order]
            del order
            m = len(rec)
            ids = (np.arange(row, row + m, dtype=np.uint64) if k == 0 else rec["id"])
            cols = {a: rec[a] for a in self.axes}
            cols[self.feature_key] = rec["f"]
            for name, _, _ in self.columns:
                cols[name] = rec[name]
            if k:
                cols["id"] = ids
            tw.append(cols)
            counts += np.bincount(tile * V + rec["f"], minlength=nx * ny * V).astype(np.uint64)
            # feature-major spill: rows of this level, by feature
            fo = np.argsort(rec["f"], kind="stable")
            fr = np.empty(m, frec)
            for a in self.axes:
                fr[a] = rec[a][fo]
            fr["row"] = np.arange(row, row + m, dtype=np.uint64)[fo]
            fr["id"] = ids[fo]
            fs = rec["f"][fo]
            cuts = np.flatnonzero(np.diff(fs)) + 1
            for lo, hi in zip(np.r_[0, cuts], np.r_[cuts, m]):
                with open(fdir / f"feat_{fs[lo]:06d}.bin", "ab") as fh:
                    fr[lo:hi].tofile(fh)
            del fr, fo, fs
            # next level (a sparse band may send none on)
            sel = rec["r"] < keep_next if keep_next else None
            if sel is not None and sel.any():
                nr = np.empty(int(sel.sum()), nrec)
                for name in rec_t.names:
                    if name != "id":
                        nr[name] = rec[name][sel]
                nr["id"] = ids[sel]
                _, nty = self.tiles(nr["x"], nr["y"], k + 1)
                nb = nty // band_rows_next
                o = np.argsort(nb, kind="stable")
                nr, nb = nr[o], nb[o]
                cuts = np.flatnonzero(np.diff(nb)) + 1
                for lo, hi in zip(np.r_[0, cuts], np.r_[cuts, len(nr)]):
                    with open(ndir / f"band_{nb[lo]:05d}.b{bi:05d}.bin", "ab") as fh:
                        nr[lo:hi].tofile(fh)
            row += m
            el = time.time() - t0
            log(f"  tile-major band {bi + 1}/{len(bands)}: {m:,} points, {row:,}/{n:,} ({el:.0f}s)")
        assert tw.close() == n == row
        counts = counts.reshape(nx * ny, V)
        # feature-major
        fm_cols = [(a, "f4", None, False) for a in self.axes]
        fm_cols.append((self.feature_key, self.fdtype, self.features, False))
        fm_cols.append(("id", "u8", None, True))
        fdf, farr = self._dataframe(level, "feature_major", n, fm_cols)
        fdf.attrs["points:copied_columns"] = []
        fw = _Writer(farr)
        per_feature = counts.sum(0)
        t0 = time.time()
        for fi in range(V):
            p = fdir / f"feat_{fi:06d}.bin"
            if not p.exists():
                assert per_feature[fi] == 0
                continue
            fr = np.fromfile(p, frec)
            p.unlink()
            assert len(fr) == per_feature[fi], (fi, len(fr), per_feature[fi])
            fr = fr[np.argsort(fr["row"])]
            cols = {a: fr[a] for a in self.axes}
            cols[self.feature_key] = np.full(len(fr), fi, self.fdtype)
            cols["id"] = fr["id"]
            fw.append(cols)
            if fi % 50 == 0:
                log(f"  feature-major {fi + 1}/{V} ({time.time() - t0:.0f}s)")
        assert fw.close() == n
        shutil.rmtree(fdir, ignore_errors=True)
        # index
        idx = level.create_group("index")
        tile_rows = np.r_[0, np.cumsum(counts.sum(1))].astype(np.uint64)
        feature_rows = np.r_[0, np.cumsum(per_feature)].astype(np.uint64)
        ft = counts.T
        nzf, nzt = np.nonzero(ft)
        ftc = ft[nzf, nzt]
        ptr = np.r_[0, np.cumsum(np.bincount(nzf, minlength=V))].astype(np.uint64)
        for name, data, dt in (
            ("tile_rows", tile_rows, "u8"),
            ("feature_rows", feature_rows, "u8"),
            ("feature_tile_ptr", ptr, "u8"),
            ("feature_tile", nzt.astype(np.uint32), "u4"),
            ("feature_tile_count", ftc.astype(np.uint32), "u4"),
        ):
            a = _array(idx, name, len(data), dt, shard=None)
            if len(data):
                a[:] = data
        log(f"level {k} done: {len(nzt):,} (feature, tile) pairs")
        return n


# ---- reader ---------------------------------------------------------------
def _col(df, name):
    names = list(df["columns"][:])
    i = names.index(name)
    node = df[f"column_{i}"]
    return node


def _read_col(df, name, lo=None, hi=None):
    node = _col(df, name)
    if isinstance(node, zarr.Group):
        codes = node["codes"][lo:hi]
        values = np.asarray(node["enum_values"][:], dtype=object)
        import pandas as pd

        return pd.Categorical.from_codes(codes.astype(np.int64), values)
    return node[lo:hi]


def open_element(path):
    g = zarr.open_group(path, mode="r")
    return g, dict(g.attrs)


def read_dataframe(path, level=0, order="tile_major", lo=None, hi=None, columns=None):
    """One order of one level (or a row range of it) as a pandas DataFrame."""
    import pandas as pd

    g, _ = open_element(path)
    df = g[str(level)][order]
    names = columns or list(df["columns"][:])
    return pd.DataFrame({c: _read_col(df, c, lo, hi) for c in names})


def feature_points(path, feature, level=0, bbox=None):
    """A feature's points (x, y(, z), id), optionally inside bbox=(x0, y0, x1, y1)."""
    import pandas as pd

    g, attrs = open_element(path)
    L = g[str(level)]
    fm = L["feature_major"]
    fcol = _col(fm, attrs["spatialdata_attrs"]["feature_key"])
    values = list(fcol["enum_values"][:])
    f = values.index(feature)
    idx = L["index"]
    lo, hi = (int(v) for v in idx["feature_rows"][f : f + 2])
    axes = attrs["axes"]
    if bbox is None:
        return pd.DataFrame({c: _read_col(fm, c, lo, hi) for c in axes + ["id"]})
    grid = attrs["multiscales"]["layout"][level]["points:grid"]
    (x0, y0), size, (ny, nx) = grid["origin"], grid["tile_size"], grid["shape"]
    p0, p1 = (int(v) for v in idx["feature_tile_ptr"][f : f + 2])
    tiles = idx["feature_tile"][p0:p1].astype(np.int64)
    cnt = idx["feature_tile_count"][p0:p1].astype(np.int64)
    starts = lo + np.r_[0, np.cumsum(cnt)[:-1]]
    tx0, ty0 = int((bbox[0] - x0) // size), int((bbox[1] - y0) // size)
    tx1, ty1 = int((bbox[2] - x0) // size), int((bbox[3] - y0) // size)
    tx, ty = tiles % nx, tiles // nx
    sel = (tx >= tx0) & (tx <= tx1) & (ty >= ty0) & (ty <= ty1)
    parts = []
    for s, c in zip(starts[sel], cnt[sel]):
        parts.append(pd.DataFrame({cn: _read_col(fm, cn, int(s), int(s + c)) for cn in axes + ["id"]}))
    if not parts:
        return pd.DataFrame(columns=axes + ["id"])
    out = pd.concat(parts, ignore_index=True)
    m = (out.x >= bbox[0]) & (out.x < bbox[2]) & (out.y >= bbox[1]) & (out.y < bbox[3])
    return out[m].reset_index(drop=True)


def validate(path, sample=200_000, seed=0) -> dict:
    """Checks the spec's invariants; returns a report, raises on failure."""
    g, attrs = open_element(path)
    rng = np.random.default_rng(seed)
    key = attrs["spatialdata_attrs"]["feature_key"]
    axes = attrs["axes"]
    report = {}
    prev_ids = None
    tm0 = g["0"]["tile_major"]
    for e in attrs["multiscales"]["layout"]:
        k = int(e["asset"])
        L = g[str(k)]
        n = e["points:n_points"]
        tm, fm, idx = L["tile_major"], L["feature_major"], L["index"]
        assert tm.attrs["df:shape"][0] == n == fm.attrs["df:shape"][0]
        tile_rows = idx["tile_rows"][:]
        feature_rows = idx["feature_rows"][:]
        assert tile_rows[-1] == n and feature_rows[-1] == n
        ptr = idx["feature_tile_ptr"][:]
        cnt = idx["feature_tile_count"][:]
        assert ptr[-1] == len(cnt)
        has = np.diff(ptr) > 0
        sums = np.add.reduceat(cnt.astype(np.int64), ptr[:-1][has].astype(np.int64))
        assert np.array_equal(sums, np.diff(feature_rows.astype(np.int64))[has])
        # sampled rows of feature_major point at identical points in level-0 tile_major
        rows = np.sort(rng.choice(n, min(sample, n), replace=False))
        fid = _col(fm, "id")
        ids = fid.get_coordinate_selection(rows)
        same = True
        for a in axes:
            va = _col(fm, a).get_coordinate_selection(rows)
            vb = _col(tm0, a).get_coordinate_selection(ids)
            same &= bool(np.array_equal(va, vb))
        fa = _col(fm, key)["codes"].get_coordinate_selection(rows)
        fb = _col(tm0, key)["codes"].get_coordinate_selection(ids)
        same &= bool(np.array_equal(fa, fb))
        assert same, f"level {k}: feature_major rows don't match level-0 tile_major"
        # nesting: this level's ids are all in the level below (checked when
        # the level below is small enough to hold its ids in memory)
        nested = None
        if prev_ids is not None:
            pos = np.clip(np.searchsorted(prev_ids, ids), 0, len(prev_ids) - 1)
            nested = bool(np.all(prev_ids[pos] == ids))
            assert nested, f"level {k} is not contained in level {k - 1}"
        prev_ids = np.sort(_col(tm, "id")[:]) if (k and n < 300_000_000) else None
        report[k] = {"n": n, "checked_rows": len(rows), "tiles": len(tile_rows) - 1,
                     "pairs": len(cnt), "nested_checked": nested}
    return report


if __name__ == "__main__":
    cmd, path = sys.argv[1], sys.argv[2]
    if cmd == "validate":
        print(json.dumps(validate(path), indent=1))
    elif cmd == "info":
        _, attrs = open_element(path)
        print(json.dumps(attrs, indent=1))
