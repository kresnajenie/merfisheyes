#!/usr/bin/env python3
"""Compare two processor output trees (Python vs browser JS) per docs/PIPELINE-SPEC.md.

    python3 scripts/parity/diff.py <python-out> <js-out> [--kind sc|sm]

Decompresses *.gz, compares JSON structurally after normalising the documented
expected differences (manifest name / created_at / dataset_id / version /
processing.*, JS-only manifest keys, JS per-gene stats in expr/index.json) and
binaries byte-for-byte (expr chunks: per gene block, ignoring the gene table's
uncompressedSize). Prints one line per differing file with the first JSON path /
byte offset, then a summary; exit 1 on any difference.
"""
import argparse
import gzip
import json
import struct
import sys
from pathlib import Path

# ---- normalisation ----------------------------------------------------------

def _norm_sc_manifest(m: dict) -> dict:
    files = m.get("files", {})
    stats = dict(m.get("statistics", {}))
    stats["available_embeddings"] = sorted(stats.get("available_embeddings", []))
    stats.pop("cluster_count", None)  # Tier-2: the two pipelines count different things
    return {
        "normalized": m.get("normalized"),
        "type": m.get("type"),
        "statistics": stats,
        "files": {
            "coordinates": sorted(files.get("coordinates", [])),
            "expression_chunks": files.get("expression_chunks"),
            "observation_columns": sorted(files.get("observation_columns", [])),
            "de_stats": sorted(files.get("de_stats", [])),
        },
    }


def _norm_sm_manifest(m: dict) -> dict:
    genes = m.get("genes", {})
    counts = genes.get("molecule_counts")
    return {
        "type": m.get("type"),
        "has_unassigned": m.get("has_unassigned"),
        "statistics": m.get("statistics"),
        "genes": {
            "unique_gene_names": genes.get("unique_gene_names"),
            "molecule_counts": counts,
        },
        "processing": {
            "coordinate_range": (m.get("processing") or {}).get("coordinate_range"),
            "scaling_factor": (m.get("processing") or {}).get("scaling_factor"),
        },
    }


def _norm_expr_index(ix: dict) -> dict:
    return {
        "total_genes": ix.get("total_genes"),
        "num_chunks": ix.get("num_chunks"),
        "chunk_size": ix.get("chunk_size"),
        "genes": [
            {"name": g.get("name"), "chunk_id": g.get("chunk_id"), "position_in_chunk": g.get("position_in_chunk")}
            for g in ix.get("genes", [])
        ],
    }


def normalise_json(rel: str, data, kind: str):
    if rel == "manifest.json" or rel == "manifest.json.gz":
        return _norm_sc_manifest(data) if kind == "sc" else _norm_sm_manifest(data)
    if rel == "expr/index.json":
        return _norm_expr_index(data)
    return data


# ---- binary-aware compares ----------------------------------------------------

def _u32s(b: bytes, off: int, n: int):
    return struct.unpack_from(f"<{n}I", b, off)


def compare_expr_chunk(a: bytes, b: bytes):
    """Expr chunk: header [u32 1][u32 G][u32 chunkId][u32 cells], G × 24-byte
    gene table, then per-gene sparse blocks. The gene table's uncompressedSize
    field is a known JS/Python difference; everything else must match."""
    if len(a) < 16 or len(b) < 16:
        return "truncated header" if a != b else None
    ha, hb = _u32s(a, 0, 4), _u32s(b, 0, 4)
    if ha != hb:
        return f"header {ha} != {hb}"
    g = ha[1]
    ta, tb = a[16:16 + 24 * g], b[16:16 + 24 * g]
    if ta != tb:
        # Entry = (geneIdx, dataOffset, compressedSize, uncompressedSize, nnz, 0).
        # Field 3 is the one documented JS/Python difference; tolerate only it.
        for i in range(g):
            ea = struct.unpack_from("<6I", ta, 24 * i)
            eb = struct.unpack_from("<6I", tb, 24 * i)
            if ea[:3] + ea[4:] != eb[:3] + eb[4:]:
                return f"gene table entry {i}: {ea} != {eb}"
    pa, pb = a[16 + 24 * g:], b[16 + 24 * g:]
    if pa != pb:
        return f"gene data differs at byte {first_diff(pa, pb) + 16 + 24 * g}"
    return None


def compare_de(a: bytes, b: bytes):
    """de/*.bin: [u32 1][u32 G][u32 C][u32 0], C × (u32 len + utf-8 name),
    u32 counts[C], f32 means[G·C], f32 pct[G·C]. Report the first differing
    section with its index and values (float32 compared bit-exactly)."""
    if a[:16] != b[:16]:
        return f"header {_u32s(a, 0, 4)} != {_u32s(b, 0, 4)}"
    _, g, c, _ = _u32s(a, 0, 4)

    def parse(buf):
        off = 16
        names = []
        for _ in range(c):
            (ln,) = struct.unpack_from("<I", buf, off)
            off += 4
            names.append(buf[off:off + ln].decode("utf-8", "replace"))
            off += ln
        counts = struct.unpack_from(f"<{c}I", buf, off)
        off += 4 * c
        means = struct.unpack_from(f"<{g * c}f", buf, off)
        off += 4 * g * c
        pct = struct.unpack_from(f"<{g * c}f", buf, off)
        return names, counts, means, pct

    try:
        na, ca, ma, pa = parse(a)
        nb, cb, mb, pb = parse(b)
    except struct.error as e:
        return f"unparseable ({e}); bytes differ at offset {first_diff(a, b)}"
    if na != nb:
        return f"celltype names differ: {na[:5]} vs {nb[:5]}"
    if ca != cb:
        i = next(i for i in range(c) if ca[i] != cb[i])
        return f"counts[{na[i]!r}]: {ca[i]} vs {cb[i]}"
    for label, xa, xb in (("means", ma, mb), ("pct", pa, pb)):
        for i in range(g * c):
            if xa[i] != xb[i] and not (xa[i] != xa[i] and xb[i] != xb[i]):
                return f"{label}[gene {i // c}, {na[i % c]!r}]: {xa[i]!r} vs {xb[i]!r}"
    return None


def first_diff(a: bytes, b: bytes) -> int:
    n = min(len(a), len(b))
    for i in range(n):
        if a[i] != b[i]:
            return i
    return n


def json_first_diff(a, b, path="$"):
    if type(a) != type(b):
        return f"{path}: {type(a).__name__} vs {type(b).__name__} ({str(a)[:60]!r} vs {str(b)[:60]!r})"
    if isinstance(a, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a:
                return f"{path}.{k}: only in JS"
            if k not in b:
                return f"{path}.{k}: only in Python"
            d = json_first_diff(a[k], b[k], f"{path}.{k}")
            if d:
                return d
        return None
    if isinstance(a, list):
        if len(a) != len(b):
            return f"{path}: length {len(a)} vs {len(b)}"
        for i, (x, y) in enumerate(zip(a, b)):
            d = json_first_diff(x, y, f"{path}[{i}]")
            if d:
                return d
        return None
    if isinstance(a, float) and isinstance(b, float):
        if a != b and not (a != a and b != b):  # NaN == NaN
            return f"{path}: {a!r} vs {b!r}"
        return None
    if a != b:
        return f"{path}: {a!r} vs {b!r}"
    return None


def read(p: Path) -> bytes:
    return gzip.open(p, "rb").read() if p.suffix == ".gz" else p.read_bytes()


def compare_trees(py: Path, js: Path, kind: str):
    fa = {str(p.relative_to(py)).replace("\\", "/") for p in py.rglob("*") if p.is_file()}
    fb = {str(p.relative_to(js)).replace("\\", "/") for p in js.rglob("*") if p.is_file()}
    problems = []
    for rel in sorted(fa - fb):
        problems.append((rel, "only in Python output"))
    for rel in sorted(fb - fa):
        problems.append((rel, "only in JS output"))
    for rel in sorted(fa & fb):
        a, b = read(py / rel), read(js / rel)
        name = rel[:-3] if rel.endswith(".gz") else rel
        if name.endswith(".json"):
            try:
                ja, jb = json.loads(a), json.loads(b)
            except Exception as e:  # noqa: BLE001
                problems.append((rel, f"invalid JSON: {e}"))
                continue
            d = json_first_diff(normalise_json(rel, ja, kind), normalise_json(rel, jb, kind))
            if d:
                problems.append((rel, d))
        elif rel.startswith("expr/chunk_"):
            d = compare_expr_chunk(a, b)
            if d:
                problems.append((rel, d))
        elif rel.startswith("de/") and a != b:
            problems.append((rel, compare_de(a, b) or "bytes differ"))
        elif a != b:
            problems.append((rel, f"bytes differ at offset {first_diff(a, b)} (sizes {len(a)} vs {len(b)})"))
    return sorted(fa | fb), problems


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("python_out", type=Path)
    ap.add_argument("js_out", type=Path)
    ap.add_argument("--kind", choices=["sc", "sm"], default="sc")
    args = ap.parse_args()
    files, problems = compare_trees(args.python_out, args.js_out, args.kind)
    for rel, why in problems:
        print(f"  FAIL {rel}: {why}")
    print(f"  {len(files) - len(problems)}/{len(files)} files identical; {len(problems)} differ")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
