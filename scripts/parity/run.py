#!/usr/bin/env python3
"""Pipeline parity harness: run BOTH pipelines on the same inputs and diff them.

    python3 scripts/parity/run.py [--work .parity] [--only CASE ...] [--skip-generate]

1. generate-inputs.py   → <work>/inputs (+ cases.json)
2. Python pipeline      → <work>/py/<case>   (scripts/process_spatial_data.py /
                          process_single_molecule.py, invoked like worker/entrypoint.py)
3. browser pipeline     → <work>/js/<case>   (scripts/parity/run-js.mts via tsx)
4. diff.py per case; summary; exit 1 if any case differs or fails.
See docs/PARITY-TESTING.md and docs/PIPELINE-SPEC.md.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.stdout.reconfigure(line_buffering=True)  # keep our lines in order with subprocess output


def run(cmd, log: Path, cwd=ROOT) -> int:
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "w") as f:
        return subprocess.run(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT).returncode


def python_cmd(case: dict, inp: Path, out: Path) -> list:
    if case["kind"] == "sc":
        return [sys.executable, str(ROOT / "scripts/process_spatial_data.py"), str(inp), str(out)]
    cmd = [sys.executable, str(ROOT / "scripts/process_single_molecule.py"), str(inp), str(out),
           "--dataset-type", case["datasetType"]]
    if case.get("cellIdCol"):
        cmd += ["--cell-id-col", case["cellIdCol"]]
    return cmd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", type=Path, default=ROOT / ".parity")
    ap.add_argument("--only", nargs="*", default=[])
    ap.add_argument("--skip-generate", action="store_true")
    args = ap.parse_args()
    work = args.work.resolve()
    inputs, py_out, js_out, logs = work / "inputs", work / "py", work / "js", work / "logs"

    if not args.skip_generate or not (inputs / "cases.json").exists():
        shutil.rmtree(inputs, ignore_errors=True)
        print("== generating inputs")
        subprocess.run([sys.executable, str(HERE / "generate-inputs.py"), str(inputs)], check=True)
    cases = json.loads((inputs / "cases.json").read_text())
    if args.only:
        cases = [c for c in cases if c["name"] in args.only]

    print("== python pipeline")
    py_fail = set()
    for c in cases:
        out = py_out / c["name"]
        shutil.rmtree(out, ignore_errors=True)
        t0 = time.time()
        rc = run(python_cmd(c, inputs / c["input"], out), logs / f"py-{c['name']}.log")
        status = "ok" if rc == 0 else f"FAILED rc={rc} (see {logs / f'py-{c["name"]}.log'})"
        print(f"  py  {c['name']:<24} {status}  {time.time() - t0:.1f}s")
        if rc:
            py_fail.add(c["name"])

    print("== browser pipeline (node)")
    # vite-node (vitest's runner) rather than tsx: the lib modules are ESM-only
    # consumers (h5wasm, hyparquet) and vite applies the same "@/" alias the unit
    # tests use (vitest.config.ts).
    js_rc = subprocess.run(
        ["npx", "vite-node", "-c", "vitest.config.ts", str(HERE / "run-js.mts"),
         str(inputs), str(js_out), *[c["name"] for c in cases]],
        cwd=ROOT, env={**os.environ, "NODE_NO_WARNINGS": "1"},
    ).returncode

    print("== diff (python vs js)")
    results = {}
    for c in cases:
        a, b = py_out / c["name"], js_out / c["name"]
        if c["name"] in py_fail or not b.exists():
            results[c["name"]] = "PIPELINE FAILED"
            print(f"  {c['name']:<24} PIPELINE FAILED")
            continue
        print(f"  {c['name']}:")
        rc = subprocess.run([sys.executable, str(HERE / "diff.py"), str(a), str(b), "--kind", c["kind"]], cwd=ROOT).returncode
        results[c["name"]] = "PASS" if rc == 0 else "DIFF"

    print("== summary")
    for name, r in results.items():
        print(f"  {r:<16} {name}")
    bad = [n for n, r in results.items() if r != "PASS"]
    print(f"  {len(results) - len(bad)}/{len(results)} cases identical" + (f"; js runner rc={js_rc}" if js_rc else ""))
    sys.exit(1 if bad or js_rc else 0)


if __name__ == "__main__":
    main()
