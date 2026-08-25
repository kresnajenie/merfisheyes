#!/usr/bin/env python3
"""Attach real infrastructure costs to server bench results.

For each server-*.json with a datasetRowId:
  - every Batch job named ingest-<rowId> (first attempt, retries, escalations):
    billed runtime = sum of attempt started->stopped, priced at Fargate
    us-west-2 Linux/x86 rates by the job's actual vCPU/memory.
  - the final chunked output under s3://<bucket>/datasets/<rowId>/ priced at
    S3 Standard per GB-month.

Writes a `cost` object into each result file:
  { runs: [{tier, vcpu, memGB, runtimeMs, usd}], totalUsd,
    outputBytes, outputUsdMonth }

    python3 scripts/bench/harvest-costs.py [--env linux-x64]
"""
import json, glob, subprocess, sys, os

VCPU_HR = 0.04048     # Fargate us-west-2, Linux/x86, per vCPU-hour
GB_HR = 0.004445      # per GB-hour
S3_GB_MONTH = 0.023   # S3 Standard, first 50 TB
REGION = os.environ.get("AWS_REGION", "us-west-2")
ENV = sys.argv[sys.argv.index("--env") + 1] if "--env" in sys.argv else "linux-x64"
BUCKET = None
for line in open(".env"):
    if line.startswith("AWS_S3_BUCKET="):
        BUCKET = line.strip().split("=", 1)[1]
assert BUCKET, "no AWS_S3_BUCKET in .env"

def aws(argv):
    r = subprocess.run(["aws"] + argv + ["--region", REGION, "--output", "json"],
                       capture_output=True, text=True)
    return json.loads(r.stdout) if r.returncode == 0 and r.stdout.strip() else None

def jobs_for(row_id):
    # NOTE: --job-status cannot be combined with --filters (API error); a
    # filters-only call returns matching jobs across every status.
    page = aws(["batch", "list-jobs", "--job-queue", "merfisheyes-ingest",
                "--filters", f"name=JOB_NAME,values=ingest-{row_id}"])
    ids = [j["jobId"] for j in (page or {}).get("jobSummaryList", [])]
    if not ids:
        return []
    detail = aws(["batch", "describe-jobs", "--jobs"] + ids) or {}
    runs = []
    for j in sorted(detail.get("jobs", []), key=lambda x: x.get("createdAt", 0)):
        res = {r["type"]: r["value"] for r in j["container"].get("resourceRequirements", [])}
        vcpu, mem_gb = float(res.get("VCPU", 0)), int(res.get("MEMORY", 0)) / 1024
        ms = 0
        for a in j.get("attempts", []):
            if a.get("startedAt") and a.get("stoppedAt"):
                ms += a["stoppedAt"] - a["startedAt"]
        if not ms and j.get("startedAt") and j.get("stoppedAt"):
            ms = j["stoppedAt"] - j["startedAt"]
        hours = ms / 3_600_000
        runs.append({"tier": f"{int(mem_gb)}GB", "vcpu": vcpu, "memGB": mem_gb,
                     "runtimeMs": ms, "status": j["status"],
                     "usd": round(hours * (vcpu * VCPU_HR + mem_gb * GB_HR), 4)})
    return runs

def output_bytes(row_id):
    r = subprocess.run(["aws", "s3", "ls", f"s3://{BUCKET}/datasets/{row_id}/",
                        "--recursive", "--summarize"], capture_output=True, text=True)
    for line in r.stdout.splitlines():
        if "Total Size:" in line:
            return int(line.split(":")[1].strip())
    return 0

for f in sorted(glob.glob(f"bench/results/{ENV}/server-*.json")):
    row = json.load(open(f))
    rid = row.get("datasetRowId")
    if not rid:
        continue
    runs = jobs_for(rid)
    ob = output_bytes(rid)
    row["cost"] = {
        "runs": runs,
        "totalUsd": round(sum(r["usd"] for r in runs), 4),
        "outputBytes": ob,
        "outputUsdMonth": round(ob / 1024**3 * S3_GB_MONTH, 4),
    }
    json.dump(row, open(f, "w"), indent=1)
    print(f"{rid:<16} runs={len(runs)} ${row['cost']['totalUsd']:<8} "
          f"out={ob/1024**2:,.0f} MB (${row['cost']['outputUsdMonth']}/mo)  {row['datasetId'][:50]}")
