# Phase 3 — AWS provisioning checklist

← [Back to design doc](README.md) · [Phase 3 doc](phase-3-trigger-callback-progress.md)

Hand-provisioning steps for the v1 spike (README §9: get it working by hand, convert
to CDK afterwards). Everything below is **your** side — the app + worker code is done
and verified locally; this is what makes `SubmitJob` → Fargate → callback real.

**Prereqs:** Docker running, and the AWS CLI configured with an **admin profile**.

The existing profiles (`default`, `staging`, `production`, `old-prod`) are all narrow
S3-only service users — none can create IAM roles, ECR repos, or Batch resources.
Create a dedicated admin user (IAM → Users → Create user → attach `AdministratorAccess`
→ Security credentials → Create access key → CLI), then:

```bash
aws configure --profile mfe-admin     # region us-west-2, output json
aws sts get-caller-identity --profile mfe-admin   # expect .../user/<your-admin-user>
export AWS_PROFILE=mfe-admin          # every command below uses this
```

> Enable MFA on that user, and **delete the access key once provisioning is done** —
> it's a long-lived admin credential on your laptop. If the account has IAM Identity
> Center, prefer `aws configure sso` (short-lived credentials) instead.

---

## 0. Decisions (fill these in first)

| Variable | Value | Notes |
|---|---|---|
| `REGION` | `us-west-2` | **Must match the S3 bucket** (`merfisheyes-staging` is us-west-2) or you pay cross-region transfer on every job. |
| `ACCOUNT_ID` | `533267148861` | From `sts get-caller-identity`. |
| `BUCKET` | `merfisheyes-staging` | Existing. |
| `ECR_REPO` | `merfisheyes-worker` | |
| Default tier | 4 vCPU / 16 GB | README §5.6 `standard`. |
| Ephemeral storage | 100 GiB | Must hold raw + output. Range 21–200 GiB (default 20). |

```bash
export REGION=us-west-2
export ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
export BUCKET=merfisheyes-staging
export ECR_REPO=merfisheyes-worker
export ECR_URI=$ACCOUNT_ID.dkr.ecr.$REGION.amazonaws.com/$ECR_REPO
```

---

## 1. ECR — repository + push the image

- [ ] Create the repo:
```bash
aws ecr create-repository --repository-name $ECR_REPO --region $REGION
```
- [ ] Log in, build, push:
```bash
aws ecr get-login-password --region $REGION \
  | docker login --username AWS --password-stdin $ACCOUNT_ID.dkr.ecr.$REGION.amazonaws.com

# ⚠️ --platform linux/amd64 is REQUIRED. Fargate runs x86_64; a default build on
# Apple Silicon produces arm64 and the job will fail to start.
docker build --platform linux/amd64 -f worker/Dockerfile -t $ECR_REPO .
docker tag $ECR_REPO:latest $ECR_URI:latest
docker push $ECR_URI:latest
```

> The amd64 build is slower locally (emulated) but only affects build time, not the
> Fargate run.

---

## 2. IAM — three roles/permissions

### 2a. Batch **execution** role (pull image, write logs)
- [ ] Create with trust policy for `ecs-tasks.amazonaws.com`, attach the managed
      `AmazonECSTaskExecutionRolePolicy`:
```bash
cat > /tmp/trust-ecs.json <<'JSON'
{"Version":"2012-10-17","Statement":[{"Effect":"Allow",
 "Principal":{"Service":"ecs-tasks.amazonaws.com"},"Action":"sts:AssumeRole"}]}
JSON

aws iam create-role --role-name merfisheyes-batch-execution \
  --assume-role-policy-document file:///tmp/trust-ecs.json
aws iam attach-role-policy --role-name merfisheyes-batch-execution \
  --policy-arn arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy
```

### 2b. Batch **task** role (what the worker itself may do)
The worker needs S3 on the bucket only. **No AWS keys are baked in** — `boto3` picks
up the task role automatically, so drop the `AWS_ACCESS_KEY_ID`/`SECRET` env vars we
used for local runs.
- [ ] Create the role (same trust policy) and attach:
```bash
cat > /tmp/worker-s3.json <<JSON
{"Version":"2012-10-17","Statement":[
 {"Effect":"Allow","Action":["s3:GetObject","s3:PutObject","s3:DeleteObject","s3:AbortMultipartUpload"],
  "Resource":"arn:aws:s3:::$BUCKET/*"},
 {"Effect":"Allow","Action":["s3:ListBucket","s3:ListBucketMultipartUploads"],
  "Resource":"arn:aws:s3:::$BUCKET"}]}
JSON

aws iam create-role --role-name merfisheyes-batch-task \
  --assume-role-policy-document file:///tmp/trust-ecs.json
aws iam put-role-policy --role-name merfisheyes-batch-task \
  --policy-name worker-s3 --policy-document file:///tmp/worker-s3.json
```

### 2c. App credentials (the Vercel keys) — allow submitting jobs
- [ ] Add to the existing `merfisheyes-staging-vercel` user:
```json
{"Version":"2012-10-17","Statement":[
 {"Effect":"Allow","Action":["batch:SubmitJob","batch:DescribeJobs","batch:TerminateJob"],
  "Resource":"*"}]}
```
- [ ] **Also add `s3:ListBucketMultipartUploads`** to that user (on
      `arn:aws:s3:::merfisheyes-staging`). This is the gap found in Phase 1 — without
      it `/api/ingest/[id]/abort` can't enumerate orphaned multipart uploads. It
      degrades gracefully today, but you want it for immediate cleanup.

---

## 3. Networking

Fargate needs **outbound internet** for S3 and the HTTPS callback.

- [ ] Simplest for a spike: use the **default VPC's public subnets** and set
      `assignPublicIp: ENABLED` (done in the job definition below). No NAT gateway,
      **no ~$32/month NAT cost**.
- [ ] Get subnet + SG ids:
```bash
VPC=$(aws ec2 describe-vpcs --filters Name=isDefault,Values=true \
      --query 'Vpcs[0].VpcId' --output text --region $REGION)
aws ec2 describe-subnets --filters Name=vpc-id,Values=$VPC \
  --query 'Subnets[].SubnetId' --output text --region $REGION
aws ec2 describe-security-groups --filters Name=vpc-id,Values=$VPC Name=group-name,Values=default \
  --query 'SecurityGroups[0].GroupId' --output text --region $REGION
```
- [ ] The default SG allows all outbound — that's all we need (no inbound).

> Production alternative: private subnets + NAT gateway, or private subnets + an **S3
> gateway VPC endpoint** (free) — but the callback still needs egress, so you'd need
> NAT regardless unless you front the callback with a VPC endpoint service.

---

## 4. AWS Batch — compute environment, queue, job definition

- [ ] **Compute environment** (Fargate):
```bash
aws batch create-compute-environment --region $REGION \
  --compute-environment-name merfisheyes-fargate \
  --type MANAGED --state ENABLED \
  --compute-resources "type=FARGATE,maxvCpus=64,subnets=[SUBNET_IDS],securityGroupIds=[SG_ID]"
```
- [ ] **Job queue**:
```bash
aws batch create-job-queue --region $REGION \
  --job-queue-name merfisheyes-ingest \
  --state ENABLED --priority 1 \
  --compute-environment-order order=1,computeEnvironment=merfisheyes-fargate
```
- [ ] **Job definition** — save as `/tmp/jobdef.json` (substitute `$ECR_URI`,
      `$ACCOUNT_ID`, `$BUCKET`) and register:
```json
{
  "jobDefinitionName": "merfisheyes-ingest",
  "type": "container",
  "platformCapabilities": ["FARGATE"],
  "containerProperties": {
    "image": "ECR_URI:latest",
    "executionRoleArn": "arn:aws:iam::ACCOUNT_ID:role/merfisheyes-batch-execution",
    "jobRoleArn": "arn:aws:iam::ACCOUNT_ID:role/merfisheyes-batch-task",
    "resourceRequirements": [
      { "type": "VCPU", "value": "4" },
      { "type": "MEMORY", "value": "16384" }
    ],
    "ephemeralStorage": { "sizeInGiB": 100 },
    "networkConfiguration": { "assignPublicIp": "ENABLED" },
    "fargatePlatformConfiguration": { "platformVersion": "LATEST" },
    "environment": [
      { "name": "AWS_S3_BUCKET", "value": "BUCKET" },
      { "name": "AWS_REGION", "value": "us-west-2" },
      { "name": "UPLOAD_CONCURRENCY", "value": "16" },
      { "name": "DELETE_RAW", "value": "true" }
    ]
  },
  "retryStrategy": { "attempts": 2 }
}
```
```bash
aws batch register-job-definition --region $REGION --cli-input-json file:///tmp/jobdef.json
```

**Notes**
- Valid Fargate vCPU↔RAM pairs: 4 vCPU → 8–30 GB, 8 → 16–60 GB, 16 → 32–120 GB.
  Tiers must snap to these.
- `DELETE_RAW=true` here is the real flow (raw deleted after success). We ran locally
  with it off to preserve the test file.
- Per-job values (`DATASET_ID`, `DATASET_KIND`, `PROCESSING_PARAMS`, `CALLBACK_URL`,
  `CALLBACK_SECRET`) are supplied at `SubmitJob` time via container overrides — the
  app does that; don't hardcode them here.
- `retryStrategy.attempts: 2` — note the callback's `COMPLETE` is idempotent, so a
  retry re-registering output is safe.

---

## 5. Secrets

- [ ] `CALLBACK_SECRET` — a 64-char hex value already exists in your local `.env`.
      Set the **same value** in:
  - [ ] Vercel env (so the app can verify signatures)
  - [ ] passed per job by the app at `SubmitJob` (it reads it from its own env)
- [ ] Prefer **AWS Secrets Manager** later and reference it via the job definition's
      `secrets` block rather than a plaintext env var.

### ⚠️ Vercel Deployment Protection will block the callback

The codebase already notes that server-to-server calls get blocked by Vercel
deployment protection. Fargate POSTing to `https://dev.merfisheyes.com/api/ingest/...`
hits the same wall. Either:
- point `CALLBACK_URL` at an **unprotected production domain**, or
- disable protection for that path, or
- have the worker send Vercel's bypass header
  (`x-vercel-protection-bypass: $VERCEL_AUTOMATION_BYPASS_SECRET`).

Decide this before the first real job, or the job will succeed while the app never
learns and the dataset stays stuck at `PROCESSING`.

---

## 6. S3 lifecycle rule on `raw/`

Backstop for cancelled/abandoned uploads (design §8.6, §10).
- [ ] Apply:
```bash
cat > /tmp/lifecycle.json <<'JSON'
{"Rules":[
 {"ID":"expire-raw-ingest","Status":"Enabled","Filter":{"Prefix":"raw/"},
  "Expiration":{"Days":7},
  "AbortIncompleteMultipartUpload":{"DaysAfterInitiation":1}}]}
JSON
aws s3api put-bucket-lifecycle-configuration --bucket $BUCKET \
  --lifecycle-configuration file:///tmp/lifecycle.json
```
> `AbortIncompleteMultipartUpload` is the important half — it reclaims orphaned parts
> from cancelled browser uploads.

---

## 7. Verify before wiring the app

- [ ] Submit a job by hand against a `QUEUED` dataset:
```bash
aws batch submit-job --region $REGION \
  --job-name test-ingest --job-queue merfisheyes-ingest \
  --job-definition merfisheyes-ingest \
  --container-overrides '{"environment":[
     {"name":"DATASET_ID","value":"ds_xxxxx"},
     {"name":"DATASET_KIND","value":"single_cell"},
     {"name":"PROCESSING_PARAMS","value":"{\"kind\":\"single_cell\",\"stages\":{\"chunk\":{\"chunkSize\":1}}}"},
     {"name":"CALLBACK_URL","value":"https://YOUR_APP/api/ingest/ds_xxxxx/callback"},
     {"name":"CALLBACK_SECRET","value":"THE_SECRET"}]}'
```
- [ ] Watch it: `aws batch describe-jobs --jobs <jobId> --region $REGION`
- [ ] Logs land in CloudWatch under `/aws/batch/job`.
- [ ] Success criteria: status goes `PROCESSING` → `COMPLETE` in the DB, output
      appears under `datasets/{id}/`, `raw/{id}/` is deleted, the email arrives, and
      `/viewer/{id}` renders.

---

## 8. Cost sanity

- Fargate scales to zero; a 4 vCPU/16 GB job for ~2 min ≈ **$0.01–0.02**. Our real
  415 MB run took ~90 s locally.
- ECR storage: cents/GB/month. CloudWatch logs: negligible at this volume.
- **No NAT gateway** if you use public subnets + `assignPublicIp` (saves ~$32/mo).
- ~500 jobs/month ≈ **$10** compute (README §12).

---

## Gotchas recap

1. **Build `--platform linux/amd64`** or the job won't start.
2. **Vercel deployment protection** will silently break the callback.
3. **Region must match the bucket** (us-west-2).
4. **Drop static AWS keys** from the job env — the task role supplies them.
5. **Fargate CPU/RAM must be a valid pair**; tiers snap to allowed combos.
