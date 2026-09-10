# Registering the spiralia datasets into the admin project

Run this to make the 45 spiralia embryos appear in MERFISHeyes as one browsable
project. It writes database rows only — the chunked data is already in S3.

## What it does

`register_datasets.ts` reads `build_report.csv` (committed here, 45 rows, all
`ok`) and writes:

| What | How many | Detail |
| --- | --- | --- |
| `Dataset` | 45 | `adminOwned: true`, `datasetType: "labelled_single_molecule"`, owned by the first admin user |
| `Project` | 1 | "Spiralia embryo atlas", created only if absent |
| `ProjectDataset` | 45 | `sortOrder` = developmental stage order, which is what the viewer's stage rail scrubs through |

Each dataset's `s3BaseUrl` is
`https://merfisheyes-bil.s3.us-west-2.amazonaws.com/yiqun-spiralia/{embryo}_lm`.
That bucket serves every environment, so nothing needs uploading.

**It is idempotent.** `s3BaseUrl` is unique, so a re-run updates rather than
duplicates. Safe to run twice.

## Prerequisites

1. This branch, which contains the portability fixes that let the script run
   off a plain clone:

   ```bash
   git checkout feat/lm-viewer-followups && git pull
   ```

2. `DATABASE_URL` pointing at **production** Supabase, and `npx prisma generate`
   having been run for that schema.

3. A user with role `ADMIN` or `SUPER_ADMIN` in that database — they become the
   owner of all 45 datasets and of the project.

No AWS credentials are needed. The script never touches S3.

## Steps

### 1. Dry run first

```bash
npx tsx scripts/spiralia/register_datasets.ts --dry-run
```

It prints the owner, the database host, and the first five rows, and writes
nothing. **Check all three before continuing:**

- **Owner email and role** — this account will own all 45 datasets. If it is not
  the intended one, stop; the script takes the oldest admin account and has no
  flag to override it.
- **Database host** — confirm it is production and not staging. The two are
  separate Supabase projects and the host is the only thing distinguishing them
  in this output.
- **Count is 45.**

### 2. Register

```bash
npx tsx scripts/spiralia/register_datasets.ts
```

Prints the number registered and the project id.

### 3. Verify

Open any embryo — it should load, and the stage rail should appear along the
bottom with one tile per stage:

```
/lm-viewer/from-s3?url=https://merfisheyes-bil.s3.us-west-2.amazonaws.com/yiqun-spiralia/MER6-2_E3_1_lm
```

The rail only appears when the dataset resolves to a project, so seeing it is
also proof the `ProjectDataset` rows landed.

## Things that will bite

- **`role` is not `isAdmin`.** The database stores `ADMIN` / `SUPER_ADMIN` on
  `User.role`. The script accepts either; an equality check on `"ADMIN"` alone
  silently fails on a database whose only administrator is `SUPER_ADMIN`.
- **Public visibility is all-or-nothing.** `/api/projects/[id]/public` serves a
  project only when *every* dataset in it is `adminOwned`. Adding one personal
  dataset to this project makes the whole rail 404 for signed-out visitors.
- **Wrong-database runs are not auto-reversible.** There is no undo command. To
  back one out, delete the `ProjectDataset` rows, then the `Dataset` rows whose
  `s3BaseUrl` starts with the `yiqun-spiralia` prefix, then the project.
- **Stage order is a hardcoded list.** `STAGE_ORDER` in the script defines the
  rail's left-to-right order. Any stage not in it sorts last; the script warns
  when that happens rather than failing.

## If the datasets already exist

Re-running is still the right move — it refreshes titles, descriptions, counts
and `sortOrder`, and adds any missing `ProjectDataset` rows. Nothing is
duplicated.
