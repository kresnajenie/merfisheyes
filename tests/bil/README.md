# tests/bil — the cluster pipeline, tested on this machine

The BIL pipeline runs on the PSC cluster under SLURM, on datasets too large
to keep here. This space lets the parts of it that **are** testable be tested
without a cluster and without the data: a synthetic multi-slice dataset that
runs in seconds, plus a runner that drives the same stages over a real
dataset when one is on disk.

Read [docs/plans/BIL-CLUSTER.md](../../docs/plans/BIL-CLUSTER.md) first for
what the pipeline is and why it is changing.

## The two halves

| | What it needs | What it proves |
|---|---|---|
| `test_combine_slices.py`, `test_slices_to_points.py` | nothing — fixtures are generated | the stages, their contracts, and the regressions below |
| `run_dry_run.sh` | a real dataset on disk | that it works at size, and how long it takes |

```bash
python -m pytest tests/bil -q              # ~30 s, no data needed
```

Both test files share one `dataset` and one `combined` fixture from
`conftest.py`, so the combine runs once per session.

## What the fixture reproduces

`conftest.py` builds a sample folder shaped like the real BIL ones, because
every one of these details has cost a debugging session:

- **non-canonical file names** — `cellpose_metadata.csv`,
  `cellpose-cell-by-gene.csv`, `cellpose-detected_transcripts.csv`, not the
  names the MERSCOPE instrument writes. The fuzzy matching in
  `combine_slices_v3.py` is what finds them.
- **slices several directories down** — `<slice>/merfish_output/cellpose_cyto2_nuclei/`.
- **slices that disagree** — one names its id column `id`, the other leaves
  it unnamed.
- **21-digit cell ids**, too wide for an int64, so they are text everywhere.
- **numeric sample ids**, which pandas reads as integers unless told not to.
- **`cell_id = -1`** on a molecule in no cell, as MERSCOPE writes it.

## The regressions these tests hold down

- **A rerun must not read its own output.** Combine walks the folder it is
  given looking for slices. Run it twice with the output inside that folder —
  which is what a downloaded copy of a finished sample looks like — and it
  used to find `combined_output/` and combine it in as a fourth slice,
  **doubling every cell**, exit 0, no warning. Measured on the real
  `ace-dip-tug`: 815,718 rows instead of 407,859.
- **The combine is the only record of where the molecules are.** Cells are
  combined into one table; molecules are far too large to be, so they stay
  where they are and the combine writes each slice's path and shift into
  `sample_offsets.json`. If it stops doing that, the single-molecule half has
  to re-discover the slices on its own and the two halves can disagree about
  which slice a file belongs to.
- **`original = combined − (dx, dy)` must hold.** Everything downstream
  places molecules by that one line.
- **The slices must not overlap** once combined — the point of the grid.
- **Each molecule carries its slice.** The combined cells carry `_sample_id`;
  the molecules carry the same tag in a `sample` column, so a molecule can be
  told back to its slice once they share a frame.

## The real-data runner

```bash
MEYES_PYTHON=/path/to/python tests/bil/run_dry_run.sh INPUT OUT NAME \
    [--labels mapping_output.csv] [--reference DIR --species mouse] [--skip-combine]
```

It runs combine → mask → annotate → prepare → points → validate: the stages
`scripts/bil-scripts/launch_combine_mmc.sh` and `launch_process_sync.sh`
chain with `sbatch --dependency=afterok`, chained here by `set -e`. `INPUT`
is never written to, and `OUT` is never inside it — the same rule the
launchers follow by writing to `/bil/data/meyes/<sample>`.

`--reference` runs MapMyCells in the pipeline; without it, `--labels` reuses
a `mapping_output.csv` that already exists, which is what the dry runs here
do (the references are 14 GB).

## The environment

`scripts/requirements.txt`, pinned, on python 3.12 — the worker's base image:

```bash
conda create -n meyes-prepare python=3.12 && conda activate meyes-prepare
pip install -r scripts/requirements.txt pytest
```

`spatialdata` is required: `prepare.py` cannot write a store without it.
