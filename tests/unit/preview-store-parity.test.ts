import { execFileSync } from "node:child_process";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";

import { afterAll, beforeAll, describe, expect, it } from "vitest";

import { keepAwake, previewOf, type Seen, storeOf } from "./helpers/seen";

/**
 * What you previewed is what you saved (docs/plans/PREVIEW-SAVE.md, D5):
 * each source opened as the browser preview opens it, against the store
 * scripts/prepare.py wrote from it opened as the SpatialData viewer does.
 *
 * The tiny sources' stores are committed
 * (scripts/testdata/make-tiny-parity-stores.py) and kept fresh by
 * tests/python/test_parity_stores.py, so this runs in CI without Python.
 * A real Xenium export (gigabytes, none committed) runs when you name one and
 * a Python with prepare.py's packages:
 *
 *   MERFISHEYES_XENIUM=/path/to/outs MERFISHEYES_PREPARE_PYTHON=/path/to/python \
 *     npx vitest run tests/unit/preview-store-parity.test.ts
 */

const ROOT = path.join(__dirname, "..", "..");
const TINY = path.join(ROOT, "tests", "data", "tiny");
const XENIUM = process.env.MERFISHEYES_XENIUM;
const PYTHON = process.env.MERFISHEYES_PREPARE_PYTHON;

const committed = (name: string) =>
  path.join(TINY, "parity", name, `${name}.zarr`);

/** prepare.py run now, for a source too big to commit a store for. */
function prepared(name: string, source: string): string {
  const out = path.join(
    fs.mkdtempSync(path.join(os.tmpdir(), "parity-")),
    name,
  );

  execFileSync(
    PYTHON!,
    [
      path.join(ROOT, "scripts", "prepare.py"),
      source,
      out,
      "--name",
      name,
      "--no-images",
    ],
    { stdio: "ignore" },
  );

  return path.join(out, `${name}.zarr`);
}

const CASES: { name: string; source: string; store: () => string }[] = [
  { name: "tiny_h5ad", source: path.join(TINY, "h5ad", "tiny.h5ad") },
  { name: "nullable", source: path.join(TINY, "h5ad", "nullable.h5ad") },
  { name: "named_index", source: path.join(TINY, "h5ad", "named_index.h5ad") },
  { name: "missing_xy", source: path.join(TINY, "h5ad", "missing_xy.h5ad") },
  { name: "merscope", source: path.join(TINY, "merscope", "tiny") },
]
  .map((c) => ({ ...c, store: () => committed(c.name) }))
  .concat(
    XENIUM && PYTHON
      ? [
          {
            name: "xenium",
            source: XENIUM,
            store: () => prepared("xenium", XENIUM),
          },
        ]
      : [],
  );

/**
 * The differences known today, each with its cause and the step of the plan
 * that removes it. A known difference is an expected failure (`it.fails`):
 * once fixed, its test passes, which fails the suite until the entry here is
 * removed — so a fix cannot go unrecorded, nor a new difference unnoticed.
 */
const KNOWN: Record<string, Partial<Record<Property, string>>> = {
  tiny_h5ad: {
    colours:
      "step 3: the store sorts labels numerically (cluster_2 before cluster_10), the preview as text, so colours are dealt in a different order",
  },
  nullable: {
    colours:
      "step 3: the store labels booleans true/false, the preview True/False",
  },
  missing_xy: {
    columns:
      "step 3: the store offers the coordinate columns center_x / center_y",
  },
  merscope: {
    columns:
      "step 3: the store offers every cell_categories column (celltype too) and the coordinate columns",
    colours: "step 3: the label sort, as tiny_h5ad",
  },
  // Only with MERFISHEYES_XENIUM (measured on the 10x lung 2-FOV export)
  xenium: {
    columns:
      "step 3: the store has no Cluster (10x's own clustering, in analysis.tar.gz, is not read) and offers z_level and nucleus_count",
    defaultColumn:
      "step 3: without Cluster the store opens on cell_area, a gradient",
    colours:
      "step 3: Cluster is missing, and transcript_counts' labels sort differently, as tiny_h5ad",
  },
};

type Property =
  | "cells"
  | "genes"
  | "expression"
  | "positions"
  | "columns"
  | "defaultColumn"
  | "colours"
  | "embeddings";

let letGo: () => void;

beforeAll(() => {
  letGo = keepAwake();
});
afterAll(() => letGo());

/**
 * What a person sees of positions: whether each cell is drawn (both
 * coordinates finite; a store keeps a lone valid y beside a missing x, the
 * preview drops both, and neither draws it) and, if so, where — equal within
 * the preview's two-decimal rounding plus the store's float32 precision.
 */
function samePositions(a: number[][], b: number[][]): string | null {
  if (a.length !== b.length) return `${a.length} vs ${b.length} positions`;
  const drawn = (p: number[]) => Number.isFinite(p[0]) && Number.isFinite(p[1]);

  for (let i = 0; i < a.length; i++) {
    if (drawn(a[i]) !== drawn(b[i])) return `cell ${i}: ${a[i]} vs ${b[i]}`;
    if (!drawn(a[i])) continue;
    for (const k of [0, 1]) {
      const tolerance = 0.005 + Math.abs(b[i][k]) * 2 ** -23;

      if (Math.abs(a[i][k] - b[i][k]) > tolerance) {
        return `cell ${i}: ${a[i]} vs ${b[i]}`;
      }
    }
  }

  return null;
}

describe.each(CASES)(
  "$name: preview vs saved store",
  ({ name, source, store: storeDir }) => {
    let preview: Seen;
    let store: Seen;

    beforeAll(async () => {
      preview = await previewOf(source);
      store = await storeOf(storeDir());
    }, 600_000);

    /** A check, or an expected failure when the difference is known. */
    const check = (
      property: Property,
      title: string,
      body: () => void | Promise<void>,
    ) => {
      const known = KNOWN[name]?.[property];

      return known
        ? it.fails(`${title} — known: ${known}`, body)
        : it(title, body);
    };

    check("cells", "the same cells", () => {
      expect(store.cells).toBe(preview.cells);
    });

    check("genes", "the same genes, in order", () => {
      expect(store.genes).toEqual(preview.genes);
    });

    check("expression", "the same expression", async () => {
      for (const gene of preview.genes.slice(0, 5)) {
        const [a, b] = [
          await preview.expression(gene),
          await store.expression(gene),
        ];

        expect(b.length, gene).toBe(a.length);
        b.forEach((v, i) =>
          expect(v, `${gene} cell ${i}`).toBeCloseTo(a[i], 4),
        );
      }
    });

    check("positions", "the same cells drawn, in the same places", () => {
      expect(samePositions(preview.positions, store.positions)).toBeNull();
    });

    check("columns", "the same colour columns, of the same kinds", () => {
      expect(store.columns).toEqual(preview.columns);
    });

    check("defaultColumn", "opens on the same column", () => {
      expect(store.defaultColumn).toBe(preview.defaultColumn);
    });

    check("colours", "the same colour for every label", () => {
      expect(store.colours).toEqual(preview.colours);
    });

    check("embeddings", "the same embeddings", () => {
      expect(store.embeddings).toEqual(preview.embeddings);
    });
  },
);
