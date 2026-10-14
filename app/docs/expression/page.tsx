import type { Metadata } from "next";

import {
  Callout,
  CompareTable,
  DocHeader,
  DocLink,
  DocSection,
  Flow,
  Prose,
} from "@/components/docs/doc-ui";

export const metadata: Metadata = {
  title: "Expression & DEGs — MERFISHEYES Docs",
  description:
    "Where single cell expression comes from, how raw counts and log values are shown, and how differentially expressed genes are calculated.",
};

export default function ExpressionPage() {
  return (
    <>
      <DocHeader
        heading="Expression & DEGs"
        lead="Where the expression values of your cells come from, how the viewer shows raw counts and log values, and how the DEG panel ranks genes."
      />

      <DocSection heading="Where the expression comes from" id="source">
        <Prose>
          <p>
            MERFISHEYES never normalizes, scales or log-transforms your
            expression. Each cell&apos;s values are read from your data as they
            are, and stay that way when a dataset is processed on our servers or
            saved. What you see depends on what you brought:
          </p>
        </Prose>
        <CompareTable
          headers={["Values read from", "Usually"]}
          labelHeader="You bring"
          rows={[
            [
              "H5AD",
              <>
                <code>X</code> — not <code>.raw</code>, <code>layers</code> or{" "}
                <code>obsm</code>
              </>,
              "Whatever X holds: often log-normalized, sometimes raw counts",
            ],
            [
              "Xenium",
              <>
                <code>cell_feature_matrix</code> (.h5 or the folder), or{" "}
                <code>cell_by_gene.csv</code>
              </>,
              "Raw counts",
            ],
            ["MERSCOPE", <code key="m">cell_by_gene.csv</code>, "Raw counts"],
            [
              "SpatialData store, or an AnnData table in Zarr",
              <>
                The table&apos;s <code>X</code>, or a CSC copy of it in{" "}
                <code>layers/</code> (e.g. <code>X_csc</code>) when X is CSR
              </>,
              "Whatever X holds",
            ],
            [
              "A processed folder (Python scripts or server processing)",
              "The chunks written from one of the above",
              "Unchanged from the source",
            ],
          ]}
        />
        <Callout>
          If your H5AD keeps raw counts somewhere other than <code>X</code> and
          you want to see counts, put them in <code>X</code> before uploading.
          The cell type CSV you bring yourself (e.g. from MapMyCells) adds
          columns; it doesn&apos;t change expression.
        </Callout>
      </DocSection>

      <DocSection heading="Raw counts or log values?" id="raw-log">
        <Prose>
          <p>
            When you select a gene, the viewer looks at its values to tell which
            kind they are: raw counts are whole numbers, normalized or log
            values are mostly fractions. If more than 5% of the nonzero values
            it samples (up to 5,000) are fractions, the data is treated as
            log-normalized; otherwise as counts. The scale bar shows what it
            decided (<b>data: counts</b> or <b>data: log-norm</b>).
          </p>
        </Prose>
        <Flow
          steps={[
            {
              title: "Your values",
              detail: "Exactly as stored, for the selected gene",
            },
            {
              title: "Counts or log?",
              detail:
                "Whole numbers → counts; mostly fractions → log-normalized",
            },
            {
              title: "Raw / Log toggle",
              detail:
                "Starts on the kind your data is; switch it under the scale bar",
            },
            {
              title: "Colours",
              detail:
                "Scaled from 0 to the 95th percentile; drag the numbers to change it",
            },
          ]}
        />
        <CompareTable
          headers={["Shown as Raw", "Shown as Log"]}
          labelHeader="Your data is"
          rows={[
            [
              "Raw counts",
              "As stored",
              <>
                <code>log1p(x)</code> = ln(1 + x) — exact
              </>,
            ],
            [
              "Log-normalized",
              <>
                <code>expm1(x)</code> = eˣ − 1 — back toward counts, but not the
                true integers (normalization can&apos;t be undone)
              </>,
              "As stored",
            ],
          ]}
        />
        <Prose>
          <p>
            The toggle changes only how values are shown — the colours, the
            scale bar and the{" "}
            <DocLink href="/docs/viewer#single-cell">plots</DocLink> — never the
            data. The log is the natural log, matching what{" "}
            <code>scanpy.pp.log1p</code> writes.
          </p>
        </Prose>
      </DocSection>

      <DocSection heading="How DEGs are calculated" id="deg">
        <Prose>
          <p>
            The DEG panel compares one group of cells (a cell type, a cluster,
            or any value of a column shown with one colour per value) against
            either all other cells or one other group you pick.
          </p>
        </Prose>
        <Flow
          steps={[
            {
              title: "Per-group summary",
              detail:
                "For every gene and every group: the mean expression, and the fraction of cells where it is above 0",
            },
            {
              title: "Choose the comparison",
              detail:
                "Group 1 against the rest (their means combined, weighted by cell count) or against group 2",
            },
            {
              title: "Back to linear units",
              detail:
                "If the data is log-normalized, both means go through expm1 first",
            },
            {
              title: "Fold change",
              detail: "mean 1 ÷ mean 2, genes ranked from highest",
            },
          ]}
        />
        <CompareTable
          headers={["What it is"]}
          labelHeader="In the panel"
          rows={[
            [
              "mean1 / mean2",
              "Mean expression of group 1 and of group 2 (or the rest), in linear units",
            ],
            [
              "FC",
              "mean1 ÷ mean2 — 5.8× means almost six times higher in group 1. The list is sorted by it",
            ],
            [
              "Hover a row",
              "The share of each group's cells where the gene is above 0",
            ],
            [
              "Above the list",
              "How many cells are in each group, e.g. 24,317 in target · 123,085 in rest",
            ],
          ]}
        />
        <Prose>
          <p>
            <b>Counts or log, for DEGs.</b> The panel decides from the per-group
            means: if the largest is above 30 it treats the data as raw counts,
            otherwise as log-normalized (log1p of size-normalized counts rarely
            goes above about 15). Its <b>Expression</b> menu shows the guess and
            lets you override it. The setting only affects the means and fold
            changes shown.
          </p>
          <p>
            <b>Means of log data.</b> For log-normalized data the mean is taken
            on the log values and then turned back with expm1, so mean1 and
            mean2 are closer to a typical cell&apos;s count than to the
            arithmetic average of counts. Fold changes between groups compare
            like with like.
          </p>
          <p>
            <b>Where the summaries come from.</b> For data opened in your
            browser they are computed there when you choose a column. For
            datasets processed on our servers (or with the Python scripts) they
            were computed once, from the same values, for every column with
            categories; SpatialData stores read them from an{" "}
            <code>annotations/</code> folder beside the store. (The Leiden
            clusters that <code>scripts/spatialdata_cluster.py</code> adds to a
            store are the one exception to &quot;as stored&quot;: their
            summaries are taken from a log-normalized copy of the table.){" "}
            <b>Copy CSV</b> downloads the table, with a note of which mode was
            used.
          </p>
        </Prose>
        <Callout kind="note">
          This is a fold-change ranking, not a statistical test: there are no
          p-values or corrections for multiple testing. Use it to find candidate
          markers quickly, and confirm them with a proper test (for example{" "}
          <code>scanpy.tl.rank_genes_groups</code>) before reporting them.
        </Callout>
      </DocSection>
    </>
  );
}
