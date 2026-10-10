"use client";

import { useState } from "react";
import NextLink from "next/link";
import {
  Modal,
  ModalContent,
  ModalHeader,
  ModalBody,
  ModalFooter,
} from "@heroui/modal";
import { Button } from "@heroui/button";
import { Tab, Tabs } from "@heroui/tabs";

type Requirement = "required" | "optional";

interface TreeNode {
  name: string;
  req?: Requirement;
  note?: string;
  children?: TreeNode[];
}

interface ColumnTable {
  caption: string;
  headers: string[];
  rows: string[][];
}

interface Format {
  key: string;
  group: "cell" | "molecule";
  title: string;
  /** What to drop on the homepage, and where. */
  drop: string;
  summary: string;
  tree?: TreeNode;
  table?: ColumnTable;
  notes: string[];
}

// Mirrors what the adapters accept (lib/adapters/*, components/file-upload.tsx,
// lib/config/moleculeColumnMappings.ts); docs/data-requirements.md has the
// exhaustive column-name lists.
const FORMATS: Format[] = [
  {
    key: "h5ad",
    group: "cell",
    title: "H5AD",
    drop: "One .h5ad file → H5AD File",
    summary:
      "An AnnData file. Cells need coordinates; everything else is optional.",
    tree: {
      name: "data.h5ad",
      children: [
        {
          name: "obsm/X_spatial",
          req: "required",
          note: "Cell coordinates, [cells × 2] or [cells × 3]. obsm/spatial also works.",
        },
        {
          name: "X",
          req: "optional",
          note: "Expression matrix, [cells × genes], dense or sparse. Needed to colour by gene.",
        },
        {
          name: "var/_index",
          req: "optional",
          note: "Gene names (var/gene or var/genes also work).",
        },
        {
          name: "obs/<column>",
          req: "optional",
          note: "Cell types, clusters and other per-cell values — each becomes a colour-by column.",
        },
        {
          name: "obsm/X_umap",
          req: "optional",
          note: "Embeddings: any obsm key starting with X_.",
        },
        {
          name: "uns/<column>_colors",
          req: "optional",
          note: "Your own colours for a categorical column.",
        },
      ],
    },
    notes: [
      "No coordinates in obsm? obs columns named center_x / center_y, x / y or centroid_x / centroid_y are used instead.",
      "A column with few distinct values is shown as categories; one with many numeric values as a gradient. Columns named leiden or louvain are always categories.",
    ],
  },
  {
    key: "xenium",
    group: "cell",
    title: "Xenium",
    drop: "The Xenium output folder → Folder",
    summary: "The folder Xenium writes. Only cells.csv is required.",
    tree: {
      name: "xenium_output/",
      children: [
        {
          name: "cells.csv",
          req: "required",
          note: "One row per cell with its centroid (x_centroid, y_centroid). .csv.gz also works.",
        },
        {
          name: "cell_feature_matrix/",
          req: "optional",
          note: "Expression. Needed to colour by gene.",
          children: [
            { name: "matrix.mtx.gz" },
            { name: "features.tsv.gz" },
            { name: "barcodes.tsv.gz" },
          ],
        },
        {
          name: "analysis/clustering/",
          req: "optional",
          note: "Cluster labels (e.g. gene_expression_graphclust/clusters.csv).",
        },
      ],
    },
    notes: [
      "A cell type or cluster column in cells.csv (cell_type, cluster, leiden, annotation, …) is picked up automatically.",
      "Xenium cells are shown in 2D.",
    ],
  },
  {
    key: "merscope",
    group: "cell",
    title: "MERSCOPE",
    drop: "The MERSCOPE output folder → Folder",
    summary: "The folder MERSCOPE writes. Only cell_metadata.csv is required.",
    tree: {
      name: "merscope_output/",
      children: [
        {
          name: "cell_metadata.csv",
          req: "required",
          note: "One row per cell with its centre (center_x, center_y).",
        },
        {
          name: "cell_by_gene.csv",
          req: "optional",
          note: "Expression: one row per cell, one column per gene. Needed to colour by gene.",
        },
        {
          name: "cell_categories.csv",
          req: "optional",
          note: "Cluster labels (leiden, cell_type, …).",
        },
        {
          name: "cell_numeric_categories.csv",
          req: "optional",
          note: "UMAP (umap_X, umap_Y).",
        },
      ],
    },
    notes: ["MERSCOPE cells are shown in 2D."],
  },
  {
    key: "zarr",
    group: "cell",
    title: "Zarr",
    drop: "The .zarr folder → Folder",
    summary:
      "An AnnData saved with write_zarr. Same contents as an H5AD, as a folder.",
    tree: {
      name: "data.zarr/",
      children: [
        {
          name: ".zgroup",
          req: "required",
          note: "The Zarr marker at the top of the folder — how it is recognised.",
        },
        {
          name: "obsm/X_spatial/",
          req: "required",
          note: "Cell coordinates. obsm/spatial also works.",
        },
        {
          name: "X/",
          req: "optional",
          note: "Expression matrix, dense or sparse. Needed to colour by gene.",
        },
        { name: "var/_index/", req: "optional", note: "Gene names." },
        {
          name: "obs/<column>/",
          req: "optional",
          note: "Cell types, clusters and other per-cell values.",
        },
      ],
    },
    notes: [
      "A dense or CSC matrix loads one gene at a time. A CSR matrix is expanded in the browser first, which is slow for large datasets — save X as CSC if you can.",
      "Zarr v2 stores are supported; Zarr v3 stores are not read yet.",
    ],
  },
  {
    key: "chunked",
    group: "cell",
    title: "Pre-chunked",
    drop: "The output folder of our Python script → Folder",
    summary:
      "For datasets too large to process in the browser: convert once in Python, then upload the result as it is.",
    tree: {
      name: "output_folder/",
      children: [
        { name: "manifest.json", req: "required" },
        {
          name: "coords/",
          children: [
            { name: "spatial.bin.gz", req: "required" },
            { name: "umap.bin.gz", req: "optional" },
          ],
        },
        {
          name: "expr/",
          children: [
            { name: "index.json", req: "required" },
            {
              name: "chunk_00000.bin.gz",
              req: "required",
              note: "At least one.",
            },
          ],
        },
        {
          name: "obs/",
          children: [
            { name: "metadata.json", req: "optional" },
            { name: "<column>.json.gz", req: "optional" },
          ],
        },
        {
          name: "palettes/",
          children: [{ name: "<column>.json", req: "optional" }],
        },
      ],
    },
    notes: [
      "Made with: python scripts/process_spatial_data.py <h5ad | Xenium folder | MERSCOPE folder> output_folder",
      "The script is in the project's GitHub repository, under scripts/.",
    ],
  },
  {
    key: "molecules",
    group: "molecule",
    title: "Parquet / CSV",
    drop: "One .parquet or .csv file → File",
    summary: "One row per detected molecule: which gene it is and where it is.",
    table: {
      caption: "The columns we look for",
      headers: ["", "Gene", "X", "Y", "Z (optional)"],
      rows: [
        ["Xenium", "feature_name", "x_location", "y_location", "z_location"],
        ["MERSCOPE", "gene", "global_x", "global_y", "global_z"],
        [
          "Anything else",
          "you choose",
          "you choose",
          "you choose",
          "you choose",
        ],
      ],
    },
    notes: [
      "Xenium (transcripts.parquet) and MERSCOPE (detected_transcripts.csv) files are recognised by their column names. For other column names, process the file in the browser and you are asked which columns to use.",
      "Without a Z column the molecules are shown in 2D.",
      "Control probes, blank codewords and unassigned entries are left out.",
    ],
  },
  {
    key: "molecules-chunked",
    group: "molecule",
    title: "Pre-chunked",
    drop: "The output folder of our Python script → Chunked Folder",
    summary:
      "For very large molecule files: convert once in Python, then upload the result as it is.",
    tree: {
      name: "output_folder/",
      children: [
        { name: "manifest.json.gz", req: "required" },
        {
          name: "genes/",
          req: "required",
          children: [{ name: "<gene>.bin.gz", note: "One file per gene." }],
        },
      ],
    },
    notes: [
      "Made with: python scripts/process_single_molecule.py transcripts.parquet output_folder",
      "The script is in the project's GitHub repository, under scripts/.",
    ],
  },
];

function RequirementTag({ req }: { req: Requirement }) {
  return (
    <span
      className={
        req === "required"
          ? "rounded-full bg-primary/20 px-2 py-0.5 text-[11px] text-primary"
          : "rounded-full bg-default-200/60 px-2 py-0.5 text-[11px] text-default-500"
      }
    >
      {req}
    </span>
  );
}

/** A file / folder structure, drawn as an indented tree. */
function Tree({ node, root = true }: { node: TreeNode; root?: boolean }) {
  return (
    <li className="list-none">
      <div className="flex flex-wrap items-baseline gap-x-2 gap-y-0.5 py-1">
        <span className={`font-mono text-sm ${root ? "font-semibold" : ""}`}>
          {node.name}
        </span>
        {node.req && <RequirementTag req={node.req} />}
        {node.note && (
          <span className="text-xs text-default-500">{node.note}</span>
        )}
      </div>
      {node.children && (
        <ul className="ml-2 border-l border-default-300 pl-4">
          {node.children.map((child) => (
            <Tree key={child.name} node={child} root={false} />
          ))}
        </ul>
      )}
    </li>
  );
}

function FormatDetail({ format }: { format: Format }) {
  return (
    <div className="flex flex-col gap-4 text-left">
      <div>
        <p className="text-sm">{format.summary}</p>
        <p className="mt-1 text-xs text-default-500">Upload: {format.drop}</p>
      </div>

      {format.tree && (
        <ul className="rounded-xl border border-default-200 bg-default-50/40 px-4 py-3">
          <Tree node={format.tree} />
        </ul>
      )}

      {format.table && (
        <div className="overflow-x-auto rounded-xl border border-default-200 bg-default-50/40 px-4 py-3">
          <table className="w-full text-sm">
            <caption className="pb-2 text-left text-xs text-default-500">
              {format.table.caption}
            </caption>
            <thead>
              <tr>
                {format.table.headers.map((h) => (
                  <th
                    key={h}
                    className="border-b border-default-200 py-1 pr-4 text-left text-xs font-medium text-default-500"
                  >
                    {h}
                  </th>
                ))}
              </tr>
            </thead>
            <tbody>
              {format.table.rows.map(([label, ...cells]) => (
                <tr key={label}>
                  <td className="py-1 pr-4 text-xs text-default-500">
                    {label}
                  </td>
                  {cells.map((cell, i) => (
                    <td key={i} className="py-1 pr-4 font-mono">
                      {cell}
                    </td>
                  ))}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}

      <ul className="flex list-disc flex-col gap-1 pl-5 text-xs text-default-500">
        {format.notes.map((note) => (
          <li key={note}>{note}</li>
        ))}
      </ul>
    </div>
  );
}

/**
 * What can be uploaded, format by format, with the expected structure of
 * each. `group` limits it to single cell or single molecule formats.
 */
export function FormatGuide({ group }: { group?: Format["group"] }) {
  const formats = group ? FORMATS.filter((f) => f.group === group) : FORMATS;

  return (
    <Tabs aria-label="Data formats" variant="underlined">
      {formats.map((format) => (
        <Tab
          key={format.key}
          title={
            group
              ? format.title
              : `${format.title}${format.group === "molecule" ? " (molecules)" : ""}`
          }
        >
          <FormatDetail format={format} />
        </Tab>
      ))}
    </Tabs>
  );
}

/** Homepage button that opens the guide for the active upload mode. */
export function FormatGuideButton({ group }: { group: Format["group"] }) {
  const [open, setOpen] = useState(false);

  return (
    <>
      <button
        className="inline-flex items-center gap-2 px-4 py-2 rounded-full text-sm border border-default-300 text-default-700 bg-default-100/40 hover:bg-default-200/60 hover:border-primary/40 focus:outline-none focus-visible:ring-2 focus-visible:ring-primary/50 transition-colors"
        type="button"
        onClick={() => setOpen(true)}
      >
        <svg
          aria-hidden
          className="h-4 w-4"
          fill="none"
          stroke="currentColor"
          strokeWidth={1.8}
          viewBox="0 0 24 24"
        >
          <circle cx="12" cy="12" r="9" />
          <path d="M12 11v5M12 8h.01" strokeLinecap="round" />
        </svg>
        <span>What can I upload?</span>
      </button>
      <Modal
        isOpen={open}
        scrollBehavior="inside"
        size="3xl"
        onClose={() => setOpen(false)}
      >
        <ModalContent>
          <ModalHeader>
            {group === "molecule"
              ? "Single molecule formats"
              : "Single cell formats"}
          </ModalHeader>
          <ModalBody>
            <FormatGuide group={group} />
          </ModalBody>
          <ModalFooter>
            <Button as={NextLink} href="/docs" variant="light">
              All formats
            </Button>
            <Button color="primary" onPress={() => setOpen(false)}>
              Close
            </Button>
          </ModalFooter>
        </ModalContent>
      </Modal>
    </>
  );
}
