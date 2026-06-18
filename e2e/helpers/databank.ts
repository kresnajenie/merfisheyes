import fs from "node:fs";
import path from "node:path";

/** Repo root (e2e/helpers -> ../../). */
export const REPO_ROOT = path.resolve(__dirname, "..", "..");
const MANIFEST_PATH = path.join(REPO_ROOT, "tests", "data", "manifest.json");

export type DatasetFormat =
  | "h5ad"
  | "xenium"
  | "merscope"
  | "chunked"
  | "single-molecule";
export type DataType = "single_cell" | "single_molecule";
export type Tier = "tiny" | "medium" | "large";

export interface ExpectedCounts {
  cells?: number;
  molecules?: number;
  genes: number;
  dimensions: number;
  clusterColumns?: string[];
  obsFields?: string[];
  knownGenes: string[];
}

export interface DatasetEntry {
  id: string;
  format: DatasetFormat;
  dataType: DataType;
  tier: Tier;
  path: string; // repo-relative
  committed: boolean;
  source: { kind: "generate" | "download"; [k: string]: unknown } | null;
  expected: ExpectedCounts;
}

interface Manifest {
  version: number;
  datasets: DatasetEntry[];
}

let cache: Manifest | null = null;

function load(): Manifest {
  if (!cache) cache = JSON.parse(fs.readFileSync(MANIFEST_PATH, "utf8"));
  return cache!;
}

export function allDatasets(): DatasetEntry[] {
  return load().datasets;
}

export function datasetsByTier(...tiers: Tier[]): DatasetEntry[] {
  return load().datasets.filter((d) => tiers.includes(d.tier));
}

export function getDataset(id: string): DatasetEntry {
  const ds = load().datasets.find((d) => d.id === id);
  if (!ds) throw new Error(`Dataset '${id}' not found in manifest`);
  return ds;
}

/** Absolute filesystem path to a dataset's file/folder. */
export function absPath(ds: DatasetEntry): string {
  return path.join(REPO_ROOT, ds.path);
}

/** True if the dataset's files are present on disk. */
export function isMaterialized(ds: DatasetEntry): boolean {
  return fs.existsSync(absPath(ds));
}
