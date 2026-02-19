import type { VisualizationMode } from "@/lib/stores/createVisualizationStore";
import type { ViewMode } from "@/lib/stores/createSingleMoleculeVisualizationStore";

import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";

// --- Compact JSON shapes for URL encoding ---

export interface CellVizUrlState {
  g?: string; // selectedGene
  c?: string; // selectedColumn
  ct?: string[]; // selectedCelltypes
  m?: VisualizationMode[]; // mode array
  gs?: [number, number]; // [geneScaleMin, geneScaleMax]
  ns?: [number, number]; // [numericalScaleMin, numericalScaleMax]
  sz?: number; // sizeScale
  e?: string; // selectedEmbedding
}

// Gene tuple: [name, color, localScale, isVisible]
export type SMGeneTuple = [string, string, number, boolean];

export interface SMVizUrlState {
  genes?: SMGeneTuple[];
  gs?: number; // globalScale
  vm?: ViewMode; // viewMode
}

// --- Base64URL helpers (Unicode-safe) ---

function toBase64Url(obj: object): string {
  const json = JSON.stringify(obj);
  const bytes = new TextEncoder().encode(json);
  const binary = String.fromCharCode(...bytes);
  const b64 = btoa(binary);

  return b64.replace(/\+/g, "-").replace(/\//g, "_").replace(/=+$/, "");
}

function fromBase64Url<T>(encoded: string): T | null {
  try {
    let b64 = encoded.replace(/-/g, "+").replace(/_/g, "/");
    const pad = b64.length % 4;

    if (pad) b64 += "=".repeat(4 - pad);

    const binary = atob(b64);
    const bytes = new Uint8Array(binary.length);

    for (let i = 0; i < binary.length; i++) {
      bytes[i] = binary.charCodeAt(i);
    }

    const json = new TextDecoder().decode(bytes);

    return JSON.parse(json) as T;
  } catch {
    return null;
  }
}

// --- Single Cell Codec ---

const DEFAULT_SCALE_MIN = VISUALIZATION_CONFIG.SCALE_BAR_DEFAULT_MIN;
const DEFAULT_SCALE_MAX = VISUALIZATION_CONFIG.SCALE_BAR_DEFAULT_MAX;
const DEFAULT_MODE: VisualizationMode[] = ["celltype"];

export function encodeCellVizState(state: {
  selectedGene: string | null;
  selectedColumn: string | null;
  selectedCelltypes: Set<string>;
  mode: VisualizationMode[];
  geneScaleMin: number;
  geneScaleMax: number;
  numericalScaleMin: number;
  numericalScaleMax: number;
  sizeScale: number;
  selectedEmbedding: string | null;
}): string | null {
  const obj: CellVizUrlState = {};

  if (state.selectedGene) obj.g = state.selectedGene;
  if (state.selectedColumn) obj.c = state.selectedColumn;
  if (state.selectedCelltypes.size > 0)
    obj.ct = Array.from(state.selectedCelltypes);

  // Only include mode if not default ["celltype"]
  const modeStr = JSON.stringify(state.mode);
  const defaultModeStr = JSON.stringify(DEFAULT_MODE);

  if (modeStr !== defaultModeStr) obj.m = state.mode;

  if (
    state.geneScaleMin !== DEFAULT_SCALE_MIN ||
    state.geneScaleMax !== DEFAULT_SCALE_MAX
  )
    obj.gs = [state.geneScaleMin, state.geneScaleMax];

  if (
    state.numericalScaleMin !== DEFAULT_SCALE_MIN ||
    state.numericalScaleMax !== DEFAULT_SCALE_MAX
  )
    obj.ns = [state.numericalScaleMin, state.numericalScaleMax];

  if (state.sizeScale !== 1.0) obj.sz = state.sizeScale;
  if (state.selectedEmbedding) obj.e = state.selectedEmbedding;

  // Don't encode if nothing interesting
  if (Object.keys(obj).length === 0) return null;

  return toBase64Url(obj);
}

export function decodeCellVizState(encoded: string): CellVizUrlState | null {
  const obj = fromBase64Url<CellVizUrlState>(encoded);

  if (!obj || typeof obj !== "object") return null;

  // Basic validation
  if (obj.g !== undefined && typeof obj.g !== "string") return null;
  if (obj.c !== undefined && typeof obj.c !== "string") return null;
  if (obj.ct !== undefined && !Array.isArray(obj.ct)) return null;
  if (obj.m !== undefined && !Array.isArray(obj.m)) return null;
  if (obj.gs !== undefined && (!Array.isArray(obj.gs) || obj.gs.length !== 2))
    return null;
  if (obj.ns !== undefined && (!Array.isArray(obj.ns) || obj.ns.length !== 2))
    return null;
  if (obj.sz !== undefined && typeof obj.sz !== "number") return null;
  if (obj.e !== undefined && typeof obj.e !== "string") return null;

  return obj;
}

// --- Single Molecule Codec ---

export function encodeSMVizState(state: {
  selectedGenes: Map<
    string,
    { gene: string; color: string; localScale: number }
  >;
  geneDataCache: Map<
    string,
    { gene: string; color: string; localScale: number }
  >;
  globalScale: number;
  viewMode: ViewMode;
}): string | null {
  const obj: SMVizUrlState = {};

  // Build gene tuples from geneDataCache (all legend genes)
  if (state.geneDataCache.size > 0) {
    const genes: SMGeneTuple[] = [];

    state.geneDataCache.forEach((viz, gene) => {
      const isVisible = state.selectedGenes.has(gene);

      genes.push([gene, viz.color, viz.localScale, isVisible]);
    });

    if (genes.length > 0) obj.genes = genes;
  }

  if (state.globalScale !== 1.0) obj.gs = state.globalScale;
  if (state.viewMode !== "2D") obj.vm = state.viewMode;

  if (Object.keys(obj).length === 0) return null;

  return toBase64Url(obj);
}

export function decodeSMVizState(encoded: string): SMVizUrlState | null {
  const obj = fromBase64Url<SMVizUrlState>(encoded);

  if (!obj || typeof obj !== "object") return null;

  if (obj.genes !== undefined && !Array.isArray(obj.genes)) return null;
  if (obj.gs !== undefined && typeof obj.gs !== "number") return null;
  if (obj.vm !== undefined && obj.vm !== "2D" && obj.vm !== "3D") return null;

  return obj;
}
