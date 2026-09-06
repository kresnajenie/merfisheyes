import type {
  VisualizationMode,
  CellViewMode,
} from "@/lib/stores/createVisualizationStore";
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
  to?: Record<string, "categorical" | "numerical">; // columnTypeOverrides
  vm?: CellViewMode; // viewMode (2D/3D)
  av?: Record<string, number>; // advanced viz settings (only non-default values)
  rot?: number; // sceneRotation (degrees)
  fx?: boolean; // flipX
  fy?: boolean; // flipY
  cm?: string; // colormap
  // Two-gene coexpression
  g2?: string; // selectedGene2
  gs2?: [number, number]; // [gene2ScaleMin, gene2ScaleMax]
  ce?: boolean; // coexpressEnabled (omitted when false)
  cw?: boolean; // coexpressSwapped (omitted when false)
}

// Gene tuple: [name, color, localScale, isVisible, showAssigned, showUnassigned, unassignedColor, unassignedLocalScale, colorSynced]
// Fields 4-8 are optional for backwards compat (default: true, true, same as color, same as localScale, true)
export type SMGeneTuple = [
  string,
  string,
  number,
  boolean,
  boolean?,
  boolean?,
  string?,
  number?,
  boolean?,
];

export interface SMVizUrlState {
  genes?: SMGeneTuple[];
  gs?: number; // globalScale
  vm?: ViewMode; // viewMode
  sa?: boolean; // global showAssigned
  su?: boolean; // global showUnassigned
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
  columnTypeOverrides: Record<string, "categorical" | "numerical">;
  viewMode: CellViewMode;
  selectedSizeMultiplier: number;
  greyedOutSizeMultiplier: number;
  greyedOutAlpha: number;
  expressionAlphaMin: number;
  expressionAlphaMax: number;
  pointSizeMultiplierMin: number;
  pointSizeMultiplierMax: number;
  targetPx: number;
  sceneRotation: number;
  flipX: boolean;
  flipY: boolean;
  colormap: string;
  selectedGene2: string | null;
  gene2ScaleMin: number;
  gene2ScaleMax: number;
  coexpressEnabled: boolean;
  coexpressSwapped: boolean;
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
  if (Object.keys(state.columnTypeOverrides).length > 0)
    obj.to = state.columnTypeOverrides;
  if (state.viewMode === "3D") obj.vm = "3D";

  // Encode advanced settings (only non-default values)
  const avDefaults: Record<string, number> = {
    selectedSizeMultiplier:
      VISUALIZATION_CONFIG.SELECTED_SIZE_MULTIPLIER as number,
    greyedOutSizeMultiplier:
      VISUALIZATION_CONFIG.GREYED_OUT_SIZE_MULTIPLIER as number,
    greyedOutAlpha: VISUALIZATION_CONFIG.GREYED_OUT_ALPHA as number,
    expressionAlphaMin: VISUALIZATION_CONFIG.EXPRESSION_ALPHA_MIN as number,
    expressionAlphaMax: VISUALIZATION_CONFIG.EXPRESSION_ALPHA_MAX as number,
    pointSizeMultiplierMin:
      VISUALIZATION_CONFIG.POINT_SIZE_MULTIPLIER_MIN as number,
    pointSizeMultiplierMax:
      VISUALIZATION_CONFIG.POINT_SIZE_MULTIPLIER_MAX as number,
    targetPx: VISUALIZATION_CONFIG.TARGET_PX_DEFAULT as number,
  };
  const av: Record<string, number> = {};

  for (const [key, defaultVal] of Object.entries(avDefaults)) {
    const val = (state as any)[key] as number;

    if (val !== defaultVal) av[key] = val;
  }
  if (Object.keys(av).length > 0) obj.av = av;
  if (state.sceneRotation !== 0) obj.rot = state.sceneRotation;
  if (state.flipX) obj.fx = true;
  if (state.flipY) obj.fy = true;
  if (state.colormap && state.colormap !== "bwr") obj.cm = state.colormap;

  // Two-gene coexpression
  if (state.selectedGene2) obj.g2 = state.selectedGene2;
  if (
    state.gene2ScaleMin !== DEFAULT_SCALE_MIN ||
    state.gene2ScaleMax !== DEFAULT_SCALE_MAX
  )
    obj.gs2 = [state.gene2ScaleMin, state.gene2ScaleMax];
  if (state.coexpressEnabled) obj.ce = true;
  if (state.coexpressSwapped) obj.cw = true;

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
  if (obj.to !== undefined && (typeof obj.to !== "object" || obj.to === null))
    return null;
  if (obj.vm !== undefined && obj.vm !== "2D" && obj.vm !== "3D") return null;
  if (obj.g2 !== undefined && typeof obj.g2 !== "string") return null;
  if (
    obj.gs2 !== undefined &&
    (!Array.isArray(obj.gs2) || obj.gs2.length !== 2)
  )
    return null;
  if (obj.ce !== undefined && typeof obj.ce !== "boolean") return null;
  if (obj.cw !== undefined && typeof obj.cw !== "boolean") return null;

  return obj;
}

// --- Single Molecule Codec ---

export function encodeSMVizState(state: {
  selectedGenes: Map<
    string,
    {
      gene: string;
      color: string;
      localScale: number;
      showAssigned: boolean;
      showUnassigned: boolean;
      unassignedColor: string;
      unassignedLocalScale: number;
      colorSynced: boolean;
    }
  >;
  geneDataCache: Map<
    string,
    {
      gene: string;
      color: string;
      localScale: number;
      showAssigned: boolean;
      showUnassigned: boolean;
      unassignedColor: string;
      unassignedLocalScale: number;
      colorSynced: boolean;
    }
  >;
  globalScale: number;
  viewMode: ViewMode;
  showAssigned: boolean;
  showUnassigned: boolean;
}): string | null {
  const obj: SMVizUrlState = {};

  // Build gene tuples from geneDataCache (all legend genes)
  if (state.geneDataCache.size > 0) {
    const genes: SMGeneTuple[] = [];

    state.geneDataCache.forEach((viz, gene) => {
      const isVisible = state.selectedGenes.has(gene);

      // Only include extended fields when they differ from defaults
      const hasExtended =
        !viz.showAssigned ||
        !viz.showUnassigned ||
        viz.unassignedColor !== viz.color ||
        viz.unassignedLocalScale !== viz.localScale ||
        !viz.colorSynced;

      if (hasExtended) {
        genes.push([
          gene,
          viz.color,
          viz.localScale,
          isVisible,
          viz.showAssigned,
          viz.showUnassigned,
          viz.unassignedColor,
          viz.unassignedLocalScale,
          viz.colorSynced,
        ]);
      } else {
        genes.push([gene, viz.color, viz.localScale, isVisible]);
      }
    });

    if (genes.length > 0) obj.genes = genes;
  }

  if (state.globalScale !== 1.0) obj.gs = state.globalScale;
  if (state.viewMode !== "2D") obj.vm = state.viewMode;
  if (!state.showAssigned) obj.sa = false;
  if (!state.showUnassigned) obj.su = false;

  if (Object.keys(obj).length === 0) return null;

  return toBase64Url(obj);
}

export function decodeSMVizState(encoded: string): SMVizUrlState | null {
  const obj = fromBase64Url<SMVizUrlState>(encoded);

  if (!obj || typeof obj !== "object") return null;

  if (obj.genes !== undefined && !Array.isArray(obj.genes)) return null;
  if (obj.gs !== undefined && typeof obj.gs !== "number") return null;
  if (obj.vm !== undefined && obj.vm !== "2D" && obj.vm !== "3D") return null;
  if (obj.sa !== undefined && typeof obj.sa !== "boolean") return null;
  if (obj.su !== undefined && typeof obj.su !== "boolean") return null;

  return obj;
}

// --- Labelled single-molecule viewer ---------------------------------------

/**
 * Compact URL shape for /lm-viewer. Genes carry their palette slot so a shared
 * link reproduces the exact colours the sender saw — slots are handed out in
 * selection order, which a bare list of names would not reconstruct.
 */
export interface LmVizUrlState {
  cb?: "gene" | "domain" | "cell"; // colorBy
  dv?: "domain_anno" | "domain_id"; // which column the domain menu reads
  g?: [string, number][]; // [gene, colourSlot]
  d?: string[]; // domain selection
  c?: string[]; // cell selection
  um?: "grey" | "hidden"; // unselectedMode
  ua?: number; // unselectedAlpha
  gs?: number; // globalScale
  ga?: number; // globalAlpha
  vm?: "2D" | "3D"; // viewMode
}

export function encodeLmVizState(state: {
  colorBy: "gene" | "domain" | "cell";
  domainVariant: "domain_anno" | "domain_id";
  selections: Record<"gene" | "domain" | "cell", Set<string>>;
  geneColorSlots: Map<string, number>;
  unselectedMode: "grey" | "hidden";
  unselectedAlpha: number;
  globalScale: number;
  globalAlpha: number;
  viewMode: "2D" | "3D";
}): string | null {
  const s: LmVizUrlState = {};

  // Defaults are omitted so an untouched viewer keeps a clean URL.
  if (state.colorBy !== "cell") s.cb = state.colorBy;
  if (state.domainVariant !== "domain_anno") s.dv = state.domainVariant;
  if (state.selections.gene.size > 0) {
    s.g = [...state.selections.gene].map((gene) => [
      gene,
      state.geneColorSlots.get(gene) ?? 0,
    ]);
  }
  if (state.selections.domain.size > 0) s.d = [...state.selections.domain];
  if (state.selections.cell.size > 0) s.c = [...state.selections.cell];
  if (state.unselectedMode !== "grey") s.um = state.unselectedMode;
  if (state.unselectedAlpha !== 0.2) s.ua = state.unselectedAlpha;
  if (state.globalScale !== 1) s.gs = state.globalScale;
  if (state.globalAlpha !== 1) s.ga = state.globalAlpha;
  if (state.viewMode !== "3D") s.vm = state.viewMode;

  if (Object.keys(s).length === 0) return null;

  return toBase64Url(s);
}

export function decodeLmVizState(encoded: string): LmVizUrlState | null {
  return fromBase64Url<LmVizUrlState>(encoded);
}
