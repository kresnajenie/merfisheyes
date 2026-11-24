import type { CellVizUrlState } from "@/lib/utils/url-visualization-state";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { sizeLocked } from "@/lib/config/technology-presets";
import {
  moleculesTarget,
  primaryTarget,
  sceneFor,
  type CellsTarget,
  type MoleculesTarget,
} from "@/lib/agent/view-targets";
import { applyCellVizState } from "@/lib/hooks/useUrlVizSync";
import { colormapRgb } from "@/lib/utils/colormaps";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";
import { getClusterValue } from "@/lib/StandardizedDataset";
import { clusterColumnFrom } from "@/lib/spatialdata/SpatialDataTableAdapter";

/**
 * What an agent can do to an open view (docs/plans/MCP.md pillar D3): read
 * it, change it, and take a picture of it. Run in the page by
 * components/agent-session.tsx for the commands the package sends; every
 * change goes through the viewer's own store actions and `applyCellVizState`
 * (the code that applies a link's settings), never a second implementation.
 */

export class ViewCommandError extends Error {}

const frames = (n = 2) =>
  new Promise<void>((done) => {
    const step = (left: number) =>
      left <= 0 ? done() : requestAnimationFrame(() => step(left - 1));

    step(n);
  });

function target(): CellsTarget & { data: StandardizedDataset } {
  const t = primaryTarget();
  const data = t?.dataset();

  if (!t) throw new ViewCommandError("No dataset is open in this tab yet.");
  if (!data)
    throw new ViewCommandError(
      "The data is still loading; try again in a moment.",
    );

  return { ...t, data };
}

const cellCount = (d: StandardizedDataset) => {
  const c = d.spatial.coordinates as Float32Array | number[][];

  return c instanceof Float32Array ? c.length / d.spatial.dimensions : c.length;
};

function columnNames(d: StandardizedDataset): string[] {
  return d.allClusterColumnNames?.length
    ? d.allClusterColumnNames
    : (d.clusters ?? []).map((c) => c.column);
}

type VizState = ReturnType<CellsTarget["store"]>;

/** A column's type as the viewer treats it: the person's override ("treat
 * as categorical"), else the loaded column's, else the store's listing. */
function columnType(
  d: StandardizedDataset,
  s: VizState,
  column: string,
): string {
  return (
    s.columnTypeOverrides?.[column] ??
    d.clusters?.find((k) => k.column === column)?.type ??
    d.allClusterColumnTypes?.[column] ??
    "unknown"
  );
}

/** The values of a column the viewer treats as categorical (loaded). */
function columnValues(
  d: StandardizedDataset,
  s: VizState,
  column: string | null,
) {
  const c = column ? d.clusters?.find((k) => k.column === column) : null;

  if (!c || columnType(d, s, column!) === "numerical") return null;

  return c.uniqueValues ?? [...new Set(c.values.map(String))].sort();
}

/**
 * Which values of the colour column are drawn, as the scene decides
 * (components/three-scene.tsx): nothing picked out — all of them; values
 * picked out — those minus the ones whose eye is off, the rest drawn at
 * `greyedOutAlpha` (0 = not drawn at all).
 */
function drawn(d: StandardizedDataset, s: VizState) {
  const values = columnValues(d, s, s.selectedColumn);

  if (!values) return null;
  if (s.selectedCelltypes.size === 0) return { values, others: null };
  const shown = values.filter(
    (v) => s.selectedCelltypes.has(v) && !s.hiddenCelltypes.has(v),
  );
  const alpha = (s as unknown as Record<string, number>).greyedOutAlpha ?? 1;

  return { values: shown, others: alpha === 0 ? "hidden" : "grey" };
}

/** How many cells are drawn in colour (all of them unless values are
 * picked out). Counted from the cells' own values. */
function visibleCells(d: StandardizedDataset, s: VizState): number {
  const n = cellCount(d);
  const pick = drawn(d, s);
  const c = s.selectedColumn
    ? d.clusters?.find((k) => k.column === s.selectedColumn)
    : null;

  if (!pick || !pick.others || !c) return n;
  const keep = new Set(pick.values);
  let count = 0;

  for (let i = 0; i < n; i++)
    if (keep.has(String(getClusterValue(c as never, i)))) count++;

  return count;
}

/**
 * Move the camera so the cells drawn in colour fill the view (all cells when
 * nothing is picked out), keeping the viewing direction. Returns how to put
 * the camera back.
 */
function fitToVisible(t: CellsTarget, data: StandardizedDataset): () => void {
  const scene = sceneFor(t.panel);
  const camera = scene?.camera();
  const controls = scene?.controls();
  const pos = scene?.positions();

  if (!scene || !camera || !controls || !pos)
    throw new ViewCommandError("The scene is not ready to be framed yet.");
  const s = t.store();
  const pick = drawn(data, s);
  const c = s.selectedColumn
    ? data.clusters?.find((k) => k.column === s.selectedColumn)
    : null;
  const keep = pick?.others && c ? new Set(pick.values) : null;
  const n = Math.min(cellCount(data), pos.length / 3);
  const lo = [Infinity, Infinity, Infinity];
  const hi = [-Infinity, -Infinity, -Infinity];

  for (let i = 0; i < n; i++) {
    if (keep && !keep.has(String(getClusterValue(c as never, i)))) continue;
    for (let k = 0; k < 3; k++) {
      const v = pos[i * 3 + k];

      if (v < lo[k]) lo[k] = v;
      if (v > hi[k]) hi[k] = v;
    }
  }
  if (!Number.isFinite(lo[0]))
    throw new ViewCommandError("No cells are drawn to frame.");
  const center = lo.map((v, k) => (v + hi[k]) / 2);
  const before = {
    position: camera.position.clone(),
    target: controls.target.clone(),
  };
  const fov = ((camera.fov ?? 50) * Math.PI) / 180;
  const aspect = camera.aspect ?? 1;
  const w = Math.max(hi[0] - lo[0], 1e-6);
  const h = Math.max(hi[1] - lo[1], 1e-6);
  // the height that fits both extents, with a margin
  const fitH = Math.max(h, w / aspect) * 1.12;
  const distance = fitH / 2 / Math.tan(fov / 2) + (hi[2] - lo[2]) / 2;
  const dir = {
    x: camera.position.x - controls.target.x,
    y: camera.position.y - controls.target.y,
    z: camera.position.z - controls.target.z,
  };
  const len = Math.hypot(dir.x, dir.y, dir.z) || 1;

  controls.target.set(center[0], center[1], center[2]);
  camera.position.copy({
    x: center[0] + (dir.x / len) * distance,
    y: center[1] + (dir.y / len) * distance,
    z: center[2] + (dir.z / len) * distance,
  });
  controls.update();

  return () => {
    camera.position.copy(before.position);
    controls.target.copy(before.target);
    controls.update();
  };
}

/** The molecules on show: genes drawn (with colours), opacity, size. */
function moleculesState(t: MoleculesTarget) {
  const d = t.dataset();
  const s = t.store();

  if (!d) return { loading: true };

  return {
    over_cells: t.overlay,
    genes_available: d.uniqueGenes.length,
    shown: [...s.selectedGenes.values()].map((g) => ({
      gene: g.gene,
      color: g.color,
    })),
    in_legend_but_hidden: [...s.selectedGenesLegend].filter(
      (g) => !s.selectedGenes.has(g),
    ),
    opacity: s.moleculeOpacity,
    size: s.globalScale,
    view: s.viewMode,
    layer_visible: s.smLayerVisible,
    examples: d.uniqueGenes.slice(0, 30),
  };
}

export interface MoleculesArgs {
  /** Exactly these genes drawn (others removed). */
  show?: string[];
  add?: string[];
  remove?: string[];
  colors?: Record<string, string>;
  opacity?: number;
  size?: number;
}

/** Change the molecules on show through the molecule store's own actions;
 * gene names matched to the dataset's without regard to case. */
async function moleculesSet(t: MoleculesTarget, a: MoleculesArgs) {
  // Molecules over cells load after the cells: wait for them a while.
  for (let i = 0; i < 120 && !t.dataset(); i++)
    await new Promise((r) => setTimeout(r, 500));
  const d = t.dataset();

  if (!d)
    throw new ViewCommandError(
      "The molecules are still loading; try again in a moment.",
    );
  const fix = (g: string) => {
    const hit = d.uniqueGenes.find((x) => x.toLowerCase() === g.toLowerCase());

    if (!hit) {
      const near = d.uniqueGenes
        .filter((x) => x.toLowerCase().includes(g.toLowerCase().slice(0, 3)))
        .slice(0, 8);

      throw new ViewCommandError(
        `No molecules of "${g}" in this dataset.${near.length ? ` Near: ${near.join(", ")}.` : ""}`,
      );
    }

    return hit;
  };
  const show = a.show?.map(fix);
  const add = a.add?.map(fix) ?? [];
  const remove = a.remove?.map(fix) ?? [];
  const s = t.store();

  if (show) {
    for (const g of [...s.selectedGenesLegend])
      if (!show.includes(g)) s.removeGene(g);
    await frames();
    for (const g of show) {
      const now = t.store();

      if (now.selectedGenes.has(g)) continue;
      if (now.selectedGenesLegend.has(g)) now.toggleGeneVisibility(g);
      else now.addGene(g, a.colors?.[g]);
    }
  }
  for (const g of add)
    if (!t.store().selectedGenes.has(g)) t.store().addGene(g, a.colors?.[g]);
  for (const g of remove) t.store().removeGene(g);
  await frames();
  for (const [g, color] of Object.entries(a.colors ?? {}))
    t.store().setGeneColor(fix(g), color);
  if (a.opacity !== undefined) t.store().setMoleculeOpacity(a.opacity);
  if (a.size !== undefined) t.store().setGlobalScale(a.size);
}

/** What is shown: dataset, colouring, hidden values, columns and values. */
export function viewState() {
  const mol = moleculesTarget();
  const cells = primaryTarget();

  if (!cells && mol) {
    const d = mol.dataset();

    return {
      viewer: "molecules",
      dataset: d?.name ?? null,
      opened: new URLSearchParams(window.location.search).get("url"),
      molecules: moleculesState(mol),
    };
  }
  const { data, store } = target();
  const s = store();

  return {
    ...(mol ? { molecules: moleculesState(mol) } : {}),
    viewer: "cells",
    dataset: data.name,
    // the data this tab shows (the address the package serves it at)
    opened: new URLSearchParams(window.location.search).get("url"),
    cells: cellCount(data),
    genes: data.genes.length,
    color_by: s.mode.includes("gene")
      ? { gene: s.selectedGene }
      : { column: s.selectedColumn },
    // what the scene draws, counted from the cells (not what was asked)
    shown_values: drawn(data, s)?.others ? drawn(data, s)!.values : "all",
    other_values: drawn(data, s)?.others ?? null,
    visible_cells: visibleCells(data, s),
    dot_size: s.sizeScale,
    // spots and bins drawn at their real size (the technology's preset)
    ...(sizeLocked(data) ? { dot_size_locked: true } : {}),
    view: s.viewMode,
    columns: columnNames(data).map((name) => {
      const values = columnValues(data, s, name);

      return {
        name,
        type: columnType(data, s, name),
        loaded: !!data.clusters?.some((c) => c.column === name),
        ...(values
          ? { n_values: values.length, values: values.slice(0, 60) }
          : {}),
      };
    }),
  };
}

export interface ViewSetArgs {
  gene?: string;
  column?: string;
  show_only?: string[];
  hide?: string[];
  show_all?: boolean;
  /** Treat the column as categories or as numbers (the viewer's switch). */
  treat_as?: "categorical" | "numerical";
  dot_size?: number;
  view?: "2D" | "3D";
  reset_camera?: boolean;
  /** Zoom the view to the cells drawn in colour. */
  fit?: boolean;
  /** The molecules on show (the molecule viewer, or molecules over cells). */
  molecules?: MoleculesArgs;
}

/** Change the view; answers with the state afterwards. Each change goes
 * through the store's own actions; the colour column is applied only when
 * it changes (changing it resets what is picked out). */
export async function viewSet(a: ViewSetArgs) {
  const mol = moleculesTarget();

  if (a.molecules) {
    if (!mol)
      throw new ViewCommandError("No molecules are drawn in this view.");
    await moleculesSet(mol, a.molecules);
  }
  if (!primaryTarget()) {
    // The molecule viewer: no cells to change.
    if (a.view) mol?.store().setViewMode(a.view);
    if (a.dot_size !== undefined) mol?.store().setGlobalScale(a.dot_size);
    await frames(3);

    return viewState();
  }
  const { data, store } = target();
  const s = store();
  let gene: string | null = null;

  if (a.gene) {
    gene =
      data.genes.find((x) => x.toLowerCase() === a.gene!.toLowerCase()) ?? null;
    if (!gene) {
      const near = data.genes
        .filter((x) =>
          x.toLowerCase().includes(a.gene!.toLowerCase().slice(0, 3)),
        )
        .slice(0, 8);

      throw new ViewCommandError(
        `This dataset does not measure "${a.gene}".${near.length ? ` Near: ${near.join(", ")}.` : ""}`,
      );
    }
  }
  if (a.column && !columnNames(data).includes(a.column))
    throw new ViewCommandError(
      `No column "${a.column}". Columns: ${columnNames(data).join(", ")}.`,
    );

  // A new colour column: through applyCellVizState, which loads it.
  if (a.column && a.column !== s.selectedColumn) {
    const decoded: CellVizUrlState = {
      c: a.column,
      m: gene ? ["gene"] : ["celltype"],
    };

    if (gene) decoded.g = gene;
    await applyCellVizState(decoded, s, data);
    await frames(3);
  } else if (gene) {
    s.setSelectedGene(gene);
    s.setMode(["gene"]);
  } else if (a.column) {
    s.setMode(["celltype"]);
  }
  // Spots and bins are drawn at their real size when the technology's
  // preset locks it: the size slider changes nothing, so neither does this.
  const sizeRefused = a.dot_size !== undefined && sizeLocked(data);

  if (a.dot_size !== undefined && !sizeRefused)
    store().setSizeScale(a.dot_size);
  if (a.view) store().setViewMode(a.view);

  const wantsValues = !!(a.show_only || a.hide || a.show_all);

  // Values are shown and hidden by category: a numeric column (fov ids, …)
  // is switched to "treat as categorical" first, as the viewer's switch does.
  const want =
    a.treat_as ?? (a.show_only || a.hide ? "categorical" : undefined);

  if (want) {
    await frames();
    const now = store();
    const column = now.selectedColumn;
    const has = column ? columnType(data, now, column) : null;

    if (column && has !== want) {
      now.toggleColumnType(
        column,
        has === "categorical" ? "categorical" : "numerical",
      );
      await frames(3);
    }
  }

  if (wantsValues) {
    await frames();
    const now = store();
    const column = now.selectedColumn;
    const values = columnValues(data, now, column);

    if (!values)
      throw new ViewCommandError(
        `Showing or hiding values needs a categorical colour column; "${column}" is not one (or is not loaded).`,
      );
    const asked = [...(a.show_only ?? []), ...(a.hide ?? [])];
    const unknown = asked.filter((v) => !values.includes(v));

    if (unknown.length)
      throw new ViewCommandError(
        `"${column}" has no value ${unknown.map((v) => `"${v}"`).join(", ")}. Its values: ${values.slice(0, 40).join(", ")}${values.length > 40 ? ", …" : ""}.`,
      );
    const set = (key: string, value: number) => now.setAdvancedViz(key, value);

    if (a.show_all) {
      // nothing picked out: every value drawn as before
      now.setCelltypes(new Set());
      set("greyedOutAlpha", VISUALIZATION_CONFIG.GREYED_OUT_ALPHA as number);
      set(
        "selectedSizeMultiplier",
        VISUALIZATION_CONFIG.SELECTED_SIZE_MULTIPLIER as number,
      );
    } else {
      // Pick out the values to keep (the panel's checkboxes) and draw the
      // rest not at all; picked-out cells keep their size.
      const current = drawn(data, now);
      const keep = new Set(
        a.show_only ?? (current?.others ? current.values : values),
      );

      for (const v of a.hide ?? []) keep.delete(v);
      if (keep.size === 0)
        throw new ViewCommandError(
          "That would hide every value; use show_all to show them again.",
        );
      now.setCelltypes(keep);
      set("greyedOutAlpha", 0);
      set("selectedSizeMultiplier", 1);
    }
  }
  if (a.reset_camera) store().resetCamera();
  if (a.fit) {
    await frames(2);
    fitToVisible(target(), data);
  }
  await frames(3);

  const after = viewState();

  return sizeRefused
    ? {
        ...after,
        note: "dot_size not applied: this dataset's spots or bins are drawn at their real size (its technology's preset locks the size).",
      }
    : after;
}

/**
 * A picture of the view: the scene as drawn, with a legend panel beside it
 * (the visible categories and their colours, or the gene and its colour
 * scale). PNG, base64.
 */
export async function viewSnapshot(
  a: { width?: number; frame?: "visible" | "current" } = {},
) {
  const cellsT = primaryTarget();
  const mol = moleculesTarget();

  if (!cellsT && !mol)
    throw new ViewCommandError("No dataset is open in this tab yet.");
  const data = cellsT?.dataset() ?? null;
  const s = cellsT?.store() ?? null;
  const canvas = document.querySelector(
    cellsT
      ? '[data-testid="sc-scene-canvas"] canvas'
      : '[data-testid="sm-scene-canvas"] canvas',
  ) as HTMLCanvasElement | null;

  if (!canvas || (cellsT && !data))
    throw new ViewCommandError("No scene is drawn in this tab yet.");
  // Framed on the cells drawn (the default), then the person's camera back.
  const restore =
    cellsT && data && a.frame !== "current" ? fitToVisible(cellsT, data) : null;

  await frames(4);
  const scene = await createImageBitmap(
    await new Promise<Blob>((done, fail) =>
      canvas.toBlob(
        (b) => (b ? done(b) : fail(new Error("empty canvas"))),
        "image/png",
      ),
    ),
  );

  restore?.();
  const width = Math.min(Math.max(a.width ?? scene.width, 400), 3000);
  const k = width / scene.width;
  const sceneH = Math.round(scene.height * k);
  const u = Math.max(1, width / 1200); // legend scale
  const legendW = Math.round(260 * u);
  const gene = s && s.mode.includes("gene") ? s.selectedGene : null;
  const column = s?.selectedColumn ?? null;
  const values = gene || !s || !data ? null : columnValues(data, s, column);
  // a gene, or a numeric column: a colour scale; a categorical one: swatches
  const scale = s ? (gene ?? (values ? null : column)) : null;
  const cluster =
    column && data ? data.clusters?.find((c) => c.column === column) : null;
  const pick = gene || !s || !data ? null : drawn(data, s);
  const shown = pick ? pick.values : [];
  const rows = Math.min(shown.length, 40);
  const molGenes = mol ? [...mol.store().selectedGenes.values()] : [];
  const legendH = Math.round(
    (60 +
      (scale ? 70 : rows * 22 + (shown.length > rows ? 22 : 0)) +
      (molGenes.length ? 40 + Math.min(molGenes.length, 20) * 22 : 0)) *
      u,
  );
  const out = document.createElement("canvas");

  out.width = width + legendW;
  out.height = Math.max(sceneH, legendH);
  const ctx = out.getContext("2d")!;

  ctx.fillStyle = "#000";
  ctx.fillRect(0, 0, out.width, out.height);
  ctx.drawImage(scene, 0, 0, width, sceneH);
  const x = width + 16 * u;
  let y = 28 * u;
  const swatch = (label: string, color: string) => {
    y += 22 * u;
    ctx.fillStyle = color;
    ctx.fillRect(x, y - 11 * u, 12 * u, 12 * u);
    ctx.fillStyle = "#ddd";
    ctx.fillText(label.slice(0, 28), x + 20 * u, y);
  };

  ctx.fillStyle = "#fff";
  ctx.font = `600 ${15 * u}px sans-serif`;
  ctx.fillText((data?.name ?? mol?.dataset()?.name ?? "").slice(0, 32), x, y);
  ctx.font = `${13 * u}px sans-serif`;
  const legend: { value: string; color: string }[] = [];

  if (s) {
    y += 22 * u;
    ctx.fillStyle = "#bbb";
    ctx.fillText(gene ? `gene: ${gene}` : `colour: ${column ?? "—"}`, x, y);
    y += 16 * u;
    if (scale) {
      const bw = legendW - 40 * u;

      for (let i = 0; i < bw; i++) {
        const [r, g, b] = colormapRgb(s.colormap, i / bw);

        ctx.fillStyle = `rgb(${r * 255},${g * 255},${b * 255})`;
        ctx.fillRect(x + i, y, 1, 14 * u);
      }
      ctx.fillStyle = "#bbb";
      ctx.fillText("low", x, y + 32 * u);
      ctx.fillText("high", x + bw - 28 * u, y + 32 * u);
      y += 40 * u;
    } else {
      for (const v of shown.slice(0, rows)) {
        const color = s.colorPalette?.[v] ?? cluster?.palette?.[v] ?? "#999";

        legend.push({ value: v, color });
        swatch(v, color);
      }
      if (shown.length > rows) {
        y += 22 * u;
        ctx.fillStyle = "#888";
        ctx.fillText(`+ ${shown.length - rows} more`, x, y);
      }
    }
  }
  if (molGenes.length) {
    y += 30 * u;
    ctx.fillStyle = "#bbb";
    ctx.fillText("molecules", x, y);
    for (const g of molGenes.slice(0, 20)) swatch(g.gene, g.color);
  }
  const dataUrl = out.toDataURL("image/png");

  return {
    png_base64: dataUrl.slice(dataUrl.indexOf(",") + 1),
    width: out.width,
    height: out.height,
    ...(gene ? { gene, colormap: s!.colormap } : s ? { column, legend } : {}),
    ...(molGenes.length
      ? { molecules: molGenes.map((g) => ({ gene: g.gene, color: g.color })) }
      : {}),
  };
}

export interface AddColumnArgs {
  /** Per-cell ids, matched to the dataset's own; absent: values in cell order. */
  ids?: string[];
  columns: { name: string; values: string[] }[];
  /** Replace a column of the same name instead of refusing. */
  replace?: boolean;
  /** Colour by the first added column. */
  color_by?: boolean;
}

/**
 * Add columns computed elsewhere (an agent's pseudotime, labels, scores —
 * pillar D4) to the open view, without touching the store: joined to the
 * cells by id (the obs index) or taken in cell order; a cell with no value
 * gets "" (a missing label). Built like a column read from a store
 * (`clusterColumnFrom`): the same type detection and palette.
 */
export async function viewAddColumn(a: AddColumnArgs) {
  const { data, store } = target();
  const n = cellCount(data);
  const existing = new Set(columnNames(data));

  for (const c of a.columns)
    if (existing.has(c.name) && !a.replace)
      throw new ViewCommandError(
        `The view already has a column "${c.name}"; pass replace=true to replace it, or use another name.`,
      );
  let order: number[];
  let unmatched = 0;

  if (a.ids) {
    const adapter = data.adapter as
      | { loadCellIds?: () => Promise<string[] | null> }
      | undefined;
    const cellIds = await adapter?.loadCellIds?.();

    if (!cellIds)
      throw new ViewCommandError(
        "This dataset's cells have no ids the viewer can read: give the values in cell order instead (one per cell, no id column).",
      );
    const at = new Map(cellIds.map((id, i) => [id, i]));

    order = a.ids.map((id) => {
      const i = at.get(id);

      if (i === undefined) unmatched++;

      return i ?? -1;
    });
    if (unmatched === a.ids.length)
      throw new ViewCommandError(
        `None of the ids match this dataset's cells. Yours look like ${a.ids
          .slice(0, 3)
          .map((v) => `"${v}"`)
          .join(", ")}; the dataset's look like ${cellIds
          .slice(0, 3)
          .map((v) => `"${v}"`)
          .join(", ")}.`,
      );
  } else {
    for (const c of a.columns)
      if (c.values.length !== n)
        throw new ViewCommandError(
          `"${c.name}" has ${c.values.length} values but the dataset has ${n} cells; give one value per cell in cell order, or add an id column.`,
        );
    order = Array.from({ length: n }, (_, i) => i);
  }

  const built = a.columns.map((c) => {
    const full = new Array<string>(n).fill("");

    for (let r = 0; r < c.values.length; r++)
      if (order[r] >= 0) full[order[r]] = String(c.values[r] ?? "");

    return clusterColumnFrom(c.name, full);
  });

  if (a.replace && data.clusters)
    data.clusters = data.clusters.filter(
      (k) => !a.columns.some((c) => c.name === k.column),
    );
  data.addClusters(built as never);
  for (const b of built) {
    if (!data.allClusterColumnNames.includes(b.column))
      data.allClusterColumnNames.push(b.column);
    data.allClusterColumnTypes[b.column] = b.type;
  }
  store().incrementClusterVersion();
  await frames(2);
  if (a.color_by) await viewSet({ column: a.columns[0].name });

  return {
    added: built.map((b) => ({
      name: b.column,
      type: b.type,
      n_values: b.uniqueValues?.length ?? 0,
    })),
    cells: n,
    matched: a.ids ? a.ids.length - unmatched : n,
    unmatched,
    cells_without_a_value: n - (a.ids ? a.ids.length - unmatched : n),
  };
}

/** Run one command the package sent. */
export async function runViewCommand(command: string, args: unknown) {
  switch (command) {
    case "state":
      return viewState();
    case "set":
      return viewSet((args ?? {}) as ViewSetArgs);
    case "snapshot":
      return viewSnapshot(
        (args ?? {}) as { width?: number; frame?: "visible" | "current" },
      );
    case "add_column":
      return viewAddColumn(args as AddColumnArgs);
    default:
      throw new ViewCommandError(`Unknown command "${command}".`);
  }
}
