import type { CellVizUrlState } from "@/lib/utils/url-visualization-state";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { primaryTarget, type CellsTarget } from "@/lib/agent/view-targets";
import { applyCellVizState } from "@/lib/hooks/useUrlVizSync";
import { colormapRgb } from "@/lib/utils/colormaps";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";
import { getClusterValue } from "@/lib/StandardizedDataset";

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

/** What is shown: dataset, colouring, hidden values, columns and values. */
export function viewState() {
  const { data, store } = target();
  const s = store();

  return {
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
}

/** Change the view; answers with the state afterwards. Each change goes
 * through the store's own actions; the colour column is applied only when
 * it changes (changing it resets what is picked out). */
export async function viewSet(a: ViewSetArgs) {
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
  if (a.dot_size !== undefined) store().setSizeScale(a.dot_size);
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
  await frames(3);

  return viewState();
}

/**
 * A picture of the view: the scene as drawn, with a legend panel beside it
 * (the visible categories and their colours, or the gene and its colour
 * scale). PNG, base64.
 */
export async function viewSnapshot(a: { width?: number } = {}) {
  const { data, store } = target();
  const s = store();
  const canvas = document.querySelector(
    '[data-testid="sc-scene-canvas"] canvas',
  ) as HTMLCanvasElement | null;

  if (!canvas) throw new ViewCommandError("No scene is drawn in this tab yet.");
  await frames(2);
  const scene = await createImageBitmap(
    await new Promise<Blob>((done, fail) =>
      canvas.toBlob(
        (b) => (b ? done(b) : fail(new Error("empty canvas"))),
        "image/png",
      ),
    ),
  );
  const width = Math.min(Math.max(a.width ?? scene.width, 400), 3000);
  const k = width / scene.width;
  const sceneH = Math.round(scene.height * k);
  const u = Math.max(1, width / 1200); // legend scale
  const legendW = Math.round(260 * u);
  const gene = s.mode.includes("gene") ? s.selectedGene : null;
  const column = s.selectedColumn;
  const values = gene ? null : columnValues(data, s, column);
  // a gene, or a numeric column: a colour scale; a categorical one: swatches
  const scale = gene ?? (values ? null : column);
  const cluster = column
    ? data.clusters?.find((c) => c.column === column)
    : null;
  const pick = gene ? null : drawn(data, s);
  const shown = pick ? pick.values : [];
  const rows = Math.min(shown.length, 40);
  const legendH = Math.round(
    (60 + (scale ? 70 : rows * 22 + (shown.length > rows ? 22 : 0))) * u,
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

  ctx.fillStyle = "#fff";
  ctx.font = `600 ${15 * u}px sans-serif`;
  ctx.fillText(data.name.slice(0, 32), x, y);
  y += 22 * u;
  ctx.font = `${13 * u}px sans-serif`;
  ctx.fillStyle = "#bbb";
  ctx.fillText(gene ? `gene: ${gene}` : `colour: ${column ?? "—"}`, x, y);
  y += 16 * u;
  const legend: { value: string; color: string }[] = [];

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
  } else {
    for (const v of shown.slice(0, rows)) {
      const color = s.colorPalette?.[v] ?? cluster?.palette?.[v] ?? "#999";

      legend.push({ value: v, color });
      y += 22 * u;
      ctx.fillStyle = color;
      ctx.fillRect(x, y - 11 * u, 12 * u, 12 * u);
      ctx.fillStyle = "#ddd";
      ctx.fillText(v.slice(0, 28), x + 20 * u, y);
    }
    if (shown.length > rows) {
      y += 22 * u;
      ctx.fillStyle = "#888";
      ctx.fillText(`+ ${shown.length - rows} more`, x, y);
    }
  }
  const dataUrl = out.toDataURL("image/png");

  return {
    png_base64: dataUrl.slice(dataUrl.indexOf(",") + 1),
    width: out.width,
    height: out.height,
    ...(gene
      ? { gene, colormap: s.colormap }
      : { column, legend, hidden: [...s.hiddenCelltypes] }),
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
      return viewSnapshot((args ?? {}) as { width?: number });
    default:
      throw new ViewCommandError(`Unknown command "${command}".`);
  }
}
