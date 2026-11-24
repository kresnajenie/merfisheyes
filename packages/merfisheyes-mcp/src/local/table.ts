/**
 * A CSV or TSV the agent wrote (pillar D4: results to add to the open
 * view): header row, quoted fields with "" escapes, comma or tab.
 */
export function parseTable(text: string): {
  header: string[];
  rows: string[][];
} {
  const firstLine = text.slice(0, text.indexOf("\n") >>> 0);
  const sep = firstLine.includes("\t") && !firstLine.includes(",") ? "\t" : ",";
  const rows: string[][] = [];
  let row: string[] = [];
  let field = "";
  let quoted = false;

  for (let i = 0; i < text.length; i++) {
    const ch = text[i];

    if (quoted) {
      if (ch === '"') {
        if (text[i + 1] === '"') {
          field += '"';
          i++;
        } else quoted = false;
      } else field += ch;
    } else if (ch === '"' && field === "") quoted = true;
    else if (ch === sep) {
      row.push(field);
      field = "";
    } else if (ch === "\n" || ch === "\r") {
      if (ch === "\r" && text[i + 1] === "\n") i++;
      row.push(field);
      field = "";
      if (row.length > 1 || row[0] !== "") rows.push(row);
      row = [];
    } else field += ch;
  }
  if (field !== "" || row.length) {
    row.push(field);
    rows.push(row);
  }
  const [header = [], ...body] = rows;

  return { header, rows: body };
}

const ID_NAMES = [
  "cell_id",
  "cell",
  "cell_ids",
  "obs_names",
  "barcode",
  "index",
  "_index",
  "cells",
];

/** The column that names cells: the one asked for; else an unnamed first
 * column (pandas' to_csv writes the index so) or one with an id-like name;
 * else none (values are then in cell order). */
export function idColumn(header: string[], asked?: string): number | null {
  if (asked) {
    const i = header.indexOf(asked);

    if (i < 0)
      throw new Error(
        `No column "${asked}" in the file; its columns: ${header.join(", ")}.`,
      );

    return i;
  }
  if (header[0] === "" || /^Unnamed: 0$/.test(header[0] ?? "")) return 0;
  const i = header.findIndex((h) => ID_NAMES.includes(h.toLowerCase()));

  return i >= 0 ? i : null;
}
