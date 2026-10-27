/**
 * The producers' names come from one file, docs/specs/conventions.json,
 * read by the viewer and by scripts/validate_store.py. These tests check
 * the rules read it, and that no reader writes one of its names in code
 * again (the drift the file exists to stop).
 */
import { readFileSync, readdirSync, statSync } from "node:fs";
import path from "node:path";

import { describe, expect, it } from "vitest";

import data from "@/docs/specs/conventions.json";
import {
  CONVENTIONS,
  NOT_DEFAULT_COLUMNS,
  isSquareBinRegion,
  pointSiblingOf,
} from "@/lib/spatialdata/conventions";
import { defaultTableName } from "@/lib/spatialdata/store";
import { selectBestClusterColumnByName } from "@/lib/utils/dataset-utils";

const ROOT = path.resolve(__dirname, "../..");

describe("conventions.json", () => {
  it("is what the module exports", () => {
    expect(CONVENTIONS).toBe(data);
    expect(CONVENTIONS.version).toBe(1);
  });

  it("drives the default table", () => {
    const { attribute, key } = CONVENTIONS.tables.curatedDefaultTable;

    expect(
      defaultTableName({ [attribute]: { [key]: "square_008um" } }, [
        "square_002um",
        "square_008um",
      ]),
    ).toBe("square_008um");
    expect(defaultTableName({}, ["other", CONVENTIONS.tables.fallbackTableName])).toBe(
      CONVENTIONS.tables.fallbackTableName,
    );
  });

  it("drives the default colouring: priority first, positional never", () => {
    const [first] = CONVENTIONS.tables.defaultColumnPriority;

    expect(
      selectBestClusterColumnByName(["region", first], { region: "categorical" }),
    ).toBe(first);
    for (const positional of CONVENTIONS.tables.notDefaultColumns) {
      expect(NOT_DEFAULT_COLUMNS.has(positional)).toBe(true);
      expect(
        selectBestClusterColumnByName([positional], { [positional]: "numerical" }),
      ).toBeNull();
    }
  });

  it("names the Point sibling and square bins", () => {
    expect(pointSiblingOf("shapes/cells_polygons")).toBe(
      `shapes/cells_polygons${CONVENTIONS.shapes.pointSiblingSuffix}`,
    );
    expect(isSquareBinRegion("shapes/sample_SQUARE_008um")).toBe(true);
    expect(isSquareBinRegion("shapes/epi_spots")).toBe(false);
  });
});

/** Every .ts / .tsx file under a folder. */
function sources(dir: string): string[] {
  return readdirSync(dir).flatMap((name) => {
    const p = path.join(dir, name);

    if (statSync(p).isDirectory()) return sources(p);

    return /\.tsx?$/.test(name) ? [p] : [];
  });
}

describe("no reader writes a producer's name in code", () => {
  // Quoted names that belong in conventions.json, not in a reader. Comments
  // may mention them; string literals may not. Plain words ("cells") are
  // used everywhere in the UI, so only the distinctive names are guarded.
  const names = [
    CONVENTIONS.shapes.pointSiblingSuffix,
    CONVENTIONS.tables.curatedDefaultTable.key,
    ...CONVENTIONS.shapes.cellOutlineElements,
    ...CONVENTIONS.tables.notDefaultColumns,
  ].filter((n) => /[_ ]/.test(n));
  const files = [
    ...sources(path.join(ROOT, "lib")),
    ...sources(path.join(ROOT, "components")),
  ].filter((f) => !f.endsWith(path.join("spatialdata", "conventions.ts")));

  it.each(names)("%s", (name) => {
    const quoted = new RegExp(`["'\`]${name.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}["'\`]`);
    const offenders = files.filter((f) => {
      const code = readFileSync(f, "utf8")
        .replace(/\/\*[\s\S]*?\*\//g, "")
        .replace(/^\s*\/\/.*$/gm, "");

      return quoted.test(code);
    });

    expect(offenders.map((f) => path.relative(ROOT, f))).toEqual([]);
  });
});
