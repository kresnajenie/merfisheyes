/**
 * The names and spellings particular producers use, from one shared file:
 * docs/specs/conventions.json (schema in docs/specs/schemas/, read by
 * scripts/validate_store.py too). Rules of the format itself stay in code;
 * see docs/specs/conventions.md for which is which and how to add one.
 */
import data from "@/docs/specs/conventions.json";

export interface Conventions {
  version: 1;
  tables: {
    curatedDefaultTable: { attribute: string; key: string };
    fallbackTableName: string;
    featureIndexColumns: string[];
    defaultColumnPriority: string[];
    notDefaultColumns: string[];
    positionColumns: string[];
  };
  shapes: {
    cellOutlineElements: string[];
    pointSiblingSuffix: string;
    squareBinPattern: string;
  };
  catalog: {
    bucketRegions: Record<string, string>;
    species: Record<string, string>;
    platforms: Record<string, string>;
  };
}

export const CONVENTIONS = data as Conventions;

/** Columns never chosen as the default colouring. */
export const NOT_DEFAULT_COLUMNS = new Set(
  CONVENTIONS.tables.notDefaultColumns,
);

/**
 * obs columns holding a cell's position: never offered as a colour, by the
 * browser preview or the SpatialData viewer.
 */
export const POSITION_COLUMNS = new Set(CONVENTIONS.tables.positionColumns);

/** The Point element beside a region that holds each cell's centre. */
export function pointSiblingOf(region: string): string {
  return `${region}${CONVENTIONS.shapes.pointSiblingSuffix}`;
}

/** Whether a polygon region's name says it is a grid of square bins. */
export function isSquareBinRegion(region: string): boolean {
  return region
    .toLowerCase()
    .includes(CONVENTIONS.shapes.squareBinPattern.toLowerCase());
}
