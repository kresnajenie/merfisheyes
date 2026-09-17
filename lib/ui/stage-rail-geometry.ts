/**
 * Geometry of the labelled-molecule viewer's stage rail.
 *
 * Its own module, separate from the component, because the global footer needs
 * the height to sit above the rail — and `components/stage-rail.tsx` imports
 * three.js, which would otherwise be pulled into every page's bundle by a
 * footer that renders on all of them.
 */

/** Side of one square preview tile. */
export const TILE = 96;

/** The stage caption under each tile: a 12px line plus its 2px top margin. */
export const LABEL_H = 14;

/** Vertical padding of the bar, top + bottom. */
export const PAD_Y = 12;

/** The stage heading over each run of tiles, plus the gap under it. */
export const STAGE_LABEL_H = 16;

/** Height of the rail; chrome above it is offset by this. */
export const STAGE_RAIL_HEIGHT = STAGE_LABEL_H + TILE + LABEL_H + PAD_Y;
