/**
 * A uniform grid over point positions, for hover picking in constant time
 * instead of a pass over every point (docs/specs/RENDERING.md). Built once
 * per set of positions; a query looks at the cells around the pointer only.
 */
/** Cells across the grid, at most. */
const GRID_MAX = 256;

export class PickingGrid {
  private cells = new Map<number, number[]>();
  private readonly minX: number;
  private readonly minY: number;
  private readonly nx: number;
  private readonly ny: number;
  readonly cellSize: number;

  /**
   * `positions` is interleaved xyz; the grid is over x and y, with cells
   * sized so that about `perCell` points share one on average.
   */
  constructor(
    private positions: ArrayLike<number>,
    private count: number,
    perCell = 16,
  ) {
    let minX = Infinity;
    let minY = Infinity;
    let maxX = -Infinity;
    let maxY = -Infinity;

    for (let i = 0; i < count; i++) {
      const x = positions[i * 3];
      const y = positions[i * 3 + 1];

      if (x < minX) minX = x;
      if (x > maxX) maxX = x;
      if (y < minY) minY = y;
      if (y > maxY) maxY = y;
    }
    if (count === 0) {
      minX = minY = 0;
      maxX = maxY = 1;
    }
    const extent = Math.max(maxX - minX, maxY - minY);
    const area = Math.max(maxX - minX, 1e-9) * Math.max(maxY - minY, 1e-9);

    // At most GRID_MAX cells across, so a query walks a bounded number of
    // cells however the points lie (on a line, on one spot)
    this.cellSize = Math.max(
      Math.sqrt((area * perCell) / Math.max(count, 1)),
      extent / GRID_MAX,
      1e-9,
    );
    if (extent === 0) this.cellSize = 1;
    this.minX = minX;
    this.minY = minY;
    this.nx = Math.floor((maxX - minX) / this.cellSize) + 1;
    this.ny = Math.floor((maxY - minY) / this.cellSize) + 1;
    for (let i = 0; i < count; i++) {
      const key = this.keyOf(positions[i * 3], positions[i * 3 + 1]);
      const cell = this.cells.get(key);

      if (cell) cell.push(i);
      else this.cells.set(key, [i]);
    }
  }

  private keyOf(x: number, y: number): number {
    const cx = Math.floor((x - this.minX) / this.cellSize);
    const cy = Math.floor((y - this.minY) / this.cellSize);

    return cy * this.nx + cx;
  }

  /**
   * The index of the point nearest to (x, y) within `radius`, or null.
   * `distance`, when given, measures between a point and the query (to pick
   * in screen space rather than in the plane); it defaults to the planar
   * distance.
   */
  nearest(
    x: number,
    y: number,
    radius: number,
    distance?: (index: number) => number,
  ): number | null {
    // No further than the grid itself
    const reach = Math.min(Math.ceil(radius / this.cellSize), GRID_MAX);
    const cx = Math.floor((x - this.minX) / this.cellSize);
    const cy = Math.floor((y - this.minY) / this.cellSize);
    let best: number | null = null;
    let bestD = radius;

    for (let dy = -reach; dy <= reach; dy++) {
      for (let dx = -reach; dx <= reach; dx++) {
        const cell = this.cells.get((cy + dy) * this.nx + cx + dx);

        if (!cell) continue;
        for (const i of cell) {
          const d = distance
            ? distance(i)
            : Math.hypot(
                this.positions[i * 3] - x,
                this.positions[i * 3 + 1] - y,
              );

          if (d <= bestD) {
            bestD = d;
            best = i;
          }
        }
      }
    }

    return best;
  }
}
