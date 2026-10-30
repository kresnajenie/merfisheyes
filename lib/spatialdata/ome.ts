/**
 * Where an OME-Zarr group keeps its NGFF metadata
 * (docs/specs/image-profile.md).
 *
 * OME-NGFF 0.5 nests everything under `attributes.ome`, and that is the only
 * layout we write. Some writers put 0.4-style attributes (`multiscales` and
 * `omero` at the top of `attributes`) inside a Zarr v3 container — ome-zarr-py
 * did so before it grew a guard, and still reads the result, warning
 * `version mismatch: detected: FormatV04`. Reading both costs one fallback.
 *
 * Every reader goes through here, so a store can never load halfway: before
 * this, the image reader refused such a group while the coordinate reader
 * accepted it. `tests/unit/ome-layout.test.ts` holds the readers to it.
 */

/** NGFF attributes of a group, whichever layout the writer used. */
export function omeAttrs(
  attrs: Record<string, any> | null | undefined,
): Record<string, any> {
  return (attrs?.ome ?? attrs ?? {}) as Record<string, any>;
}

/** An image group's `multiscales[0]`, or undefined when it has none. */
export function imageMultiscale(
  attrs: Record<string, any> | null | undefined,
): Record<string, any> | undefined {
  return omeAttrs(attrs).multiscales?.[0];
}

/** The axis names of a multiscale, e.g. `"cyx"`. */
export function axisNames(multiscale: Record<string, any> | undefined): string {
  return (multiscale?.axes ?? [])
    .map((a: { name?: string } | string) =>
      typeof a === "string" ? a : a.name,
    )
    .join("");
}
