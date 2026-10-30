import { readdirSync, readFileSync, statSync } from "node:fs";
import path from "node:path";

import { afterEach, describe, expect, it, vi } from "vitest";

import { elementTransforms } from "@/lib/spatialdata/coordinates";
import { axisNames, imageMultiscale, omeAttrs } from "@/lib/spatialdata/ome";
import { resolveZarrUrl } from "@/lib/spatialdata/resolveUrl";
import { isOmeZarrImage } from "@/lib/upload/uploadImageToS3";

/**
 * An image's NGFF metadata is read wherever the writer put it
 * (docs/specs/image-profile.md): under `attributes.ome` as 0.5 says, or at
 * the top of `attributes` as ome-zarr-py wrote before it guarded against
 * 0.4-in-v3 — the layout the CoQ1s stitcher produced.
 *
 * Every reader must agree. Before this, `image.ts` refused such a group
 * while `coordinates.ts` accepted it, so a store could load halfway: the
 * image failed to draw but its placement was read. These tests run each
 * reader over both layouts, and the last one stops a new reader reaching
 * past the helper.
 */

const ROOT = path.resolve(__dirname, "../..");

const multiscales = [
  {
    axes: [
      { name: "c", type: "channel" },
      { name: "y", type: "space" },
      { name: "x", type: "space" },
    ],
    datasets: [
      { path: "0", coordinateTransformations: [{ type: "scale", scale: [1, 0.65, 0.65] }] },
    ],
    coordinateTransformations: [{ type: "scale", scale: [1, 2, 2] }],
  },
];
const omero = { channels: [{ label: "DAPI" }] };

/** The same image, written the two ways. */
const LAYOUTS = {
  "0.5, under `ome`": { ome: { version: "0.5", multiscales, omero } },
  "0.4-style, at the top (CoQ1s)": { multiscales, omero },
} as const;
const names = Object.keys(LAYOUTS) as (keyof typeof LAYOUTS)[];

describe("both layouts read the same", () => {
  it.each(names)("%s: the helper finds multiscales and omero", (layout) => {
    const attrs = LAYOUTS[layout];

    expect(imageMultiscale(attrs)).toBe(multiscales[0]);
    expect(axisNames(imageMultiscale(attrs))).toBe("cyx");
    expect(omeAttrs(attrs).omero).toBe(omero);
  });

  it.each(names)("%s: it is detected as an image to draw", (layout) => {
    expect(isOmeZarrImage({ attributes: LAYOUTS[layout] })).toBe(true);
  });

  it.each(names)("%s: its transforms and axes are read", (layout) => {
    expect(elementTransforms("images", LAYOUTS[layout])).toEqual({
      transforms: multiscales[0].coordinateTransformations,
      axes: ["c", "y", "x"],
    });
  });

  it.each(names)("%s: a URL resolves to an image", async (layout) => {
    vi.stubGlobal(
      "fetch",
      vi.fn(async (url: string) =>
        url.endsWith("/zarr.json")
          ? new Response(
              JSON.stringify({ zarr_format: 3, attributes: LAYOUTS[layout] }),
            )
          : new Response("", { status: 404 }),
      ),
    );

    expect(await resolveZarrUrl("https://bucket.example/img.zarr")).toMatchObject(
      { kind: "image", url: "https://bucket.example/img.zarr" },
    );
  });
});

afterEach(() => vi.unstubAllGlobals());

describe("what is not an image stays not an image", () => {
  it.each([
    ["a table", { spatialdata_attrs: { version: "0.2" } }],
    ["axes z, y, x", { ome: { multiscales: [{ axes: [{ name: "z" }, { name: "y" }, { name: "x" }] }] } }],
    ["nothing", {}],
  ])("%s", (_what, attrs) => {
    expect(isOmeZarrImage({ attributes: attrs })).toBe(false);
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

describe("no reader reaches past the helper", () => {
  // `attrs.ome` / `attributes.ome` in a reader means it sees only the 0.5
  // layout — the split this helper exists to end. Comments may mention it.
  const files = [
    ...sources(path.join(ROOT, "lib")),
    ...sources(path.join(ROOT, "components")),
  ].filter((f) => !f.endsWith(path.join("spatialdata", "ome.ts")));

  it("reads ome metadata only through lib/spatialdata/ome.ts", () => {
    const reach = /\b(attrs|attributes)\s*(\?\.|\.)\s*ome\b/;
    const offenders = files.filter((f) => {
      const code = readFileSync(f, "utf8")
        .replace(/\/\*[\s\S]*?\*\//g, "")
        .replace(/^\s*\/\/.*$/gm, "");

      return reach.test(code);
    });

    expect(offenders.map((f) => path.relative(ROOT, f))).toEqual([]);
  });
});
