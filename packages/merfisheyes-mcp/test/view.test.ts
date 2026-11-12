import { expect, it } from "vitest";

import { compactCard } from "../src/tools/catalog.js";
import { toBase64Url, viewState, viewerKind } from "../src/tools/view.js";

/** Decode as the site does (fromBase64Url, lib/utils/url-visualization-state.ts). */
const decode = (v: string) =>
  JSON.parse(Buffer.from(v, "base64url").toString("utf8"));

it("knows which viewers read cell and molecule settings", () => {
  const at = (p: string) => viewerKind(new URL(`https://x.test${p}?url=a`));

  expect(at("/viewer/from-s3")).toBe("cells");
  expect(at("/viewer/ds_abc")).toBe("cells");
  expect(at("/spatialdata-viewer/from-s3")).toBe("cells");
  expect(at("/sm-viewer/from-s3")).toBe("molecules");
  expect(at("/image-viewer/from-s3")).toBe("none");
});

it("writes the cells viewer's state: one gene, co-expression, column, cell types, 3D", () => {
  expect(viewState("cells", { genes: ["Slc17a7"] })).toEqual({
    g: "Slc17a7",
    m: ["gene"],
  });
  expect(
    viewState("cells", { genes: ["Gad1", "Slc17a7"], mode3d: true }),
  ).toEqual({
    g: "Gad1",
    m: ["gene"],
    g2: "Slc17a7",
    ce: true,
    vm: "3D",
  });
  expect(
    viewState("cells", {
      genes: [],
      column: "subclass",
      celltypes: ["L2/3 IT"],
    }),
  ).toEqual({ c: "subclass", ct: ["L2/3 IT"], m: ["celltype"] });
  expect(viewState("cells", { genes: [] })).toBeNull();
});

it("gives each molecule gene its own colour", () => {
  const s = viewState("molecules", { genes: ["Gad1", "Slc17a7"] }) as {
    genes: [string, string, number, boolean][];
  };

  expect(s.genes.map((g) => g[0])).toEqual(["Gad1", "Slc17a7"]);
  expect(s.genes[0][1]).not.toBe(s.genes[1][1]);
  expect(s.genes.every((g) => g[2] === 1 && g[3] === true)).toBe(true);
});

it("encodes the state as the site decodes it, Unicode included", () => {
  const state = { c: "cell type", ct: ["Astro–Epen", "L2/3 IT"] };

  expect(decode(toBase64Url(state))).toEqual(state);
  expect(toBase64Url(state)).not.toMatch(/[+/=]/);
});

it("cuts a search card to one viewer link per kind, with the counts", () => {
  const entry = (label: string, kind: string) => ({
    label,
    kind,
    datasetType: "",
    datasetId: null,
    url: null,
    viewerUrl: `https://x.test/${label}`,
  });
  const card = {
    id: "c1",
    title: "Atlas",
    species: null,
    tissue: null,
    platform: null,
    bilCode: null,
    cells: null,
    genes: null,
    source: "bil",
    license: null,
    pageUrl: "",
    entries: [
      entry("All cells", "single_cell"),
      ...Array.from({ length: 58 }, (_, i) =>
        entry(`S${i}`, "single_molecule"),
      ),
    ],
  };
  const c = compactCard(card);

  expect(c.entries.map((e) => e.label)).toEqual(["All cells", "S0"]);
  expect(c.entryCounts).toEqual({ single_cell: 1, single_molecule: 58 });
  // a short card is left alone
  expect(
    compactCard({ ...card, entries: card.entries.slice(0, 3) }).entryCounts,
  ).toBeUndefined();
});
