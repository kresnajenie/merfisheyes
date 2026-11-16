import { readFileSync } from "node:fs";
import path from "node:path";

import { afterEach, describe, expect, it, vi } from "vitest";

import schema from "@/docs/specs/schemas/preset.schema.json";
import { TECHNOLOGY_PRESETS, readPreset } from "@/lib/config/technology-presets";
import { presetConfigError, presetHash } from "@/lib/presets/validate";
import { COLORMAPS } from "@/lib/utils/colormaps";

const file = (name: string) =>
  JSON.parse(
    readFileSync(path.join(__dirname, "..", "data", "presets", name), "utf8"),
  );
const bins = { cellFootprint: { shape: "square" as const, size: 8 } };

describe("the preset contract (docs/specs/schemas/preset.schema.json)", () => {
  it("names exactly the colormaps the viewer has", () => {
    expect(schema.properties.colormap.enum).toEqual(Object.keys(COLORMAPS));
  });

  it("accepts the built-in presets and a hand-written file", () => {
    for (const preset of Object.values(TECHNOLOGY_PRESETS))
      expect(presetConfigError(preset)).toBeNull();
    expect(presetConfigError(file("valid.json"))).toBeNull();
    expect(presetConfigError({})).toBeNull();
  });

  it("refuses a malformed one, saying what is wrong", () => {
    const problem = presetConfigError(file("invalid.json"));

    expect(problem).toContain("/colormap");
    expect(problem).toContain("/cellOpacity");
    expect(problem).toContain("additional properties");
    expect(presetConfigError("viridis")).not.toBeNull();
    expect(presetConfigError({ hideZeroExpression: "yes" })).not.toBeNull();
  });

  it("gives a config one hash whatever its field order, per technology", () => {
    const a = presetHash("Visium HD", { colormap: "viridis", cellOpacity: 0.8 });
    const b = presetHash("Visium HD", { cellOpacity: 0.8, colormap: "viridis" });

    expect(a).toBe(b);
    expect(a).toMatch(/^[0-9a-f]{16}$/);
    expect(presetHash("Xenium", { colormap: "viridis", cellOpacity: 0.8 })).not.toBe(a);
    expect(presetHash("Visium HD", { colormap: "viridis" })).not.toBe(a);
  });

  it("a viewer keeps the fields it knows and drops the rest", () => {
    expect(
      readPreset({ colormap: "plasma", cellOpacity: 0.5, fromTheFuture: 1 }),
    ).toEqual({ colormap: "plasma", cellOpacity: 0.5 });
    expect(readPreset({ colormap: "rainbow", cellOpacity: 2 })).toEqual({});
    expect(readPreset(null)).toEqual({});
  });
});

describe("which preset a dataset opens with", () => {
  afterEach(() => {
    vi.unstubAllGlobals();
    vi.resetModules();
  });

  /** A fresh copy of the module, with the API answering as given. */
  async function load(opts: {
    search?: string;
    presets?: Record<string, unknown> | "down";
    version?: unknown;
  }) {
    vi.stubGlobal("window", { location: { search: opts.search ?? "" } });
    vi.stubGlobal(
      "fetch",
      vi.fn(async (url: string) => {
        if (url === "/api/presets") {
          if (opts.presets === "down") return { ok: false, status: 500 };

          return { ok: true, json: async () => ({ presets: opts.presets ?? {} }) };
        }

        return opts.version
          ? { ok: true, json: async () => opts.version }
          : { ok: false, status: 404 };
      }),
    );
    const m = await import("@/lib/config/technology-presets");

    await m.loadTechnologyPresets();

    return m;
  }

  it("the built-in one when nothing is published or the API is down", async () => {
    for (const presets of [{}, "down"] as const) {
      const m = await load({ presets });

      expect(m.technologyPresetOf(bins)).toEqual(TECHNOLOGY_PRESETS["Visium HD"]);
      vi.resetModules();
    }
  });

  it("the published one, and none when it is turned off", async () => {
    let m = await load({
      presets: { "Visium HD": { hash: "a", config: { colormap: "plasma" } } },
    });

    expect(m.technologyPresetOf(bins)).toEqual({ colormap: "plasma" });
    expect(m.technologyPresetOf({ cellFootprint: null })).toEqual({});

    vi.resetModules();
    m = await load({ presets: { "Visium HD": null } });
    expect(m.technologyPresetOf(bins)).toEqual({});
  });

  it("the version the link names, over what is published", async () => {
    const m = await load({
      search: "?preset=abc",
      presets: { "Visium HD": { hash: "a", config: { colormap: "plasma" } } },
      version: { hash: "abc", technology: "Visium HD", config: { cellOpacity: 0.3 } },
    });

    expect(m.technologyPresetOf(bins)).toEqual({ cellOpacity: 0.3 });
    expect(fetch).toHaveBeenCalledWith(
      "/api/presets/versions/abc",
      expect.anything(),
    );
  });

  it("the published one when the link names a version that does not exist", async () => {
    const m = await load({
      search: "?preset=nope",
      presets: { "Visium HD": { hash: "a", config: { colormap: "plasma" } } },
    });

    expect(m.technologyPresetOf(bins)).toEqual({ colormap: "plasma" });
  });
});
