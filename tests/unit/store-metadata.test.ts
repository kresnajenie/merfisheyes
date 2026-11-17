/**
 * readStoreMetadata (lib/spatialdata/metadata.ts): the info panel's one
 * reader of a store's root attributes (docs/specs/metadata.md §4). Our full
 * block, a SpatialData-DB store with only their labels (synthetic, and a
 * real root from tests/data/sddb/roots), and a store with none.
 */
import { readdirSync, readFileSync } from "node:fs";
import path from "node:path";

import { describe, expect, it } from "vitest";

import { ontologyHref, readStoreMetadata } from "@/lib/spatialdata/metadata";
import type { StoreMetadata } from "@/lib/spatialdata/metadata";

/** tests/python/test_specs.py metadata_example(), which the schema validates. */
function metadataExample() {
  return {
    spatialdata_attrs: { version: "0.2" },
    metadata_general: {
      Product: "MERSCOPE", Assay: "Vizgen MERSCOPE", "Biomaterial Type": "Specimen from Organism",
      Organism: "mouse", Tissue: "brain", Modality: "RNA", Disease: "unknown",
      "Default Table": "table", "Dataset Url": "https://api.brainimagelibrary.org/web/view?bildid=ace-dip-van",
      "Sample ID at Source": "ace-dip-van", License: "CC BY-NC 4.0", "Development Stage": "13-week-old stage",
      "HnE Image Available": false, "IF Image Available": false, "FTU Annotation Available": false,
      "Number of Transcripts": 472686422, "General Metadata Version": "0.3.2",
    },
    metadata_merscope: { "Merscope Metadata Version": "0.1.0" },
    meyes_metadata: {
      version: "1", canonical_key: "bil:ace-dip-van", study_id: "bil:ace-dip-van",
      title: "MERFISH Imaging of mesocorticolimbic regions in the mouse brain (infected with EcoHIV virus; systemic)",
      organism: "mouse", organism_id: "NCBITaxon:10090",
      tissue: "brain", tissue_id: "UBERON:0000955",
      disease: "unknown", disease_id: null,
      development_stage: "13-week-old stage", meyes_development_stage_id: "MmusDv:0000159",
      technology: "vizgen_merscope", vendor: "vizgen", assay: "Vizgen MERSCOPE",
      n_obs: 243224, n_features: 300, n_transcripts: 472686422,
      timepoint: "91.0 Days", timepoint_kind: "age", timepoint_days: 91.0,
      genotype: "C57BL/6J", condition: "EcoHIV infection, systemic",
      meyes_sex: "male", meyes_sex_id: "PATO:0000384", license: "CC-BY-NC-4.0",
      meyes_kind: "single_cell", meyes_n_sections: 3,
      meyes_dims: 2, meyes_sections: [{ id: "s1" }, { id: "s2" }, { id: "s3" }],
      meyes_matrix_hash: "mh1:d0ff37d6072dd0d8ef11bae317b5fa6e573f59364d40c61962dd183432dfb025",
    },
    meyes_provenance: {
      canonical_key: "bil:ace-dip-van", pipeline: "93ab1c0", input: "merscope",
      written_at: "2026-11-16T00:00:00Z",
      validation: {
        status: "pass", at: "2026-11-16T00:00:00Z",
        sources: { organism: "ncbitaxon 2023-06-20", tissue: "uberon 2024-08-07" },
      },
    },
  };
}

const fact = (m: StoreMetadata | null, key: string) =>
  m?.groups.flatMap((g) => g.facts).find((f) => f.key === key);

describe("readStoreMetadata", () => {
  it("reads our block: ids beside labels, ours over theirs", () => {
    const m = readStoreMetadata(metadataExample());

    expect(m).not.toBeNull();
    expect(m!.title).toMatch(/^MERFISH Imaging of mesocorticolimbic/);
    expect(fact(m, "tissue")).toMatchObject({
      value: "brain",
      id: "UBERON:0000955",
      href: ontologyHref("UBERON:0000955"),
      from: "meyes",
    });
    expect(fact(m, "organism")).toMatchObject({ id: "NCBITaxon:10090" });
    expect(fact(m, "development_stage")).toMatchObject({ id: "MmusDv:0000159" });
    expect(fact(m, "meyes_sex")).toMatchObject({ value: "male", id: "PATO:0000384" });
    // unknown is a value; its id is null
    expect(fact(m, "disease")).toMatchObject({ value: "unknown", id: null });
    expect(fact(m, "timepoint")?.value).toBe("91.0 Days (age)");
    expect(fact(m, "n_obs")?.value).toBe(243224);
    // A licence id is shown by its label, with its page
    expect(fact(m, "license")).toMatchObject({
      value: "CC BY-NC 4.0",
      id: "CC-BY-NC-4.0",
      href: "https://creativecommons.org/licenses/by-nc/4.0/",
    });
    // Only in theirs: filled from metadata_general
    expect(fact(m, "product")).toMatchObject({ value: "MERSCOPE", from: "general" });
    expect(fact(m, "dataset_url")?.href).toMatch(/^https:\/\/api\.brainimagelibrary/);

    expect(m!.matrixHash).toMatch(/^mh1:[0-9a-f]{64}$/);
    expect(m!.sections.map((s) => s.id)).toEqual(["s1", "s2", "s3"]);
    expect(m!.validation).toMatchObject({ status: "pass" });
    expect(m!.validation!.sources.tissue).toBe("uberon 2024-08-07");
    expect(m!.provenance).toMatchObject({
      writer: "merfisheyes",
      version: "93ab1c0",
      input: "merscope",
    });
    // Every block, every key: nothing dropped
    expect(m!.blocks.map((b) => b.name)).toEqual([
      "meyes_metadata",
      "metadata_general",
      "metadata_merscope",
      "meyes_provenance",
    ]);
    expect(m!.blocks[0].entries).toHaveLength(
      Object.keys(metadataExample().meyes_metadata).length,
    );
  });

  it("says not validated until the check has passed", () => {
    const attrs = metadataExample();

    (attrs.meyes_provenance.validation as { status: string }).status =
      "not_validated";
    expect(readStoreMetadata(attrs)!.validation!.status).toBe("not_validated");
    delete (attrs.meyes_provenance as Partial<typeof attrs.meyes_provenance>)
      .validation;
    expect(readStoreMetadata(attrs)!.validation).toBeNull();
  });

  it("reads a SpatialData-DB store: their labels, no ids", () => {
    const m = readStoreMetadata({
      spatialdata_attrs: { version: "0.2" },
      metadata_general: {
        Organism: "human", Tissue: "tonsil", Assay: "CosMx", Disease: "unknown",
        License: "CC BY 4.0", DOI: "10.1000/xyz", "Number of Transcripts": null,
        "Count Matrix Hash": "cmh1:abc",
      },
      sddb_provenance: {
        sddb_id: "sddb_x", canonical_key: "bruker:protein_tonsil_25",
        container_version: "unknown", converted_at: "2026-10-20T01:38:28Z",
      },
    });

    expect(m).not.toBeNull();
    expect(m!.title).toBeNull();
    expect(fact(m, "tissue")).toMatchObject({
      value: "tonsil",
      id: null,
      href: null,
      from: "general",
    });
    expect(fact(m, "license")).toMatchObject({ id: "CC-BY-4.0" });
    expect(fact(m, "doi")?.href).toBe("https://doi.org/10.1000/xyz");
    // null is absent, not shown
    expect(fact(m, "n_transcripts")).toBeUndefined();
    expect(m!.validation).toBeNull();
    expect(m!.matrixHash).toBeNull();
    expect(m!.provenance).toMatchObject({
      writer: "spatialdata-db",
      canonicalKey: "bruker:protein_tonsil_25",
      writtenAt: "2026-10-20T01:38:28Z",
    });
    // Their hash is theirs (cmh1): kept under all fields, not taken for ours
    const general = m!.blocks.find((b) => b.name === "metadata_general")!;

    expect(general.entries).toContainEqual(["Count Matrix Hash", "cmh1:abc"]);
  });

  it("reads every real SpatialData-DB root", () => {
    const dir = path.resolve(__dirname, "../data/sddb/roots");
    const attrs = (name: string) =>
      JSON.parse(readFileSync(path.join(dir, name), "utf8")).attributes;

    for (const name of readdirSync(dir).filter((n) => n.endsWith(".json"))) {
      const each = readStoreMetadata(attrs(name));

      expect(fact(each, "organism"), name).toBeDefined();
      expect(fact(each, "tissue"), name).toBeDefined();
    }
    const m = readStoreMetadata(attrs("571yTEm4hqlkXRPf.json"));

    expect(fact(m, "organism")?.value).toBe("human");
    expect(fact(m, "tissue")?.value).toBe("tonsil");
    expect(fact(m, "modality")?.value).toBe("protein");
    expect(m!.blocks.map((b) => b.name)).toEqual([
      "metadata_general",
      "metadata_cosmx",
      "sddb_provenance",
    ]);
    expect(m!.provenance?.writer).toBe("spatialdata-db");
  });

  it("is null for a store with none, and never throws", () => {
    expect(readStoreMetadata({ spatialdata_attrs: { version: "0.2" } })).toBeNull();
    expect(readStoreMetadata({})).toBeNull();
    expect(readStoreMetadata(undefined)).toBeNull();
    expect(readStoreMetadata(null)).toBeNull();
    expect(readStoreMetadata("zarr")).toBeNull();
    // Blocks of the wrong shape are skipped, not trusted
    expect(
      readStoreMetadata({ metadata_general: "x", meyes_metadata: [1, 2] }),
    ).toBeNull();
    const m = readStoreMetadata({
      meyes_metadata: { organism: null, meyes_sections: "nope", n_obs: "many" },
    });

    expect(m).not.toBeNull();
    expect(m!.groups.flatMap((g) => g.facts).map((f) => f.key)).toEqual(["n_obs"]);
    expect(m!.sections).toEqual([]);
  });
});
