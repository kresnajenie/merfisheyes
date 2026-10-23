import { beforeEach, describe, expect, it, vi } from "vitest";

/**
 * Every link the app sends out opens the viewer of what the dataset is
 * (lib/datasets/kind.ts `viewerHref`): a SpatialData store — including one
 * that keeps the sm_ id of the single molecule card it came from — an image,
 * cells, molecules. Checked on the helper and on the emails as sent.
 */

const sent: Array<{ Message: { Subject: { Data: string }; Body: { Text: { Data: string }; Html: { Data: string } } } }> = [];

vi.mock("@aws-sdk/client-ses", () => ({
  SESClient: class {
    send = vi.fn(async (command: { input: (typeof sent)[number] }) => {
      sent.push(command.input);
    });
  },
  SendEmailCommand: class {
    constructor(public input: unknown) {}
  },
}));

const { datasetTypeOfViewer, viewerHref } = await import("@/lib/datasets/kind");
const { sendDatasetReadyEmail, sendDuplicateDatasetEmail } = await import("@/lib/ses");

const BASE = "http://localhost:3000";

/** What a link ought to be, for each kind of dataset. */
const CASES = [
  { id: "ds_cells", datasetType: "single_cell", path: "/viewer/ds_cells" },
  { id: "ds_h5ad", datasetType: "h5ad", path: "/viewer/ds_h5ad" },
  { id: "sm_mols", datasetType: "single_molecule", path: "/sm-viewer/sm_mols" },
  { id: "ds_store", datasetType: "spatialdata", path: "/spatialdata-viewer/ds_store" },
  // A molecule table from the single molecule card, prepared into a store
  { id: "sm_store", datasetType: "spatialdata", path: "/spatialdata-viewer/sm_store" },
  { id: "ds_img", datasetType: "image", path: "/image-viewer/ds_img" },
];

describe("viewerHref", () => {
  it.each(CASES)("$datasetType $id → $path", ({ id, datasetType, path }) => {
    expect(viewerHref({ id, datasetType })).toBe(path);
  });

  it("reads a viewer link back as the type it opens", () => {
    expect(datasetTypeOfViewer("viewer")).toBe("single_cell");
    expect(datasetTypeOfViewer("sm-viewer")).toBe("single_molecule");
    expect(datasetTypeOfViewer("spatialdata-viewer")).toBe("spatialdata");
    expect(datasetTypeOfViewer("image-viewer")).toBe("image");
    expect(datasetTypeOfViewer("explore")).toBeNull();
  });
});

describe("the emails", () => {
  beforeEach(() => {
    sent.length = 0;
  });

  it.each(CASES)("dataset ready: $datasetType $id links $path", async ({ id, datasetType, path }) => {
    await sendDatasetReadyEmail({ email: "a@example.com", datasetId: id, datasetType });

    expect(sent[0].Message.Body.Text.Data).toContain(`${BASE}${path}`);
    expect(sent[0].Message.Body.Html.Data).toContain(`${BASE}${path}`);
  });

  it("counts molecules for a molecule store, cells for a cells store", async () => {
    await sendDatasetReadyEmail({
      email: "a@example.com",
      datasetId: "sm_store",
      datasetType: "spatialdata",
      molecules: true,
      metadata: { numCells: 851_483_449 },
    });
    await sendDatasetReadyEmail({
      email: "a@example.com",
      datasetId: "ds_store",
      datasetType: "spatialdata",
      molecules: false,
      metadata: { numCells: 124_364 },
    });

    expect(sent[0].Message.Body.Text.Data).toContain("Molecules: 851,483,449");
    expect(sent[1].Message.Body.Text.Data).toContain("Cells: 124,364");
  });

  it.each(CASES)("duplicate: the existing $datasetType $id links $path", async ({ id, datasetType, path }) => {
    await sendDuplicateDatasetEmail({
      email: "a@example.com",
      existingId: id,
      existingType: datasetType,
    });

    expect(JSON.stringify(sent[0].Message.Body)).toContain(`${BASE}${path}`);
  });

  it("the single molecule upload's email links what the dataset is", async () => {
    vi.doMock("@/lib/prisma", () => ({
      prisma: {
        dataset: {
          findUnique: vi.fn(async () => ({
            title: "set1",
            numCells: 10,
            numGenes: 2,
            datasetType: "spatialdata",
          })),
        },
      },
    }));
    const { POST } = await import("@/app/api/send-email-single-molecule/route");

    await POST(
      new Request("http://localhost/api/send-email-single-molecule", {
        method: "POST",
        body: JSON.stringify({ email: "a@example.com", datasetId: "sm_store" }),
      }) as never,
    );

    expect(JSON.stringify(sent[0].Message.Body)).toContain(
      "/spatialdata-viewer/sm_store",
    );
  });
});
