import { beforeEach, describe, expect, it, vi } from "vitest";

/**
 * The "ready" email for a save of several datasets links each of them, and a
 * name typed by the user cannot add markup to it.
 */

const send = vi.fn();

vi.mock("@aws-sdk/client-ses", () => ({
  SESClient: class {
    send = send;
  },
  SendEmailCommand: class {
    constructor(public input: any) {}
  },
}));

const { sendDatasetReadyEmail } = await import("@/lib/ses");

const sent = () => {
  const { input } = send.mock.calls[0][0];

  return {
    to: input.Destination.ToAddresses,
    subject: input.Message.Subject.Data as string,
    text: input.Message.Body.Text.Data as string,
    html: input.Message.Body.Html.Data as string,
  };
};

describe("the ready email", () => {
  beforeEach(() => {
    send.mockReset();
    process.env.NEXT_PUBLIC_APP_URL = "https://example.org";
  });

  it("links the other datasets of the save, each in its own viewer", async () => {
    await sendDatasetReadyEmail({
      email: "owner@example.com",
      datasetId: "ds_cells",
      datasetType: "spatialdata",
      datasetName: "set1 (cells)",
      alsoSaved: [
        { datasetId: "sm_mol", datasetType: "single_molecule", datasetName: "set1 (molecules)" },
        { datasetId: "ds_img", datasetType: "image", datasetName: "set1 (image)" },
      ],
    });

    const mail = sent();

    expect(send).toHaveBeenCalledTimes(1);
    expect(mail.to).toEqual(["owner@example.com"]);
    expect(mail.subject).toBe("set1 (cells) - Dataset Ready - MERFISHEYES");
    for (const body of [mail.text, mail.html]) {
      expect(body).toContain("https://example.org/spatialdata-viewer/ds_cells");
      expect(body).toContain("https://example.org/sm-viewer/sm_mol");
      expect(body).toContain("https://example.org/image-viewer/ds_img");
      expect(body).toContain("Also in this save");
    }
  });

  it("is as before for one dataset", async () => {
    await sendDatasetReadyEmail({
      email: "owner@example.com",
      datasetId: "ds_img",
      datasetType: "image",
      datasetName: "slide",
    });

    const mail = sent();

    expect(mail.text).not.toContain("Also in this save");
    expect(mail.html).not.toContain("Also in this save");
  });

  it("escapes a name in the list", async () => {
    await sendDatasetReadyEmail({
      email: "owner@example.com",
      datasetId: "ds_cells",
      datasetType: "spatialdata",
      alsoSaved: [{ datasetId: "ds_img", datasetType: "image", datasetName: "<b>x</b>" }],
    });

    expect(sent().html).toContain("&lt;b&gt;x&lt;/b&gt;");
    expect(sent().html).not.toContain("<b>x</b>");
  });
});
