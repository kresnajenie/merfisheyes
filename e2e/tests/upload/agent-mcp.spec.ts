import { mkdtempSync } from "node:fs";
import { tmpdir } from "node:os";
import path from "node:path";

import { PrismaClient } from "@prisma/client";
import { expect, test } from "@playwright/test";

import { signInDev } from "../../helpers/auth";
import { McpStdio } from "../../helpers/mcp-stdio";

/**
 * Step 4's gate (docs/plans/MCP.md): the merfisheyes-mcp package, run as an
 * agent runs it (stdio), against this site — finds a public dataset without
 * signing in; is told to sign in for the person's own data; signs in through
 * the browser, the person allowing one project; then reaches that project
 * and its dataset and nothing else; signs out.
 *
 * Needs the package built (`npm ci && npm run build` in
 * packages/merfisheyes-mcp). Tagged @upload: it needs the suite's database;
 * E2E_AGENT=1 runs it against any migrated database.
 */

const ENABLED = process.env.E2E_UPLOAD === "1" || process.env.E2E_AGENT === "1";

test.describe("@upload agent via MCP", () => {
  test.skip(!ENABLED, "Set E2E_UPLOAD=1 (docker infra) or E2E_AGENT=1 (any migrated database).");

  test("public search, sign in for one project, sign out", async ({ page, baseURL }) => {
    test.setTimeout(120_000);
    const origin = new URL(baseURL!).origin;
    const run = Date.now();
    const email = `e2e-mcp-${run}@example.com`;
    const prisma = new PrismaClient();

    // The person, two projects, a dataset in the first; a public catalog row.
    await signInDev(page, email);
    const project = async (title: string) => (await (await page.request.post("/api/projects", { data: { title } })).json()).id as string;
    const chosen = await project(`MCP chosen ${run}`);
    const other = await project(`MCP other ${run}`);
    const user = await prisma.user.findUniqueOrThrow({ where: { email } });
    const datasetId = `ds_mcp${run}`.slice(0, 50);

    await prisma.dataset.create({
      data: {
        id: datasetId,
        fingerprint: datasetId,
        title: `MCP dataset ${run}`,
        numCells: 12,
        numGenes: 3,
        genes: ["Gad1", "Gad2", "Slc17a7"],
        status: "COMPLETE",
        datasetType: "spatialdata",
        ownerId: user.id,
        projects: { create: { projectId: chosen } },
      },
    });
    await prisma.catalogDataset.create({
      data: {
        title: `Public MCP atlas ${run}`,
        species: "Mouse",
        isPublished: true,
        genes: ["Gad1"],
        entries: { create: { label: "Cells", datasetType: "spatialdata", s3BaseUrl: `https://example.org/mcp-${run}.zarr` } },
      },
    });

    const mcp = await McpStdio.start({
      MERFISHEYES_URL: origin,
      MERFISHEYES_NO_BROWSER: "1",
      MERFISHEYES_CREDENTIALS: path.join(mkdtempSync(path.join(tmpdir(), "meyes-e2e-")), "credentials.json"),
    });

    try {
      const { tools } = await mcp.request<{ tools: { name: string }[] }>("tools/list");

      expect(tools.map((t) => t.name).sort()).toEqual([
        "meyes_get_dataset",
        "meyes_get_processing_status",
        "meyes_get_project",
        "meyes_list_facets",
        "meyes_list_genes",
        "meyes_list_my_datasets",
        "meyes_list_projects",
        "meyes_search_datasets",
        "meyes_sign_in",
        "meyes_sign_out",
        "meyes_whoami",
      ]);

      // Public: no sign-in.
      const found = await mcp.call("meyes_search_datasets", { query: `atlas ${run}` });

      expect(found.isError).toBeFalsy();
      expect(found.content[0].text).toContain(`Public MCP atlas ${run}`);

      // Own data: told how to sign in.
      const before = await mcp.call("meyes_list_my_datasets");

      expect(before.isError).toBe(true);
      expect(before.content[0].text).toContain("meyes_sign_in");

      // Sign in; the person opens the link and allows one project.
      const signIn = await mcp.call("meyes_sign_in");
      const link = /(http\S+\/oauth\/authorize\S+)/.exec(signIn.content[0].text)![1];

      await page.goto(link);
      await expect(page.getByRole("heading", { name: /Allow MERFISHEYES MCP/ })).toBeVisible();
      await page.getByLabel("Only these projects").check();
      await page.getByLabel(`MCP chosen ${run}`).check();
      await page.getByRole("button", { name: "Allow" }).click();
      await expect(page.getByText("Signed in to MERFISHEYES")).toBeVisible();

      const me = await mcp.call("meyes_whoami");

      expect(me.content[0].text).toContain(email);

      const projects = await mcp.call("meyes_list_projects", { response_format: "json" });

      expect((projects.structuredContent!.items as { id: string }[]).map((p) => p.id)).toEqual([chosen]);
      expect((await mcp.call("meyes_get_project", { id: other })).isError).toBe(true);

      const mine = await mcp.call("meyes_get_dataset", { id: datasetId });

      expect(mine.isError).toBeFalsy();
      expect(mine.content[0].text).toContain(`MCP dataset ${run}`);
      const genes = await mcp.call("meyes_list_genes", { id: datasetId, query: "gad" });

      expect(genes.content[0].text).toContain("Gad1, Gad2");

      // Sign out: the token is revoked.
      expect((await mcp.call("meyes_sign_out")).content[0].text).toContain("Signed out");
      expect((await mcp.call("meyes_list_my_datasets")).isError).toBe(true);
    } finally {
      mcp.close();
      await prisma.$disconnect();
    }
  });
});
