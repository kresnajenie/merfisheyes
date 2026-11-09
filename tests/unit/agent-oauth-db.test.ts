import { afterAll, beforeAll, beforeEach, describe, expect, it, vi } from "vitest";

/**
 * The agent API against a real PostgreSQL (docs/specs/MCP-PROTOCOL.md §2) —
 * the gate of step 2 of docs/plans/MCP.md: nothing unpublished, unapproved,
 * internal, or someone else's private comes back, and Explore's own answers
 * are those of the shared rule.
 *
 * Needs a migrated database; skipped without one (CI's unit-tests job
 * provides one):
 *
 *   TEST_DATABASE_URL=postgresql://postgres@127.0.0.1:55432/meyes_mcp_test \
 *     npx vitest run tests/unit/agent-oauth-db.test.ts
 *
 * It empties the tables it uses — never point it at a database you keep.
 */

const DB = process.env.TEST_DATABASE_URL;
const ISSUER = "http://localhost:3000";

const session = vi.fn(async (): Promise<unknown> => null);

vi.mock("@/lib/auth", () => ({ auth: () => session() }));

const as = (id: string, role = "USER") => session.mockResolvedValue({ user: { id, role } });

describe.skipIf(!DB)("agent API on a database", async () => {
  process.env.DATABASE_URL = DB;
  process.env.OAUTH_ISSUER = ISSUER;

  const { prisma } = await import("@/lib/prisma");

  // The route handlers under test, by their path below app/api.
  const ROUTES = {
    "agent/v1/catalog/search": () => import("@/app/api/agent/v1/catalog/search/route"),
    "agent/v1/catalog/[id]": () => import("@/app/api/agent/v1/catalog/[id]/route"),
    "agent/v1/catalog/[id]/genes": () => import("@/app/api/agent/v1/catalog/[id]/genes/route"),
    "agent/v1/datasets/[id]": () => import("@/app/api/agent/v1/datasets/[id]/route"),
    "agent/v1/me/datasets": () => import("@/app/api/agent/v1/me/datasets/route"),
    "agent/v1/projects": () => import("@/app/api/agent/v1/projects/route"),
    explore: () => import("@/app/api/explore/route"),
    "explore/[id]": () => import("@/app/api/explore/[id]/route"),
    "catalog/facets": () => import("@/app/api/catalog/facets/route"),
  };
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const route = (p: keyof typeof ROUTES): Promise<any> => ROUTES[p]();

  const call = async (
    path: string,
    handler: Parameters<typeof route>[0],
    opts: { token?: string; params?: Record<string, string>; query?: string } = {},
  ) => {
    const mod = await route(handler);
    const req = new Request(`${ISSUER}${path}${opts.query ? `?${opts.query}` : ""}`, {
      headers: opts.token ? { authorization: `Bearer ${opts.token}` } : {},
    });

    return (mod.GET as (r: Request, c?: unknown) => Promise<Response>)(
      req,
      opts.params ? { params: Promise.resolve(opts.params) } : undefined,
    );
  };

  async function reset() {
    await prisma.projectDataset.deleteMany();
    await prisma.project.deleteMany();
    await prisma.catalogDatasetEntry.deleteMany();
    await prisma.catalogDataset.deleteMany();
    await prisma.dataset.deleteMany();
    await prisma.user.deleteMany();
  }

  const ds = (id: string, extra: Record<string, unknown> = {}) =>
    prisma.dataset.create({
      data: { id, fingerprint: id, numCells: 10, numGenes: 2, genes: ["Gad1", "Slc17a7"], status: "COMPLETE", datasetType: "spatialdata", ...extra },
    });
  const row = (id: string, extra: Record<string, unknown>, entry?: Record<string, unknown>) =>
    prisma.catalogDataset.create({
      data: {
        id,
        title: id,
        species: `species-${id}`,
        isPublished: true,
        genes: ["Gad1"],
        ...extra,
        ...(entry ? { entries: { create: { label: "Cells", datasetType: "spatialdata", ...entry } } } : {}),
      },
    });

  beforeAll(async () => {
    await reset();
    await prisma.user.createMany({
      data: [
        { id: "uA", email: "a@test" },
        { id: "uB", email: "b@test" },
        { id: "uC", email: "c@test", role: "ADMIN" },
      ],
    });
    await ds("ds_a1", { ownerId: "uA" });
    await ds("ds_a2", { ownerId: "uA" });
    await ds("ds_b", { ownerId: "uB" });
    await ds("ds_admin", { adminOwned: true });
    await ds("ds_comm", { ownerId: "uB" });
    await ds("ds_pending", { ownerId: "uB" });
    await ds("ds_curated", { s3BaseUrl: "https://example.org/curated.zarr" });
    await prisma.project.create({ data: { id: "pA1", title: "A one", ownerId: "uA", datasets: { create: { datasetId: "ds_a1" } } } });
    await prisma.project.create({ data: { id: "pB", title: "B", ownerId: "uB", datasets: { create: { datasetId: "ds_b" } } } });

    await row("curated", {}, { s3BaseUrl: "https://example.org/curated.zarr" });
    await row("unpublished", { isPublished: false }, { s3BaseUrl: "https://example.org/u.zarr" });
    await row("internal", { isInternal: true }, { s3BaseUrl: "https://example.org/i.zarr" });
    await row("approved", { isCommunity: true, reviewStatus: "APPROVED" }, { datasetId: "ds_comm" });
    await row("pending", { isCommunity: true, reviewStatus: "PENDING", isPublished: true }, { datasetId: "ds_pending" });
    await row("noentry", {});
  });

  afterAll(async () => {
    await reset();
    await prisma.$disconnect();
  });

  beforeEach(() => session.mockResolvedValue(null));

  it("searches only what Explore shows", async () => {
    const res = await call("/api/agent/v1/catalog/search", "agent/v1/catalog/search");
    const ids = (await res.json()).items.map((i: { id: string }) => i.id).sort();

    expect(ids).toEqual(["approved", "curated"]);
  });

  it("Explore's list, detail and facets give the same answers", async () => {
    const list = await (await (await route("explore")).GET(new Request(`${ISSUER}/api/explore?include=meta`))).json();

    expect(list.items.map((i: { id: string }) => i.id).sort()).toEqual(["approved", "curated"]);
    expect(list.filters.species).not.toContain("species-internal");
    expect(list.filters.species).not.toContain("species-unpublished");

    const facets = await (await (await route("catalog/facets")).GET()).json();

    expect(facets.species.sort()).toEqual(["species-approved", "species-curated", "species-noentry"]);
    as("uC", "ADMIN");
    expect((await (await (await route("catalog/facets")).GET()).json()).species).toContain("species-internal");
    session.mockResolvedValue(null);

    const detail = await route("explore/[id]");
    const status = async (id: string) => (await detail.GET(new Request(`${ISSUER}/x`), { params: Promise.resolve({ id }) })).status;

    expect(await status("curated")).toBe(200);
    expect(await status("internal")).toBe(404);
    expect(await status("unpublished")).toBe(404);
    expect(await status("pending")).toBe(404);
    as("uC", "ADMIN");
    expect(await status("internal")).toBe(200);
  });

  it("opens visible catalog rows by id or code, and nothing else", async () => {
    const get = (id: string) => call(`/api/agent/v1/catalog/${id}`, "agent/v1/catalog/[id]", { params: { id } });

    expect((await get("curated")).status).toBe(200);
    expect((await get("noentry")).status).toBe(200);
    for (const id of ["unpublished", "internal", "pending", "nope"]) expect((await get(id)).status, id).toBe(404);
    as("uC", "ADMIN");
    expect((await get("internal")).status).toBe(200);
  });

  it("opens a dataset only when a visible row names it, or it is the caller's own", async () => {
    const get = async (id: string) =>
      (await call(`/api/agent/v1/datasets/${id}`, "agent/v1/datasets/[id]", { params: { id } })).status;

    expect(await get("ds_comm")).toBe(200); // named by an approved community row
    expect(await get("ds_curated")).toBe(200); // named by URL
    expect(await get("ds_pending")).toBe(404); // its row is not approved
    expect(await get("ds_b")).toBe(404); // someone's private dataset, id or not
    expect(await get("ds_a1")).toBe(404);
    as("uA");
    expect(await get("ds_a1")).toBe(200);
    expect(await get("ds_a2")).toBe(200);
    expect(await get("ds_b")).toBe(404);
    expect(await get("ds_admin")).toBe(404);
    as("uC", "ADMIN");
    expect(await get("ds_admin")).toBe(200);
  });

  it("lists the caller's own datasets and projects only", async () => {
    as("uA");
    const ids = async (p: string, h: Parameters<typeof route>[0]) =>
      (await (await call(p, h)).json()).items.map((i: { id: string }) => i.id).sort();

    expect(await ids("/api/agent/v1/me/datasets", "agent/v1/me/datasets")).toEqual(["ds_a1", "ds_a2"]);
    expect(await ids("/api/agent/v1/projects", "agent/v1/projects")).toEqual(["pA1"]);
  });

  it("pages genes", async () => {
    const res = await call("/api/agent/v1/catalog/curated/genes", "agent/v1/catalog/[id]/genes", { params: { id: "curated" }, query: "q=gad" });

    expect(await res.json()).toMatchObject({ total: 1, items: ["Gad1"], has_more: false });
  });

  it("answers 401 without sign-in, naming the metadata and the scopes", async () => {
    const res = await call("/api/agent/v1/me/datasets", "agent/v1/me/datasets");

    expect(res.status).toBe(401);
    const h = res.headers.get("www-authenticate")!;

    expect(h).toContain(`resource_metadata="${ISSUER}/.well-known/oauth-protected-resource/api/agent/v1"`);
    expect(h).toContain('scope="projects:read"');
  });

  it("refuses a bearer token rather than treating it as anonymous", async () => {
    const res = await call("/api/agent/v1/catalog/search", "agent/v1/catalog/search", { token: "mea_nope" });

    expect(res.status).toBe(401);
    expect(res.headers.get("www-authenticate")).toContain('error="invalid_token"');
  });
});
