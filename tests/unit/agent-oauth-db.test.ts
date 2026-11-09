import { createHash, randomBytes } from "crypto";

import { afterAll, beforeAll, beforeEach, describe, expect, it, vi } from "vitest";

/**
 * The agent API and the authorization server against a real PostgreSQL
 * (docs/specs/MCP-PROTOCOL.md §2, §6) — the gates of steps 2 and 3 of
 * docs/plans/MCP.md:
 *
 *   - nothing unpublished, unapproved, internal, or someone else's private
 *     comes back; a grant limited to projects reaches only those;
 *   - Explore's own answers are those of the shared rule;
 *   - the code flow: PKCE, single-use codes, resource binding, refresh
 *     rotation with theft detection, revocation, expiry, scope challenges.
 *
 * Needs a migrated database; skipped without one:
 *
 *   TEST_DATABASE_URL=postgresql://postgres@127.0.0.1:55432/meyes_mcp_test \
 *     npx vitest run tests/unit/agent-oauth-db.test.ts
 *
 * It empties the tables it uses — never point it at a database you keep.
 */

const DB = process.env.TEST_DATABASE_URL;
const ISSUER = "http://localhost:3000";
const RESOURCE = `${ISSUER}/api/agent/v1`;

const session = vi.fn(async (): Promise<unknown> => null);

vi.mock("@/lib/auth", () => ({ auth: () => session() }));

describe.skipIf(!DB)("agent API and OAuth on a database", async () => {
  process.env.DATABASE_URL = DB;
  process.env.OAUTH_ISSUER = ISSUER;
  process.env.AUTH_SECRET = "test-secret";

  const { prisma } = await import("@/lib/prisma");
  const flow = await import("@/lib/oauth/flow");
  const { registerClient } = await import("@/lib/oauth/clients");
  // The route handlers under test, by their path below app/api.
  const ROUTES = {
    "agent/v1/catalog/search": () => import("@/app/api/agent/v1/catalog/search/route"),
    "agent/v1/catalog/facets": () => import("@/app/api/agent/v1/catalog/facets/route"),
    "oauth/register": () => import("@/app/api/oauth/register/route"),
    "agent/v1/catalog/[id]": () => import("@/app/api/agent/v1/catalog/[id]/route"),
    "agent/v1/catalog/[id]/genes": () => import("@/app/api/agent/v1/catalog/[id]/genes/route"),
    "agent/v1/datasets/[id]": () => import("@/app/api/agent/v1/datasets/[id]/route"),
    "agent/v1/me/datasets": () => import("@/app/api/agent/v1/me/datasets/route"),
    "agent/v1/projects": () => import("@/app/api/agent/v1/projects/route"),
    "agent/v1/projects/[id]": () => import("@/app/api/agent/v1/projects/[id]/route"),
    explore: () => import("@/app/api/explore/route"),
    "explore/[id]": () => import("@/app/api/explore/[id]/route"),
    "catalog/facets": () => import("@/app/api/catalog/facets/route"),
    "account/agents": () => import("@/app/api/account/agents/route"),
    "account/agents/[id]": () => import("@/app/api/account/agents/[id]/route"),
    "oauth/token": () => import("@/app/api/oauth/token/route"),
    "oauth/revoke": () => import("@/app/api/oauth/revoke/route"),
    "oauth/authorize": () => import("@/app/api/oauth/authorize/route"),
    "oauth/metadata": () => import("@/app/api/oauth/metadata/route"),
    "oauth/protected-resource": () => import("@/app/api/oauth/protected-resource/route"),
  };
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const route = (p: keyof typeof ROUTES): Promise<any> => ROUTES[p]();

  const call = async (
    path: string,
    handler: Parameters<typeof route>[0],
    opts: { token?: string; params?: Record<string, string>; query?: string; ip?: string } = {},
  ) => {
    const mod = await route(handler);
    const req = new Request(`${ISSUER}${path}${opts.query ? `?${opts.query}` : ""}`, {
      headers: {
        ...(opts.token ? { authorization: `Bearer ${opts.token}` } : {}),
        ...(opts.ip ? { "x-forwarded-for": opts.ip } : {}),
      },
    });

    return (mod.GET as (r: Request, c?: unknown) => Promise<Response>)(
      req,
      opts.params ? { params: Promise.resolve(opts.params) } : undefined,
    );
  };

  async function reset() {
    await prisma.rateLimitCounter.deleteMany();
    await prisma.oAuthToken.deleteMany();
    await prisma.oAuthCode.deleteMany();
    await prisma.oAuthGrant.deleteMany();
    await prisma.oAuthClient.deleteMany();
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
    await prisma.project.create({ data: { id: "pA2", title: "A two", ownerId: "uA" } });
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

  /** A grant and an access token, straight from the flow (the HTTP
   * round-trip is its own test below). */
  async function tokenFor(userId: string, scopes: string[], reach: { all: true } | { all: false; projectIds: string[] }) {
    const client = await registerClient({ redirect_uris: ["http://127.0.0.1/cb"], client_name: `t-${userId}-${randomBytes(3).toString("hex")}` });
    const verifier = randomBytes(32).toString("base64url");
    const q = new URLSearchParams({
      response_type: "code",
      client_id: client.client_id,
      redirect_uri: "http://127.0.0.1:5000/cb",
      code_challenge: createHash("sha256").update(verifier).digest("base64url"),
      code_challenge_method: "S256",
      scope: scopes.join(" "),
      resource: RESOURCE,
      state: "s",
    });
    const v = await flow.validateAuthorize(q);

    if (!("request" in v)) throw new Error(`authorize refused: ${v.redirect}`);
    const back = new URL(await flow.approve(userId, v.request, reach));
    const tokens = await flow.exchangeCode({
      code: back.searchParams.get("code")!,
      clientId: client.client_id,
      redirectUri: "http://127.0.0.1:5000/cb",
      codeVerifier: verifier,
      resource: RESOURCE,
    });

    return { ...tokens, clientId: client.client_id, back };
  }

  describe("what an agent sees (step 2)", () => {
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
      session.mockResolvedValue({ user: { id: "uC", role: "ADMIN" } });
      expect((await (await (await route("catalog/facets")).GET()).json()).species).toContain("species-internal");
      session.mockResolvedValue(null);

      const detail = await route("explore/[id]");
      const status = async (id: string) => (await detail.GET(new Request(`${ISSUER}/x`), { params: Promise.resolve({ id }) })).status;

      expect(await status("curated")).toBe(200);
      expect(await status("internal")).toBe(404);
      expect(await status("unpublished")).toBe(404);
      expect(await status("pending")).toBe(404);
      session.mockResolvedValue({ user: { id: "uC", role: "ADMIN" } });
      expect(await status("internal")).toBe(200);
    });

    it("opens visible catalog rows by id or code, and nothing else", async () => {
      const get = (id: string, token?: string) =>
        call(`/api/agent/v1/catalog/${id}`, "agent/v1/catalog/[id]", { params: { id }, token });

      expect((await get("curated")).status).toBe(200);
      expect((await get("noentry")).status).toBe(200);
      for (const id of ["unpublished", "internal", "pending", "nope"]) expect((await get(id)).status, id).toBe(404);
      const admin = await tokenFor("uC", ["catalog:read"], { all: true });

      expect((await get("internal", admin.access_token)).status).toBe(200);
    });

    it("opens a dataset only when a visible row names it, or it is the caller's own", async () => {
      const a = await tokenFor("uA", ["projects:read"], { all: true });
      const admin = await tokenFor("uC", ["projects:read"], { all: true });
      const get = async (id: string, token?: string) =>
        (await call(`/api/agent/v1/datasets/${id}`, "agent/v1/datasets/[id]", { params: { id }, token })).status;

      expect(await get("ds_comm")).toBe(200); // named by an approved community row
      expect(await get("ds_curated")).toBe(200); // named by URL
      expect(await get("ds_pending")).toBe(404); // its row is not approved
      expect(await get("ds_b")).toBe(404); // someone's private dataset, id or not
      expect(await get("ds_a1")).toBe(404);
      expect(await get("ds_a1", a.access_token)).toBe(200);
      expect(await get("ds_a2", a.access_token)).toBe(200);
      expect(await get("ds_b", a.access_token)).toBe(404);
      expect(await get("ds_admin", a.access_token)).toBe(404);
      expect(await get("ds_admin", admin.access_token)).toBe(200);
    });

    it("a grant limited to projects reaches only those", async () => {
      const a = await tokenFor("uA", ["projects:read"], { all: false, projectIds: ["pA1"] });
      const ids = async (p: string, h: string) =>
        (await (await call(p, h, { token: a.access_token })).json()).items.map((i: { id: string }) => i.id);

      expect(await ids("/api/agent/v1/projects", "agent/v1/projects")).toEqual(["pA1"]);
      expect(await ids("/api/agent/v1/me/datasets", "agent/v1/me/datasets")).toEqual(["ds_a1"]);
      const get = async (id: string) =>
        (await call(`/api/agent/v1/datasets/${id}`, "agent/v1/datasets/[id]", { params: { id }, token: a.access_token })).status;

      expect(await get("ds_a1")).toBe(200);
      expect(await get("ds_a2")).toBe(404);
      const proj = async (id: string) =>
        (await call(`/api/agent/v1/projects/${id}`, "agent/v1/projects/[id]", { params: { id }, token: a.access_token })).status;

      expect(await proj("pA2")).toBe(404);
      expect(await proj("pB")).toBe(404);
    });

    it("refuses a grant reaching someone else's project", async () => {
      await expect(tokenFor("uA", ["projects:read"], { all: false, projectIds: ["pB"] })).rejects.toThrow(/your projects/);
    });

    it("pages genes", async () => {
      const res = await call("/api/agent/v1/catalog/curated/genes", "agent/v1/catalog/[id]/genes", { params: { id: "curated" }, query: "q=gad" });

      expect(await res.json()).toMatchObject({ total: 1, items: ["Gad1"], has_more: false });
    });
  });

  describe("challenges (step 3, §6.5)", () => {
    it("401 without a token, naming the metadata and the scopes", async () => {
      const res = await call("/api/agent/v1/me/datasets", "agent/v1/me/datasets");

      expect(res.status).toBe(401);
      const h = res.headers.get("www-authenticate")!;

      expect(h).toContain(`resource_metadata="${ISSUER}/.well-known/oauth-protected-resource/api/agent/v1"`);
      expect(h).toContain('scope="projects:read"');
    });

    it("401 for a bad token — never anonymous", async () => {
      const res = await call("/api/agent/v1/catalog/search", "agent/v1/catalog/search", { token: "mea_nope" });

      expect(res.status).toBe(401);
      expect(res.headers.get("www-authenticate")).toContain('error="invalid_token"');
    });

    it("403 insufficient_scope for a valid token without the scope", async () => {
      const t = await tokenFor("uA", ["catalog:read"], { all: true });
      const res = await call("/api/agent/v1/me/datasets", "agent/v1/me/datasets", { token: t.access_token });

      expect(res.status).toBe(403);
      expect(res.headers.get("www-authenticate")).toContain('error="insufficient_scope", scope="projects:read"');
    });
  });

  describe("the code flow (step 3)", () => {
    it("redirects back with code, state and iss", async () => {
      const t = await tokenFor("uA", ["projects:read"], { all: true });

      expect(t.back.searchParams.get("state")).toBe("s");
      expect(t.back.searchParams.get("iss")).toBe(ISSUER);
      expect(t.back.port).toBe("5000");
      expect(t.scope.split(" ").sort()).toEqual(["catalog:read", "projects:read"]);
      expect(t.expires_in).toBe(3600);
    });

    it("refuses untrusted clients and redirect URIs without redirecting", async () => {
      const bad = new URLSearchParams({ client_id: "mcl_unknown", redirect_uri: "http://127.0.0.1/cb" });

      await expect(flow.validateAuthorize(bad)).rejects.toBeInstanceOf(flow.UntrustedRequestError);
      const c = await registerClient({ redirect_uris: ["https://app.example/cb"] });
      const wrong = new URLSearchParams({ client_id: c.client_id, redirect_uri: "https://evil.example/cb" });

      await expect(flow.validateAuthorize(wrong)).rejects.toBeInstanceOf(flow.UntrustedRequestError);
    });

    it("sends errors back once the redirect is trusted", async () => {
      const c = await registerClient({ redirect_uris: ["https://app.example/cb"] });
      const base = { client_id: c.client_id, redirect_uri: "https://app.example/cb", response_type: "code", state: "z" };
      const err = async (extra: Record<string, string>) => {
        const v = await flow.validateAuthorize(new URLSearchParams({ ...base, ...extra }));

        return "redirect" in v ? new URL(v.redirect).searchParams.get("error") : null;
      };
      const challenge = createHash("sha256").update("v".repeat(43)).digest("base64url");

      expect(await err({})).toBe("invalid_request"); // no PKCE
      expect(await err({ code_challenge: challenge, code_challenge_method: "plain" })).toBe("invalid_request");
      expect(await err({ code_challenge: challenge, code_challenge_method: "S256", scope: "admin" })).toBe("invalid_scope");
      expect(await err({ code_challenge: challenge, code_challenge_method: "S256", resource: "https://other.example" })).toBe("invalid_target");
      expect(await err({ code_challenge: challenge, code_challenge_method: "S256", response_type: "token" })).toBe("unsupported_response_type");
    });

    it("a code is single use; replaying it revokes what it gave", async () => {
      const client = await registerClient({ redirect_uris: ["http://127.0.0.1/cb"] });
      const verifier = randomBytes(32).toString("base64url");
      const v = await flow.validateAuthorize(
        new URLSearchParams({
          response_type: "code", client_id: client.client_id, redirect_uri: "http://127.0.0.1:7/cb",
          code_challenge: createHash("sha256").update(verifier).digest("base64url"), code_challenge_method: "S256",
        }),
      );

      if (!("request" in v)) throw new Error("refused");
      const code = new URL(await flow.approve("uA", v.request, { all: true })).searchParams.get("code")!;
      const args = { code, clientId: client.client_id, redirectUri: "http://127.0.0.1:7/cb", codeVerifier: verifier, resource: null };

      await expect(flow.exchangeCode({ ...args, codeVerifier: verifier.slice(1) + "A" })).rejects.toThrow(/code_verifier/);
      await expect(flow.exchangeCode({ ...args, redirectUri: "http://127.0.0.1:8/cb" })).rejects.toThrow(/redirect_uri/);
      await expect(flow.exchangeCode({ ...args, clientId: "mcl_other" })).rejects.toThrow(/another client/);
      const first = await flow.exchangeCode(args);

      expect(await flow.verifyAccessToken(first.access_token, RESOURCE)).not.toBeNull();
      await expect(flow.exchangeCode(args)).rejects.toThrow(/already used/);
      expect(await flow.verifyAccessToken(first.access_token, RESOURCE)).toBeNull();
    });

    it("binds a token to its resource", async () => {
      const t = await tokenFor("uA", ["projects:read"], { all: true });

      expect(await flow.verifyAccessToken(t.access_token, RESOURCE)).not.toBeNull();
      expect(await flow.verifyAccessToken(t.access_token, `${ISSUER}/api/mcp`)).toBeNull();
    });

    it("rotates refresh tokens and treats reuse as theft", async () => {
      const t = await tokenFor("uA", ["projects:read"], { all: true });
      const next = await flow.refresh({ refreshToken: t.refresh_token, clientId: t.clientId, scope: null, resource: RESOURCE });

      expect(next.refresh_token).not.toBe(t.refresh_token);
      expect(await flow.verifyAccessToken(next.access_token, RESOURCE)).not.toBeNull();
      await expect(
        flow.refresh({ refreshToken: t.refresh_token, clientId: t.clientId, scope: null, resource: null }),
      ).rejects.toThrow(/already used/);
      expect(await flow.verifyAccessToken(next.access_token, RESOURCE)).toBeNull();
      await expect(
        flow.refresh({ refreshToken: next.refresh_token, clientId: t.clientId, scope: null, resource: null }),
      ).rejects.toThrow(/revoked/);
    });

    it("a refresh may narrow scopes, never widen them", async () => {
      const t = await tokenFor("uA", ["projects:read"], { all: true });

      await expect(
        flow.refresh({ refreshToken: t.refresh_token, clientId: t.clientId, scope: "datasets:write", resource: null }),
      ).rejects.toThrow(/narrow/);
      const narrowed = await flow.refresh({ refreshToken: t.refresh_token, clientId: t.clientId, scope: "catalog:read", resource: null });

      expect(narrowed.scope).toBe("catalog:read");
    });

    it("expired and revoked tokens stop working", async () => {
      const t = await tokenFor("uA", ["projects:read"], { all: true });
      const hash = createHash("sha256").update(t.access_token).digest("hex");

      await prisma.oAuthToken.update({ where: { tokenHash: hash }, data: { expiresAt: new Date(Date.now() - 1000) } });
      expect(await flow.verifyAccessToken(t.access_token, RESOURCE)).toBeNull();

      const u = await tokenFor("uA", ["projects:read"], { all: true });

      await (await route("oauth/revoke")).POST(
        new Request(`${ISSUER}/api/oauth/revoke`, {
          method: "POST",
          headers: { "content-type": "application/x-www-form-urlencoded" },
          body: new URLSearchParams({ token: u.refresh_token }),
        }),
      );
      expect(await flow.verifyAccessToken(u.access_token, RESOURCE)).toBeNull();
    });

    it("disconnecting in /account revokes the grant", async () => {
      const t = await tokenFor("uB", ["projects:read"], { all: true });
      const grant = await prisma.oAuthGrant.findFirstOrThrow({ where: { clientId: t.clientId } });
      const del = (await route("account/agents/[id]")).DELETE;

      session.mockResolvedValue({ user: { id: "uA", role: "USER" } });
      expect((await del(new Request(`${ISSUER}/x`), { params: Promise.resolve({ id: grant.id }) })).status).toBe(404);
      session.mockResolvedValue({ user: { id: "uB", role: "USER" } });
      const list = await (await (await route("account/agents")).GET()).json();

      expect(list.map((a: { id: string }) => a.id)).toContain(grant.id);
      expect((await del(new Request(`${ISSUER}/x`), { params: Promise.resolve({ id: grant.id }) })).status).toBe(200);
      expect(await flow.verifyAccessToken(t.access_token, RESOURCE)).toBeNull();
    });

    it("the token endpoint speaks RFC 6749 over HTTP", async () => {
      const post = (body: Record<string, string>) =>
        route("oauth/token").then((m) =>
          m.POST(
            new Request(`${ISSUER}/api/oauth/token`, {
              method: "POST",
              headers: { "content-type": "application/x-www-form-urlencoded" },
              body: new URLSearchParams(body),
            }),
          ),
        );
      const t = await tokenFor("uA", ["projects:read"], { all: true });
      const ok = await post({ grant_type: "refresh_token", refresh_token: t.refresh_token, client_id: t.clientId });

      expect(ok.status).toBe(200);
      expect(ok.headers.get("cache-control")).toBe("no-store");
      expect(await ok.json()).toMatchObject({ token_type: "Bearer", expires_in: 3600 });
      const bad = await post({ grant_type: "password", client_id: t.clientId });

      expect(bad.status).toBe(400);
      expect(await bad.json()).toMatchObject({ error: "unsupported_grant_type" });
      expect((await post({ grant_type: "refresh_token", refresh_token: "x" })).status).toBe(401);
    });
  });

  describe("consent decision (step 3)", () => {
    const decide = async (form: Record<string, string>, origin = ISSUER) =>
      (await route("oauth/authorize")).POST(
        new Request(`${ISSUER}/api/oauth/authorize`, {
          method: "POST",
          headers: { origin, "content-type": "application/x-www-form-urlencoded" },
          body: new URLSearchParams(form),
        }),
      );

    it("accepts only the signed request, from the site, for the person it was shown to", async () => {
      const { sign } = await import("@/lib/oauth/crypto");
      const client = await registerClient({ redirect_uris: ["http://127.0.0.1/cb"] });
      const v = await flow.validateAuthorize(
        new URLSearchParams({
          response_type: "code", client_id: client.client_id, redirect_uri: "http://127.0.0.1:9/cb", state: "st",
          code_challenge: createHash("sha256").update("w".repeat(43)).digest("base64url"), code_challenge_method: "S256",
        }),
      );

      if (!("request" in v)) throw new Error("refused");
      const signed = sign({ req: v.request, userId: "uA", exp: Date.now() + 60_000 });

      session.mockResolvedValue({ user: { id: "uA", role: "USER" } });
      expect((await decide({ request: signed, decision: "approve" }, "https://evil.example")).status).toBe(403);
      expect((await decide({ request: signed.slice(0, -2) + "xx", decision: "approve" })).status).toBe(400);
      expect((await decide({ request: sign({ req: v.request, userId: "uA", exp: Date.now() - 1 }), decision: "approve" })).status).toBe(400);
      session.mockResolvedValue({ user: { id: "uB", role: "USER" } });
      expect((await decide({ request: signed, decision: "approve" })).status).toBe(403);
      session.mockResolvedValue({ user: { id: "uA", role: "USER" } });

      const denied = await decide({ request: signed, decision: "deny" });

      expect(denied.status).toBe(303);
      expect(new URL(denied.headers.get("location")!).searchParams.get("error")).toBe("access_denied");
      const ok = await decide({ request: signed, decision: "approve", reach: "some", projectIds: "pA1" });
      const loc = new URL(ok.headers.get("location")!);

      expect(ok.status).toBe(303);
      expect(loc.searchParams.get("code")).toMatch(/^mec_/);
      expect(loc.searchParams.get("state")).toBe("st");
      const grant = await prisma.oAuthGrant.findFirstOrThrow({ where: { clientId: client.client_id } });

      expect(grant).toMatchObject({ allProjects: false, projectIds: ["pA1"] });
    });
  });

  describe("rate limits (D14)", () => {
    it("count per key and window, and say how long to wait", async () => {
      const { checkRateLimits } = await import("@/lib/rate-limit");
      const t = new Date("2026-11-09T12:00:30Z");
      const limit = [{ key: "test:a", limit: 3, windowSeconds: 60 }];

      for (let i = 0; i < 3; i++) expect((await checkRateLimits(limit, t)).ok).toBe(true);
      const over = await checkRateLimits(limit, t);

      expect(over).toEqual({ ok: false, key: "test:a", retryAfter: 30 });
      expect((await checkRateLimits([{ ...limit[0], key: "test:b" }], t)).ok).toBe(true);
      expect((await checkRateLimits(limit, new Date("2026-11-09T12:01:00Z"))).ok).toBe(true);
    });

    it("the agent API refuses the 61st anonymous call a minute from one address", async () => {
      const facets = () => call("/api/agent/v1/catalog/facets", "agent/v1/catalog/facets", { ip: "203.0.113.7" });

      for (let i = 0; i < 60; i++) expect((await facets()).status).toBe(200);
      const res = await facets();

      expect(res.status).toBe(429);
      expect(Number(res.headers.get("retry-after"))).toBeGreaterThan(0);
      expect((await res.json()).error).toBe("rate_limited");
      // another address is not affected
      expect((await call("/api/agent/v1/catalog/facets", "agent/v1/catalog/facets", { ip: "203.0.113.8" })).status).toBe(200);
    });

    it("registration allows ten clients an hour from one address", async () => {
      const register = async () =>
        (await route("oauth/register")).POST(
          new Request(`${ISSUER}/api/oauth/register`, {
            method: "POST",
            headers: { "content-type": "application/json", "x-forwarded-for": "203.0.113.9" },
            body: JSON.stringify({ redirect_uris: ["http://127.0.0.1/cb"] }),
          }),
        );

      for (let i = 0; i < 10; i++) expect((await register()).status).toBe(201);
      const res = await register();

      expect(res.status).toBe(429);
      expect(res.headers.get("retry-after")).toBeTruthy();
    });
  });

  describe("discovery (step 3, §6.1)", () => {
    it("advertises what MCP clients check", async () => {
      const as = await (await (await route("oauth/metadata")).GET()).json();

      expect(as).toMatchObject({
        issuer: ISSUER,
        code_challenge_methods_supported: ["S256"],
        client_id_metadata_document_supported: true,
        authorization_response_iss_parameter_supported: true,
        token_endpoint_auth_methods_supported: ["none"],
      });
      const prm = (await route("oauth/protected-resource")).GET;
      const at = async (q: string) => (await prm(new Request(`${ISSUER}/api/oauth/protected-resource${q}`))).json();

      expect(await at("?path=api/agent/v1")).toMatchObject({ resource: RESOURCE, authorization_servers: [ISSUER] });
      expect((await at("")).resource).toBe(RESOURCE);
      expect((await prm(new Request(`${ISSUER}/x?path=nope`))).status).toBe(404);
    });
  });
});
