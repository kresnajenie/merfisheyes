import {
  test,
  expect,
  startSpatialDataServer,
  storeUrl,
  cellsUrl,
  largeMoleculesUrl,
  SELECTORS,
} from "../../fixtures/merfish-fixtures";
import { type FixtureServer } from "../../helpers/fixture-server";

/**
 * `/open?url=` sends a URL to the viewer that reads what it holds
 * (lib/viewer/viewer-for-url.ts) and passes the other parameters on. The
 * stand-in bucket runs on 127.0.0.1, as `merfisheyes-mcp view` serves a
 * person's own folder: so these also check that a local address is never
 * sent to the site (no registration lookup, no claim banner).
 */

let server: FixtureServer;

test.beforeAll(async () => {
  server = await startSpatialDataServer();
});
test.afterAll(() => server.close());

const cases = () => [
  {
    name: "a SpatialData store",
    url: storeUrl(server, { cells: true }),
    path: "/spatialdata-viewer/from-s3",
    canvas: SELECTORS.scCanvas,
  },
  {
    name: "a MERFISHEYES cells folder",
    url: cellsUrl(server, "o"),
    path: "/viewer/from-s3",
    canvas: SELECTORS.scCanvas,
  },
  {
    name: "a MERFISHEYES molecules folder",
    url: largeMoleculesUrl(server, "o"),
    path: "/sm-viewer/from-s3",
    canvas: SELECTORS.smCanvas,
  },
];

for (const i of [0, 1, 2]) {
  test(`@quick /open sends ${["a SpatialData store", "a cells folder", "a molecules folder"][i]} to its viewer, local address kept private`, async ({
    page,
  }) => {
    const c = cases()[i];
    const lookups: string[] = [];

    page.on("request", (r) => {
      if (r.url().includes("/api/datasets/by-url")) lookups.push(r.url());
    });
    await page.goto(`/open?url=${encodeURIComponent(c.url)}&v=e30`);
    await page.waitForURL((u) => u.pathname === c.path, { timeout: 60_000 });
    expect(new URL(page.url()).searchParams.get("v")).not.toBeNull();
    await expect(page.locator(c.canvas).first()).toBeVisible({
      timeout: 120_000,
    });
    await expect(page.getByText("isn't owned yet")).toHaveCount(0);
    expect(lookups).toEqual([]);
  });
}

test("@quick /open says what to check when nothing is there", async ({
  page,
}) => {
  await page.goto(
    `/open?url=${encodeURIComponent(`${server.url}/nothing/here`)}`,
  );
  await expect(page.getByText("Couldn't open")).toBeVisible({
    timeout: 60_000,
  });
  await expect(page.getByText("merfisheyes-mcp view")).toBeVisible();
});
