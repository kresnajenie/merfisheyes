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
 * (lib/viewer/viewer-for-url.ts) and passes the other parameters on. (Whether
 * a local address is kept from the site depends on where the page itself
 * comes from — tests/unit/local-url.test.ts; here both are on this computer.)
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
  test(`@quick /open sends ${["a SpatialData store", "a cells folder", "a molecules folder"][i]} to its viewer, and passes its settings on`, async ({
    page,
  }) => {
    const c = cases()[i];

    await page.goto(`/open?url=${encodeURIComponent(c.url)}&v=e30`);
    await page.waitForURL((u) => u.pathname === c.path, { timeout: 60_000 });
    expect(new URL(page.url()).searchParams.get("v")).not.toBeNull();
    await expect(page.locator(c.canvas).first()).toBeVisible({
      timeout: 120_000,
    });
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
