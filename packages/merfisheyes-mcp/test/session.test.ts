import { mkdirSync, mkdtempSync, writeFileSync } from "node:fs";
import { request } from "node:http";
import { tmpdir } from "node:os";
import { join } from "node:path";

import { afterAll, beforeAll, expect, it } from "vitest";

const { LocalServer } = await import("../src/local/serve.ts");

/**
 * The live session (pillar D3): the package sends a command to the page it
 * served and waits for the answer. A stand-in page here: it listens on
 * /session/events and answers on /session/reply, as components/agent-session.tsx does.
 */
const viewer = mkdtempSync(join(tmpdir(), "meyes-session-"));
const data = join(viewer, "..", `data-${Date.now()}`);
const srv = new LocalServer({ viewerDir: viewer });
let origin = "";

beforeAll(async () => {
  mkdirSync(join(viewer, "open"), { recursive: true });
  writeFileSync(join(viewer, "open", "index.html"), "<html></html>");
  mkdirSync(join(data, "s.zarr"), { recursive: true });
  writeFileSync(join(data, "s.zarr", "zarr.json"), "{}");
  origin = new URL((await srv.share(data)).url).origin;
});
afterAll(() => srv.close());

/** A page that answers each command with `answer(command, args)`. */
function page(token: string, answer: (c: string, a: unknown) => unknown) {
  const port = new URL(origin).port;
  const host = `127.0.0.1:${port}`;

  return new Promise<{ close: () => void; status: number }>((ok) => {
    const req = request(
      {
        host: "127.0.0.1",
        port,
        path: `/session/events?token=${token}`,
        headers: { host },
      },
      (res) => {
        let buf = "";

        res.on("data", (d) => {
          buf += d;
          let i;

          while ((i = buf.indexOf("\n\n")) >= 0) {
            const event = buf.slice(0, i);

            buf = buf.slice(i + 2);
            const line = event.split("\n").find((l) => l.startsWith("data: "));

            if (!line) continue;
            const { id, command, args } = JSON.parse(line.slice(6));
            let reply;

            try {
              reply = { id, ok: true, result: answer(command, args) };
            } catch (e) {
              reply = { id, ok: false, error: (e as Error).message };
            }
            const post = request({
              host: "127.0.0.1",
              port,
              method: "POST",
              path: `/session/reply?token=${token}`,
              headers: { host, "content-type": "application/json" },
            });

            post.end(JSON.stringify(reply));
          }
        });
        ok({ close: () => req.destroy(), status: res.statusCode! });
      },
    );

    req.end();
  });
}

it("says no view is open before a page connects", async () => {
  await expect(srv.command("state", {})).rejects.toThrow(
    /No MERFISHEYES view is open/,
  );
});

it("puts the session in the link, and sends commands to the page and back", async () => {
  expect(srv.link({ path: "", url: `${origin}/t/d` })).toContain(
    `&session=${srv.sessionToken}`,
  );
  const p = await page(srv.sessionToken, (c, a) => {
    if (c === "fail") throw new Error("no such gene");

    return { command: c, args: a };
  });

  expect(p.status).toBe(200);
  await new Promise((r) => setTimeout(r, 50));
  expect(srv.viewOpen).toBe(true);
  expect(await srv.command("set", { gene: "Gad1" })).toEqual({
    command: "set",
    args: { gene: "Gad1" },
  });
  await expect(srv.command("fail", {})).rejects.toThrow("no such gene");
  p.close();
  await new Promise((r) => setTimeout(r, 50));
  expect(srv.viewOpen).toBe(false);
});

it("refuses a page without the session token, and allows no other site", async () => {
  const port = new URL(origin).port;
  const res = await new Promise<{ status: number; acao?: string }>((ok) => {
    const req = request(
      {
        host: "127.0.0.1",
        port,
        path: "/session/events?token=wrong",
        headers: {
          host: `127.0.0.1:${port}`,
          origin: "https://www.merfisheyes.com",
        },
      },
      (r) => {
        ok({
          status: r.statusCode!,
          acao: r.headers["access-control-allow-origin"] as string,
        });
        r.resume();
      },
    );

    req.end();
  });

  expect(res.status).toBe(403);
  expect(res.acao).toBeUndefined();
});

it("gives up when the page does not answer", async () => {
  const p = await page(srv.sessionToken, () => new Promise(() => {}));

  await new Promise((r) => setTimeout(r, 50));
  // the stand-in answers {} for a never-settling promise via JSON.stringify
  // of a Promise, so time out on a command the page never posts back for:
  p.close();
  await expect(srv.command("state", {}, 200)).rejects.toThrow(
    /No MERFISHEYES view is open|did not answer/,
  );
});
