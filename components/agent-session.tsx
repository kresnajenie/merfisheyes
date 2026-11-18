"use client";

import { useEffect } from "react";

import { runViewCommand } from "@/lib/agent/view-commands";

/**
 * The page's end of a live session with an agent (docs/plans/MCP.md pillar
 * D3). When the viewer was opened with `?session=<token>` by the package
 * that serves it (`merfisheyes view`, meyes_open_local), it listens on that
 * same address for commands (`/session/events`, server-sent events), runs
 * them against the open view (lib/agent/view-commands.ts) and posts each
 * answer back (`/session/reply`). Nothing without a token; same origin only.
 */
export function AgentSession() {
  useEffect(() => {
    const token = new URLSearchParams(window.location.search).get("session");

    if (!token) return;
    const q = `token=${encodeURIComponent(token)}`;
    let events: EventSource;
    const onMessage = async (event: MessageEvent) => {
      const { id, command, args } = JSON.parse(event.data);
      let reply: unknown;

      try {
        reply = { id, ok: true, result: await runViewCommand(command, args) };
      } catch (e) {
        reply = { id, ok: false, error: (e as Error).message };
      }
      await fetch(`/session/reply?${q}`, {
        method: "POST",
        headers: { "content-type": "application/json" },
        body: JSON.stringify(reply),
      }).catch(() => {});
    };

    // The package sends commands to the tab that connected last: a tab
    // brought to the front reconnects, so the one the person is looking at
    // is the one the agent drives.
    const connect = () => {
      events?.close();
      events = new EventSource(`/session/events?${q}`);
      events.onmessage = onMessage;
    };
    const front = () => {
      if (document.visibilityState === "visible") connect();
    };

    connect();
    window.addEventListener("focus", front);
    document.addEventListener("visibilitychange", front);

    return () => {
      window.removeEventListener("focus", front);
      document.removeEventListener("visibilitychange", front);
      events.close();
    };
  }, []);

  return null;
}
