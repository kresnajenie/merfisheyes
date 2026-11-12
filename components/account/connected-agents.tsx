"use client";

import { Button } from "@heroui/button";
import { useCallback, useEffect, useState } from "react";

import { SCOPES } from "@/lib/oauth/config";

interface Agent {
  id: string;
  name: string;
  verified: boolean;
  publisher: string | null;
  scopes: (keyof typeof SCOPES)[];
  projects: "all" | string[];
  connectedAt: string;
  lastUsedAt: string | null;
}

const when = (iso: string | null) =>
  iso ? new Date(iso).toLocaleDateString() : "never";

/** The agents this person allowed (docs/specs/MCP-PROTOCOL.md §6.4), each
 * with what it may do and a Disconnect button. Renders nothing when there
 * are none. Kept short: name, what it can do, which projects, last used. */
export function ConnectedAgents() {
  const [agents, setAgents] = useState<Agent[] | null>(null);
  const [busy, setBusy] = useState<string | null>(null);

  const load = useCallback(() => {
    fetch("/api/account/agents")
      .then((r) => (r.ok ? r.json() : []))
      .then(setAgents)
      .catch(() => setAgents([]));
  }, []);

  useEffect(load, [load]);

  if (!agents || agents.length === 0) return null;

  const disconnect = async (id: string) => {
    setBusy(id);
    await fetch(`/api/account/agents/${id}`, { method: "DELETE" });
    setBusy(null);
    load();
  };

  return (
    <section className="mt-12">
      <h2 className="text-lg font-semibold">
        AI agents connected to your account
      </h2>
      <ul className="mt-4 divide-y divide-default-200 rounded-xl border border-default-200">
        {agents.map((a) => (
          <li key={a.id} className="flex items-start justify-between gap-4 p-4">
            <div className="min-w-0 text-sm">
              <p className="font-medium">
                {a.name}
                {a.publisher && (
                  <span className="ml-2 font-normal text-default-500">
                    {a.publisher}
                  </span>
                )}
                {!a.verified && (
                  <span className="ml-2 font-normal text-warning-600">
                    unverified
                  </span>
                )}
              </p>
              <p className="mt-1 text-default-600">
                Permissions:{" "}
                {a.scopes
                  .filter((s) => s !== "catalog:read")
                  .map((s) => SCOPES[s] ?? s)
                  .join(" · ") || "Search the public catalog"}
              </p>
              {a.scopes.includes("projects:read") && (
                <p className="mt-1 text-default-600">
                  Projects:{" "}
                  {a.projects === "all"
                    ? "All your projects"
                    : a.projects.join(", ")}
                </p>
              )}
              <p className="mt-1 text-default-500">
                {a.lastUsedAt
                  ? `Last used ${when(a.lastUsedAt)}`
                  : "Not used yet"}
              </p>
            </div>
            <Button
              color="danger"
              isLoading={busy === a.id}
              size="sm"
              variant="flat"
              onPress={() => disconnect(a.id)}
            >
              Disconnect
            </Button>
          </li>
        ))}
      </ul>
    </section>
  );
}
