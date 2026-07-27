"use client";

import { useCallback, useEffect, useState } from "react";
import Link from "next/link";
import { CheckCircle2, Loader2, XCircle } from "lucide-react";

interface Upload {
  id: string;
  title: string | null;
  status: string;
  datasetType: string | null;
  progress: { stage?: string; percent?: number } | null;
  errorMessage: string | null;
  numCells: number;
  numGenes: number;
  createdAt: string;
  completedAt: string | null;
  viewerUrl: string | null;
}

const IN_PROGRESS = new Set(["UPLOADING", "QUEUED", "PROCESSING"]);

/**
 * "Your uploads" strip at the top of /explore. Shows the signed-in user's own
 * datasets — the destination a server-side upload redirects to, so the in-flight
 * dataset is visible processing here instead of behind a full-screen overlay on
 * the homepage. Polls while anything is still processing; self-hides when the
 * user has no uploads.
 */
export function YourUploads({ isSignedIn }: { isSignedIn: boolean }) {
  const [uploads, setUploads] = useState<Upload[] | null>(null);

  const load = useCallback(async () => {
    try {
      const r = await fetch("/api/ingest/mine", { cache: "no-store" });

      if (!r.ok) {
        setUploads([]);

        return;
      }
      const j = await r.json();

      setUploads(j.datasets ?? []);
    } catch {
      setUploads([]);
    }
  }, []);

  useEffect(() => {
    if (isSignedIn) load();
  }, [isSignedIn, load]);

  // Keep polling while any upload is still moving through the pipeline.
  useEffect(() => {
    if (!uploads || !uploads.some((u) => IN_PROGRESS.has(u.status))) return;
    const t = setInterval(load, 4000);

    return () => clearInterval(t);
  }, [uploads, load]);

  if (!isSignedIn || !uploads || uploads.length === 0) return null;

  return (
    <section className="mb-10">
      <h2 className="mb-3 text-lg font-semibold text-foreground">
        Your uploads
      </h2>
      <div className="grid grid-cols-1 gap-3 sm:grid-cols-2 lg:grid-cols-3">
        {uploads.map((u) => (
          <UploadCard key={u.id} upload={u} />
        ))}
      </div>
    </section>
  );
}

function UploadCard({ upload: u }: { upload: Upload }) {
  const inProgress = IN_PROGRESS.has(u.status);
  const name = u.title || "Untitled dataset";
  const isMolecule = u.datasetType === "single_molecule";
  const countLabel = isMolecule ? "molecules" : "cells";

  const body = (
    <div
      className={`rounded-2xl border p-4 transition-colors ${
        u.status === "FAILED"
          ? "border-danger/40 bg-danger/5"
          : "border-default-200 bg-default-100/40 hover:border-primary/40"
      }`}
    >
      <div className="flex items-start gap-2">
        {u.status === "COMPLETE" ? (
          <CheckCircle2 className="mt-0.5 h-5 w-5 shrink-0 text-emerald-400" />
        ) : u.status === "FAILED" ? (
          <XCircle className="mt-0.5 h-5 w-5 shrink-0 text-danger" />
        ) : (
          <Loader2 className="mt-0.5 h-5 w-5 shrink-0 animate-spin text-blue-400" />
        )}
        <div className="min-w-0 flex-1">
          <p className="truncate text-base font-semibold text-foreground">
            {name}
          </p>
          <p className="mt-0.5 text-sm text-default-500">
            {u.status === "COMPLETE"
              ? `${u.numCells.toLocaleString()} ${countLabel} · ${u.numGenes.toLocaleString()} genes`
              : u.status === "FAILED"
                ? "Failed"
                : u.progress?.stage || statusLabel(u.status)}
          </p>
        </div>
      </div>

      {inProgress && (
        <div className="mt-3">
          <div className="h-2 w-full overflow-hidden rounded-full bg-default-200">
            {typeof u.progress?.percent === "number" ? (
              <div
                className="h-full rounded-full bg-blue-400 transition-all duration-300"
                style={{ width: `${Math.min(100, u.progress.percent)}%` }}
              />
            ) : (
              <div
                className="h-full w-1/4 rounded-full bg-blue-400"
                style={{
                  animation: "glassIndeterminate 1.2s ease-in-out infinite",
                }}
              />
            )}
          </div>
        </div>
      )}

      {u.status === "FAILED" && u.errorMessage && (
        <p className="mt-2 break-words text-sm text-default-500">
          {u.errorMessage}
        </p>
      )}

      {u.status === "COMPLETE" && (
        <p className="mt-2 text-sm font-medium text-primary">Open dataset →</p>
      )}
    </div>
  );

  // A finished dataset is a link into its viewer; everything else is static.
  return u.status === "COMPLETE" && u.viewerUrl ? (
    <Link className="block" href={u.viewerUrl}>
      {body}
    </Link>
  ) : (
    body
  );
}

function statusLabel(status: string): string {
  switch (status) {
    case "UPLOADING":
      return "Uploading…";
    case "QUEUED":
      return "Queued…";
    case "PROCESSING":
      return "Processing…";
    default:
      return status;
  }
}
