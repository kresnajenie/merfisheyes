"use client";

import { Button } from "@heroui/button";
import { useEffect, useState } from "react";
import { toast } from "react-toastify";

interface SmOption {
  id: string;
  title: string | null;
}

/**
 * Attach / change / remove a single-molecule overlay on a single-cell dataset
 * the user owns. Writes datasets/{scId}/mapping.json via the overlay API and
 * groups the two in a Project. Self-contained: fetches the user's own
 * single-molecule datasets and the current overlay on mount.
 */
export function OverlayManager({
  scDatasetId,
  onChange,
}: {
  scDatasetId: string;
  onChange?: (smId: string | null) => void;
}) {
  const [current, setCurrent] = useState<string | null>(null);
  const [options, setOptions] = useState<SmOption[]>([]);
  const [selected, setSelected] = useState("");
  const [busy, setBusy] = useState(false);
  const [loaded, setLoaded] = useState(false);

  useEffect(() => {
    let cancelled = false;

    (async () => {
      try {
        const [ov, mine] = await Promise.all([
          fetch(`/api/datasets/${scDatasetId}/overlay`).then((r) =>
            r.ok ? r.json() : null,
          ),
          fetch("/api/ingest/mine?type=single_molecule&limit=100").then((r) =>
            r.ok ? r.json() : null,
          ),
        ]);

        if (cancelled) return;
        setCurrent(ov?.smDatasetId ?? null);
        setOptions(
          (mine?.datasets ?? []).map((d: any) => ({ id: d.id, title: d.title })),
        );
      } finally {
        if (!cancelled) setLoaded(true);
      }
    })();

    return () => {
      cancelled = true;
    };
  }, [scDatasetId]);

  const attach = async (smId: string, force = false) => {
    setBusy(true);
    try {
      const res = await fetch(`/api/datasets/${scDatasetId}/overlay`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ smDatasetId: smId, force }),
      });

      if (res.status === 409) {
        if (confirm("This dataset already has an overlay. Replace it?")) {
          return attach(smId, true);
        }

        return;
      }
      if (!res.ok) {
        const j = await res.json().catch(() => ({}));

        toast.error(j.error ?? "Couldn't attach overlay.");

        return;
      }
      setCurrent(smId);
      setSelected("");
      onChange?.(smId);
      toast.success("Single-molecule overlay attached.");
    } finally {
      setBusy(false);
    }
  };

  const remove = async () => {
    if (!confirm("Remove the single-molecule overlay from this dataset?"))
      return;
    setBusy(true);
    try {
      const res = await fetch(`/api/datasets/${scDatasetId}/overlay`, {
        method: "DELETE",
      });

      if (!res.ok) {
        toast.error("Couldn't remove overlay.");

        return;
      }
      setCurrent(null);
      onChange?.(null);
      toast.success("Overlay removed.");
    } finally {
      setBusy(false);
    }
  };

  if (!loaded) return null;

  const currentTitle =
    options.find((o) => o.id === current)?.title ?? current ?? "";

  return (
    <div className="w-full max-w-md space-y-3">
      <div>
        <h3 className="text-sm font-semibold text-foreground">
          Single-molecule overlay
        </h3>
        <p className="mt-0.5 text-xs text-default-500">
          Render one of your single-molecule datasets on top of these cells.
        </p>
      </div>

      {current ? (
        <div className="flex items-center justify-between gap-3 rounded-lg bg-content2 px-3 py-2">
          <div className="min-w-0">
            <span className="text-[11px] uppercase tracking-wide text-default-400">
              Overlaid
            </span>
            <p className="truncate text-sm font-medium text-foreground">
              {currentTitle}
            </p>
          </div>
          <Button
            className="shrink-0"
            color="danger"
            isDisabled={busy}
            size="sm"
            variant="flat"
            onPress={remove}
          >
            Remove
          </Button>
        </div>
      ) : null}

      {options.length === 0 ? (
        <p className="text-xs text-default-400">
          You have no single-molecule datasets to overlay yet.
        </p>
      ) : (
        <div className="space-y-2">
          <select
            className="w-full rounded-medium border border-default-200 bg-content2 px-3 py-2 text-sm text-foreground focus:border-primary focus:outline-none"
            disabled={busy}
            value={selected}
            onChange={(e) => setSelected(e.target.value)}
          >
            <option value="">
              {current
                ? "Change to a different dataset…"
                : "Choose a single-molecule dataset…"}
            </option>
            {options
              .filter((o) => o.id !== current)
              .map((o) => (
                <option key={o.id} value={o.id}>
                  {o.title ?? o.id}
                </option>
              ))}
          </select>
          <div className="flex justify-end">
            <Button
              color="primary"
              isDisabled={busy || !selected}
              isLoading={busy}
              size="sm"
              onPress={() => selected && attach(selected)}
            >
              {current ? "Change overlay" : "Attach overlay"}
            </Button>
          </div>
        </div>
      )}
    </div>
  );
}
