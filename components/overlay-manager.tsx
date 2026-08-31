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
    <div className="rounded-xl border border-default-200 p-4">
      <h3 className="text-sm font-semibold">Single-molecule overlay</h3>
      <p className="mt-1 text-xs text-default-500">
        Render one of your single-molecule datasets on top of these cells.
      </p>

      {current ? (
        <div className="mt-3 flex items-center justify-between gap-3">
          <span className="text-sm">
            Overlaid:{" "}
            <strong className="font-medium">{currentTitle}</strong>
          </span>
          <Button
            color="danger"
            isDisabled={busy}
            size="sm"
            variant="light"
            onPress={remove}
          >
            Remove
          </Button>
        </div>
      ) : null}

      {options.length === 0 ? (
        <p className="mt-3 text-xs text-default-400">
          You have no single-molecule datasets to overlay yet.
        </p>
      ) : (
        <div className="mt-3 flex items-center gap-2">
          <select
            className="flex-1 rounded-medium border border-default-200 bg-content1 px-3 py-2 text-sm"
            disabled={busy}
            value={selected}
            onChange={(e) => setSelected(e.target.value)}
          >
            <option value="">
              {current ? "Change to…" : "Choose a single-molecule dataset…"}
            </option>
            {options
              .filter((o) => o.id !== current)
              .map((o) => (
                <option key={o.id} value={o.id}>
                  {o.title ?? o.id}
                </option>
              ))}
          </select>
          <Button
            color="primary"
            isDisabled={busy || !selected}
            size="sm"
            onPress={() => selected && attach(selected)}
          >
            {current ? "Change" : "Attach"}
          </Button>
        </div>
      )}
    </div>
  );
}
