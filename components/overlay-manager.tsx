"use client";

import { Button } from "@heroui/button";
import { Input } from "@heroui/input";
import { Spinner } from "@heroui/spinner";
import { useCallback, useEffect, useState } from "react";
import { toast } from "react-toastify";

interface SmItem {
  id: string;
  title: string | null;
  numCells: number | null; // molecule count for SM datasets
}

const PAGE_SIZE = 8;

function fmtCount(n: number | null): string {
  if (!n) return "";
  if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(1)}M molecules`;
  if (n >= 1_000) return `${Math.round(n / 1_000)}K molecules`;
  return `${n} molecules`;
}

/**
 * Attach / change / remove a single-molecule overlay on a single-cell dataset
 * the user owns. Searchable + paginated so it scales to large libraries (it
 * fetches one page of the user's own single-molecule datasets at a time).
 */
export function OverlayManager({
  scDatasetId,
  onChange,
}: {
  scDatasetId: string;
  onChange?: (smId: string | null) => void;
}) {
  const [current, setCurrent] = useState<string | null>(null);
  const [currentTitle, setCurrentTitle] = useState<string | null>(null);
  const [search, setSearch] = useState("");
  const [items, setItems] = useState<SmItem[]>([]);
  const [total, setTotal] = useState(0);
  const [page, setPage] = useState(1);
  const [loading, setLoading] = useState(true);
  const [busy, setBusy] = useState(false);

  // Current overlay (once).
  useEffect(() => {
    let cancelled = false;

    fetch(`/api/datasets/${scDatasetId}/overlay`)
      .then((r) => (r.ok ? r.json() : null))
      .then((d) => !cancelled && setCurrent(d?.smDatasetId ?? null))
      .catch(() => {});

    return () => {
      cancelled = true;
    };
  }, [scDatasetId]);

  // Paginated + searched list of the user's single-molecule datasets.
  const fetchPage = useCallback(async () => {
    setLoading(true);
    try {
      const params = new URLSearchParams({
        type: "single_molecule",
        page: String(page),
        limit: String(PAGE_SIZE),
      });

      if (search.trim()) params.set("search", search.trim());
      const res = await fetch(`/api/ingest/mine?${params}`);
      const data = res.ok ? await res.json() : null;

      setItems(data?.datasets ?? []);
      setTotal(data?.total ?? 0);
    } finally {
      setLoading(false);
    }
  }, [page, search]);

  useEffect(() => {
    const t = setTimeout(fetchPage, 250);

    return () => clearTimeout(t);
  }, [fetchPage]);

  // Reset to page 1 when the search changes.
  useEffect(() => {
    setPage(1);
  }, [search]);

  // Keep the current overlay's title in sync as pages load.
  useEffect(() => {
    if (current) {
      const hit = items.find((i) => i.id === current);

      if (hit) setCurrentTitle(hit.title);
    }
  }, [current, items]);

  const attach = async (smId: string, title: string | null, force = false) => {
    setBusy(true);
    try {
      const res = await fetch(`/api/datasets/${scDatasetId}/overlay`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ smDatasetId: smId, force }),
      });

      if (res.status === 409) {
        if (confirm("This dataset already has an overlay. Replace it?")) {
          return attach(smId, title, true);
        }

        return;
      }
      if (!res.ok) {
        const j = await res.json().catch(() => ({}));

        toast.error(j.error ?? "Couldn't attach overlay.");

        return;
      }
      setCurrent(smId);
      setCurrentTitle(title);
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
      setCurrentTitle(null);
      onChange?.(null);
      toast.success("Overlay removed.");
    } finally {
      setBusy(false);
    }
  };

  const totalPages = Math.max(1, Math.ceil(total / PAGE_SIZE));

  return (
    <div className="flex w-full flex-col gap-3">
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
              Currently overlaid
            </span>
            <p className="truncate text-sm font-medium text-foreground">
              {currentTitle ?? current}
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

      <Input
        isClearable
        placeholder="Search your single-molecule datasets…"
        size="sm"
        value={search}
        onClear={() => setSearch("")}
        onValueChange={setSearch}
      />

      <div className="min-h-[220px]">
        {loading ? (
          <div className="flex h-[220px] items-center justify-center">
            <Spinner size="sm" />
          </div>
        ) : items.length === 0 ? (
          <p className="py-8 text-center text-xs text-default-400">
            {search
              ? "No single-molecule datasets match."
              : "You have no single-molecule datasets to overlay yet."}
          </p>
        ) : (
          <div className="flex flex-col gap-1.5">
            {items.map((it) => {
              const isCurrent = it.id === current;

              return (
                <button
                  key={it.id}
                  className={`flex items-center justify-between gap-3 rounded-lg border px-3 py-2 text-left transition-colors ${
                    isCurrent
                      ? "border-primary bg-primary/10"
                      : "border-default-200 hover:bg-content2"
                  }`}
                  disabled={busy}
                  type="button"
                  onClick={() => !isCurrent && attach(it.id, it.title)}
                >
                  <div className="min-w-0">
                    <p className="truncate text-sm font-medium text-foreground">
                      {it.title ?? it.id}
                    </p>
                    {fmtCount(it.numCells) ? (
                      <p className="text-[11px] text-default-400">
                        {fmtCount(it.numCells)}
                      </p>
                    ) : null}
                  </div>
                  {isCurrent ? (
                    <span className="shrink-0 text-[11px] font-medium text-primary">
                      Overlaid
                    </span>
                  ) : (
                    <span className="shrink-0 text-[11px] text-default-400">
                      Attach →
                    </span>
                  )}
                </button>
              );
            })}
          </div>
        )}
      </div>

      {total > PAGE_SIZE ? (
        <div className="flex items-center justify-between text-xs text-default-400">
          <Button
            isDisabled={page <= 1 || loading}
            size="sm"
            variant="flat"
            onPress={() => setPage((p) => Math.max(1, p - 1))}
          >
            Previous
          </Button>
          <span>
            Page {page} of {totalPages} · {total} datasets
          </span>
          <Button
            isDisabled={page >= totalPages || loading}
            size="sm"
            variant="flat"
            onPress={() => setPage((p) => Math.min(totalPages, p + 1))}
          >
            Next
          </Button>
        </div>
      ) : null}
    </div>
  );
}
