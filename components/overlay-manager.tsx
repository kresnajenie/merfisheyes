"use client";

import { Button } from "@heroui/button";
import { Input } from "@heroui/input";
import { Spinner } from "@heroui/spinner";
import { useCallback, useEffect, useState } from "react";
import { toast } from "react-toastify";

import { ImageOverlayPicker } from "@/components/image-overlay-picker";
import { useLocalBundleStore } from "@/lib/stores/localBundleStore";
import { dataUrlOf, moleculeUrlFormat } from "@/lib/overlays/moleculeSource";

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
 * Attach / change / remove the overlays on a single-cell dataset the user
 * owns: a single-molecule dataset drawn over the cells (searchable +
 * paginated so it scales to large libraries — it fetches one page of the
 * user's own single-molecule datasets at a time) or public molecules by
 * their link, and one of their images drawn under them (see
 * ImageOverlayPicker).
 *
 * `local`: a preview not saved yet. Nothing is written; the molecules picked
 * are drawn over it to check they line up (localBundleStore `check`), the
 * way a saved overlay would draw them.
 */
export function OverlayManager({
  scDatasetId,
  local = false,
  onChange,
}: {
  scDatasetId: string | null;
  local?: boolean;
  onChange?: (smId: string | null) => void;
}) {
  const check = useLocalBundleStore((st) => st.check);
  const [saved, setSaved] = useState<string | null>(null);
  const [savedTitle, setSavedTitle] = useState<string | null>(null);
  const current = local
    ? (check?.datasetId ?? check?.links.__all__ ?? null)
    : saved;
  const currentTitle = local ? (checkTitle(check) ?? null) : savedTitle;
  const setCurrent = setSaved;
  const setCurrentTitle = setSavedTitle;
  const [linkDraft, setLinkDraft] = useState("");
  const [search, setSearch] = useState("");
  const [items, setItems] = useState<SmItem[]>([]);
  const [total, setTotal] = useState(0);
  const [page, setPage] = useState(1);
  const [loading, setLoading] = useState(true);
  const [busy, setBusy] = useState(false);

  // Current overlays (once): a dataset, or public molecules by their link
  useEffect(() => {
    if (local || !scDatasetId) return;
    let cancelled = false;

    fetch(`/api/datasets/${scDatasetId}/overlay`)
      .then((r) => (r.ok ? r.json() : null))
      .then((d) => {
        if (cancelled) return;
        setCurrent(d?.smDatasetId ?? d?.smUrl ?? null);
        if (!d?.smDatasetId && d?.smUrl) setCurrentTitle(d.smUrl);
      })
      .catch(() => {});

    return () => {
      cancelled = true;
    };
  }, [scDatasetId]);

  // Paginated + searched list of the user's molecule datasets: single
  // molecule ones, and SpatialData stores (what the server makes of a molecule
  // table; the route checks the one picked has molecules).
  const fetchPage = useCallback(async () => {
    setLoading(true);
    try {
      const params = new URLSearchParams({
        type: "single_molecule,spatialdata",
        page: String(page),
        limit: String(PAGE_SIZE),
      });

      if (search.trim()) params.set("search", search.trim());
      const res = await fetch(`/api/ingest/mine?${params}`);
      const data = res.ok ? await res.json() : null;

      setItems(
        ((data?.datasets ?? []) as typeof items).filter(
          (d) => d.id !== scDatasetId,
        ),
      );
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
      if (local) {
        // How a saved overlay would link it, drawn over the preview only
        const res = await fetch(`/api/datasets/${smId}/molecule-link`);

        if (!res.ok) {
          toast.error("Couldn't open those molecules.");

          return;
        }
        useLocalBundleStore
          .getState()
          .set({ check: { ...(await res.json()), title: title ?? undefined } });
        toast.success("Drawn over this preview to check. Not saved.");

        return;
      }
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

  // Public molecules by their link: a store (or the folder holding it) with
  // a points element, a merfisheyes molecule folder, or a viewer link that
  // opened one from S3
  const addByLink = async (force = false) => {
    const url = dataUrlOf(linkDraft);

    if (!url) {
      toast.error("Paste a link starting with https://");

      return;
    }
    setBusy(true);
    try {
      if (local) {
        // Checked here the way the viewer will open it
        const linkFormat = await moleculeUrlFormat(url);

        if (!linkFormat) {
          toast.error(NO_MOLECULES);

          return;
        }
        useLocalBundleStore.getState().set({
          check: {
            linkColumn: "__all__",
            links: { __all__: url },
            source: "s3",
            linkFormat,
            title: url,
          },
        });
        setLinkDraft("");
        toast.success("Drawn over this preview to check. Not saved.");

        return;
      }
      const res = await fetch(`/api/datasets/${scDatasetId}/overlay`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ url, force }),
      });

      if (res.status === 409) {
        if (confirm("This dataset already has an overlay. Replace it?")) {
          return addByLink(true);
        }

        return;
      }
      if (!res.ok) {
        const j = await res.json().catch(() => ({}));

        toast.error(j.error ?? NO_MOLECULES);

        return;
      }
      setCurrent(url);
      setCurrentTitle(url);
      setLinkDraft("");
      onChange?.(null);
      toast.success("Molecules overlaid by their link.");
    } finally {
      setBusy(false);
    }
  };

  const remove = async () => {
    if (local) {
      useLocalBundleStore.getState().set({ check: null });

      return;
    }
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
          {local
            ? "Check how molecules line up with this preview before you save it. Nothing here is saved."
            : "Render one of your single-molecule datasets, or public molecules by their link, on top of these cells."}
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
        placeholder="Search your molecule datasets…"
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
              ? "No molecule datasets match."
              : "You have no molecule datasets to overlay yet."}
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

      <form
        className="flex gap-2"
        onSubmit={(e) => {
          e.preventDefault();
          void addByLink();
        }}
      >
        <Input
          aria-label="Molecules by link"
          isDisabled={busy}
          placeholder="Add by link: https://bucket.s3.region.amazonaws.com/folder or …/store.zarr"
          size="sm"
          value={linkDraft}
          onValueChange={setLinkDraft}
        />
        <Button
          isDisabled={busy || !linkDraft.trim()}
          size="sm"
          type="submit"
          variant="flat"
        >
          Add
        </Button>
      </form>

      {!local && scDatasetId ? (
        <div className="mt-2 border-t border-default-200 pt-4">
          <ImageOverlayPicker scDatasetId={scDatasetId} />
        </div>
      ) : null}
    </div>
  );
}

const NO_MOLECULES =
  "No molecules at that link. Paste a public SpatialData store with a points element (or the folder holding it), or a merfisheyes molecule folder, that allows this site in CORS.";

/** What to call the molecules a preview is checked against. */
function checkTitle(
  check: {
    title?: string;
    datasetId?: string;
    links: Record<string, string>;
  } | null,
): string | null {
  if (!check) return null;

  return check.title ?? check.datasetId ?? check.links.__all__ ?? null;
}
