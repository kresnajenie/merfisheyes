"use client";

import { Button } from "@heroui/button";
import { Input } from "@heroui/input";
import { Spinner } from "@heroui/spinner";
import { useCallback, useEffect, useState } from "react";
import { toast } from "react-toastify";

interface ImageItem {
  id: string;
  title: string | null;
}

const PAGE_SIZE = 6;

/**
 * Choose the image drawn under a single cell dataset's cells: browse the
 * user's own image datasets (searchable, a page at a time), or add one by
 * URL — which first adds it to their account, like any dataset opened by URL.
 */
export function ImageOverlayPicker({ scDatasetId }: { scDatasetId: string }) {
  // What the dataset's mapping.json names now: a URL or a dataset id, and
  // the dataset id when it is a registered image.
  const [current, setCurrent] = useState<{
    image: string;
    datasetId: string | null;
  } | null>(null);
  const [currentTitle, setCurrentTitle] = useState<string | null>(null);
  const [search, setSearch] = useState("");
  const [items, setItems] = useState<ImageItem[]>([]);
  const [total, setTotal] = useState(0);
  const [page, setPage] = useState(1);
  const [loading, setLoading] = useState(true);
  const [busy, setBusy] = useState(false);
  const [urlDraft, setUrlDraft] = useState("");

  useEffect(() => {
    let cancelled = false;

    fetch(`/api/datasets/${scDatasetId}/overlay`)
      .then((r) => (r.ok ? r.json() : null))
      .then((d) => {
        if (cancelled || !d?.image) return;
        setCurrent({ image: d.image, datasetId: d.imageDatasetId ?? null });
      })
      .catch(() => {});

    return () => {
      cancelled = true;
    };
  }, [scDatasetId]);

  const fetchPage = useCallback(async () => {
    setLoading(true);
    try {
      const params = new URLSearchParams({
        type: "image",
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

  useEffect(() => {
    setPage(1);
  }, [search]);

  // Name the current image once a page that holds it has loaded.
  useEffect(() => {
    const hit = items.find((i) => i.id === current?.datasetId);

    if (hit) setCurrentTitle(hit.title);
  }, [current, items]);

  /** Save the overlay; resolves to what the server stored, or null on failure. */
  const save = async (
    body: { imageDatasetId: string } | { image: string | null },
  ) => {
    const res = await fetch(`/api/datasets/${scDatasetId}/overlay`, {
      method: "PATCH",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
    const j = await res.json().catch(() => ({}));

    if (!res.ok) {
      toast.error(j.error ?? "Couldn't update the image overlay.");

      return null;
    }

    return j as { image: string | null; imageDatasetId: string | null };
  };

  const attach = async (item: ImageItem) => {
    setBusy(true);
    try {
      const saved = await save({ imageDatasetId: item.id });

      if (!saved?.image) return;
      setCurrent({ image: saved.image, datasetId: item.id });
      setCurrentTitle(item.title);
      toast.success("Image overlay attached. Reload the viewer to see it.");
    } finally {
      setBusy(false);
    }
  };

  const addByUrl = async () => {
    const url = urlDraft.trim().replace(/\/+$/, "");

    if (!url) return;
    setBusy(true);
    try {
      // Open it the way the viewer will, so a wrong URL, a private bucket or
      // an unsupported image is caught here rather than as a missing layer.
      try {
        const { openOmeZarrImage } = await import("@/lib/spatialdata/image");

        await openOmeZarrImage(url);
      } catch {
        toast.error(
          "Couldn't read an image at that URL. It must be a public OME-Zarr image (Zarr v3, axes c/y/x) that allows this site in CORS.",
        );

        return;
      }

      // Add it to the account, so it is listed here from now on. An image
      // someone else already owns can still be drawn, by its URL.
      const reg = await fetch("/api/datasets/register-s3", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ url, kind: "image" }),
      })
        .then((r) => (r.ok ? r.json() : null))
        .catch(() => null);
      const mine = reg?.dataset?.datasetType === "image" && !reg.ownedByOther;
      const saved = await save(
        mine ? { imageDatasetId: reg.dataset.id } : { image: url },
      );

      if (!saved?.image) return;
      setCurrent({ image: saved.image, datasetId: saved.imageDatasetId });
      setCurrentTitle(mine ? reg.dataset.title : null);
      setUrlDraft("");
      void fetchPage();
      toast.success(
        mine
          ? "Image added to your account and attached. Reload the viewer to see it."
          : "Image overlay attached. Reload the viewer to see it.",
      );
    } finally {
      setBusy(false);
    }
  };

  const remove = async () => {
    if (!confirm("Remove the image overlay from this dataset?")) return;
    setBusy(true);
    try {
      if (!(await save({ image: null }))) return;
      setCurrent(null);
      setCurrentTitle(null);
      toast.success("Image overlay removed.");
    } finally {
      setBusy(false);
    }
  };

  const totalPages = Math.max(1, Math.ceil(total / PAGE_SIZE));

  return (
    <div className="flex w-full flex-col gap-3">
      <div>
        <h3 className="text-sm font-semibold text-foreground">Image overlay</h3>
        <p className="mt-0.5 text-xs text-default-500">
          Draw one of your images, such as a stitched DAPI mosaic, under these
          cells. Its coordinates must match the cells (the dataset must keep its
          original, un-normalized coordinates).
        </p>
      </div>

      {current ? (
        <div className="flex items-center justify-between gap-3 rounded-lg bg-content2 px-3 py-2">
          <div className="min-w-0">
            <span className="text-[11px] uppercase tracking-wide text-default-400">
              Current image
            </span>
            <p className="truncate text-sm font-medium text-foreground">
              {currentTitle ?? current.image}
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
        aria-label="Search your images"
        placeholder="Search your images…"
        size="sm"
        value={search}
        onClear={() => setSearch("")}
        onValueChange={setSearch}
      />

      {loading ? (
        <div className="flex h-24 items-center justify-center">
          <Spinner size="sm" />
        </div>
      ) : items.length === 0 ? (
        <p className="py-4 text-center text-xs text-default-400">
          {search
            ? "No images match."
            : "You have no images yet. Add one by URL below."}
        </p>
      ) : (
        <div className="flex flex-col gap-1.5" data-testid="image-overlay-list">
          {items.map((it) => {
            const isCurrent = it.id === current?.datasetId;

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
                onClick={() => !isCurrent && attach(it)}
              >
                <p className="min-w-0 truncate text-sm font-medium text-foreground">
                  {it.title ?? it.id}
                </p>
                <span
                  className={`shrink-0 text-[11px] ${isCurrent ? "font-medium text-primary" : "text-default-400"}`}
                >
                  {isCurrent ? "Overlaid" : "Attach →"}
                </span>
              </button>
            );
          })}
        </div>
      )}

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
            Page {page} of {totalPages} · {total} images
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

      <div className="flex items-center gap-2">
        <Input
          aria-label="Image URL"
          placeholder="Add by URL: https://bucket.s3.region.amazonaws.com/folder/image.zarr"
          size="sm"
          value={urlDraft}
          onValueChange={setUrlDraft}
        />
        <Button
          className="shrink-0"
          color="primary"
          isDisabled={busy || !urlDraft.trim()}
          size="sm"
          onPress={addByUrl}
        >
          Add
        </Button>
      </div>
    </div>
  );
}
