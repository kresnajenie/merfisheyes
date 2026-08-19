"use client";

import type { CatalogDatasetItem, CatalogDatasetEntry } from "./types";

import { useState, useRef, useEffect, useCallback } from "react";
import { createPortal } from "react-dom";
import { Card, CardBody } from "@heroui/card";
import { Chip } from "@heroui/chip";
import NextLink from "next/link";
import { useRouter } from "next/navigation";
import { motion, AnimatePresence } from "framer-motion";

/**
 * Compute a static destination URL for the card if exactly one exists.
 * Returns null when the card needs runtime branching (e.g. multi-entry
 * datasets that expand a dropdown). Used to render the card as an actual
 * <a>/<Link> so right-click and middle-click open in new tab.
 */
function getCardHref(
  dataset: CatalogDatasetItem,
  entries: CatalogDatasetEntry[],
): { href: string; external?: boolean } | null {
  if (dataset.bilCode) {
    return { href: `/explore/bil/${dataset.bilCode}` };
  }
  // Community projects open a detail page, like BIL — always when the row is
  // project-backed (sourceProjectId), even with a single dataset, and as a
  // fallback for multi-entry rows from payloads without that field.
  if (dataset.isCommunity && (dataset.sourceProjectId || entries.length > 1)) {
    return { href: `/explore/${dataset.id}` };
  }
  if (entries.length === 0) {
    return dataset.externalLink
      ? { href: dataset.externalLink, external: true }
      : null;
  }
  if (entries.length === 1) {
    const e = entries[0];
    const base = e.datasetType === "single_molecule" ? "/sm-viewer" : "/viewer";

    if (e.s3BaseUrl) {
      const url = e.s3BaseUrl.replace(/\/+$/, "");

      return { href: `${base}/from-s3?url=${encodeURIComponent(url)}` };
    }
    if (e.datasetId) {
      return { href: `${base}/${e.datasetId}` };
    }
  }

  return null;
}

interface ExploreDatasetCardProps {
  dataset: CatalogDatasetItem;
  /** Use portal-based dropdown instead of inline expansion (for horizontal scroll rows) */
  usePopover?: boolean;
  /** When provided, calls this instead of navigating via router */
  onSelect?: (dataset: CatalogDatasetItem, entry: CatalogDatasetEntry) => void;
  /**
   * When provided, intercepts the card click — skips link navigation,
   * inline dropdown expansion, and `onSelect` shortcut. The host (e.g. the
   * explore modal) decides what to do next (load entry directly, route to
   * an in-host detail view, etc.).
   */
  onCardClick?: (dataset: CatalogDatasetItem) => void;
  /** Highlight matching genes from search */
  geneHighlight?: string;
  /**
   * Owner controls rendered in the top-right corner (e.g. a ⋯ menu). When
   * present, the built-in entry-count badge is suppressed so the two don't
   * overlap. Used by the account cards.
   */
  ownerOverlay?: React.ReactNode;
  /** Badges overlaid on the header's bottom-left (e.g. project membership). */
  headerBadges?: React.ReactNode;
}

function navigateToEntry(
  entry: CatalogDatasetEntry,
  router: ReturnType<typeof useRouter>,
) {
  const base =
    entry.datasetType === "single_molecule" ? "/sm-viewer" : "/viewer";

  if (entry.s3BaseUrl) {
    const url = entry.s3BaseUrl.replace(/\/+$/, "");

    router.push(`${base}/from-s3?url=${encodeURIComponent(url)}`);
  } else if (entry.datasetId) {
    router.push(`${base}/${entry.datasetId}`);
  }
}

function EntryList({
  dataset,
  entries,
  router,
  onSelect,
}: {
  dataset: CatalogDatasetItem;
  entries: CatalogDatasetEntry[];
  router: ReturnType<typeof useRouter>;
  onSelect?: (dataset: CatalogDatasetItem, entry: CatalogDatasetEntry) => void;
}) {
  return (
    <div className="bg-content1 dark:bg-default-100 rounded-xl shadow-lg border border-default-200 p-2 flex flex-col gap-1 max-h-[220px] overflow-y-auto">
      {entries.map((entry) => {
        const hasLink = entry.s3BaseUrl || entry.datasetId;

        return (
          <button
            key={entry.id}
            className={`flex items-center gap-2 px-3 py-2 rounded-lg text-left text-sm transition-colors ${
              hasLink
                ? "hover:bg-default-100 dark:hover:bg-default-200/50 cursor-pointer"
                : "opacity-50 cursor-default"
            }`}
            disabled={!hasLink}
            onClick={() => {
              if (!hasLink) return;
              if (onSelect) {
                onSelect(dataset, entry);
              } else {
                navigateToEntry(entry, router);
              }
            }}
          >
            <Chip
              color={
                entry.datasetType === "single_cell" ? "primary" : "secondary"
              }
              size="sm"
              variant="flat"
            >
              {entry.datasetType === "single_cell" ? "SC" : "SM"}
            </Chip>
            <span className="truncate">{entry.label}</span>
          </button>
        );
      })}
    </div>
  );
}

export function ExploreDatasetCard({
  dataset,
  usePopover,
  onSelect,
  onCardClick,
  geneHighlight,
  ownerOverlay,
  headerBadges,
}: ExploreDatasetCardProps) {
  const router = useRouter();
  const [isExpanded, setIsExpanded] = useState(false);
  const cardRef = useRef<HTMLDivElement>(null);
  const dropdownRef = useRef<HTMLDivElement>(null);
  const [dropdownStyle, setDropdownStyle] = useState<React.CSSProperties>({});
  const [expandUp, setExpandUp] = useState(false);

  const entries = dataset.entries ?? [];
  const hasMultipleEntries = entries.length > 1;
  const uniqueTypes = [...new Set(entries.map((e) => e.datasetType))];
  // The card is a fixed size everywhere (Explore + account): a header area is
  // always shown (a placeholder when there's no thumbnail) so grids stay
  // uniform regardless of how much metadata each item has.

  // Recompute dropdown position from current card rect
  const updatePosition = useCallback(() => {
    if (!cardRef.current) return;
    const rect = cardRef.current.getBoundingClientRect();
    const spaceBelow = window.innerHeight - rect.bottom;
    const shouldExpandUp = spaceBelow < 250;

    setExpandUp(shouldExpandUp);

    if (usePopover) {
      if (shouldExpandUp) {
        setDropdownStyle({
          position: "fixed",
          bottom: window.innerHeight - rect.top + 4,
          left: rect.left,
          width: rect.width,
          zIndex: 9999,
        });
      } else {
        setDropdownStyle({
          position: "fixed",
          top: rect.bottom + 4,
          left: rect.left,
          width: rect.width,
          zIndex: 9999,
        });
      }
    }
  }, [usePopover]);

  // Position on open + reposition on scroll
  useEffect(() => {
    if (!isExpanded) return;
    updatePosition();
    const onScroll = () => updatePosition();

    window.addEventListener("scroll", onScroll, true);

    return () => window.removeEventListener("scroll", onScroll, true);
  }, [isExpanded, updatePosition]);

  // Close dropdown when clicking outside
  const handleClickOutside = useCallback(
    (e: MouseEvent) => {
      if (!isExpanded) return;
      const target = e.target as Node;

      if (cardRef.current?.contains(target)) return;
      if (dropdownRef.current?.contains(target)) return;
      setIsExpanded(false);
    },
    [isExpanded],
  );

  useEffect(() => {
    if (!isExpanded) return;
    document.addEventListener("mousedown", handleClickOutside);

    return () => {
      document.removeEventListener("mousedown", handleClickOutside);
    };
  }, [isExpanded, handleClickOutside]);

  const handlePress = () => {
    // Host (e.g. explore modal) takes over — no navigation, no inline dropdown.
    if (onCardClick) {
      onCardClick(dataset);

      return;
    }

    // Navigate to detail page
    if (dataset.bilCode) {
      router.push(`/explore/bil/${dataset.bilCode}`);

      return;
    }

    // Project-backed community rows open their detail page, even with a
    // single entry (mirrors getCardHref). onSelect hosts keep their own flow.
    if (dataset.isCommunity && dataset.sourceProjectId && !onSelect) {
      router.push(`/explore/${dataset.id}`);

      return;
    }

    // Non-BIL datasets: keep existing behavior
    // Zero entries: open external link or do nothing
    if (entries.length === 0) {
      if (dataset.externalLink) {
        window.open(dataset.externalLink, "_blank", "noopener,noreferrer");
      }

      return;
    }

    // Single entry: navigate directly (or call onSelect)
    if (entries.length === 1) {
      if (onSelect) {
        onSelect(dataset, entries[0]);
      } else {
        navigateToEntry(entries[0], router);
      }

      return;
    }

    // Multiple entries: toggle expansion
    setIsExpanded((prev) => !prev);
  };

  const formatCount = (n: number | null) => {
    if (n == null) return null;
    if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(1)}M`;
    if (n >= 1_000) return `${(n / 1_000).toFixed(0)}K`;

    return String(n);
  };

  const hasSomeEntry = entries.length > 0 || dataset.externalLink;

  // If the card has a single static destination AND no onSelect handler,
  // render the whole thing as a real link (right-click → open in new tab,
  // middle-click, cmd-click, etc. all work). Multi-entry cards or onSelect
  // overrides keep the click-driven button behavior so the dropdown / handler
  // can run.
  // onCardClick fully overrides default behavior — never linkify in that mode.
  const cardHref =
    onSelect || onCardClick ? null : getCardHref(dataset, entries);
  const linkified = cardHref !== null;

  // Portal-based dropdown for usePopover mode
  const portalDropdown =
    usePopover &&
    isExpanded &&
    hasMultipleEntries &&
    typeof document !== "undefined"
      ? createPortal(
          <div ref={dropdownRef} style={dropdownStyle}>
            <motion.div
              animate={{ opacity: 1, y: 0 }}
              exit={{ opacity: 0, y: expandUp ? 4 : -4 }}
              initial={{ opacity: 0, y: expandUp ? 4 : -4 }}
              transition={{ duration: 0.15 }}
            >
              <EntryList
                dataset={dataset}
                entries={entries}
                router={router}
                onSelect={onSelect}
              />
            </motion.div>
          </div>,
          document.body,
        )
      : null;

  const cardElement = (
    <Card
      isBlurred
      className="border-none bg-background/60 dark:bg-default-100/50 hover:bg-default-200/50 hover:scale-[1.02] transition-all duration-200 cursor-pointer h-[400px] w-full"
      isPressable={!linkified && (!!hasSomeEntry || !!onCardClick)}
      shadow="sm"
      onPress={linkified ? undefined : handlePress}
    >
      <CardBody className="p-0 overflow-hidden flex h-full flex-col">
        {/* Header — thumbnail, or a gradient placeholder when there's none.
            The gradient always renders underneath; a broken thumbnail hides
            itself (otherwise the browser paints its alt text across the
            header). FIXED tall height (~half the card) so every card's text
            starts at the same y — rows stay aligned regardless of how much
            metadata each card has. Sized so a full card (3-line title, byline,
            2-line description, chips, stats) keeps its stats visible. */}
        <div className="relative h-48 w-full shrink-0 bg-default-200/30">
          <div className="h-full w-full bg-gradient-to-br from-default-100 to-default-200/40" />
          {dataset.thumbnailUrl && (
            <img
              alt={dataset.title}
              className="absolute inset-0 w-full h-full object-cover"
              src={dataset.thumbnailUrl}
              onError={(e) => {
                e.currentTarget.style.display = "none";
              }}
            />
          )}
          {/* Type badge overlay */}
          <div className="absolute top-2 left-2 flex gap-1">
            {uniqueTypes.map((type) => (
              <Chip
                key={type}
                color={type === "single_cell" ? "primary" : "secondary"}
                size="sm"
                variant="solid"
              >
                {type === "single_cell" ? "SC" : "SM"}
              </Chip>
            ))}
          </div>
          {/* Entry count badge */}
          {hasMultipleEntries && !ownerOverlay && (
            <div className="absolute top-2 right-2">
              <Chip
                className="bg-default-800/70 text-white"
                size="sm"
                variant="solid"
              >
                {entries.length} entries
              </Chip>
            </div>
          )}
          {/* Link icons */}
          {(dataset.externalLink || dataset.publicationLink) && (
            <div className="absolute bottom-2 right-2 flex gap-1.5">
              {dataset.publicationLink && (
                <div
                  className="w-7 h-7 rounded-full bg-default-800/60 flex items-center justify-center hover:bg-default-800/80 transition-colors cursor-pointer"
                  role="link"
                  title="Open publication"
                  onClick={(e) => {
                    e.preventDefault();
                    e.stopPropagation();
                    window.open(
                      dataset.publicationLink!,
                      "_blank",
                      "noopener,noreferrer",
                    );
                  }}
                >
                  <svg
                    className="w-3.5 h-3.5 text-white"
                    fill="none"
                    stroke="currentColor"
                    strokeWidth={2}
                    viewBox="0 0 24 24"
                  >
                    <path
                      d="M12 6.253v13m0-13C10.832 5.477 9.246 5 7.5 5S4.168 5.477 3 6.253v13C4.168 18.477 5.754 18 7.5 18s3.332.477 4.5 1.253m0-13C13.168 5.477 14.754 5 16.5 5c1.747 0 3.332.477 4.5 1.253v13C19.832 18.477 18.247 18 16.5 18c-1.746 0-3.332.477-4.5 1.253"
                      strokeLinecap="round"
                      strokeLinejoin="round"
                    />
                  </svg>
                </div>
              )}
              {dataset.externalLink && (
                <div
                  className="w-7 h-7 rounded-full bg-default-800/60 flex items-center justify-center hover:bg-default-800/80 transition-colors cursor-pointer"
                  role="link"
                  title="Open source page"
                  onClick={(e) => {
                    e.preventDefault();
                    e.stopPropagation();
                    window.open(
                      dataset.externalLink!,
                      "_blank",
                      "noopener,noreferrer",
                    );
                  }}
                >
                  <svg
                    className="w-3.5 h-3.5 text-white"
                    fill="none"
                    stroke="currentColor"
                    strokeWidth={2}
                    viewBox="0 0 24 24"
                  >
                    <path
                      d="M18 13v6a2 2 0 01-2 2H5a2 2 0 01-2-2V8a2 2 0 012-2h6M15 3h6v6M10 14L21 3"
                      strokeLinecap="round"
                      strokeLinejoin="round"
                    />
                  </svg>
                </div>
              )}
            </div>
          )}
          {/* Header badges (e.g. project membership) — bottom-left */}
          {headerBadges && (
            <div className="absolute bottom-2 left-2 z-10 flex max-w-[70%] flex-wrap gap-1">
              {headerBadges}
            </div>
          )}
        </div>

        {/* Content. Every row that shows text is shrink-0: the column has
            overflow-hidden, and without it flexbox compresses the rows to fit,
            slicing lines mid-glyph. With shrink-0 each row keeps its full
            height and any overflow crops whole rows at the bottom instead. */}
        <div className="p-4 flex flex-col gap-2 min-h-0 flex-1 overflow-hidden">
          <h3
            className="font-semibold text-foreground text-sm line-clamp-3 min-w-0 [overflow-wrap:anywhere] shrink-0"
            title={dataset.title}
          >
            {dataset.title}
          </h3>

          {/* Byline: investigator · age · publication year — the metadata a
              first-time browser needs to judge a dataset at a glance. */}
          {(() => {
            const meta = (dataset.metadata ?? {}) as Record<string, unknown>;
            const parts = [meta.investigator, meta.age, meta.publicationYear]
              .filter(Boolean)
              .map(String);

            if (parts.length === 0) return null;
            const byline = parts.join(" · ");

            return (
              <p
                className="text-xs text-default-400 truncate shrink-0"
                title={byline}
              >
                {byline}
              </p>
            );
          })()}

          {/* Description */}
          {dataset.description && (
            <p
              className="text-xs text-default-500 line-clamp-2 min-w-0 [overflow-wrap:anywhere] shrink-0"
              title={dataset.description}
            >
              {dataset.description}
            </p>
          )}

          {/* Metadata row */}
          {(dataset.species || dataset.tissue) && (
            <div className="flex flex-wrap gap-1 shrink-0">
              {dataset.species && (
                <Chip
                  className="max-w-full"
                  classNames={{ content: "truncate" }}
                  size="sm"
                  title={dataset.species}
                  variant="flat"
                >
                  {dataset.species}
                </Chip>
              )}
              {dataset.tissue && (
                <Chip
                  className="max-w-full"
                  classNames={{ content: "truncate" }}
                  size="sm"
                  title={dataset.tissue}
                  variant="flat"
                >
                  {dataset.tissue}
                </Chip>
              )}
            </div>
          )}

          {/* Matching genes — the API computes these per card when a gene
              search is active (list payloads no longer carry the full genes
              array). Falls back to client-side filtering when `genes` is
              present (e.g. the detail page's full payload). Rendered above the
              bottom-pinned stats so it survives cropping before tags do. */}
          {geneHighlight &&
            (dataset.matchedGenes || dataset.genes) &&
            (() => {
              const matches =
                dataset.matchedGenes ??
                (() => {
                  const terms = geneHighlight
                    .toLowerCase()
                    .split(/[,\s]+/)
                    .filter(Boolean);

                  return (dataset.genes ?? []).filter((g) =>
                    terms.some((t) => g.toLowerCase().includes(t)),
                  );
                })();

              if (matches.length === 0) return null;

              return (
                <div className="flex flex-wrap gap-1 shrink-0">
                  {matches.slice(0, 8).map((gene) => (
                    <Chip
                      key={gene}
                      className="text-[10px]"
                      color="secondary"
                      size="sm"
                      variant="flat"
                    >
                      {gene}
                    </Chip>
                  ))}
                  {matches.length > 8 && (
                    <span className="text-[10px] text-default-400 self-center">
                      +{matches.length - 8} more
                    </span>
                  )}
                </div>
              );
            })()}

          {/* Platform & stats. Flows right after the metadata above (no
              mt-auto bottom-pinning) so sparse cards — title + stats only —
              don't open a blank gap in the middle; unused space collects at
              the bottom of the fixed-height card instead. */}
          <div className="flex flex-col gap-0.5 text-xs text-default-500 shrink-0">
            {dataset.platform && <span>Platform: {dataset.platform}</span>}
            <div className="flex gap-3">
              {dataset.numCells != null && (
                <span>{formatCount(dataset.numCells)} cells/molecules</span>
              )}
              {dataset.numGenes != null && (
                <span>{formatCount(dataset.numGenes)} genes</span>
              )}
            </div>
          </div>

          {/* Tags */}
          {dataset.tags.length > 0 && (
            <div className="flex flex-wrap gap-1 mt-1 min-w-0">
              {dataset.tags.slice(0, 3).map((tag) => (
                <Chip
                  key={tag}
                  className="text-[10px] max-w-full"
                  classNames={{ content: "truncate" }}
                  size="sm"
                  variant="dot"
                >
                  {tag}
                </Chip>
              ))}
            </div>
          )}
        </div>
      </CardBody>
    </Card>
  );

  // Wrap in a real anchor when a static href exists. Internal links use
  // NextLink (prefetch + client nav). External (e.g. dataset.externalLink)
  // gets a plain <a target="_blank">.
  const wrappedCard =
    cardHref === null ? (
      cardElement
    ) : cardHref.external ? (
      <a
        className="block"
        href={cardHref.href}
        rel="noopener noreferrer"
        target="_blank"
      >
        {cardElement}
      </a>
    ) : (
      <NextLink className="block" href={cardHref.href}>
        {cardElement}
      </NextLink>
    );

  return (
    <div ref={cardRef} className="relative">
      {wrappedCard}

      {/* Owner controls (top-right) — sibling of the card so clicks here don't
          trigger the card's navigation. */}
      {ownerOverlay && (
        <div className="absolute top-2 right-2 z-20 flex items-center gap-1">
          {ownerOverlay}
        </div>
      )}

      {/* Expanded entry list — inline (grid cards) */}
      {!usePopover && (
        <AnimatePresence>
          {isExpanded && hasMultipleEntries && (
            <motion.div
              animate={{ opacity: 1, y: 0, height: "auto" }}
              className={expandUp ? "mb-2 order-first" : "mt-2"}
              exit={{ opacity: 0, y: expandUp ? 4 : -4, height: 0 }}
              initial={{ opacity: 0, y: expandUp ? 4 : -4, height: 0 }}
              style={
                expandUp
                  ? {
                      position: "absolute",
                      bottom: "100%",
                      left: 0,
                      right: 0,
                      zIndex: 50,
                    }
                  : undefined
              }
              transition={{ duration: 0.2 }}
            >
              <EntryList
                dataset={dataset}
                entries={entries}
                router={router}
                onSelect={onSelect}
              />
            </motion.div>
          )}
        </AnimatePresence>
      )}

      {/* Portal-based dropdown (horizontal scroll rows) */}
      {portalDropdown}
    </div>
  );
}
