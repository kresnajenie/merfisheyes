"use client";

import { useState, useRef, useEffect, useCallback } from "react";
import { createPortal } from "react-dom";
import { Card, CardBody } from "@heroui/card";
import { Chip } from "@heroui/chip";
import { useRouter } from "next/navigation";
import { motion, AnimatePresence } from "framer-motion";

import type { CatalogDatasetItem, CatalogDatasetEntry } from "./types";

interface ExploreDatasetCardProps {
  dataset: CatalogDatasetItem;
  /** Use portal-based dropdown instead of inline expansion (for horizontal scroll rows) */
  usePopover?: boolean;
  /** When provided, calls this instead of navigating via router */
  onSelect?: (dataset: CatalogDatasetItem, entry: CatalogDatasetEntry) => void;
}

function navigateToEntry(entry: CatalogDatasetEntry, router: ReturnType<typeof useRouter>) {
  const base = entry.datasetType === "single_molecule" ? "/sm-viewer" : "/viewer";
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
    <div className="bg-content1 dark:bg-default-100 rounded-xl shadow-lg border border-default-200 p-2 flex flex-col gap-1">
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
              color={entry.datasetType === "single_cell" ? "primary" : "secondary"}
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

export function ExploreDatasetCard({ dataset, usePopover, onSelect }: ExploreDatasetCardProps) {
  const router = useRouter();
  const [isExpanded, setIsExpanded] = useState(false);
  const cardRef = useRef<HTMLDivElement>(null);
  const dropdownRef = useRef<HTMLDivElement>(null);
  const [dropdownStyle, setDropdownStyle] = useState<React.CSSProperties>({});

  const entries = dataset.entries ?? [];
  const hasMultipleEntries = entries.length > 1;
  const uniqueTypes = [...new Set(entries.map((e) => e.datasetType))];

  // Position the portal dropdown below the card
  useEffect(() => {
    if (!isExpanded || !usePopover || !cardRef.current) return;
    const rect = cardRef.current.getBoundingClientRect();
    setDropdownStyle({
      position: "fixed",
      top: rect.bottom + 4,
      left: rect.left,
      width: rect.width,
      zIndex: 9999,
    });
  }, [isExpanded, usePopover]);

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
    const closeOnScroll = () => setIsExpanded(false);
    window.addEventListener("scroll", closeOnScroll, true);
    return () => {
      document.removeEventListener("mousedown", handleClickOutside);
      window.removeEventListener("scroll", closeOnScroll, true);
    };
  }, [isExpanded, handleClickOutside]);

  const handlePress = () => {
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

  // Portal-based dropdown for usePopover mode
  const portalDropdown =
    usePopover && isExpanded && hasMultipleEntries && typeof document !== "undefined"
      ? createPortal(
          <div ref={dropdownRef} style={dropdownStyle}>
            <motion.div
              animate={{ opacity: 1, y: 0 }}
              exit={{ opacity: 0, y: -4 }}
              initial={{ opacity: 0, y: -4 }}
              transition={{ duration: 0.15 }}
            >
              <EntryList dataset={dataset} entries={entries} router={router} onSelect={onSelect} />
            </motion.div>
          </div>,
          document.body,
        )
      : null;

  return (
    <div ref={cardRef} className="relative">
      <Card
        isBlurred
        isPressable={!!hasSomeEntry}
        className="border-none bg-background/60 dark:bg-default-100/50 hover:bg-default-200/50 hover:scale-[1.02] transition-all duration-200 cursor-pointer"
        shadow="sm"
        onPress={handlePress}
      >
        <CardBody className="p-0 overflow-hidden">
          {/* Thumbnail — only shown when image exists */}
          {dataset.thumbnailUrl && (
            <div className="relative h-40 w-full bg-default-200/30">
              <img
                alt={dataset.title}
                className="w-full h-full object-cover"
                src={dataset.thumbnailUrl}
              />
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
              {hasMultipleEntries && (
                <div className="absolute top-2 right-2">
                  <Chip size="sm" variant="solid" className="bg-default-800/70 text-white">
                    {entries.length} entries
                  </Chip>
                </div>
              )}
              {/* External link icon */}
              {dataset.externalLink && (
                <div
                  className="absolute bottom-2 right-2 w-7 h-7 rounded-full bg-default-800/60 flex items-center justify-center hover:bg-default-800/80 transition-colors cursor-pointer"
                  role="link"
                  title="Open source page"
                  onClick={(e) => {
                    e.stopPropagation();
                    window.open(dataset.externalLink!, "_blank", "noopener,noreferrer");
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

          {/* Content */}
          <div className="p-4 flex flex-col gap-2">
            {/* Inline badges when no thumbnail */}
            {!dataset.thumbnailUrl && (
              <div className="flex items-center gap-1 flex-wrap">
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
                {hasMultipleEntries && (
                  <Chip size="sm" variant="flat" className="text-default-500">
                    {entries.length} entries
                  </Chip>
                )}
                {dataset.externalLink && (
                  <div
                    className="ml-auto w-6 h-6 rounded-full bg-default-200/60 flex items-center justify-center hover:bg-default-300/80 transition-colors cursor-pointer"
                    role="link"
                    title="Open source page"
                    onClick={(e) => {
                      e.stopPropagation();
                      window.open(dataset.externalLink!, "_blank", "noopener,noreferrer");
                    }}
                  >
                    <svg
                      className="w-3 h-3 text-default-600"
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

            <h3 className="font-semibold text-foreground text-sm line-clamp-2">
              {dataset.title}
            </h3>

            {/* Description — show more lines when no thumbnail */}
            {dataset.description && (
              <p className={`text-xs text-default-500 ${dataset.thumbnailUrl ? "line-clamp-2" : "line-clamp-3"}`}>
                {dataset.description}
              </p>
            )}

            {/* Metadata row */}
            <div className="flex flex-wrap gap-1">
              {dataset.species && (
                <Chip size="sm" variant="flat">{dataset.species}</Chip>
              )}
              {dataset.tissue && (
                <Chip size="sm" variant="flat">{dataset.tissue}</Chip>
              )}
            </div>

            {/* Platform & stats */}
            <div className="flex flex-col gap-0.5 text-xs text-default-500">
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
              <div className="flex flex-wrap gap-1 mt-1">
                {dataset.tags.slice(0, 3).map((tag) => (
                  <Chip key={tag} className="text-[10px]" size="sm" variant="dot">
                    {tag}
                  </Chip>
                ))}
              </div>
            )}
          </div>
        </CardBody>
      </Card>

      {/* Expanded entry list — inline (grid cards) */}
      {!usePopover && (
        <AnimatePresence>
          {isExpanded && hasMultipleEntries && (
            <motion.div
              animate={{ opacity: 1, y: 0, height: "auto" }}
              className="mt-2"
              exit={{ opacity: 0, y: -4, height: 0 }}
              initial={{ opacity: 0, y: -4, height: 0 }}
              transition={{ duration: 0.2 }}
            >
              <EntryList dataset={dataset} entries={entries} router={router} onSelect={onSelect} />
            </motion.div>
          )}
        </AnimatePresence>
      )}

      {/* Portal-based dropdown (horizontal scroll rows) */}
      {portalDropdown}
    </div>
  );
}
