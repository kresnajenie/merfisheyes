"use client";

import type { PanelType } from "@/lib/stores/splitScreenStore";

import { useState } from "react";
import { Button, Input } from "@heroui/react";
import { usePathname, useSearchParams } from "next/navigation";

import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";

export function SplitPanelPicker() {
  const { setRightPanel, setRightPanelS3 } = useSplitScreenStore();
  const pathname = usePathname();
  const searchParams = useSearchParams();
  const [linkInput, setLinkInput] = useState("");
  const [showLinkInput, setShowLinkInput] = useState(false);

  const isFromS3Page = pathname.includes("/from-s3");

  // Extract current dataset info from URL
  const getCurrentDatasetInfo = (): {
    id: string | null;
    s3Url: string | null;
    type: PanelType;
  } => {
    const isSm = pathname.startsWith("/sm-viewer");
    const type: PanelType = isSm ? "sm" : "cell";

    // For from-s3 pages, return the S3 URL
    if (isFromS3Page) {
      const s3Url = searchParams.get("url");

      return {
        id: null,
        s3Url: s3Url ? decodeURIComponent(s3Url) : null,
        type,
      };
    }

    // Extract ID from /viewer/[id] or /sm-viewer/[id]
    const parts = pathname.split("/");
    const id = parts.length >= 3 ? parts[parts.length - 1] : null;

    // Don't return ID for non-ID pages like /viewer or /sm-viewer
    if (id === "viewer" || id === "sm-viewer") {
      return { id: null, s3Url: null, type };
    }

    return { id, s3Url: null, type };
  };

  const handleSameDataset = () => {
    const { id, s3Url, type } = getCurrentDatasetInfo();

    if (s3Url) {
      setRightPanelS3(s3Url, type);
    } else if (id) {
      setRightPanel(id, type);
    }
  };

  // Detect if a URL is a raw S3 URL (e.g., https://bucket.s3.region.amazonaws.com/...)
  const isRawS3Url = (url: string): boolean => {
    try {
      const parsed = new URL(url);

      return (
        parsed.hostname.includes("s3") &&
        parsed.hostname.includes("amazonaws.com")
      );
    } catch {
      return false;
    }
  };

  const handlePasteLink = () => {
    if (!linkInput.trim()) return;

    const input = linkInput.trim();

    // Check if it's a raw S3 URL
    if (isRawS3Url(input)) {
      const isSm = pathname.startsWith("/sm-viewer");

      setRightPanelS3(input, isSm ? "sm" : "cell");

      return;
    }

    try {
      const url = new URL(input);
      const pathParts = url.pathname.split("/").filter(Boolean);

      let type: PanelType = "cell";
      let id: string | null = null;

      if (pathParts[0] === "sm-viewer") {
        type = "sm";
        if (pathParts[1] === "from-s3") {
          // Handle from-s3 URLs: /sm-viewer/from-s3?url=...
          const s3UrlParam = url.searchParams.get("url");

          if (s3UrlParam) {
            setRightPanelS3(decodeURIComponent(s3UrlParam), type);

            return;
          }
        } else if (pathParts.length >= 2) {
          id = pathParts[1];
        }
      } else if (pathParts[0] === "viewer") {
        type = "cell";
        if (pathParts[1] === "from-s3") {
          // Handle from-s3 URLs: /viewer/from-s3?url=...
          const s3UrlParam = url.searchParams.get("url");

          if (s3UrlParam) {
            setRightPanelS3(decodeURIComponent(s3UrlParam), type);

            return;
          }
        } else if (pathParts.length >= 2) {
          id = pathParts[1];
        }
      }

      if (id) {
        setRightPanel(id, type);
      }
    } catch {
      // Invalid URL — try treating as dataset ID
      if (input.length > 0) {
        const isSm = pathname.startsWith("/sm-viewer");

        setRightPanel(input, isSm ? "sm" : "cell");
      }
    }
  };

  const { id: currentId, s3Url: currentS3Url } = getCurrentDatasetInfo();
  const hasCurrentDataset = !!(currentId || currentS3Url);

  return (
    <div className="absolute inset-0 flex items-center justify-center bg-black">
      <div className="flex flex-col items-center gap-6 max-w-sm w-full px-8">
        <div className="text-center mb-2">
          <h3 className="text-lg font-semibold text-white">Choose a Dataset</h3>
          <p className="text-sm text-white/50 mt-1">
            Select what to display in this panel
          </p>
        </div>

        <div className="flex flex-col gap-3 w-full">
          {/* Same Dataset */}
          {hasCurrentDataset && (
            <Button
              className="w-full justify-start"
              color="primary"
              variant="bordered"
              onPress={handleSameDataset}
            >
              <svg
                className="w-5 h-5 mr-2 flex-shrink-0"
                fill="none"
                stroke="currentColor"
                strokeWidth={1.5}
                viewBox="0 0 24 24"
              >
                <path
                  d="M15.75 17.25v3.375c0 .621-.504 1.125-1.125 1.125h-9.75a1.125 1.125 0 01-1.125-1.125V7.875c0-.621.504-1.125 1.125-1.125H6.75a9.06 9.06 0 011.5.124m7.5 10.376h3.375c.621 0 1.125-.504 1.125-1.125V11.25c0-4.46-3.243-8.161-7.5-8.876a9.06 9.06 0 00-1.5-.124H9.375c-.621 0-1.125.504-1.125 1.125v3.5m7.5 10.375H9.375a1.125 1.125 0 01-1.125-1.125v-9.25m12 6.625v-1.875a3.375 3.375 0 00-3.375-3.375h-1.5a1.125 1.125 0 01-1.125-1.125v-1.5a3.375 3.375 0 00-3.375-3.375H9.75"
                  strokeLinecap="round"
                  strokeLinejoin="round"
                />
              </svg>
              Same Dataset
            </Button>
          )}

          {/* Paste Link */}
          {showLinkInput ? (
            <div className="flex flex-col gap-2">
              <Input
                classNames={{
                  input: "text-white",
                  inputWrapper: "bg-white/5 border-white/20 hover:bg-white/10",
                }}
                placeholder="Paste dataset URL or ID..."
                size="sm"
                value={linkInput}
                onChange={(e) => setLinkInput(e.target.value)}
                onKeyDown={(e) => {
                  if (e.key === "Enter") handlePasteLink();
                }}
              />
              <div className="flex gap-2">
                <Button
                  className="flex-1"
                  color="primary"
                  isDisabled={!linkInput.trim()}
                  size="sm"
                  onPress={handlePasteLink}
                >
                  Load
                </Button>
                <Button
                  className="flex-1"
                  size="sm"
                  variant="bordered"
                  onPress={() => {
                    setShowLinkInput(false);
                    setLinkInput("");
                  }}
                >
                  Cancel
                </Button>
              </div>
            </div>
          ) : (
            <Button
              className="w-full justify-start"
              variant="bordered"
              onPress={() => setShowLinkInput(true)}
            >
              <svg
                className="w-5 h-5 mr-2 flex-shrink-0"
                fill="none"
                stroke="currentColor"
                strokeWidth={1.5}
                viewBox="0 0 24 24"
              >
                <path
                  d="M13.19 8.688a4.5 4.5 0 011.242 7.244l-4.5 4.5a4.5 4.5 0 01-6.364-6.364l1.757-1.757m9.86-3.06a4.5 4.5 0 00-1.242-7.244l-4.5-4.5a4.5 4.5 0 00-6.364 6.364l1.757 1.757"
                  strokeLinecap="round"
                  strokeLinejoin="round"
                />
              </svg>
              Paste Link
            </Button>
          )}
        </div>
      </div>
    </div>
  );
}
