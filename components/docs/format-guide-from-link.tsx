"use client";

import { useSearchParams } from "next/navigation";

import { FormatGuide } from "@/components/format-guide";

/** The format guide, opened on the tab a link names (`?format=image`). */
export function FormatGuideFromLink() {
  return <FormatGuide initialKey={useSearchParams().get("format")} />;
}
