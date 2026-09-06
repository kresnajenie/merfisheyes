"use client";

import LabelledMoleculeViewerPage from "@/components/labelled-molecule-viewer";

/**
 * Alpha-blended counterpart to /lm-viewer: no depth writes, so nothing occludes
 * anything and every fragment shades. Kept only for comparison.
 */
export default function Page() {
  return <LabelledMoleculeViewerPage renderMode="blended" />;
}
