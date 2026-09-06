"use client";

import LabelledMoleculeViewerPage from "@/components/labelled-molecule-viewer";

/**
 * Depth-correct opaque rasterisation — the default. /lm2-viewer is kept as the
 * blended counterpart for comparison.
 */
export default function Page() {
  return <LabelledMoleculeViewerPage renderMode="opaque" />;
}
