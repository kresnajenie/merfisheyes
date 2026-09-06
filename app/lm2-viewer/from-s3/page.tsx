"use client";

import LabelledMoleculeViewerPage from "@/components/labelled-molecule-viewer";

/**
 * Same viewer, opaque rasterisation: points write depth and discard below
 * alphaTest, so the GPU rejects buried points early. Side-by-side counterpart
 * to /lm-viewer for judging occlusion and frame rate.
 */
export default function Page() {
  return <LabelledMoleculeViewerPage renderMode="opaque" />;
}
