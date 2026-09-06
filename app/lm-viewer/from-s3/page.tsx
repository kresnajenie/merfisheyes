"use client";

import LabelledMoleculeViewerPage from "@/components/labelled-molecule-viewer";

/** Alpha-blended points: partial opacity, but no depth occlusion. */
export default function Page() {
  return <LabelledMoleculeViewerPage renderMode="blended" />;
}
