"use client";

import LabelledMoleculeViewerPage from "@/components/labelled-molecule-viewer";
import { SplitScreenContainer } from "@/components/split-screen-container";

export default function Page() {
  // The container always renders the same tree, so the left panel's WebGL
  // context survives opening and closing the split.
  return (
    <SplitScreenContainer>
      <LabelledMoleculeViewerPage />
    </SplitScreenContainer>
  );
}
