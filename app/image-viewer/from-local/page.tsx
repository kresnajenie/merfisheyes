"use client";

import { useEffect } from "react";
import { useRouter } from "next/navigation";

import { ImageScene } from "@/components/image-viewer";
import { ImageViewerControls } from "@/components/image-viewer-controls";
import { SplitScreenContainer } from "@/components/split-screen-container";
import { useLocalBundleStore } from "@/lib/stores/localBundleStore";

/**
 * `/image-viewer/from-local` — an OME-Zarr image folder dropped on the
 * homepage, read from disk. Nothing is uploaded.
 */
export default function ImageViewerFromLocalPage() {
  const router = useRouter();
  const image = useLocalBundleStore((s) => s.image);

  // Nothing survives a reload: the files were only ever in this tab
  useEffect(() => {
    if (!image) router.replace("/");
  }, [image, router]);

  // As the other previews: a split beside it, but no link to share (the
  // files are only in this tab)
  return (
    image && (
      <SplitScreenContainer>
        <ImageViewerControls share={false} />
        <ImageScene image={image} />
      </SplitScreenContainer>
    )
  );
}
