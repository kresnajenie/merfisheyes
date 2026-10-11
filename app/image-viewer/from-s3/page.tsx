"use client";

import { Suspense } from "react";
import { useSearchParams } from "next/navigation";

import { ImageViewer } from "@/components/image-viewer";

/** `/image-viewer/from-s3?url=<public URL of an OME-Zarr image>` */
function FromS3Content() {
  const url = useSearchParams().get("url");

  if (!url) {
    return (
      <div className="flex items-center justify-center h-full">
        Add ?url= with the public URL of an OME-Zarr image.
      </div>
    );
  }

  return <ImageViewer url={url} />;
}

export default function ImageViewerFromS3Page() {
  return (
    <Suspense fallback={null}>
      <FromS3Content />
    </Suspense>
  );
}
