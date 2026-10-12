"use client";

import { Suspense } from "react";
import { useSearchParams } from "next/navigation";

import { SpatialDataViewer } from "@/components/spatialdata-viewer";

/** `/spatialdata-viewer/from-s3?url=<public URL of a .zarr store>` */
function FromS3Content() {
  const url = useSearchParams().get("url");

  if (!url) {
    return (
      <div className="flex items-center justify-center h-full">
        Add ?url= with the public URL of a SpatialData .zarr store.
      </div>
    );
  }

  return <SpatialDataViewer url={url} />;
}

export default function SpatialDataViewerFromS3Page() {
  return (
    <Suspense fallback={null}>
      <FromS3Content />
    </Suspense>
  );
}
