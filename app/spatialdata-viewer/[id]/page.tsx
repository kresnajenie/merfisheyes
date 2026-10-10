"use client";

import { useParams } from "next/navigation";

import { SpatialDataViewer } from "@/components/spatialdata-viewer";

export default function SpatialDataViewerPage() {
  return <SpatialDataViewer datasetId={useParams().id as string} />;
}
