"use client";

import { use } from "react";

import { ImageViewer } from "@/components/image-viewer";

/** `/image-viewer/{id}` — an image dataset uploaded through the app. */
export default function ImageViewerByIdPage({
  params,
}: {
  params: Promise<{ id: string }>;
}) {
  const { id } = use(params);

  return <ImageViewer datasetId={id} />;
}
