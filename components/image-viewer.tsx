"use client";

import { useEffect, useMemo, useRef, useState } from "react";
import { Spinner } from "@heroui/react";
import * as THREE from "three";

import { subtitle, title } from "@/components/primitives";
import { LayerControls } from "@/components/layer-controls";
import { SpatialScaleBar } from "@/components/spatial-scale-bar";
import { WebGLErrorOverlay } from "@/components/webgl-error-overlay";
import { LayerHost } from "@/lib/scene/LayerHost";
import { useLayers } from "@/lib/scene/useLayers";
import { ClaimDatasetBanner } from "@/components/claim-dataset-banner";
import { ImageViewerControls } from "@/components/image-viewer-controls";
import { SplitScreenPage } from "@/components/split-screen-container";
import {
  imageBounds,
  openOmeZarrImage,
  openOmeZarrImageDataset,
  type SpatialDataImage,
} from "@/lib/spatialdata/image";
import { useViewerRegistrationStore } from "@/lib/stores/viewerRegistrationStore";
import { normalizeS3Url } from "@/lib/utils/viewer-config";
import { initializeScene } from "@/lib/webgl/scene-manager";

/**
 * Viewer for an image on its own, with no cells or molecules: an OME-Zarr
 * multiscale image at a public URL (`url`), or an image dataset in the app's
 * own bucket (`datasetId`).
 *
 * An image is a dataset like any other. One opened by URL is looked up in the
 * database, and when nobody owns it yet the banner offers to add it to the
 * viewer's account.
 *
 * `panel`: the right side of a split view. The page's registration (whose
 * dataset is open, who owns it) belongs to the left panel, so a panel neither
 * sets it nor counts a view, and its messages stay inside the panel.
 */
export function ImageViewer({
  url,
  datasetId,
  panel = false,
}: {
  url?: string;
  datasetId?: string;
  panel?: boolean;
}) {
  const [image, setImage] = useState<SpatialDataImage | null>(null);
  const [error, setError] = useState<string | null>(null);
  // Database id of the image, once known; views are counted against it.
  const [registeredId, setRegisteredId] = useState<string | null>(null);
  const setRegistration = useViewerRegistrationStore((s) => s.set);
  const resetRegistration = useViewerRegistrationStore((s) => s.reset);

  useEffect(() => {
    let cancelled = false;

    setImage(null);
    setError(null);
    (url ? openOmeZarrImage(url) : openOmeZarrImageDataset(datasetId!))
      .then((opened) => !cancelled && setImage(opened))
      .catch((e) => {
        console.error("Error opening image:", e);
        if (!cancelled) {
          setError(
            url
              ? "It must be a public OME-Zarr image (Zarr v3, axes c/y/x) that allows this site in CORS."
              : "This dataset has no OME-Zarr image (Zarr v3, axes c/y/x) the viewer can read.",
          );
        }
      });

    return () => {
      cancelled = true;
    };
  }, [url, datasetId]);

  // Who owns it. A URL nobody has registered still opens; the banner then
  // offers to add it.
  useEffect(() => {
    let cancelled = false;

    setRegisteredId(datasetId ?? null);
    if (!url || panel) return;
    const s3Url = normalizeS3Url(url);

    fetch(`/api/datasets/by-url?url=${encodeURIComponent(s3Url)}`)
      .then((r) => (r.ok ? r.json() : null))
      .catch(() => null)
      .then((reg) => {
        if (cancelled) return;
        setRegisteredId(reg?.registered ? reg.id : null);
        setRegistration({
          dbId: reg?.registered ? reg.id : null,
          ownerId: reg?.ownerId ?? null,
          adminOwned: !!reg?.adminOwned,
          registered: !!reg?.registered,
          viewerConfig: null,
          s3Url,
          kind: "image",
        });
      });

    return () => {
      cancelled = true;
      resetRegistration();
    };
  }, [url, datasetId]);

  // One view per open (the server dedups per session and day).
  useEffect(() => {
    if (!image || !registeredId || panel) return;
    fetch(`/api/datasets/${registeredId}/view`, { method: "POST" }).catch(
      () => {},
    );
  }, [image, registeredId]);

  if (error) {
    return (
      <div className="flex flex-col items-center justify-center h-full gap-4 p-8 text-center">
        <h2 className={title({ size: "md", color: "pink" })}>
          {url
            ? "Couldn't read an image at that URL"
            : "Couldn't open this image"}
        </h2>
        <p className={subtitle()}>{error}</p>
      </div>
    );
  }
  if (!image) {
    return (
      <div
        className={`${panel ? "absolute" : "fixed"} inset-0 flex flex-col items-center justify-center gap-4`}
      >
        <Spinner color="primary" size="lg" />
        <p className={subtitle({ class: "text-center" })}>Opening image...</p>
      </div>
    );
  }

  const scene = (
    <>
      <ImageViewerControls />
      <ImageScene image={image} />
    </>
  );

  // A panel is already inside a split; the page is one that can split
  return panel ? (
    scene
  ) : (
    <SplitScreenPage>
      {url && <ClaimDatasetBanner />}
      {scene}
    </SplitScreenPage>
  );
}

/**
 * An image with nothing over it. It draws the same tiled layer the other
 * viewers put under their points, in the image's own coordinates, so what is
 * seen here lines up with the image as an overlay.
 */
export function ImageScene({ image }: { image: SpatialDataImage }) {
  const containerRef = useRef<HTMLDivElement>(null);
  const cameraRef = useRef<THREE.PerspectiveCamera | null>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const controlsRef = useRef<any>(null);
  const [webglFailed, setWebglFailed] = useState(false);
  // The image as the one layer of the contract, through the registry
  const { onSceneReady, layers } = useLayers(
    useMemo(() => ({ image }), [image]),
  );

  useEffect(() => {
    const container = containerRef.current;

    if (!container) return;

    // Frame the whole image as placed in the scene (pixel i of the finest
    // level spans i ± 0.5; a rotated placement makes the box a little larger)
    const bounds = imageBounds(image);
    const width = bounds.maxX - bounds.minX;
    const height = bounds.maxY - bounds.minY;
    const centre = new THREE.Vector3(
      bounds.minX + width / 2,
      bounds.minY + height / 2,
      0,
    );
    const aspect = container.clientWidth / container.clientHeight;
    const distance =
      (1.05 * Math.max(height, width / aspect)) /
      (2 * Math.tan(THREE.MathUtils.degToRad(75) / 2));
    let init: ReturnType<typeof initializeScene>;
    let host: LayerHost | null = null;

    try {
      init = initializeScene(container, {
        is2D: true,
        cameraPosition: new THREE.Vector3(centre.x, centre.y, distance),
        lookAtPosition: centre,
        near: distance / 1e5,
        far: distance * 10,
        onFrame: (now) => host?.tick(now),
      });
    } catch (e) {
      console.error("[ImageViewer] WebGL init failed:", e);
      setWebglFailed(true);

      return;
    }
    // Past the far plane the image would vanish
    init.controls.maxDistance = distance * 5;

    const group = new THREE.Group();

    init.scene.add(group);
    host = new LayerHost(group, init.camera, init.renderer, init.invalidate);
    const removeLayers = onSceneReady(host.handle);

    cameraRef.current = init.camera;
    rendererRef.current = init.renderer;
    controlsRef.current = init.controls;
    init.animate();

    return () => {
      removeLayers();
      host?.dispose();
      init.dispose();
      cameraRef.current = null;
      rendererRef.current = null;
      controlsRef.current = null;
    };
  }, [image]);

  return (
    // Fills its side of the split it is always placed in
    <div className="absolute inset-0 overflow-hidden">
      <div
        ref={containerRef}
        className="absolute inset-0"
        data-testid="image-scene-canvas"
      />
      {webglFailed && <WebGLErrorOverlay />}
      <LayerControls slot="panel" />
      {layers.image && (
        <SpatialScaleBar
          cameraRef={cameraRef}
          controlsRef={controlsRef}
          rendererRef={rendererRef}
        />
      )}
    </div>
  );
}
