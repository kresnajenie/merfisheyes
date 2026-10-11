"use client";

import { useEffect, useRef, useState } from "react";
import { Spinner } from "@heroui/react";
import * as THREE from "three";

import { subtitle, title } from "@/components/primitives";
import { SpatialDataImageControls } from "@/components/spatialdata-image-controls";
import { SpatialScaleBar } from "@/components/spatial-scale-bar";
import { WebGLErrorOverlay } from "@/components/webgl-error-overlay";
import { ImageLayer } from "@/lib/spatialdata/ImageLayer";
import {
  openOmeZarrImage,
  type SpatialDataImage,
} from "@/lib/spatialdata/image";
import { initializeScene } from "@/lib/webgl/scene-manager";

/**
 * Viewer for an image on its own: an OME-Zarr multiscale image at a public
 * URL, with no cells or molecules.
 */
export function ImageViewer({ url }: { url: string }) {
  const [image, setImage] = useState<SpatialDataImage | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;

    setImage(null);
    setError(null);
    openOmeZarrImage(url)
      .then((opened) => !cancelled && setImage(opened))
      .catch((e) => {
        console.error("Error opening image:", e);
        if (!cancelled) {
          setError(
            "It must be a public OME-Zarr image (Zarr v3, axes c/y/x) that allows this site in CORS.",
          );
        }
      });

    return () => {
      cancelled = true;
    };
  }, [url]);

  if (error) {
    return (
      <div className="flex flex-col items-center justify-center h-full gap-4 p-8 text-center">
        <h2 className={title({ size: "md", color: "pink" })}>
          Couldn&apos;t read an image at that URL
        </h2>
        <p className={subtitle()}>{error}</p>
      </div>
    );
  }
  if (!image) {
    return (
      <div className="fixed inset-0 flex flex-col items-center justify-center gap-4">
        <Spinner color="primary" size="lg" />
        <p className={subtitle()}>Opening image...</p>
      </div>
    );
  }

  return <ImageScene image={image} />;
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
  const [layer, setLayer] = useState<ImageLayer | null>(null);
  const [webglFailed, setWebglFailed] = useState(false);

  useEffect(() => {
    const container = containerRef.current;

    if (!container) return;

    // Frame the whole image (pixel i of the finest level spans i ± 0.5)
    const finest = image.levels[0];
    const width = finest.width * finest.pixelWidth;
    const height = finest.height * finest.pixelHeight;
    const centre = new THREE.Vector3(
      finest.originX - finest.pixelWidth / 2 + width / 2,
      finest.originY - finest.pixelHeight / 2 + height / 2,
      0,
    );
    const aspect = container.clientWidth / container.clientHeight;
    const distance =
      (1.05 * Math.max(height, width / aspect)) /
      (2 * Math.tan(THREE.MathUtils.degToRad(75) / 2));
    let init: ReturnType<typeof initializeScene>;

    try {
      init = initializeScene(container, {
        is2D: true,
        cameraPosition: new THREE.Vector3(centre.x, centre.y, distance),
        lookAtPosition: centre,
        near: distance / 1e5,
        far: distance * 10,
      });
    } catch (e) {
      console.error("[ImageViewer] WebGL init failed:", e);
      setWebglFailed(true);

      return;
    }
    // Past the far plane the image would vanish
    init.controls.maxDistance = distance * 5;

    const imageLayer = new ImageLayer(
      image,
      init.scene,
      init.camera,
      init.renderer,
    );

    cameraRef.current = init.camera;
    rendererRef.current = init.renderer;
    controlsRef.current = init.controls;
    setLayer(imageLayer);
    init.animate();

    return () => {
      imageLayer.dispose();
      setLayer(null);
      init.dispose();
      cameraRef.current = null;
      rendererRef.current = null;
      controlsRef.current = null;
    };
  }, [image]);

  return (
    <div className="fixed inset-0 overflow-hidden">
      <div
        ref={containerRef}
        className="absolute inset-0"
        data-testid="image-scene-canvas"
      />
      {webglFailed && <WebGLErrorOverlay />}
      {layer && <SpatialDataImageControls layer={layer} />}
      {layer && (
        <SpatialScaleBar
          cameraRef={cameraRef}
          controlsRef={controlsRef}
          rendererRef={rendererRef}
        />
      )}
    </div>
  );
}
