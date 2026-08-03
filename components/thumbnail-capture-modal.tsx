"use client";

import {
  Modal,
  ModalContent,
  ModalHeader,
  ModalBody,
  ModalFooter,
} from "@heroui/modal";
import { Button } from "@heroui/button";
import { useCallback, useEffect, useRef, useState } from "react";
import { toast } from "react-toastify";

// Recommended thumbnail aspect ratio (width / height). Matches the explore
// detail card's large tile (aspect-[4/3]); cards elsewhere object-cover it.
const RATIO = 4 / 3;
const MIN_W = 40; // min crop width in displayed px

/** Grab the current 3D scene canvas as a JPEG blob (last rendered frame). */
function captureScene(viewerKind: "sc" | "sm"): Promise<Blob | null> {
  const selector =
    viewerKind === "sm"
      ? '[data-testid="sm-scene-canvas"] canvas'
      : '[data-testid="sc-scene-canvas"] canvas';
  const canvas = document.querySelector<HTMLCanvasElement>(selector);

  if (!canvas) return Promise.resolve(null);

  return new Promise((resolve) =>
    canvas.toBlob((b) => resolve(b), "image/jpeg", 0.9),
  );
}

interface Rect {
  x: number;
  y: number;
  w: number;
  h: number;
}

/**
 * Owner-only "set thumbnail from screenshot" flow: captures the scene, shows it
 * with an aspect-locked crop box (move + resize), and uploads the cropped region.
 */
export function ThumbnailCaptureModal({
  isOpen,
  onClose,
  dbId,
  viewerKind,
  onSaved,
}: {
  isOpen: boolean;
  onClose: () => void;
  dbId: string;
  viewerKind: "sc" | "sm";
  onSaved?: (url: string) => void;
}) {
  const [blob, setBlob] = useState<Blob | null>(null);
  const [preview, setPreview] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [imgSize, setImgSize] = useState<{ w: number; h: number } | null>(null);
  const [crop, setCrop] = useState<Rect | null>(null);

  const imgRef = useRef<HTMLImageElement | null>(null);
  const dragRef = useRef<{
    mode: "move" | "resize";
    startX: number;
    startY: number;
    orig: Rect;
  } | null>(null);

  // Capture the scene when the modal opens.
  useEffect(() => {
    if (!isOpen) return;
    let url: string | null = null;

    captureScene(viewerKind).then((b) => {
      if (!b) {
        toast.error("Couldn't capture the scene.");

        return;
      }
      url = URL.createObjectURL(b);
      setBlob(b);
      setPreview(url);
    });

    return () => {
      if (url) URL.revokeObjectURL(url);
      setBlob(null);
      setPreview(null);
      setImgSize(null);
      setCrop(null);
    };
  }, [isOpen, viewerKind]);

  // Center a default 4:3 crop once we know the displayed image size.
  const onImgLoad = () => {
    const el = imgRef.current;

    if (!el) return;
    const w = el.clientWidth;
    const h = el.clientHeight;

    setImgSize({ w, h });
    // Largest 4:3 box at 90% that fits inside the image, centered.
    let cw = w * 0.9;
    let ch = cw / RATIO;

    if (ch > h * 0.9) {
      ch = h * 0.9;
      cw = ch * RATIO;
    }
    setCrop({ x: (w - cw) / 2, y: (h - ch) / 2, w: cw, h: ch });
  };

  const clampCrop = useCallback(
    (r: Rect, size: { w: number; h: number }): Rect => {
      const w = Math.min(r.w, size.w);
      const h = w / RATIO;
      const x = Math.min(Math.max(0, r.x), size.w - w);
      const y = Math.min(Math.max(0, r.y), size.h - h);

      return { x, y, w, h };
    },
    [],
  );

  // Drag / resize handling on the window so it keeps tracking outside the box.
  useEffect(() => {
    if (!crop || !imgSize) return;
    const onMove = (e: MouseEvent) => {
      const d = dragRef.current;

      if (!d) return;
      const dx = e.clientX - d.startX;
      const dy = e.clientY - d.startY;

      if (d.mode === "move") {
        setCrop(
          clampCrop({ ...d.orig, x: d.orig.x + dx, y: d.orig.y + dy }, imgSize),
        );
      } else {
        // Resize from the bottom-right corner, aspect-locked.
        const maxW = Math.min(imgSize.w - d.orig.x, (imgSize.h - d.orig.y) * RATIO);
        const w = Math.max(MIN_W, Math.min(d.orig.w + dx, maxW));

        setCrop({ ...d.orig, w, h: w / RATIO });
      }
    };
    const onUp = () => {
      dragRef.current = null;
      document.body.style.userSelect = "";
    };

    window.addEventListener("mousemove", onMove);
    window.addEventListener("mouseup", onUp);

    return () => {
      window.removeEventListener("mousemove", onMove);
      window.removeEventListener("mouseup", onUp);
    };
  }, [crop, imgSize, clampCrop]);

  const startDrag = (mode: "move" | "resize") => (e: React.MouseEvent) => {
    if (!crop) return;
    e.preventDefault();
    e.stopPropagation();
    dragRef.current = {
      mode,
      startX: e.clientX,
      startY: e.clientY,
      orig: crop,
    };
    document.body.style.userSelect = "none";
  };

  const upload = async () => {
    if (!blob || !crop || !imgSize) return;
    setBusy(true);
    try {
      const src = await createImageBitmap(blob);
      const scaleX = src.width / imgSize.w;
      const scaleY = src.height / imgSize.h;
      const sw = Math.max(1, Math.round(crop.w * scaleX));
      const sh = Math.max(1, Math.round(crop.h * scaleY));
      const out = document.createElement("canvas");

      out.width = sw;
      out.height = sh;
      const ctx = out.getContext("2d");

      if (!ctx) {
        toast.error("Couldn't crop the image.");

        return;
      }
      ctx.drawImage(
        src,
        Math.round(crop.x * scaleX),
        Math.round(crop.y * scaleY),
        sw,
        sh,
        0,
        0,
        sw,
        sh,
      );
      const cropped = await new Promise<Blob | null>((resolve) =>
        out.toBlob((b) => resolve(b), "image/jpeg", 0.85),
      );

      if (!cropped) {
        toast.error("Couldn't crop the image.");

        return;
      }

      const res = await fetch(`/api/datasets/${dbId}/thumbnail`, {
        method: "POST",
        headers: { "Content-Type": "image/jpeg" },
        body: cropped,
      });

      if (!res.ok) {
        toast.error("Upload failed.");

        return;
      }
      const { thumbnailUrl } = await res.json();

      toast.success("Thumbnail saved.");
      onSaved?.(thumbnailUrl);
      onClose();
    } catch {
      toast.error("Upload failed.");
    } finally {
      setBusy(false);
    }
  };

  return (
    <Modal isOpen={isOpen} size="lg" onClose={onClose}>
      <ModalContent>
        <ModalHeader>Set thumbnail from screenshot</ModalHeader>
        <ModalBody>
          {preview ? (
            <div className="relative select-none">
              {/* eslint-disable-next-line @next/next/no-img-element */}
              <img
                ref={imgRef}
                alt="Scene preview"
                className="w-full rounded-lg border border-default-200 block"
                draggable={false}
                src={preview}
                onLoad={onImgLoad}
              />
              {crop && (
                <div
                  className="absolute border-2 border-primary cursor-move"
                  style={{
                    left: crop.x,
                    top: crop.y,
                    width: crop.w,
                    height: crop.h,
                    boxShadow: "0 0 0 9999px rgba(0,0,0,0.45)",
                  }}
                  onMouseDown={startDrag("move")}
                >
                  {/* Bottom-right resize handle */}
                  <div
                    className="absolute -right-1.5 -bottom-1.5 h-3.5 w-3.5 rounded-sm bg-primary border border-white cursor-nwse-resize"
                    onMouseDown={startDrag("resize")}
                  />
                </div>
              )}
            </div>
          ) : (
            <p className="text-sm text-default-500 py-8 text-center">
              Capturing…
            </p>
          )}
          <p className="text-xs text-default-400">
            Drag the box to move it, or the corner to resize (locked to a 4:3
            thumbnail ratio). To change the shot, cancel, reposition the view,
            then reopen.
          </p>
        </ModalBody>
        <ModalFooter>
          <Button isDisabled={busy} variant="light" onPress={onClose}>
            Cancel
          </Button>
          <Button
            color="primary"
            isDisabled={!blob || !crop}
            isLoading={busy}
            onPress={upload}
          >
            Upload as thumbnail
          </Button>
        </ModalFooter>
      </ModalContent>
    </Modal>
  );
}
