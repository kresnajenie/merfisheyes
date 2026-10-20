"use client";

import type { ViewerConfig } from "@/lib/utils/viewer-config";

import { Checkbox } from "@heroui/checkbox";

/**
 * Collapsed "Advanced" section of the upload dialogs. Holds the defaults a new
 * dataset starts with — today only whether it opens in the 3D camera. Off by
 * default: thin-z datasets look flat but are technically 3D, so 2D stays the
 * default unless the uploader says otherwise. The owner can change it later
 * from the camera panel ("Save current view as default").
 */
export function UploadAdvancedOptions({
  open3D,
  onOpen3DChange,
}: {
  open3D: boolean;
  onOpen3DChange: (v: boolean) => void;
}) {
  return (
    <details className="group rounded-xl border border-default-200 bg-default-100/50 px-3 py-2">
      <summary className="cursor-pointer select-none text-sm text-default-600">
        Advanced
      </summary>
      <div className="mt-2 pb-1">
        <Checkbox isSelected={open3D} size="sm" onValueChange={onOpen3DChange}>
          <span className="text-sm">Open in 3D by default</span>
        </Checkbox>
        <p className="mt-1 ml-7 text-xs text-default-500">
          For datasets with real depth. Viewers start with the rotating 3D
          camera instead of the flat 2D view.
        </p>
      </div>
    </details>
  );
}

/** The viewerConfig to send with the upload, or undefined for no defaults. */
export function uploadViewerConfig(open3D: boolean): ViewerConfig | undefined {
  return open3D ? { version: 1, viewMode: "3D" } : undefined;
}
