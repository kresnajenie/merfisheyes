"use client";

import type { LayersStore } from "../stores/createLayersStore";
import type { ImageChannel, ImageLayer } from "../spatialdata/ImageLayer";

import { useEffect } from "react";
import { useStore } from "zustand";

import { layersStore } from "../stores/layersStore";
import { useSplitScreenStore } from "../stores/splitScreenStore";
import { encodeImageState } from "../spatialdata/imageUrlState";
import { scheduleUrlUpdate } from "../utils/url-state-writer";

/**
 * What to change in `to` so its channels are on and off as `from`'s are,
 * matched by label: two different images share what they have in common
 * (DAPI is DAPI in both), and a channel only one has is left as it is.
 */
export function channelVisibilityChanges(
  from: Pick<ImageChannel, "label" | "visible">[],
  to: Pick<ImageChannel, "label" | "visible">[],
): { index: number; visible: boolean }[] {
  const on = new Map(from.map((c) => [c.label, c.visible]));

  return to.flatMap((c, index) => {
    const visible = on.get(c.label);

    return visible === undefined || visible === c.visible
      ? []
      : [{ index, visible }];
  });
}

/**
 * Bidirectional sync of which image channels are on, between the left
 * panel's image layer and the right one's. Active when syncEnabled and both
 * panels have an image (`available`, from useSyncCapabilities). Colours,
 * windows and the camera stay each panel's own.
 *
 * As the other sync hooks: turning sync on brings the right panel to the
 * left's; a link with sync=1 first lets each panel apply its own img= / rimg=
 * for 3 s. A change made here is written to the changed panel's link too, as
 * the image controls write the ones made by hand.
 */
export function useSyncImageChannels(
  rightLayersStore: LayersStore,
  available: boolean,
) {
  const syncEnabled = useSplitScreenStore((s) => s.syncEnabled);
  const left = useStore(layersStore, (s) => s.layers.image) as
    | ImageLayer
    | undefined;
  const right = useStore(rightLayersStore, (s) => s.layers.image) as
    | ImageLayer
    | undefined;

  const isActive = syncEnabled && available && !!left && !!right;

  useEffect(() => {
    if (!isActive || !left || !right) return;

    let settling = useSplitScreenStore.getState().syncFromUrl;
    let applying = false;

    const mirror = (
      from: ImageLayer,
      to: ImageLayer,
      slot: "image" | "rightImage",
    ) => {
      if (settling || applying) return;
      const changes = channelVisibilityChanges(from.channels, to.channels);

      if (!changes.length) return;
      applying = true;
      try {
        for (const { index, visible } of changes) {
          to.setChannel(index, { visible });
        }
      } finally {
        applying = false;
      }
      scheduleUrlUpdate(slot, encodeImageState(to));
    };

    // A layer pushes every change it makes (contrast found, a tile in); the
    // comparison makes the ones that aren't about visibility cost nothing,
    // and a mirrored change coming back finds nothing left to do
    const offLeft = left.subscribe(() => mirror(left, right, "rightImage"));
    const offRight = right.subscribe(() => mirror(right, left, "image"));
    let settleTimer: ReturnType<typeof setTimeout> | undefined;

    if (settling) {
      settleTimer = setTimeout(() => {
        settling = false;
        useSplitScreenStore.getState().setSyncFromUrl(false);
      }, 3000);
    } else {
      mirror(left, right, "rightImage");
    }

    return () => {
      offLeft();
      offRight();
      if (settleTimer) clearTimeout(settleTimer);
    };
  }, [isActive, left, right]);
}
