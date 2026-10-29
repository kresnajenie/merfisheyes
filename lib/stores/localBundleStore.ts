import type { DatasetLinkConfig } from "@/lib/config/dataset-links";
import type { LocalBundle } from "@/lib/local/classify-bundle";

/** Molecules a local preview is checked against, as a saved overlay links them. */
export type CheckOverlay = DatasetLinkConfig & {
  /** The dataset they are, when picked from the user's own. */
  datasetId?: string;
  /** What the dialog calls them. */
  title?: string;
};
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import type { SpatialDataImage } from "@/lib/spatialdata/image";
import type { CellPolygons } from "@/lib/spatialdata/shapes";

import { create } from "zustand";

/**
 * What was opened from a folder dropped for a browser preview, beyond the one
 * dataset the viewer page shows: the molecules drawn over the cells, the image
 * drawn under them (or under the molecules, or on its own), and the dropped
 * files themselves, which is what "Upload & Save" sends as they are.
 *
 * `anchorId` is the id of the dataset the extras belong with — the cells, or
 * the molecules when there are no cells — so a viewer showing some other
 * dataset leaves them alone.
 */
interface LocalBundleState {
  bundle: LocalBundle | null;
  anchorId: string | null;
  molecules: SingleMoleculeDataset | null;
  image: SpatialDataImage | null;
  shapes: (() => Promise<CellPolygons>) | null;
  /**
   * Molecules picked to check the preview against (the overlay dialog's
   * local mode): drawn over the cells instead of any dropped with them, and
   * not saved — Upload & Save sends only the dropped files.
   */
  check: CheckOverlay | null;
  set: (patch: Partial<Omit<LocalBundleState, "set" | "reset">>) => void;
  reset: () => void;
}

const initial = {
  bundle: null,
  anchorId: null,
  molecules: null,
  image: null,
  shapes: null,
  check: null,
};

export const useLocalBundleStore = create<LocalBundleState>((set) => ({
  ...initial,
  set: (patch) => set(patch),
  reset: () => set(initial),
}));
