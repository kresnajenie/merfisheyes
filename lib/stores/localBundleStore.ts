import type { LocalBundle } from "@/lib/local/classify-bundle";
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
  set: (patch: Partial<Omit<LocalBundleState, "set" | "reset">>) => void;
  reset: () => void;
}

const initial = {
  bundle: null,
  anchorId: null,
  molecules: null,
  image: null,
  shapes: null,
};

export const useLocalBundleStore = create<LocalBundleState>((set) => ({
  ...initial,
  set: (patch) => set(patch),
  reset: () => set(initial),
}));
