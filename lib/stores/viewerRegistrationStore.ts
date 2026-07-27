import { create } from "zustand";

import type { ViewerConfig } from "@/lib/utils/viewer-config";

/**
 * DB registration for the dataset currently open in the (left/main) viewer:
 * its database id, owner, saved viewer config, and — for from-s3 loads — the
 * S3 URL so it can be registered/claimed. Populated by the viewer pages and
 * read by the owner Save-defaults control and the "claim ownership" banner.
 *
 * Scoped to the primary viewer (not per split panel) — owner edits target the
 * main dataset.
 */
interface ViewerRegistrationState {
  dbId: string | null; // Dataset.id, or null when unregistered
  ownerId: string | null;
  adminOwned: boolean; // owned by "admin" collectively (any admin can edit)
  registered: boolean;
  viewerConfig: ViewerConfig | null;
  s3Url: string | null; // set for from-s3 viewers (needed to register/claim)
  set: (patch: Partial<Omit<ViewerRegistrationState, "set" | "reset">>) => void;
  reset: () => void;
}

const initial = {
  dbId: null,
  ownerId: null,
  adminOwned: false,
  registered: false,
  viewerConfig: null,
  s3Url: null,
};

export const useViewerRegistrationStore = create<ViewerRegistrationState>(
  (set) => ({
    ...initial,
    set: (patch) => set(patch),
    reset: () => set(initial),
  }),
);
