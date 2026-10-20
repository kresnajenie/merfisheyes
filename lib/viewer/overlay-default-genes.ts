import type { ViewerConfig } from "@/lib/utils/viewer-config";

import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";

/**
 * Which genes a cell viewer's molecule overlay shows first: the cell
 * dataset's own saved overlay defaults, else the overlay dataset's saved
 * defaults (looked up by the public URL it was loaded from), else null for
 * the heuristic. Only the main viewer's page reads the overlay's source URL,
 * from the global store.
 */
export function overlayDefaultGenes(config: ViewerConfig | null) {
  return async (smDs: { uniqueGenes: string[] }): Promise<string[] | null> => {
    const sc = (config?.defaultGenes ?? []).filter((g) =>
      smDs.uniqueGenes.includes(g),
    );

    if (sc.length > 0) return sc;
    const url = useSingleMoleculeStore.getState().overlaySourceUrl;

    if (url) {
      try {
        const res = await fetch(
          `/api/datasets/by-url?url=${encodeURIComponent(url)}`,
        );

        if (res.ok) {
          const j = await res.json();
          const sm = (
            (j?.viewerConfig?.defaultGenes as string[] | undefined) ?? []
          ).filter((g) => smDs.uniqueGenes.includes(g));

          if (sm.length > 0) return sm;
        }
      } catch {
        // ignore — fall through to the heuristic
      }
    }

    return null;
  };
}
