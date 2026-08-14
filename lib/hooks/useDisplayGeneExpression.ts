import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { useEffect, useState } from "react";

import {
  detectExpressionKind,
  transformExpression,
  type DetectedExpressionKind,
  type ExpressionScaleMode,
} from "@/lib/utils/expression-scale";

/**
 * Fetch a gene's expression and return it in the current display space
 * (raw counts vs log1p), so the plot panel matches the 3D scene's scale-bar
 * toggle. `detectedKind` normally comes from the store (set when the scene
 * renders the gene); we fall back to sniffing the fetched values so the plots
 * are still correct if the scene hasn't run yet.
 */
export function useDisplayGeneExpression(
  dataset: StandardizedDataset,
  gene: string | null,
  mode: ExpressionScaleMode,
  detectedKind: DetectedExpressionKind | null,
): { expression: number[] | null; loading: boolean } {
  const [expression, setExpression] = useState<number[] | null>(null);
  const [loading, setLoading] = useState(false);

  useEffect(() => {
    if (!gene) {
      setExpression(null);

      return;
    }
    let cancelled = false;

    setLoading(true);
    dataset
      .getGeneExpression(gene)
      .then((vals) => {
        if (cancelled) return;
        if (vals) {
          const kind = detectedKind ?? detectExpressionKind(vals);

          setExpression(transformExpression(vals, kind, mode));
        } else {
          setExpression(null);
        }
        setLoading(false);
      })
      .catch(() => {
        if (!cancelled) {
          setExpression(null);
          setLoading(false);
        }
      });

    return () => {
      cancelled = true;
    };
  }, [dataset, gene, mode, detectedKind]);

  return { expression, loading };
}
