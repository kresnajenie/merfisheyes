import { SVGProps } from "react";
import { z } from "zod";

export type IconSvgProps = SVGProps<SVGSVGElement> & {
  size?: number;
};

// Zod schemas for data validation
export const SpatialDataSchema = z.object({
  coordinates: z.array(z.array(z.number())),
  dimensions: z.union([z.literal(2), z.literal(3)]),
});

export const ClusterDataSchema = z.object({
  column: z.string(),
  values: z.array(z.any()),
  palette: z.record(z.string(), z.string()),
});

export const StandardizedDatasetSchema = z.object({
  id: z.string(),
  name: z.string(),
  type: z.string(),
  spatial: SpatialDataSchema,
  embeddings: z.record(z.string(), z.array(z.array(z.number()))),
  genes: z.array(z.string()),
  clusters: ClusterDataSchema.nullable(),
  metadata: z.record(z.string(), z.any()),
  rawData: z.any().optional(),
  adapter: z.any().optional(),
});

// TypeScript types derived from Zod schemas
export type SpatialData = z.infer<typeof SpatialDataSchema>;
export type ClusterData = z.infer<typeof ClusterDataSchema>;
export type StandardizedDatasetType = z.infer<typeof StandardizedDatasetSchema>;
