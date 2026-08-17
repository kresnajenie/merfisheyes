import { DatasetDetailSkeleton } from "@/components/explore/dataset-detail-skeleton";

// Detail-route loading UI. Overrides the parent /explore grid skeleton so a
// navigation into a dataset/project detail page shows a shape-matched
// placeholder instead of the Explore card grid.
export default function DatasetDetailLoading() {
  return <DatasetDetailSkeleton />;
}
