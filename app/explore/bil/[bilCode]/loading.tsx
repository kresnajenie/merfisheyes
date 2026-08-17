import { DatasetDetailSkeleton } from "@/components/explore/dataset-detail-skeleton";

// Detail-route loading UI for BIL datasets. Overrides the parent /explore grid
// skeleton with a shape-matched detail placeholder.
export default function BilDatasetLoading() {
  return <DatasetDetailSkeleton />;
}
