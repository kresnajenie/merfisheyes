/**
 * Loading skeleton for the dataset / project detail page
 * (components/explore/dataset-detail-page.tsx). Mirrors its layout — back
 * button, header, entries grid, details card — so a client navigation into a
 * detail route streams a shape-matched placeholder instead of the Explore
 * grid skeleton it would otherwise inherit from the parent segment.
 */
export function DatasetDetailSkeleton() {
  const bar = "animate-pulse rounded bg-default-100/50";

  return (
    <div aria-label="Loading dataset" className="space-y-6" role="status">
      {/* Back button */}
      <div className={`h-8 w-32 ${bar}`} />

      {/* Header */}
      <div className="space-y-3">
        <div className="flex gap-2">
          <div className={`h-6 w-20 rounded-full ${bar}`} />
          <div className={`h-6 w-16 rounded-full ${bar}`} />
          <div className={`h-6 w-16 rounded-full ${bar}`} />
        </div>
        <div className={`h-8 w-2/3 ${bar}`} />
        <div className={`h-4 w-1/2 ${bar}`} />
        <div className={`h-4 w-full max-w-2xl ${bar}`} />
      </div>

      {/* Entries — SC large on left, SM grid on right */}
      <div className="grid grid-cols-1 gap-4 md:grid-cols-3">
        <div className="space-y-2">
          <div className={`h-4 w-24 ${bar}`} />
          <div className={`aspect-[4/3] w-full rounded-xl ${bar}`} />
        </div>
        <div className="space-y-2 md:col-span-2">
          <div className={`h-4 w-40 ${bar}`} />
          <div className="grid grid-cols-2 gap-2 lg:grid-cols-3">
            {Array.from({ length: 6 }).map((_, i) => (
              <div
                key={i}
                className={`aspect-[16/10] w-full rounded-xl ${bar}`}
              />
            ))}
          </div>
        </div>
      </div>

      {/* Details card */}
      <div className="h-40 w-full animate-pulse rounded-2xl border border-default-200 bg-default-100/40" />
    </div>
  );
}
