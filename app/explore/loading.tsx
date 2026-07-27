import { title } from "@/components/primitives";
import { ExploreBackground } from "@/components/explore/explore-background";

/**
 * Route-level loading UI for /explore.
 *
 * Without this, a client-side navigation to /explore blocks on the page's
 * ~8 server-side DB queries before anything renders — a visible pause after
 * the click. This boundary lets Next.js paint the header + skeleton instantly
 * and stream the real content in when the queries resolve.
 */
export default function ExploreLoading() {
  return (
    <>
      <ExploreBackground />
      <div className="w-full relative z-10">
        <div className="mb-8">
          <h1 className={title({ size: "lg" })}>
            <span className={title({ color: "violet", size: "lg" })}>
              Explore
            </span>{" "}
            Datasets
          </h1>
          <p className="text-default-500 mt-4">
            Browse and visualize our curated collection of spatial
            transcriptomics datasets
          </p>
        </div>

        <div
          aria-label="Loading datasets"
          className="grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-3"
          role="status"
        >
          {Array.from({ length: 6 }).map((_, i) => (
            <div
              key={i}
              className="h-48 animate-pulse rounded-2xl border border-default-200 bg-default-100/40"
            />
          ))}
        </div>
      </div>
    </>
  );
}
