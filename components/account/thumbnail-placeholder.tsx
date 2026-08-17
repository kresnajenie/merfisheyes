"use client";

/**
 * Fallback visual for cards whose dataset/project has no saved thumbnail.
 * A subtle gradient + scattered dots (evoking spatial points) so cards read as
 * intentional and uniform — matching the Explore card look — instead of empty.
 * Tinted by kind: single-cell (blue/primary), molecule (purple/secondary),
 * project (neutral).
 */
export function ThumbnailPlaceholder({
  kind,
  className = "",
}: {
  kind: "single_cell" | "single_molecule" | "project";
  className?: string;
}) {
  const from =
    kind === "single_molecule"
      ? "from-secondary-500/25"
      : kind === "project"
        ? "from-default-400/20"
        : "from-primary-500/25";
  const dot =
    kind === "single_molecule"
      ? "text-secondary-400/50"
      : kind === "project"
        ? "text-default-400/50"
        : "text-primary-400/50";

  return (
    <div
      className={`relative w-full bg-gradient-to-br ${from} to-transparent overflow-hidden ${className}`}
    >
      <svg
        aria-hidden
        className={`absolute inset-0 w-full h-full ${dot}`}
        preserveAspectRatio="xMidYMid slice"
        viewBox="0 0 100 60"
      >
        {[
          [18, 20],
          [30, 34],
          [44, 16],
          [58, 40],
          [70, 22],
          [82, 38],
          [26, 48],
          [50, 28],
          [64, 12],
          [38, 44],
          [74, 50],
          [88, 18],
        ].map(([cx, cy], i) => (
          <circle
            key={i}
            cx={cx}
            cy={cy}
            fill="currentColor"
            r={i % 3 === 0 ? 2.4 : 1.6}
          />
        ))}
      </svg>
    </div>
  );
}
