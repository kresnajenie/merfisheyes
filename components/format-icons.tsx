/**
 * Placeholder format glyphs for the upload UI.
 *
 * These are deliberately simple inline SVGs so we don't ship third-party
 * brand marks. Swap any of these for real brand logos later by replacing
 * the contents while keeping the same component API.
 */
import type { SVGProps } from "react";

const baseIcon = "rounded-md";

/**
 * AnnData / .h5ad — uses /logos/scverse.png. The scverse mark is monochrome
 * (black/gray on transparent), so we invert it to render light-on-dark.
 *
 * `className` is accepted for API compatibility with the SVG icons but isn't
 * threaded through (the tile drives the look).
 */
export function AnnDataIcon(_props: SVGProps<SVGSVGElement>) {
  return <LogoTile alt="AnnData / scverse" invert src="/logos/scverse.png" />;
}

/**
 * Brand-logo tile. The page is dark, so by default the wrapper uses a
 * glass-on-dark style (subtle border + faint white wash), which works for
 * transparent-PNG logos. For monochrome marks (e.g. scverse) we invert
 * the image so a black glyph becomes white. For raster logos with baked-in
 * white backgrounds (JPGs) we use `tone="light"` to keep a small white pill
 * — inverting would just flip the white square to black.
 */
type LogoTone = "dark" | "light";
function LogoTile({
  src,
  alt,
  className,
  tone = "dark",
  invert = false,
}: {
  src: string;
  alt: string;
  className?: string;
  tone?: LogoTone;
  invert?: boolean;
}) {
  const tile =
    tone === "light"
      ? "bg-white/95 border border-white/30"
      : "bg-white/10 border border-white/15 backdrop-blur-sm";

  return (
    <span
      className={`inline-flex items-center justify-center rounded-md ${tile} ${className ?? ""}`}
      style={{ width: 48, height: 48, padding: 6 }}
    >
      {/* Plain <img> so we don't pull next/image in for tiny static logos */}
      {/* eslint-disable-next-line @next/next/no-img-element */}
      <img
        alt={alt}
        className="max-w-full max-h-full object-contain"
        src={src}
        style={invert ? { filter: "invert(1) brightness(1.4)" } : undefined}
      />
    </span>
  );
}

/** Zarr — uses /logos/zarr.png. */
export function ZarrIcon(props: { className?: string }) {
  return (
    <LogoTile alt="Zarr" className={props.className} src="/logos/zarr.png" />
  );
}

/** Xenium — uses /logos/10x.png (transparent PNG, sits on the dark glass). */
export function XeniumIcon(props: { className?: string }) {
  return (
    <LogoTile
      alt="Xenium (10x Genomics)"
      className={props.className}
      src="/logos/10x.png"
    />
  );
}

/**
 * MERSCOPE — Vizgen logo is shipped as a JPG with a baked-in white
 * background, so we keep a white pill rather than inverting (which would
 * flip the whole tile black).
 */
export function MerscopeIcon(props: { className?: string }) {
  return (
    <LogoTile
      alt="MERSCOPE (Vizgen)"
      className={props.className}
      src="/logos/vizgen.png"
      tone="light"
    />
  );
}

/** Chunked (our preprocessed format) — stacked-blocks glyph at tile size. */
export function ChunkedIcon(props: SVGProps<SVGSVGElement>) {
  return (
    <span
      className="inline-flex items-center justify-center rounded-md bg-white/10 border border-white/15 backdrop-blur-sm"
      style={{ width: 48, height: 48, padding: 6 }}
    >
      <svg
        className={baseIcon}
        fill="none"
        height="32"
        viewBox="0 0 32 32"
        width="32"
        xmlns="http://www.w3.org/2000/svg"
        {...props}
      >
        <rect
          fill="currentColor"
          fillOpacity="0.7"
          height="4"
          rx="1"
          width="22"
          x="5"
          y="9"
        />
        <rect
          fill="currentColor"
          fillOpacity="0.5"
          height="4"
          rx="1"
          width="22"
          x="5"
          y="15"
        />
        <rect
          fill="currentColor"
          fillOpacity="0.3"
          height="4"
          rx="1"
          width="22"
          x="5"
          y="21"
        />
      </svg>
    </span>
  );
}

/** Cloud / S3 — used by the small "Load from S3" button. */
export function CloudIcon(props: SVGProps<SVGSVGElement>) {
  return (
    <svg
      fill="none"
      height="20"
      stroke="currentColor"
      strokeLinecap="round"
      strokeLinejoin="round"
      strokeWidth="1.6"
      viewBox="0 0 24 24"
      width="20"
      xmlns="http://www.w3.org/2000/svg"
      {...props}
    >
      <path d="M16 18a4 4 0 1 0-1-7.87A6 6 0 0 0 4 12a4 4 0 0 0 4 4h8z" />
      <path d="M12 13v5" />
      <path d="M9 16l3 3 3-3" />
    </svg>
  );
}
