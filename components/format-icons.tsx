/**
 * Placeholder format glyphs for the upload UI.
 *
 * These are deliberately simple inline SVGs so we don't ship third-party
 * brand marks. Swap any of these for real brand logos later by replacing
 * the contents while keeping the same component API.
 */
import type { SVGProps } from "react";

const baseIcon = "rounded-md";

/** AnnData / .h5ad — abstract "data matrix in a cell" glyph. */
export function AnnDataIcon(props: SVGProps<SVGSVGElement>) {
  return (
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
        fillOpacity="0.1"
        height="32"
        rx="6"
        width="32"
      />
      <circle cx="16" cy="16" r="9" stroke="currentColor" strokeWidth="1.4" />
      <rect fill="currentColor" height="2" rx="1" width="2" x="11" y="11" />
      <rect fill="currentColor" height="2" rx="1" width="2" x="15" y="11" />
      <rect fill="currentColor" height="2" rx="1" width="2" x="19" y="11" />
      <rect fill="currentColor" height="2" rx="1" width="2" x="11" y="15" />
      <rect fill="currentColor" height="2" rx="1" width="2" x="15" y="15" />
      <rect fill="currentColor" height="2" rx="1" width="2" x="19" y="15" />
      <rect fill="currentColor" height="2" rx="1" width="2" x="11" y="19" />
      <rect fill="currentColor" height="2" rx="1" width="2" x="15" y="19" />
      <rect fill="currentColor" height="2" rx="1" width="2" x="19" y="19" />
    </svg>
  );
}

/**
 * Brand-logo tile. Wraps an image at /logos/{file} on a small white pill
 * so dark logos remain visible against the dark page background. Pixel
 * size matches the SVG glyphs (32px) so they sit in the same row cleanly.
 */
function LogoTile({
  src,
  alt,
  className,
}: {
  src: string;
  alt: string;
  className?: string;
}) {
  return (
    <span
      className={`inline-flex items-center justify-center bg-white/95 rounded-md ${className ?? ""}`}
      style={{ width: 32, height: 32, padding: 4 }}
    >
      {/* Plain <img> so we don't pull next/image in for tiny static logos */}
      {/* eslint-disable-next-line @next/next/no-img-element */}
      <img alt={alt} className="max-w-full max-h-full object-contain" src={src} />
    </span>
  );
}

/** Zarr — uses /logos/zarr.png. */
export function ZarrIcon(props: { className?: string }) {
  return <LogoTile alt="Zarr" className={props.className} src="/logos/zarr.png" />;
}

/** Xenium — uses /logos/10x.png. */
export function XeniumIcon(props: { className?: string }) {
  return (
    <LogoTile alt="Xenium (10x Genomics)" className={props.className} src="/logos/10x.png" />
  );
}

/** MERSCOPE — uses /logos/vizgen.jpg. */
export function MerscopeIcon(props: { className?: string }) {
  return (
    <LogoTile alt="MERSCOPE (Vizgen)" className={props.className} src="/logos/vizgen.jpg" />
  );
}

/** Chunked (our preprocessed format) — stacked-blocks glyph. */
export function ChunkedIcon(props: SVGProps<SVGSVGElement>) {
  return (
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
        fillOpacity="0.1"
        height="32"
        rx="6"
        width="32"
      />
      <rect
        fill="currentColor"
        fillOpacity="0.7"
        height="4"
        rx="1"
        width="18"
        x="7"
        y="8"
      />
      <rect
        fill="currentColor"
        fillOpacity="0.5"
        height="4"
        rx="1"
        width="18"
        x="7"
        y="14"
      />
      <rect
        fill="currentColor"
        fillOpacity="0.3"
        height="4"
        rx="1"
        width="18"
        x="7"
        y="20"
      />
    </svg>
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
