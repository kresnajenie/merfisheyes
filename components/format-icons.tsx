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

/** Zarr — chunked-grid glyph. */
export function ZarrIcon(props: SVGProps<SVGSVGElement>) {
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
        fillOpacity="0.6"
        height="8"
        rx="1.5"
        width="8"
        x="6"
        y="6"
      />
      <rect
        fill="currentColor"
        fillOpacity="0.4"
        height="8"
        rx="1.5"
        width="8"
        x="18"
        y="6"
      />
      <rect
        fill="currentColor"
        fillOpacity="0.4"
        height="8"
        rx="1.5"
        width="8"
        x="6"
        y="18"
      />
      <rect
        fill="currentColor"
        fillOpacity="0.6"
        height="8"
        rx="1.5"
        width="8"
        x="18"
        y="18"
      />
    </svg>
  );
}

/** Xenium — labeled tile (real 10x brand mark not shipped). */
export function XeniumIcon(props: SVGProps<SVGSVGElement>) {
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
      <rect fill="currentColor" fillOpacity="0.15" height="32" rx="6" width="32" />
      <text
        fill="currentColor"
        fontFamily="ui-sans-serif, system-ui"
        fontSize="13"
        fontWeight="700"
        textAnchor="middle"
        x="16"
        y="20"
      >
        Xn
      </text>
    </svg>
  );
}

/** MERSCOPE — labeled tile. */
export function MerscopeIcon(props: SVGProps<SVGSVGElement>) {
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
      <rect fill="currentColor" fillOpacity="0.15" height="32" rx="6" width="32" />
      <text
        fill="currentColor"
        fontFamily="ui-sans-serif, system-ui"
        fontSize="13"
        fontWeight="700"
        textAnchor="middle"
        x="16"
        y="20"
      >
        Ms
      </text>
    </svg>
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
