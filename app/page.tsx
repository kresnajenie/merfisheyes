"use client";

import { Button } from "@heroui/button";
import { useSession } from "next-auth/react";
import Link from "next/link";
import {
  useCallback,
  useEffect,
  useMemo,
  memo,
  useRef,
  useState,
  Suspense,
} from "react";
import clsx from "clsx";
import { usePathname, useRouter, useSearchParams } from "next/navigation";

import { title, subtitle } from "@/components/primitives";
import { FileUpload } from "@/components/file-upload";
import {
  AnnDataIcon,
  ChunkedIcon,
  CloudIcon,
  MerscopeIcon,
  XeniumIcon,
  ZarrIcon,
} from "@/components/format-icons";
import LightRays from "@/components/react-bits/LightRays";
import { BrainToggle } from "@/components/brain-toggle";
import { LoadFromS3Modal } from "@/components/load-from-s3-modal";
import { FeaturedDatasets } from "@/components/explore/featured-datasets";
import type { CatalogDatasetItem } from "@/components/explore/types";

const MemoizedLightRays = memo(LightRays);

// The two data-mode pills differ only by accent colour, so they share one
// component rather than two near-identical class-string copies.
const MODE_ACCENT = {
  blue: {
    ring: "focus-visible:ring-blue-400/80",
    active: "bg-blue-500 text-white shadow-[0_12px_25px_rgba(59,130,246,0.35)]",
    idle: "bg-transparent text-default-500 border border-blue-400/30 hover:bg-blue-500/15 hover:text-white",
  },
  purple: {
    ring: "focus-visible:ring-purple-400/80",
    active:
      "bg-purple-500 text-white shadow-[0_12px_25px_rgba(168,85,247,0.35)]",
    idle: "bg-transparent text-default-500 border border-purple-400/30 hover:bg-purple-500/15 hover:text-white",
  },
} as const;

function ModeToggleButton({
  active,
  accent,
  label,
  onClick,
}: {
  active: boolean;
  accent: keyof typeof MODE_ACCENT;
  label: string;
  onClick: () => void;
}) {
  const a = MODE_ACCENT[accent];

  return (
    <button
      aria-pressed={active}
      className={clsx(
        "px-5 py-2.5 md:px-6 md:py-3 rounded-full text-sm md:text-base font-semibold tracking-[0.2em] uppercase transition-colors duration-300 cursor-pointer focus:outline-none focus-visible:ring-2 focus-visible:ring-offset-2 focus-visible:ring-offset-black",
        a.ring,
        active ? a.active : a.idle,
      )}
      type="button"
      onClick={onClick}
    >
      {label}
    </button>
  );
}

// One option in the "where does processing happen" segmented control. Kept
// large and clearly labelled (not a 12px switch) so the choice — and which
// side is active — is obvious at a glance.
function UploadModeOption({
  active,
  label,
  onClick,
}: {
  active: boolean;
  label: string;
  onClick: () => void;
}) {
  return (
    <button
      aria-checked={active}
      className={clsx(
        "px-3.5 py-1.5 rounded-full text-xs md:text-sm font-medium transition-colors duration-200 focus:outline-none focus-visible:ring-2 focus-visible:ring-primary/60 focus-visible:ring-offset-2 focus-visible:ring-offset-black cursor-pointer",
        active
          ? "bg-primary text-white shadow-[0_6px_16px_rgba(59,130,246,0.3)]"
          : "text-default-600 hover:text-foreground",
      )}
      role="radio"
      type="button"
      onClick={onClick}
    >
      {label}
    </button>
  );
}

function LoadFromS3Button({ onClick }: { onClick: () => void }) {
  return (
    <button
      className="inline-flex items-center gap-2 px-4 py-2 rounded-full text-sm border border-default-300 text-default-700 bg-default-100/40 hover:bg-default-200/60 hover:border-primary/40 focus:outline-none focus-visible:ring-2 focus-visible:ring-primary/50 transition-colors"
      type="button"
      onClick={onClick}
    >
      <CloudIcon />
      <span>Load from S3</span>
    </button>
  );
}

function useModeToggleState() {
  const router = useRouter();
  const pathname = usePathname();
  const searchParams = useSearchParams();
  const modeParam = searchParams.get("mode");

  const initialMode = useMemo(() => {
    if (modeParam) {
      return modeParam === "sm";
    }

    return false;
  }, [modeParam]);

  const [isSingleMolecule, setIsSingleMolecule] = useState(initialMode);

  useEffect(() => {
    setIsSingleMolecule(initialMode);
  }, [initialMode]);

  const handleModeChange = useCallback(
    (selected: boolean) => {
      setIsSingleMolecule(selected);

      const params = new URLSearchParams(searchParams.toString());

      if (selected) {
        params.set("mode", "sm");
      } else {
        params.delete("mode");
      }

      const queryString = params.toString();
      const newUrl = queryString ? `${pathname}?${queryString}` : pathname;

      router.replace(newUrl, { scroll: false });
    },
    [pathname, router, searchParams],
  );

  return { isSingleMolecule, handleModeChange };
}

function HomeContent() {
  const name = "MERFISH";
  const { isSingleMolecule, handleModeChange } = useModeToggleState();
  const animationFrameRef = useRef<number | null>(null);
  const currentColorRef = useRef("#5EA2EF");
  const [animatedRaysColor, setAnimatedRaysColor] = useState("#5EA2EF");
  const [isS3ModalOpen, setIsS3ModalOpen] = useState(false);
  const [featured, setFeatured] = useState<CatalogDatasetItem[]>([]);

  useEffect(() => {
    let cancelled = false;

    fetch("/api/explore?limit=1&include=meta")
      .then((r) => (r.ok ? r.json() : null))
      .then((d) => {
        if (!cancelled && d?.featured) setFeatured(d.featured);
      })
      .catch(() => {});

    return () => {
      cancelled = true;
    };
  }, []);
  const { data: session } = useSession();
  const [serverMode, setServerMode] = useState(false);
  // Optional per-cell label CSV (e.g. a MapMyCells result the user generated
  // themselves). Uploaded with the dataset; the chunk stage merges its columns
  // into obs, with palettes and DE stats like any other cluster column.
  const [annotationCsv, setAnnotationCsv] = useState<File | null>(null);
  const [annotationDragOver, setAnnotationDragOver] = useState(false);

  const targetRaysColor = useMemo(
    () => (isSingleMolecule ? "#FF1CF7" : "#5EA2EF"),
    [isSingleMolecule],
  );

  useEffect(() => {
    const duration = 650;

    const hexToRgb = (hex: string): [number, number, number] => {
      const parsed = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);

      if (!parsed) return [0, 0, 0];

      return [
        parseInt(parsed[1], 16),
        parseInt(parsed[2], 16),
        parseInt(parsed[3], 16),
      ];
    };

    const rgbToHex = (rgb: [number, number, number]) => {
      const toHex = (value: number) => {
        const clamped = Math.round(Math.min(Math.max(value, 0), 255));
        const hex = clamped.toString(16);

        return hex.length === 1 ? `0${hex}` : hex;
      };

      return `#${rgb.map(toHex).join("")}`;
    };

    const startColor = currentColorRef.current;

    if (startColor === targetRaysColor) return;

    const start = performance.now();
    const startRgb = hexToRgb(startColor);
    const endRgb = hexToRgb(targetRaysColor);

    const animate = (now: number) => {
      const progress = Math.min((now - start) / duration, 1);

      const interpolated: [number, number, number] = [
        startRgb[0] + (endRgb[0] - startRgb[0]) * progress,
        startRgb[1] + (endRgb[1] - startRgb[1]) * progress,
        startRgb[2] + (endRgb[2] - startRgb[2]) * progress,
      ];

      const nextColor = rgbToHex(interpolated);

      currentColorRef.current = nextColor;
      setAnimatedRaysColor(nextColor);

      if (progress < 1) {
        animationFrameRef.current = requestAnimationFrame(animate);

        return;
      }

      currentColorRef.current = targetRaysColor;
      setAnimatedRaysColor(targetRaysColor);
      animationFrameRef.current = null;
    };

    animationFrameRef.current = requestAnimationFrame(animate);

    return () => {
      if (animationFrameRef.current) {
        cancelAnimationFrame(animationFrameRef.current);
        animationFrameRef.current = null;
      }
    };
  }, [targetRaysColor]);

  return (
    <>
      <div className="fixed inset-0 w-full h-full z-0">
        <MemoizedLightRays
          lightSpread={1.0}
          mouseInfluence={0.1}
          pulsating={false}
          rayLength={10}
          raysColor={animatedRaysColor}
          raysOrigin="top-center"
          raysSpeed={1.0}
        />
      </div>
      <section className="relative flex flex-col items-center gap-1 -mt-4 py-3 md:py-4 px-4 md:px-8">
        <div className="relative z-10 flex flex-col items-center gap-4 max-w-3xl w-full">
          <div className="flex flex-col items-center">
            <h1 className="flex flex-col items-center text-center">
              <span className={title({ size: "lg" })}>
                Bring{" "}
                <span
                  className={title({
                    color: isSingleMolecule ? "violet" : "blue",
                    size: "lg",
                  })}
                  style={{ transition: "color 0.3s ease" }}
                >
                  {name}
                </span>{" "}
                to Life
              </span>
            </h1>
            <div className="flex items-center gap-6 mt-6 justify-center w-full">
              <ModeToggleButton
                accent="blue"
                active={!isSingleMolecule}
                label="single cell"
                onClick={() => handleModeChange(false)}
              />
              <BrainToggle
                isActive={isSingleMolecule}
                onToggle={handleModeChange}
              />
              <ModeToggleButton
                accent="purple"
                active={isSingleMolecule}
                label="single molecule"
                onClick={() => handleModeChange(true)}
              />
            </div>
            <div className={subtitle({ class: "mt-4 text-center" })}>
              Explore your{" "}
              <span className="font-bold">
                {isSingleMolecule
                  ? ".parquet, .csv, or Xenium"
                  : ".h5ad, Merscope & Xenium"}
              </span>{" "}
              datasets
            </div>

            <div className="mt-5 flex flex-col items-center gap-2">
              <div
                aria-label="Where to process your dataset"
                className="inline-flex items-center gap-1 rounded-full border border-default-200 bg-default-100/40 p-1"
                role="radiogroup"
              >
                <UploadModeOption
                  active={!serverMode}
                  label="Preview in browser"
                  onClick={() => setServerMode(false)}
                />
                <UploadModeOption
                  active={serverMode}
                  label="Upload & process on server"
                  onClick={() => setServerMode(true)}
                />
              </div>

              <p className="max-w-2xl text-center text-sm text-default-500">
                {serverMode
                  ? session?.user
                    ? "Uploads to the server; we email you a link when it's ready. Best for large datasets."
                    : "Drop a file, verify with an email code, and we'll email you a link when it's ready."
                  : "Parsed in your browser and opens instantly. Nothing leaves your computer."}
              </p>

              {/* MapMyCells cell-type CSV is single-cell only — there are no
                  cells to annotate in single-molecule data. */}
              {serverMode && !isSingleMolecule ? (
                <div className="mt-2 flex flex-col items-center gap-2">
                  <span className="text-sm font-semibold text-default-600">
                    <a
                      className="text-primary underline"
                      href="https://knowledge.brain-map.org/mapmycells/process"
                      rel="noreferrer"
                      target="_blank"
                    >
                      MapMyCells
                    </a>{" "}
                    celltype CSV{" "}
                    <span className="group relative inline-block">
                      <span className="cursor-help font-normal text-default-500 underline decoration-dotted underline-offset-2">
                        (optional)
                      </span>
                      <span
                        className="pointer-events-none absolute bottom-full left-1/2 z-[var(--z-modal)] mb-2 hidden w-64 -translate-x-1/2 rounded-lg border border-default-200 bg-default-100 px-3 py-2 text-xs font-normal leading-snug text-default-600 shadow-lg group-hover:block"
                        role="tooltip"
                      >
                        Upload a MapMyCells result (or any per-cell label CSV)
                        to add cell type columns.
                      </span>
                    </span>
                  </span>
                  <label
                    className={clsx(
                      "flex w-full max-w-md cursor-pointer flex-col items-center justify-center gap-1 rounded-xl border-2 border-dashed px-4 py-6 text-center transition-colors",
                      annotationDragOver
                        ? "border-primary bg-primary-50 dark:bg-primary-100/10"
                        : "border-default-300 hover:border-primary/50 hover:bg-default-100/50",
                    )}
                    htmlFor="annotation-csv"
                    onDragLeave={() => setAnnotationDragOver(false)}
                    onDragOver={(e) => {
                      e.preventDefault();
                      setAnnotationDragOver(true);
                    }}
                    onDrop={(e) => {
                      e.preventDefault();
                      setAnnotationDragOver(false);
                      setAnnotationCsv(e.dataTransfer.files?.[0] ?? null);
                    }}
                  >
                    <svg
                      className="h-6 w-6 text-default-400"
                      fill="none"
                      stroke="currentColor"
                      strokeWidth={1.5}
                      viewBox="0 0 24 24"
                    >
                      <path
                        d="M12 16.5V6m0 0-3.5 3.5M12 6l3.5 3.5M4.5 15.75v1.5A2.25 2.25 0 0 0 6.75 19.5h10.5a2.25 2.25 0 0 0 2.25-2.25v-1.5"
                        strokeLinecap="round"
                        strokeLinejoin="round"
                      />
                    </svg>
                    {annotationCsv ? (
                      <span className="text-sm text-default-600">
                        {annotationCsv.name}
                        <button
                          className="ml-2 underline hover:text-foreground"
                          type="button"
                          onClick={(e) => {
                            e.preventDefault();
                            setAnnotationCsv(null);
                          }}
                        >
                          Remove
                        </button>
                      </span>
                    ) : (
                      <>
                        <span className="text-sm font-medium text-default-700">
                          Drop a .csv here, or click to browse
                        </span>
                        <span className="text-[11px] text-default-400">
                          One row per cell with a label column
                        </span>
                      </>
                    )}
                    <input
                      accept=".csv"
                      className="hidden"
                      id="annotation-csv"
                      type="file"
                      onChange={(e) =>
                        setAnnotationCsv(e.target.files?.[0] ?? null)
                      }
                    />
                  </label>
                </div>
              ) : null}
            </div>
          </div>
        </div>

        <div className="relative z-10 flex flex-col items-center gap-4 max-w-5xl w-full">
          <div className="relative w-full max-w-5xl min-h-[232px] flex flex-col gap-3 items-center justify-center lg:min-h-[232px]">
            <div
              className={clsx(
                "hidden lg:flex absolute inset-0 items-center justify-center transition-opacity duration-[1100ms] ease-out pointer-events-none",
                isSingleMolecule ? "opacity-100" : "opacity-0",
              )}
            >
              <div className="h-[120%] w-[120%] rounded-full bg-gradient-to-br from-blue-500 via-purple-500 to-fuchsia-500 opacity-30 blur-3xl" />
            </div>

            <div
              className={clsx(
                "w-full transition-all duration-[1100ms] ease-[cubic-bezier(0.22,1,0.36,1)] origin-center lg:absolute lg:left-1/2 lg:top-1/2 lg:-translate-x-1/2 lg:-translate-y-1/2",
                isSingleMolecule
                  ? "hidden lg:block lg:scale-[1.2] lg:opacity-0 lg:blur-sm lg:pointer-events-none"
                  : "block lg:block lg:scale-100 lg:opacity-100 lg:blur-0 lg:pointer-events-auto",
              )}
            >
              <div className="flex flex-col items-center gap-3 mx-auto max-w-2xl w-full">
                <div className="grid w-full grid-cols-1 gap-4 md:grid-cols-2">
                  <FileUpload
                    annotationCsv={annotationCsv}
                    description="Single .h5ad file"
                    icons={<AnnDataIcon className="text-primary" />}
                    serverUpload={serverMode}
                    title="H5AD File"
                    type="h5ad"
                  />
                  <FileUpload
                    annotationCsv={annotationCsv}
                    description="Xenium, MERSCOPE, chunked, or .zarr folder"
                    icons={
                      <>
                        <ZarrIcon />
                        <XeniumIcon />
                        <MerscopeIcon />
                        <ChunkedIcon className="text-primary" />
                      </>
                    }
                    serverUpload={serverMode}
                    title="Folder"
                    type="folder"
                  />
                </div>

                <LoadFromS3Button onClick={() => setIsS3ModalOpen(true)} />
              </div>
            </div>

            <div
              className={clsx(
                "w-full max-w-3xl transition-all duration-[1100ms] ease-[cubic-bezier(0.22,1,0.36,1)] lg:absolute lg:left-1/2 lg:top-1/2 lg:-translate-x-1/2 lg:-translate-y-1/2",
                isSingleMolecule
                  ? "block lg:block lg:opacity-100 lg:scale-100 lg:pointer-events-auto"
                  : "hidden lg:block lg:opacity-0 lg:scale-90 lg:pointer-events-none",
              )}
            >
              <div className="flex flex-col items-center gap-3 mx-auto max-w-2xl w-full">
                <div className="grid w-full grid-cols-1 gap-4 md:grid-cols-2">
                  <FileUpload
                    description=".parquet or .csv (auto-detects schema)"
                    icons={
                      <>
                        <XeniumIcon />
                        <MerscopeIcon />
                      </>
                    }
                    serverUpload={serverMode}
                    singleMolecule={true}
                    title="File"
                    type="file"
                  />
                  <FileUpload
                    description="Pre-processed chunked folder"
                    icons={<ChunkedIcon className="text-primary" />}
                    singleMolecule={true}
                    title="Chunked Folder"
                    type="chunked"
                  />
                </div>

                <LoadFromS3Button onClick={() => setIsS3ModalOpen(true)} />
              </div>
            </div>
          </div>
        </div>

        {featured.length > 0 && (
          <div className="relative z-10 w-full max-w-6xl mx-auto mt-6 px-4">
            <FeaturedDatasets datasets={featured} />
          </div>
        )}

        <div className="relative z-10 flex justify-center mt-3">
          <Button
            as={Link}
            color="primary"
            href="/explore"
            radius="full"
            size="md"
          >
            Explore Example Data
          </Button>
        </div>
      </section>

      <LoadFromS3Modal
        datasetType={isSingleMolecule ? "single_molecule" : "single_cell"}
        isOpen={isS3ModalOpen}
        onClose={() => setIsS3ModalOpen(false)}
      />
    </>
  );
}

export default function Home() {
  return (
    <Suspense
      fallback={
        <div className="flex items-center justify-center min-h-[60vh]">
          <span className="text-sm text-default-400">Loading…</span>
        </div>
      }
    >
      <HomeContent />
    </Suspense>
  );
}
