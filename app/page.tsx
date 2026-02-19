"use client";

import { Button } from "@heroui/button";
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
import LightRays from "@/components/react-bits/LightRays";
import { BrainToggle } from "@/components/brain-toggle";
import { LoadFromS3Modal } from "@/components/load-from-s3-modal";

const MemoizedLightRays = memo(LightRays);

function useModeToggleState() {
  const router = useRouter();
  const pathname = usePathname();
  const searchParams = useSearchParams();
  const modeParam = searchParams.get("mode");

  const initialMode = useMemo(() => {
    if (modeParam) {
      return modeParam === "sm";
    }

    if (typeof window !== "undefined") {
      const storedMode = window.localStorage.getItem("lastDatasetMode");

      if (storedMode) {
        return storedMode === "sm";
      }
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
      <section className="relative flex flex-col items-center gap-1 py-12 md:py-20 px-4 md:px-8">
        <div className="relative z-10 flex flex-col items-center gap-10 max-w-3xl w-full">
          <div className="flex flex-col items-center">
            <h1 className="flex flex-col items-center text-center">
              <span className={title({ size: "xl" })}>
                Bring{" "}
                <span
                  className={title({
                    color: isSingleMolecule ? "violet" : "blue",
                    size: "xl",
                  })}
                  style={{ transition: "color 0.3s ease" }}
                >
                  {name}
                </span>{" "}
                to Life
              </span>
            </h1>
            <div className="flex items-center gap-6 mt-6 justify-center w-full">
              <button
                aria-pressed={!isSingleMolecule}
                className={clsx(
                  "px-4 py-2 rounded-full text-xs md:text-sm font-semibold tracking-[0.28em] uppercase transition-colors duration-300 cursor-pointer focus:outline-none focus-visible:ring-2 focus-visible:ring-offset-2 focus-visible:ring-blue-400/80 focus-visible:ring-offset-slate-900",
                  !isSingleMolecule
                    ? "bg-blue-500 text-white shadow-[0_12px_25px_rgba(59,130,246,0.35)]"
                    : "bg-transparent text-default-500 border border-blue-400/30 hover:bg-blue-500/15 hover:text-white",
                )}
                type="button"
                onClick={() => handleModeChange(false)}
              >
                single cell
              </button>
              <BrainToggle
                isActive={isSingleMolecule}
                onToggle={handleModeChange}
              />
              <button
                aria-pressed={isSingleMolecule}
                className={clsx(
                  "px-4 py-2 rounded-full text-xs md:text-sm font-semibold tracking-[0.28em] uppercase transition-colors duration-300 cursor-pointer focus:outline-none focus-visible:ring-2 focus-visible:ring-offset-2 focus-visible:ring-purple-400/80 focus-visible:ring-offset-slate-900",
                  isSingleMolecule
                    ? "bg-purple-500 text-white shadow-[0_12px_25px_rgba(168,85,247,0.35)]"
                    : "bg-transparent text-default-500 border border-purple-400/30 hover:bg-purple-500/15 hover:text-white",
                )}
                type="button"
                onClick={() => handleModeChange(true)}
              >
                single molecule
              </button>
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
          </div>
        </div>

        <div className="relative z-10 flex flex-col items-center gap-8 max-w-5xl w-full">
          <div className="relative w-full max-w-5xl min-h-[360px] flex flex-col gap-6 items-center justify-center lg:min-h-[420px]">
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
              <div className="grid grid-cols-1 gap-4 md:grid-cols-2 lg:grid-cols-5">
                <FileUpload
                  description="Single .h5ad file"
                  title="H5AD File"
                  type="h5ad"
                />
                <FileUpload
                  description="Pre-processed chunked folder"
                  title="Chunked Folder"
                  type="chunked"
                />
                <FileUpload
                  description="Select Xenium output folder"
                  title="Xenium Folder"
                  type="xenium"
                />
                <FileUpload
                  description="Select Merscope output folder"
                  title="Merscope Folder"
                  type="merscope"
                />
                <div className="w-full">
                  <button
                    className="relative w-full border-2 border-dashed rounded-lg p-6 text-center cursor-pointer transition-all duration-200 ease-in-out aspect-square flex items-center justify-center border-default-300 hover:border-primary/50 hover:bg-default-100/50 focus:outline-none focus-visible:ring-2 focus-visible:ring-primary/50"
                    type="button"
                    onClick={() => setIsS3ModalOpen(true)}
                  >
                    <div className="flex flex-col items-center gap-2 w-full">
                      <p className="text-lg font-semibold text-foreground">
                        Load from S3
                      </p>
                      <p className="text-xs text-default-500">
                        Your own S3 bucket
                      </p>
                    </div>
                  </button>
                </div>
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
              <div className="grid grid-cols-1 gap-4 md:grid-cols-2 lg:grid-cols-4">
                <FileUpload
                  description="Pre-processed chunked folder"
                  singleMolecule={true}
                  title="Chunked Folder"
                  type="chunked"
                />
                <FileUpload
                  description="Select .parquet or .csv file"
                  singleMolecule={true}
                  title="Xenium Parquet/CSV"
                  type="xenium"
                />
                <FileUpload
                  description="Select .parquet or .csv file"
                  singleMolecule={true}
                  title="MERSCOPE Parquet/CSV"
                  type="merscope"
                />
                <div className="w-full">
                  <button
                    className="relative w-full border-2 border-dashed rounded-lg p-6 text-center cursor-pointer transition-all duration-200 ease-in-out aspect-square flex items-center justify-center border-default-300 hover:border-primary/50 hover:bg-default-100/50 focus:outline-none focus-visible:ring-2 focus-visible:ring-primary/50"
                    type="button"
                    onClick={() => setIsS3ModalOpen(true)}
                  >
                    <div className="flex flex-col items-center gap-2 w-full">
                      <p className="text-lg font-semibold text-foreground">
                        Load from S3
                      </p>
                      <p className="text-xs text-default-500">
                        Your own S3 bucket
                      </p>
                    </div>
                  </button>
                </div>
              </div>
            </div>
          </div>
        </div>

        <div className="relative z-10 flex justify-center mt-8">
          <Button
            as={Link}
            color="primary"
            href="/explore"
            radius="full"
            size="lg"
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
    <Suspense fallback={<div>Loading...</div>}>
      <HomeContent />
    </Suspense>
  );
}
