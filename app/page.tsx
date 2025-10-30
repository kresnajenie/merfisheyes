"use client";

import { Button } from "@heroui/button";
import Link from "next/link";
import { memo, useEffect, useMemo, useRef, useState } from "react";
import clsx from "clsx";

import { BrainToggle } from "@/components/brain-toggle";
import { FileUpload } from "@/components/file-upload";
import { title, subtitle } from "@/components/primitives";
import LightRays from "@/components/react-bits/LightRays";

const MemoizedLightRays = memo(LightRays);

export default function Home() {
  const name = "MERFISH";
  const [isSingleMolecule, setIsSingleMolecule] = useState(false);
  const animationFrameRef = useRef<number | null>(null);
  const currentColorRef = useRef("#5EA2EF");
  const [animatedRaysColor, setAnimatedRaysColor] = useState("#5EA2EF");
  const [viewportHeight, setViewportHeight] = useState<number | null>(
    typeof window === "undefined" ? null : window.innerHeight,
  );
  const [viewportWidth, setViewportWidth] = useState<number | null>(
    typeof window === "undefined" ? null : window.innerWidth,
  );

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

  useEffect(() => {
    const updateViewport = () => {
      if (typeof window !== "undefined") {
        setViewportHeight(window.innerHeight);
        setViewportWidth(window.innerWidth);
      }
    };

    updateViewport();
    window.addEventListener("resize", updateViewport);
    window.addEventListener("orientationchange", updateViewport);

    return () => {
      window.removeEventListener("resize", updateViewport);
      window.removeEventListener("orientationchange", updateViewport);
    };
  }, []);

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
      <section
        className="relative flex flex-col items-center gap-1 px-4 md:px-8"
        style={{ minHeight: viewportHeight ? `${viewportHeight}px` : "100vh" }}
      >
        <div className="relative z-10 flex flex-col items-center gap-10 max-w-3xl w-full pt-12 md:pt-16">
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
            <div className="flex flex-col items-center gap-3 mt-4 sm:flex-row sm:gap-6">
              <button
                type="button"
                onClick={() => setIsSingleMolecule(false)}
                className={`text-sm transition-all duration-300 px-3 py-1 rounded-full cursor-pointer hover:scale-[1.04] hover:shadow-lg hover:shadow-blue-500/20 focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-offset-2 focus-visible:ring-blue-500 ${
                  !isSingleMolecule
                    ? "font-medium text-white bg-blue-500"
                    : "font-normal text-default-500"
                }`}
              >
                single cell
              </button>
              <BrainToggle
                isActive={isSingleMolecule}
                onToggle={setIsSingleMolecule}
              />
              <button
                type="button"
                onClick={() => setIsSingleMolecule(true)}
                className={`text-sm transition-all duration-300 px-3 py-1 rounded-full cursor-pointer hover:scale-[1.04] hover:shadow-lg hover:shadow-purple-500/20 focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-offset-2 focus-visible:ring-purple-500 ${
                  isSingleMolecule
                    ? "font-medium text-white bg-purple-500"
                    : "font-normal text-default-500"
                }`}
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

        <div className="relative z-10 flex flex-col items-center gap-8 max-w-5xl w-full pb-12">
          <div className="relative w-full max-w-5xl flex items-center justify-center">
            <div
              className={clsx(
                "absolute inset-0 flex items-center justify-center transition-opacity duration-[1100ms] ease-out pointer-events-none",
                isSingleMolecule ? "opacity-100" : "opacity-0",
              )}
            >
              <div className="h-[120%] w-[120%] rounded-full bg-gradient-to-br from-blue-500 via-purple-500 to-fuchsia-500 opacity-30 blur-3xl" />
            </div>
          </div>

          <div className="relative w-full max-w-2xl min-h-[320px] md:min-h-[280px]">
            <div
              className={clsx(
                "absolute inset-0 transition-all duration-[1100ms] ease-[cubic-bezier(0.22,1,0.36,1)]",
                isSingleMolecule
                  ? "opacity-0 pointer-events-none translate-y-6 scale-95"
                  : "opacity-100 pointer-events-auto translate-y-0 scale-100",
              )}
            >
              <div className="grid grid-cols-1 gap-4 lg:grid-cols-3">
                <FileUpload
                  description="Single .h5ad file"
                  title="H5AD File"
                  type="h5ad"
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
              </div>
            </div>

            <div
              className={clsx(
                "absolute inset-0 transition-all duration-[1100ms] ease-[cubic-bezier(0.22,1,0.36,1)]",
                isSingleMolecule
                  ? "opacity-100 pointer-events-auto translate-y-0 scale-100"
                  : "opacity-0 pointer-events-none -translate-y-6 scale-95",
              )}
            >
              <div className="grid grid-cols-1 gap-4 md:grid-cols-2">
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
              </div>
            </div>
          </div>
        </div>

        <div
          className={clsx(
            "relative z-10 flex justify-center transition-all duration-500",
            viewportWidth && viewportWidth < 500
              ? "mt-[800px] mb-16"
              : "mt-8 mb-12",
          )}
        >
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
    </>
  );
}
