"use client";

import { Button } from "@heroui/button";
import Link from "next/link";
import {
  useEffect,
  useMemo,
  memo,
  useRef,
  useState,
} from "react";
import clsx from "clsx";
import {
  usePathname,
  useRouter,
  useSearchParams,
} from "next/navigation";

import { title, subtitle } from "@/components/primitives";
import { FileUpload } from "@/components/file-upload";
import LightRays from "@/components/react-bits/LightRays";
import { Switch } from "@heroui/react";

const MemoizedLightRays = memo(LightRays);

export default function Home() {
  const name = "MERFISH";
  const router = useRouter();
  const pathname = usePathname();
  const searchParams = useSearchParams();
  const [isSingleMolecule, setIsSingleMolecule] = useState(
    () => searchParams.get("mode") === "sm",
  );
  const animationFrameRef = useRef<number | null>(null);
  const currentColorRef = useRef("#5EA2EF");
  const [animatedRaysColor, setAnimatedRaysColor] = useState("#5EA2EF");
  const modeParam = searchParams.get("mode");

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
    const shouldBeSingleMolecule = modeParam === "sm";

    setIsSingleMolecule((prev) =>
      prev === shouldBeSingleMolecule ? prev : shouldBeSingleMolecule,
    );
  }, [modeParam]);

  const handleModeChange = (selected: boolean) => {
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
  };

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
            <div className="flex items-center gap-3 mt-4">
              <span
                className={`text-sm transition-all duration-300 px-3 py-1 rounded-full ${
                  !isSingleMolecule
                    ? "font-medium text-white bg-blue-500"
                    : "font-normal text-default-500"
                }`}
              >
                single cell
              </span>
              <Switch
                isSelected={isSingleMolecule}
                onValueChange={handleModeChange}
                aria-label="Toggle between single cell and single molecule"
                classNames={{
                  wrapper: isSingleMolecule
                    ? "bg-purple-500 group-data-[selected=true]:bg-purple-500"
                    : "bg-blue-500",
                  thumb: "bg-white",
                }}
              />
              <span
                className={`text-sm transition-all duration-300 px-3 py-1 rounded-full ${
                  isSingleMolecule
                    ? "font-medium text-white bg-purple-500"
                    : "font-normal text-default-500"
                }`}
              >
                single molecule
              </span>
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
          <div className="relative w-full max-w-5xl min-h-[360px] flex items-center justify-center">
            <div
              className={clsx(
                "absolute inset-0 flex items-center justify-center transition-opacity duration-[1100ms] ease-out pointer-events-none",
                isSingleMolecule ? "opacity-100" : "opacity-0",
              )}
            >
              <div className="h-[120%] w-[120%] rounded-full bg-gradient-to-br from-blue-500 via-purple-500 to-fuchsia-500 opacity-30 blur-3xl" />
            </div>

            <div
              className={clsx(
                "absolute left-1/2 top-1/2 w-full -translate-x-1/2 -translate-y-1/2 transition-all duration-[1100ms] ease-[cubic-bezier(0.22,1,0.36,1)] origin-center",
                isSingleMolecule
                  ? "scale-[1.2] opacity-0 blur-sm pointer-events-none"
                  : "scale-100 opacity-100 blur-0 pointer-events-auto",
              )}
            >
              <div className="grid grid-cols-1 gap-4 md:grid-cols-2 lg:grid-cols-4">
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
              </div>
            </div>

            <div
              className={clsx(
                "absolute left-1/2 top-1/2 w-full max-w-3xl -translate-x-1/2 -translate-y-1/2 transition-all duration-[1100ms] ease-[cubic-bezier(0.22,1,0.36,1)]",
                isSingleMolecule
                  ? "opacity-100 scale-100 pointer-events-auto"
                  : "opacity-0 scale-90 pointer-events-none",
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
    </>
  );
}
