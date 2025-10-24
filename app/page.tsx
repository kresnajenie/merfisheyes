"use client";

import { Button } from "@heroui/button";
import Link from "next/link";
import { useState, useMemo, memo } from "react";

import { title, subtitle } from "@/components/primitives";
import { FileUpload } from "@/components/file-upload";
import LightRays from "@/components/react-bits/LightRays";
import { Switch } from "@heroui/react";

const MemoizedLightRays = memo(LightRays);

export default function Home() {
  const name = "MERFISH";
  const [isSingleMolecule, setIsSingleMolecule] = useState(false);

  const lightRaysColor = useMemo(
    () => (isSingleMolecule ? "#FF1CF7" : "#5EA2EF"),
    [isSingleMolecule],
  );

  return (
    <>
      <div className="fixed inset-0 w-full h-full z-0">
        <MemoizedLightRays
          lightSpread={1.0}
          mouseInfluence={0.1}
          pulsating={false}
          rayLength={10}
          raysColor={lightRaysColor}
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
                onValueChange={setIsSingleMolecule}
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
          {/* Single Cell File Uploaders */}
          {!isSingleMolecule && (
            <div className="grid grid-cols-3 gap-4 w-full">
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
          )}

          {/* Single Molecule File Uploaders */}
          {isSingleMolecule && (
            <div className="grid grid-cols-2 gap-4 w-full">
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
          )}
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
