"use client";

import { Button } from "@heroui/button";
import Link from "next/link";

import { title, subtitle } from "@/components/primitives";
import { FileUpload } from "@/components/file-upload";
import LightRays from "@/components/react-bits/LightRays";

export default function Home() {
  const name = "MERFISH";

  return (
    <>
      <div className="fixed inset-0 w-full h-full z-0">
        <LightRays
          lightSpread={1.0}
          mouseInfluence={0.1}
          pulsating={false}
          rayLength={10}
          raysColor="#5EA2EF"
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
                <span className={title({ color: "blue", size: "xl" })}>
                  {name}
                </span>{" "}
                to Life
              </span>
            </h1>
            <div className={subtitle({ class: "mt-4 text-center" })}>
              Explore your{" "}
              <span className="font-bold">.h5ad , Merscope & Xenium</span>{" "}
              datasets
            </div>
          </div>
        </div>

        <div className="relative z-10 flex flex-col items-center gap-8 max-w-5xl w-full">
          {/* Three separate file uploaders for testing */}
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
