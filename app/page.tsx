"use client";

import { title, subtitle } from "@/components/primitives";
import { FileUpload } from "@/components/file-upload";
import LightRays from "@/components/react-bits/LightRays";
import { Button } from "@heroui/button";
import Link from "next/link";

export default function Home() {
  const name = "MERFISH";
  return (
    <>
      <div className="fixed inset-0 w-full h-full z-0">
        <LightRays
          raysOrigin="top-center"
          raysColor="#5EA2EF"
          rayLength={10}
          raysSpeed={1.0}
          lightSpread={1.0}
          pulsating={false}
          mouseInfluence={0.1}
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
              <span className="font-bold">.h5ad , MERSCOPE & XENIUM</span>{" "}
              datasets
            </div>
          </div>
        </div>

        <div className="relative z-10 flex flex-col items-center gap-8 max-w-5xl w-full">
          {/* Three separate file uploaders for testing */}
          <div className="grid grid-cols-3 gap-4 w-full">
            <FileUpload
              type="h5ad"
              title="H5AD File"
              description="Single .h5ad file"
            />
            <FileUpload
              type="xenium"
              title="Xenium Folder"
              description="Select Xenium output folder"
            />
            <FileUpload
              type="merscope"
              title="MERSCOPE Folder"
              description="Select MERSCOPE output folder"
            />
          </div>
        </div>

        <div className="relative z-10 flex justify-center mt-8">
          <Button
            as={Link}
            href="/explore"
            color="primary"
            size="lg"
            radius="full"
          >
            Explore Example Data
          </Button>
        </div>
      </section>
    </>
  );
}
