"use client";

import { title, subtitle } from "@/components/primitives";
import LightRays from "@/components/react-bits/LightRays";

export default function AboutPage() {
  return (
    <div>
      <div className="fixed inset-0 w-full h-full z-0">
        <LightRays
          raysOrigin="bottom-center"
          raysColor="#6FEE8D"
          rayLength={10}
          raysSpeed={1.0}
          lightSpread={1.0}
          pulsating={false}
          mouseInfluence={0.1}
        />
      </div>
      <div className="relative z-10 flex flex-col items-center gap-1 py-12 md:py-20 px-4 md:px-8">
        <h1 className={title({ size: "md" })}>
          <span className={title({ color: "green" })}>About</span> Us
        </h1>
        <div className={subtitle({ class: "mt-4 text-center" })}>
          We are a team of researchers and developers dedicated to advancing the
          future of spatial transcriptomics.
        </div>
      </div>
    </div>
  );
}
