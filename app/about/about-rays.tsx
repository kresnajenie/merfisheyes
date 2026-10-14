"use client";

import LightRays from "@/components/react-bits/LightRays";

/** The page's animated background; it runs in the browser only. */
export function AboutRays() {
  return (
    <LightRays
      lightSpread={1.0}
      mouseInfluence={0.1}
      pulsating={false}
      rayLength={10}
      raysColor="#6FEE8D"
      raysOrigin="bottom-center"
      raysSpeed={1.0}
    />
  );
}
