"use client";

import LightRays from "@/components/react-bits/LightRays";

export function ExploreBackground() {
  return (
    <div className="fixed inset-0 w-full h-full z-0">
      <LightRays
        lightSpread={1.2}
        mouseInfluence={0.1}
        pulsating={false}
        rayLength={10}
        raysColor="#FF1CF7"
        raysOrigin="top-left"
        raysSpeed={1.0}
      />
    </div>
  );
}
