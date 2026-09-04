"use client";

import { glassPanel } from "@/components/primitives";

/**
 * Shown when THREE.WebGLRenderer fails to create a WebGL context (GPU
 * blocklisted, stale browser update, software-only GL, …). Without this the
 * viewer loads data fine but leaves a silently blank canvas.
 */
export function WebGLErrorOverlay() {
  return (
    <div className="absolute inset-0 z-50 flex items-center justify-center p-6">
      <div className={`${glassPanel()} max-w-md p-6 text-center`}>
        <h2 className="text-lg font-semibold">3D view unavailable</h2>
        <p className="mt-2 text-sm text-default-500">
          Your browser couldn&apos;t create a WebGL context, so the scene
          can&apos;t render. The dataset itself loaded fine.
        </p>
        <ul className="mt-3 list-inside list-disc text-left text-sm text-default-500">
          <li>
            Fully quit and reopen your browser (especially after an update)
          </li>
          <li>
            Check <span className="font-mono">chrome://gpu</span> — WebGL should
            say &quot;Hardware accelerated&quot;
          </li>
          <li>
            Enable &quot;Use graphics acceleration&quot; in your browser&apos;s
            system settings
          </li>
        </ul>
      </div>
    </div>
  );
}
