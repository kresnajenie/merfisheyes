"use client";

import { useEffect, useRef } from "react";
import * as THREE from "three";
import gsap from "gsap";

import { useTourStore } from "@/lib/stores/tourStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import { resolveStopTiming } from "@/lib/tour/types";

function hslToHex(h: number, s: number, l: number): string {
  const sf = s / 100;
  const lf = l / 100;
  const k = (n: number) => (n + h / 30) % 12;
  const a = sf * Math.min(lf, 1 - lf);
  const f = (n: number) => {
    const c = lf - a * Math.max(-1, Math.min(k(n) - 3, Math.min(9 - k(n), 1)));

    return Math.round(255 * c)
      .toString(16)
      .padStart(2, "0");
  };

  return `#${f(0)}${f(8)}${f(4)}`;
}

// Per-stop genes cycle through these two colors. With most stops at 1 gene
// and the rest at 2, we never need a third — but if a stop has more it
// repeats blue, green, blue, ...
const PER_STOP_GENE_COLORS = [
  hslToHex(233, 96.8, 72.4),
  hslToHex(138.1, 99.1, 43.8),
];
const PER_STOP_GENE_SCALE = 2.0;
// Per-stop genes are only shown while orbiting (close-zoom phase). The
// fly-in lands at `zoomRadius` (default 100µm, well inside the user's
// requested 1000µm threshold), and the zoom-out tween starts immediately
// after orbit ends — so adding at orbit start and removing at orbit end
// satisfies "active when zoomed in, not while zooming out".

type OrbitTicker = (time: number, deltaTime: number) => void;

interface UseTourPlayerOpts {
  cameraRef: React.MutableRefObject<THREE.Camera | null>;
  controlsRef: React.MutableRefObject<any>;
  orbitTickerRef: React.MutableRefObject<OrbitTicker | null>;
  isNormalized: boolean;
  // Distance-from-target cell-rendering filter. The hook drives these for
  // the orbit phase only, then restores the user's pre-tour values.
  targetFilterEnabledRef: React.MutableRefObject<boolean>;
  targetFilterRadiusRef: React.MutableRefObject<number>;
  setTargetFilterEnabled: (enabled: boolean) => void;
  setTargetFilterRadius: (radius: number) => void;
  // SC point cloud — hidden during the orbit phase so only the SM molecules
  // are visible while we sweep around the target.
  scPointCloudRef: React.MutableRefObject<THREE.Points | null>;
}

export function useTourPlayer({
  cameraRef,
  controlsRef,
  orbitTickerRef,
  isNormalized,
  targetFilterEnabledRef,
  targetFilterRadiusRef,
  setTargetFilterEnabled,
  setTargetFilterRadius,
  scPointCloudRef,
}: UseTourPlayerOpts) {
  const isPlaying = useTourStore((s) => s.isPlaying);
  const config = useTourStore((s) => s.config);
  const setProgress = useTourStore((s) => s.setProgress);
  const finish = useTourStore((s) => s.finish);
  const cancelRef = useRef<(() => void) | null>(null);

  useEffect(() => {
    if (!isPlaying || !config) return;

    let cancelled = false;
    const cleanups: Array<() => void> = [];

    const cancel = () => {
      cancelled = true;
      const list = cleanups.splice(0);

      list.forEach((fn) => {
        try {
          fn();
        } catch {
          // best-effort cleanup
        }
      });
    };

    cancelRef.current = cancel;

    const cs = isNormalized ? 500 : 1;

    // Snapshot the user's pre-tour distance-filter settings so we can restore
    // them on finish/cancel. The radius is in world microns for raw-coord
    // datasets and gets multiplied by `cs` for normalized ones so it lands
    // at the same physical distance regardless of coordinate space.
    const originalFilterEnabled = targetFilterEnabledRef.current;
    const originalFilterRadius = targetFilterRadiusRef.current;

    // Snapshot SC visibility so a mid-orbit Stop always restores it.
    const originalScVisible = scPointCloudRef.current?.visible ?? true;

    cleanups.push(() => {
      setTargetFilterEnabled(originalFilterEnabled);
      setTargetFilterRadius(originalFilterRadius);
      if (scPointCloudRef.current) {
        scPointCloudRef.current.visible = originalScVisible;
      }
    });

    // Set of always-on gene names; used to skip them when diffing per-stop
    // genes so the user's manual color/scale edits stick across stops.
    const alwaysOnNames = new Set(
      (config.alwaysOnGenes ?? []).map((g) => g.name),
    );

    const smStore = useSingleMoleculeVisualizationStore.getState();

    // Apply always-on genes once at tour start with their declared colors.
    // They stay visible throughout the tour — only per-stop genes cycle on
    // and off around each orbit.
    for (const ao of config.alwaysOnGenes ?? []) {
      smStore.removeGene(ao.name);
      smStore.addGene(ao.name, ao.color);
    }

    // Tracks the per-stop genes currently added to the store (so we can
    // remove them after each orbit without touching always-on entries).
    let activeStopGenes: string[] = [];

    const addStopGenes = (stopGenes: string[]) => {
      const current =
        useSingleMoleculeVisualizationStore.getState().selectedGenes;

      // Count only genes we actually add so the first one added is always
      // blue (slot 0), the second green (slot 1) — skipped genes don't
      // burn a color slot.
      let colorIdx = 0;

      for (const g of stopGenes) {
        if (alwaysOnNames.has(g) || current.has(g)) continue;
        const color =
          PER_STOP_GENE_COLORS[colorIdx % PER_STOP_GENE_COLORS.length];

        smStore.addGene(g, color, PER_STOP_GENE_SCALE);
        colorIdx++;
      }
      activeStopGenes = stopGenes.slice();
    };

    const removeStopGenes = () => {
      for (const g of activeStopGenes) {
        if (alwaysOnNames.has(g)) continue;
        smStore.removeGene(g);
      }
      activeStopGenes = [];
    };

    // Also drop any stray per-stop genes left over if the tour is cancelled
    // mid-orbit (Stop pressed while genes are still on screen).
    cleanups.push(() => {
      removeStopGenes();
    });

    const wait = (ms: number): Promise<void> =>
      new Promise((resolve) => {
        if (ms <= 0 || cancelled) {
          resolve();

          return;
        }
        const tid = window.setTimeout(() => resolve(), ms);

        cleanups.push(() => {
          window.clearTimeout(tid);
          resolve();
        });
      });

    const tweenCamera = (
      targetEnd: THREE.Vector3,
      posEnd: THREE.Vector3,
      durationMs: number,
    ): Promise<void> =>
      new Promise((resolve) => {
        const cam = cameraRef.current;
        const ctrls = controlsRef.current;

        if (!cam || !ctrls) {
          resolve();

          return;
        }
        const targetStart = ctrls.target.clone();
        const posStart = cam.position.clone();
        const tweenObj = { t: 0 };
        const tween = gsap.to(tweenObj, {
          t: 1,
          duration: Math.max(0.001, durationMs / 1000),
          ease: "power2.inOut",
          overwrite: true,
          onUpdate: () => {
            if (cancelled) return;
            ctrls.target.lerpVectors(targetStart, targetEnd, tweenObj.t);
            cam.position.lerpVectors(posStart, posEnd, tweenObj.t);
            ctrls.update();
          },
          onComplete: resolve,
        });

        cleanups.push(() => {
          tween.kill();
          resolve();
        });
      });

    // Combined zoom-out + pan: arcs the camera from `targetStart` to
    // `targetEnd` while pushing the camera distance up to `peakRadius` at
    // the midpoint, then back down to `endRadius`. Direction (which way
    // the camera looks at the target) is captured at tween start and
    // preserved throughout. Replaces the old "zoom-out, then fly-in"
    // pair so the two motions happen simultaneously.
    const archedTransition = (
      targetStart: THREE.Vector3,
      targetEnd: THREE.Vector3,
      peakRadius: number,
      endRadius: number,
      durationMs: number,
    ): Promise<void> =>
      new Promise((resolve) => {
        const cam = cameraRef.current;
        const ctrls = controlsRef.current;

        if (!cam || !ctrls) {
          resolve();

          return;
        }
        const posStart = cam.position.clone();
        const dir = posStart.clone().sub(targetStart).normalize();

        if (!Number.isFinite(dir.x) || dir.lengthSq() < 1e-6) {
          dir.set(0, 0, 1);
        }
        const startRadius = posStart.distanceTo(targetStart);
        // Only "lift" the camera above the start/end radii — if peakRadius
        // is already lower than them, the arc is flat (no upward bump).
        const peakBoost = Math.max(
          0,
          peakRadius - Math.max(startRadius, endRadius),
        );
        const tweenObj = { t: 0 };
        const tween = gsap.to(tweenObj, {
          t: 1,
          duration: Math.max(0.001, durationMs / 1000),
          ease: "power2.inOut",
          overwrite: true,
          onUpdate: () => {
            if (cancelled) return;
            const targetNow = targetStart
              .clone()
              .lerp(targetEnd, tweenObj.t);
            const baseRadius =
              startRadius + (endRadius - startRadius) * tweenObj.t;
            const arc = peakBoost * Math.sin(tweenObj.t * Math.PI);
            const distance = baseRadius + arc;

            cam.position
              .copy(targetNow)
              .add(dir.clone().multiplyScalar(distance));
            ctrls.target.copy(targetNow);
            ctrls.update();
          },
          onComplete: resolve,
        });

        cleanups.push(() => {
          tween.kill();
          resolve();
        });
      });

    const runOrbit = (durationMs: number, revolutions: number): Promise<void> =>
      new Promise((resolve) => {
        const cam = cameraRef.current;
        const ctrls = controlsRef.current;

        if (!cam || !ctrls || durationMs <= 0 || revolutions <= 0) {
          resolve();

          return;
        }
        const target = ctrls.target.clone();
        const startRadius = cam.position.clone().sub(target);
        const radiusXZ = Math.sqrt(
          startRadius.x * startRadius.x + startRadius.z * startRadius.z,
        );
        const baseY = startRadius.y;
        const periodMs = durationMs / revolutions;
        let angle = 0;
        let elev = 0;
        let finished = false;

        const finishOrbit = () => {
          if (finished) return;
          finished = true;
          if (orbitTickerRef.current === ticker) {
            orbitTickerRef.current = null;
          }
          gsap.ticker.remove(ticker);
          resolve();
        };

        const ticker: OrbitTicker = (_time, deltaTime) => {
          if (cancelled || orbitTickerRef.current !== ticker) {
            finishOrbit();

            return;
          }
          angle += (deltaTime / periodMs) * Math.PI * 2;
          elev += (deltaTime / (periodMs * 2)) * Math.PI * 2;
          const sin = Math.sin(angle);
          const cos = Math.cos(angle);
          const yOsc = Math.sin(elev) * radiusXZ * 0.3;

          cam.position.set(
            target.x + startRadius.x * cos + startRadius.z * sin,
            target.y + baseY + yOsc,
            target.z + -startRadius.x * sin + startRadius.z * cos,
          );
          ctrls.update();
        };

        // Replace any orbit currently in flight (e.g. user had pressed 'o').
        if (orbitTickerRef.current) {
          gsap.ticker.remove(orbitTickerRef.current);
        }
        orbitTickerRef.current = ticker;
        gsap.ticker.add(ticker);

        const tid = window.setTimeout(finishOrbit, durationMs);

        cleanups.push(() => {
          window.clearTimeout(tid);
          finishOrbit();
        });
      });

    const runTour = async () => {
      const cam = cameraRef.current;
      const ctrls = controlsRef.current;

      if (!cam || !ctrls) {
        finish();

        return;
      }

      // Snapshot the camera at tour start so we can fly back to it after
      // the last stop's zoom-out completes.
      const homeTarget = ctrls.target.clone();
      const homePosition = cam.position.clone();

      for (let i = 0; i < config.stops.length; i++) {
        if (cancelled) return;

        const stop = config.stops[i];
        const timing = resolveStopTiming(stop, config);

        setProgress(i, "fly");

        const targetWorld = new THREE.Vector3(
          stop.x * cs,
          stop.y * cs,
          stop.z * cs,
        );

        // Preserve current view direction so the camera flies in from
        // wherever it currently is rather than snapping to a fixed angle.
        const dirIn = cameraRef
          .current!.position.clone()
          .sub(controlsRef.current.target)
          .normalize();

        if (!Number.isFinite(dirIn.x) || dirIn.lengthSq() < 1e-6) {
          dirIn.set(0, 0, 1);
        }
        const flyEnd = targetWorld
          .clone()
          .add(dirIn.clone().multiplyScalar(timing.zoomRadius));

        await tweenCamera(targetWorld, flyEnd, timing.flyInMs);
        if (cancelled) return;

        // Wait phase: add per-stop genes (they show immediately) and hold
        // the camera still for waitMs so the renderer can settle before
        // the orbit starts. Persistent genes are already visible.
        setProgress(i, "wait");
        addStopGenes(stop.genes);
        await wait(timing.waitMs);
        if (cancelled) return;

        setProgress(i, "orbit");
        setTargetFilterRadius(timing.filterRadius * cs);
        setTargetFilterEnabled(true);
        if (scPointCloudRef.current) {
          scPointCloudRef.current.visible = false;
        }
        await runOrbit(timing.orbitMs, timing.orbitRevolutions);
        removeStopGenes();
        if (scPointCloudRef.current) {
          scPointCloudRef.current.visible = originalScVisible;
        }
        setTargetFilterEnabled(originalFilterEnabled);
        setTargetFilterRadius(originalFilterRadius);
        if (cancelled) return;

        // Post-orbit wait: per-stop genes have just been deselected. Camera
        // holds still for waitAfterMs before the zoom-out tween begins so
        // the change reads as a deliberate beat instead of a snap.
        setProgress(i, "wait");
        await wait(timing.waitAfterMs);
        if (cancelled) return;

        // Skip the zoom-out for the final stop — the post-loop home tween
        // already flies the camera straight back to where it started, so a
        // separate zoom-out first would just pull us partway out and then
        // make a second move.
        const isLastStop = i === config.stops.length - 1;

        if (!isLastStop) {
          setProgress(i, "zoom-out");
          const dirOut = cameraRef
            .current!.position.clone()
            .sub(controlsRef.current.target)
            .normalize();

          if (!Number.isFinite(dirOut.x) || dirOut.lengthSq() < 1e-6) {
            dirOut.set(0, 0, 1);
          }
          const zoomOutEnd = targetWorld
            .clone()
            .add(dirOut.clone().multiplyScalar(timing.zoomOutRadius * cs));

          await tweenCamera(targetWorld, zoomOutEnd, timing.zoomOutMs);
          if (cancelled) return;
        }
      }

      // Tour complete — fly the camera back to where the user started.
      // Reuses the last stop's flyInMs for the return tween duration.
      if (!cancelled && config.stops.length > 0) {
        const lastTiming = resolveStopTiming(
          config.stops[config.stops.length - 1],
          config,
        );

        setProgress(config.stops.length - 1, "fly");
        await tweenCamera(homeTarget, homePosition, lastTiming.flyInMs);
        if (cancelled) return;
      }

      if (!cancelled) finish();
    };

    runTour();

    return () => {
      cancelRef.current = null;
      cancel();
    };
  }, [
    isPlaying,
    config,
    isNormalized,
    cameraRef,
    controlsRef,
    orbitTickerRef,
    targetFilterEnabledRef,
    targetFilterRadiusRef,
    setTargetFilterEnabled,
    setTargetFilterRadius,
    scPointCloudRef,
    setProgress,
    finish,
  ]);
}
