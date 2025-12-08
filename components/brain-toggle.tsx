"use client";

import { AnimatePresence, motion } from "framer-motion";
import clsx from "clsx";
import { useEffect, useMemo, useState } from "react";

interface BrainToggleProps {
  isActive: boolean;
  onToggle: (nextValue: boolean) => void;
  className?: string;
}

interface Point2D {
  x: number;
  y: number;
}

interface BrainPoint {
  id: string;
  home: Point2D;
  scatter: Point2D;
  size: number;
  delay: number;
}

interface MoleculePoint {
  id: string;
  position: Point2D;
  size: number;
  delay: number;
}

const HOME_SCALE = 46;
const SCATTER_SCALE = 220;
const INNER_SCALE = 44;
const CENTER_OFFSET_X = 0;
const CENTER_OFFSET_Y = 0;

const ease = [0.22, 1, 0.36, 1] as const;

const seeded = (seed: number) => {
  const x = Math.sin(seed) * 43758.5453123;

  return x - Math.floor(x);
};

const generateBrainPoints = (count: number): BrainPoint[] => {
  const half = Math.floor(count / 2);
  const points: BrainPoint[] = [];

  for (let i = 0; i < count; i++) {
    const hemisphere = i < half ? -1 : 1;
    const hemiIndex = i < half ? i : i - half;
    const ratio = hemiIndex / half;
    const angle = ratio * Math.PI * 2;
    const wobble = Math.sin(angle * 3 + hemisphere * 0.6) * 0.08;
    const radius = 0.72 + wobble;

    const baseX =
      hemisphere * 0.28 +
      Math.cos(angle) * radius * (hemisphere === -1 ? 0.42 : 0.46);
    const baseY = Math.sin(angle) * radius * 0.58;

    const scatterX = (seeded(i * 19 + 11) - 0.5) * 6.2;
    const scatterY = (seeded(i * 23 + 7) - 0.5) * 5.8;

    const size = 5 + seeded(i * 17 + 31) * 9;
    const delay = seeded(i * 13 + 29) * 0.22;

    points.push({
      id: `brain-${i}`,
      home: { x: baseX, y: baseY },
      scatter: { x: scatterX, y: scatterY },
      size,
      delay,
    });
  }

  return points;
};

const generateMoleculePoints = (count: number): MoleculePoint[] => {
  const points: MoleculePoint[] = [];
  const steps = Math.max(12, count);
  const turns = 3.4;
  const radius = 0.36;
  const height = 1.6;

  for (let i = 0; i < steps; i++) {
    const t = steps <= 1 ? 0 : i / (steps - 1);
    const angle = t * turns * Math.PI * 2;
    const y = (t - 0.5) * height;

    const x = Math.cos(angle) * radius;
    const emphasis = 1 - Math.abs(t - 0.5);
    const size = 4.1 + emphasis * 3.4 + seeded(i * 73 + 19) * 1.1;
    const delay = 0.18 + t * 0.48 + seeded(i * 89 + 31) * 0.15;

    points.push({
      id: `molecule-helix-${i}`,
      position: { x, y },
      size,
      delay,
    });
  }

  return points.slice(0, count);
};

export function BrainToggle({
  isActive,
  onToggle,
  className,
}: BrainToggleProps) {
  const points = useMemo(() => generateBrainPoints(120), []);
  const molecules = useMemo(() => generateMoleculePoints(42), []);
  const [isHovering, setIsHovering] = useState(false);

  useEffect(() => {
    if (!isActive) {
      setIsHovering(false);
    }
  }, [isActive]);

  return (
    <button
      aria-label={
        isActive
          ? "Switch to single cell datasets"
          : "Switch to single molecule datasets"
      }
      aria-pressed={isActive}
      className={clsx(
        "group relative flex items-center justify-center rounded-full p-2 transition-shadow cursor-pointer focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-offset-2 focus-visible:ring-slate-100/70 focus-visible:ring-offset-slate-950/40",
        className,
      )}
      type="button"
      onClick={() => onToggle(!isActive)}
      onMouseEnter={() => setIsHovering(true)}
      onMouseLeave={() => setIsHovering(false)}
    >
      <div className="absolute inset-0 rounded-full bg-gradient-to-br from-slate-900/10 via-slate-500/5 to-slate-50/5 blur-xl" />
      <motion.div
        animate={{
          boxShadow: isActive
            ? "0 30px 60px rgba(168, 85, 247, 0.35)"
            : "0 22px 50px rgba(14, 165, 233, 0.25)",
          borderColor: isActive
            ? "rgba(124, 58, 237, 0.35)"
            : "rgba(59, 130, 246, 0.3)",
        }}
        className="relative flex h-28 w-28 items-center justify-center rounded-full border border-white/10 bg-slate-900/30 backdrop-blur-xl shadow-[0_25px_45px_rgba(15,23,42,0.25)] overflow-visible"
        initial={false}
        transition={{ duration: 1.1, ease }}
        whileFocus={
          isActive
            ? {
                scale: 1.04,
                boxShadow: "0 30px 65px rgba(168, 85, 247, 0.4)",
              }
            : {
                scale: 1.04,
                boxShadow: "0 30px 65px rgba(14, 165, 233, 0.35)",
              }
        }
        whileHover={
          isActive
            ? {
                scale: 1.04,
                boxShadow: "0 30px 65px rgba(168, 85, 247, 0.4)",
              }
            : {
                scale: 1.04,
                boxShadow: "0 30px 65px rgba(14, 165, 233, 0.35)",
              }
        }
        whileTap={{ scale: 0.98 }}
      >
        <motion.div
          animate={{
            opacity: isActive ? 0.45 : 0.25,
            scale: isActive ? 1.12 : 1,
          }}
          className="absolute inset-0 rounded-full bg-gradient-to-br from-sky-500/20 via-purple-400/15 to-fuchsia-500/10"
          initial={false}
          transition={{ duration: 1.1, ease }}
        />

        <motion.div
          animate={{
            borderColor: isActive
              ? "rgba(168, 85, 247, 0.35)"
              : "rgba(59, 130, 246, 0.3)",
            boxShadow: isActive
              ? [
                  "0 0 0 0 rgba(168, 85, 247, 0.4)",
                  "0 0 0 16px rgba(168, 85, 247, 0)",
                ]
              : [
                  "0 0 0 0 rgba(59, 130, 246, 0.35)",
                  "0 0 0 16px rgba(59, 130, 246, 0)",
                ],
          }}
          className="absolute inset-0 rounded-full border-2 border-transparent"
          initial={false}
          transition={{
            repeat: Infinity,
            repeatDelay: 1.6,
            duration: 1.8,
            ease: "easeOut",
          }}
        />

        {points.map((point) => (
          <motion.span
            key={point.id}
            animate={
              isActive
                ? {
                    x: point.scatter.x * SCATTER_SCALE,
                    y: point.scatter.y * SCATTER_SCALE,
                    opacity: 0.55,
                  }
                : {
                    x: point.home.x * HOME_SCALE + CENTER_OFFSET_X,
                    y: point.home.y * HOME_SCALE + CENTER_OFFSET_Y,
                    opacity: 0.92,
                  }
            }
            className="absolute rounded-full bg-white shadow-[0_0_15px_rgba(255,255,255,0.35)]"
            initial={false}
            style={{ width: point.size, height: point.size }}
            transition={{
              duration: 1.1,
              ease,
              delay: isActive ? point.delay * 0.8 : point.delay * 0.5,
            }}
          />
        ))}

        <motion.span
          animate={
            isActive
              ? {
                  scale: isHovering ? 3.4 : 3.1,
                  opacity: isHovering ? 0.16 : 0.12,
                }
              : { scale: 1.15, opacity: 0.28 }
          }
          className="absolute left-1/2 top-1/2 -translate-x-1/2 -translate-y-1/2 rounded-full bg-gradient-to-br from-sky-300/40 via-violet-400/40 to-fuchsia-500/40"
          initial={false}
          style={{ width: 48, height: 48 }}
          transition={{
            duration: isHovering ? 0.4 : 1.1,
            ease,
          }}
        />

        <AnimatePresence>
          {isActive && (
            <motion.div
              key="inner-molecules"
              animate={{
                opacity: 1,
                scale: isHovering ? 1.08 : 1,
              }}
              className="absolute left-1/2 top-1/2 -translate-x-1/2 -translate-y-1/2"
              exit={{ opacity: 0, scale: 0.3 }}
              initial={{ opacity: 0, scale: 0.35 }}
              transition={{
                duration: isHovering ? 0.45 : 0.9,
                ease,
                delay: isHovering ? 0 : 0.2,
              }}
            >
              <div className="relative h-24 w-24 rounded-full bg-gradient-to-br from-violet-500/20 via-fuchsia-400/15 to-rose-300/20 blur-xl" />
              <div className="absolute left-1/2 top-1/2 -translate-x-1/2 -translate-y-1/2">
                {molecules.map((point) => (
                  <motion.span
                    key={point.id}
                    animate={{
                      opacity: 0.95,
                      scale: 1,
                      x: point.position.x * INNER_SCALE + CENTER_OFFSET_X,
                      y: point.position.y * INNER_SCALE + CENTER_OFFSET_Y,
                    }}
                    className="absolute rounded-full bg-white shadow-[0_0_12px_rgba(255,255,255,0.4)]"
                    initial={{ opacity: 0, scale: 0 }}
                    style={{ width: point.size, height: point.size }}
                    transition={{
                      duration: 0.7,
                      ease,
                      delay: point.delay,
                    }}
                  />
                ))}
              </div>
            </motion.div>
          )}
        </AnimatePresence>
      </motion.div>
    </button>
  );
}
