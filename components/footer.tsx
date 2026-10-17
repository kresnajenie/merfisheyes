"use client";

import { Link } from "@heroui/link";
import { usePathname } from "next/navigation";

import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import { STAGE_RAIL_HEIGHT } from "@/lib/ui/stage-rail-geometry";

export function Footer() {
  const isSplitMode = useSplitScreenStore((s) => s.isSplitMode);
  const pathname = usePathname();
  // The labelled-molecule viewer is the only page whose bottom edge is taken
  // by something — the stage rail — so it is the only one that moves the
  // credit. Every other page keeps the footer in normal flow.
  const overRail = pathname?.startsWith("/lm-viewer") ?? false;

  return (
    <footer
      className={
        overRail
          ? // Right-aligned: the bottom-left holds the molecule count and the
            // scale bar, and the stage card opens centred above the rail.
            "fixed right-6 z-[var(--z-chrome)] flex items-center gap-1 text-sm"
          : `w-full flex items-center py-3 ${
              isSplitMode ? "justify-end pr-6" : "justify-center"
            }`
      }
      data-ui-overlay={overRail ? "" : undefined}
      style={overRail ? { bottom: STAGE_RAIL_HEIGHT + 10 } : undefined}
    >
      <Link
        isExternal
        className="flex items-center gap-1 text-current"
        href="https://www.linkedin.com/in/ignatius-jenie-1023521b3/"
        title="linkedin.com kresnajenie"
      >
        <span className="text-default-600">Made by</span>
        <p className="text-primary">Ignatius Jenie</p>
      </Link>
      <Link
        isExternal
        className="flex items-center gap-1 text-current"
        href="https://b.bintulab.com"
        title="b.bintulab.com"
      >
        <span className="text-default-600" />
        <span className="text-default-600">at</span>
        <p className="text-primary">Bintu Lab</p>
      </Link>
    </footer>
  );
}
