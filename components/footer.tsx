"use client";

import { Link } from "@heroui/link";

import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";

export function Footer() {
  const isSplitMode = useSplitScreenStore((s) => s.isSplitMode);

  return (
    <footer
      className={`w-full flex items-center py-3 ${
        isSplitMode ? "justify-end pr-6" : "justify-center"
      }`}
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
