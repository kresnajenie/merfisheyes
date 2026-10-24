"use client";

import type { HTMLAttributes } from "react";

import { useEffect, useRef, useState } from "react";

const FADE_PX = 28;

/**
 * A column that scrolls once it is taller than the room it has — the
 * viewers' button rails and legends, which grow with what a dataset holds.
 * No scrollbar over the scene: an edge fades out where there is more to see,
 * and the wheel or a swipe scrolls it.
 */
export function FadeScroll({
  className = "",
  style,
  children,
  ...rest
}: HTMLAttributes<HTMLDivElement>) {
  const ref = useRef<HTMLDivElement>(null);
  const [more, setMore] = useState({ above: false, below: false });

  useEffect(() => {
    const el = ref.current;

    if (!el) return;
    const update = () => {
      const above = el.scrollTop > 1;
      const below = el.scrollTop + el.clientHeight < el.scrollHeight - 1;

      setMore((m) =>
        m.above === above && m.below === below ? m : { above, below },
      );
    };
    // Its own size, and what is in it, both change what is out of view
    const resized = new ResizeObserver(update);
    const observeChildren = () => {
      for (const child of Array.from(el.children)) resized.observe(child);
    };
    const added = new MutationObserver(() => {
      observeChildren();
      update();
    });

    resized.observe(el);
    observeChildren();
    added.observe(el, { childList: true, subtree: true });
    el.addEventListener("scroll", update, { passive: true });
    update();

    return () => {
      resized.disconnect();
      added.disconnect();
      el.removeEventListener("scroll", update);
    };
  }, []);

  const mask = `linear-gradient(to bottom, ${
    more.above ? "transparent" : "black"
  } 0, black ${FADE_PX}px, black calc(100% - ${FADE_PX}px), ${
    more.below ? "transparent" : "black"
  } 100%)`;

  return (
    <div
      ref={ref}
      className={`overflow-y-auto overflow-x-hidden overscroll-contain [scrollbar-width:none] [&>*]:shrink-0 [&::-webkit-scrollbar]:hidden ${className}`}
      data-fade-scroll={more.above || more.below ? "overflowing" : "fits"}
      style={{ maskImage: mask, WebkitMaskImage: mask, ...style }}
      {...rest}
    >
      {children}
    </div>
  );
}
