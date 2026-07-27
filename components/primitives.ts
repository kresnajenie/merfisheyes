import { tv } from "tailwind-variants";

export const title = tv({
  base: "tracking-tight inline font-semibold",
  variants: {
    color: {
      violet: "from-[#FF1CF7] to-[#b249f8]",
      yellow: "from-[#FF705B] to-[#FFB457]",
      blue: "from-[#5EA2EF] to-[#0072F5]",
      cyan: "from-[#00b7fa] to-[#01cfea]",
      green: "from-[#6FEE8D] to-[#17c964]",
      pink: "from-[#FF72E1] to-[#F54C7A]",
      foreground: "dark:from-[#FFFFFF] dark:to-[#4B4B4B]",
    },
    size: {
      sm: "text-3xl lg:text-4xl",
      md: "text-[2.3rem] lg:text-5xl",
      lg: "text-4xl lg:text-6xl",
      xl: "text-5xl lg:text-7xl",
    },
    fullWidth: {
      true: "w-full block",
    },
  },
  defaultVariants: {
    size: "md",
  },
  compoundVariants: [
    {
      color: [
        "violet",
        "yellow",
        "blue",
        "cyan",
        "green",
        "pink",
        "foreground",
      ],
      class: "bg-clip-text text-transparent bg-gradient-to-b",
    },
  ],
});

export const subtitle = tv({
  base: "w-full md:w-1/2 my-2 text-lg lg:text-xl text-default-600 block max-w-full",
  variants: {
    fullWidth: {
      true: "!w-full",
    },
  },
  defaultVariants: {
    fullWidth: true,
  },
});

/**
 * Chrome-tier glass for buttons and small capsules.
 *
 * Just the `.glass` class — the old `!bg-[…]` and `backdrop-blur-[50px]` were
 * both dead: `.glass` already sets the background, and its `blur(15px)
 * !important` overrode the 50px every time.
 */
export const glassButton = tv({
  base: "glass",
});

/**
 * Content-tier glass for flyout panels (viz / DEG / camera / legends / plot).
 * More opaque than `glassButton` so text and charts stay readable over the 3D
 * scene. Replaces the `border-2 border-white/20 rounded-3xl shadow-lg
 * ${glassButton()}` string that was copy-pasted across the panel components;
 * callers keep their own width/position. Border and shadow come from
 * `.glass-panel`, so don't re-add them.
 */
export const glassPanel = tv({
  base: "glass-panel rounded-3xl",
});
