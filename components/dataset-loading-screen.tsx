"use client";

import { Button, Progress, Spinner } from "@heroui/react";
import { useRouter } from "next/navigation";

import LightRays from "@/components/react-bits/LightRays";
import { subtitle, title } from "@/components/primitives";

/** Which viewer family the screen stands in for; picks its colours. */
export type LoadingTone = "cell" | "sm";

/**
 * The one loading screen of every viewer: a full page (light rays behind a
 * spinner, a progress bar and the loader's message) or, inside a split
 * panel, the compact black version.
 */
export function DatasetLoadingScreen({
  progress,
  message,
  label = "Loading dataset...",
  tone = "cell",
  variant = "page",
}: {
  progress: number;
  message: string;
  label?: string;
  tone?: LoadingTone;
  variant?: "page" | "panel";
}) {
  if (variant === "panel") {
    return (
      <div className="absolute inset-0 flex items-center justify-center bg-black">
        <div className="flex flex-col items-center gap-4 w-full max-w-xs px-4">
          <Spinner color="primary" size="lg" />
          <p className="text-sm text-white/70">{label}</p>
          <Progress
            aria-label="Loading progress"
            className="w-full"
            color="primary"
            size="sm"
            value={progress}
          />
          <p className="text-xs text-white/50">{message}</p>
        </div>
      </div>
    );
  }
  const colour = tone === "sm" ? "secondary" : "primary";
  const rays = tone === "sm" ? ["#667eea", "#764ba2"] : ["#FFD700", "#FFD700"];

  return (
    <>
      <div className="fixed inset-0 w-full h-full z-0">
        <LightRays
          lightSpread={1.0}
          mouseInfluence={0.1}
          pulsating={false}
          rayLength={10}
          raysColor={rays[0]}
          raysOrigin="top-left"
          raysSpeed={0.8}
        />
      </div>
      <div className="fixed inset-0 w-full h-full z-0">
        <LightRays
          lightSpread={1.0}
          mouseInfluence={0.1}
          pulsating={false}
          rayLength={10}
          raysColor={rays[1]}
          raysOrigin="top-right"
          raysSpeed={0.8}
        />
      </div>
      <div className="relative z-10 flex items-center justify-center h-full">
        <div className="flex flex-col items-center gap-4 w-full max-w-md px-4">
          <Spinner color={colour} size="lg" />
          <p className={subtitle({ class: "text-center" })}>{label}</p>
          <Progress
            aria-label="Loading progress"
            className="w-full"
            color={colour}
            size="md"
            value={progress}
          />
          <p className="text-sm text-default-500 text-center">{message}</p>
        </div>
      </div>
    </>
  );
}

/** What a viewer shows when its dataset could not be opened. */
export function DatasetErrorScreen({
  error,
  heading = "Failed to load dataset",
  homeHref,
  onRetry,
  tone = "cell",
  variant = "page",
}: {
  error: string;
  heading?: string;
  /** Where "Go to Home" leads; no button without it. */
  homeHref?: string;
  onRetry?: () => void;
  tone?: LoadingTone;
  variant?: "page" | "panel";
}) {
  const router = useRouter();

  if (variant === "panel") {
    return (
      <div className="absolute inset-0 flex items-center justify-center bg-black">
        <div className="text-center px-4">
          <p className="text-sm text-danger">{error}</p>
        </div>
      </div>
    );
  }

  return (
    <>
      <div className="fixed inset-0 w-full h-full z-0">
        <LightRays
          lightSpread={1.0}
          mouseInfluence={0.1}
          pulsating={false}
          rayLength={10}
          raysColor="#FF72E1"
          raysOrigin="top-center"
          raysSpeed={0.8}
        />
      </div>
      <div className="relative z-10 flex items-center justify-center h-full p-8">
        <div className="flex flex-col items-center gap-6 max-w-2xl w-full">
          <div className="text-center">
            <h2 className={title({ size: "md", color: "pink" })}>{heading}</h2>
            <p className={subtitle({ class: "mt-4" })}>{error}</p>
          </div>
          {(homeHref || onRetry) && (
            <div className="flex gap-4">
              {homeHref && (
                <Button
                  color={tone === "sm" ? "secondary" : "primary"}
                  onPress={() => router.push(homeHref)}
                >
                  Go to Home
                </Button>
              )}
              {onRetry && (
                <Button color="default" variant="bordered" onPress={onRetry}>
                  Retry
                </Button>
              )}
            </div>
          )}
        </div>
      </div>
    </>
  );
}
