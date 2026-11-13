"use client";

import { Suspense, useEffect, useState } from "react";
import { useRouter, useSearchParams } from "next/navigation";
import { Spinner } from "@heroui/react";

import { isLocalUrl } from "@/lib/utils/local-url";
import { viewerPathForUrl } from "@/lib/viewer/viewer-for-url";

/**
 * `/open?url=<URL>` — opens whatever the URL holds in the viewer that reads
 * it (`viewerPathForUrl`): a merfisheyes folder, a SpatialData store or a
 * folder holding one, an OME-Zarr image. Every other parameter (`v`, `split`,
 * `splitS3Url`, `rv`, …) is passed on. It is the link `merfisheyes-mcp view`
 * opens for a folder it serves from this computer.
 */
function OpenContent() {
  const params = useSearchParams();
  const router = useRouter();
  const url = params.get("url");
  const [step, setStep] = useState("Opening…");
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!url) return;
    let cancelled = false;

    viewerPathForUrl(url, { onStep: setStep }).then(
      (path) => {
        if (cancelled) return;
        const target = new URL(path, window.location.origin);

        params.forEach((value, key) => {
          if (key !== "url") target.searchParams.set(key, value);
        });
        router.replace(`${target.pathname}${target.search}`);
      },
      (e: Error) => !cancelled && setError(e.message),
    );

    return () => {
      cancelled = true;
    };
    // `params` is read once, with the URL it came with.
  }, [url, router]);

  if (!url) {
    return (
      <div className="flex h-full items-center justify-center p-6">
        Add ?url= with the address of a dataset, store or image.
      </div>
    );
  }
  if (error) {
    return (
      <div className="mx-auto flex h-full max-w-xl flex-col justify-center gap-3 p-6">
        <p className="font-semibold">Couldn&apos;t open {url}</p>
        <p className="whitespace-pre-line text-sm text-default-500">{error}</p>
        {isLocalUrl(url) && (
          <p className="text-sm text-default-500">
            This is a folder on your computer: is `merfisheyes-mcp view` still
            running? In Chrome, allow www.merfisheyes.com to reach apps on this
            device when it asks.
          </p>
        )}
      </div>
    );
  }

  return (
    <div className="flex h-full flex-col items-center justify-center gap-3">
      <Spinner />
      <p className="text-sm text-default-500">{step}</p>
    </div>
  );
}

export default function OpenPage() {
  return (
    <Suspense fallback={null}>
      <OpenContent />
    </Suspense>
  );
}
