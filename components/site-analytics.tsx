"use client";

import { Analytics, type BeforeSendEvent } from "@vercel/analytics/next";

import { withoutLocalUrls } from "@/lib/utils/local-url";

/** Vercel Analytics, for every page. */
export function SiteAnalytics() {
  return (
    <Analytics
      beforeSend={(event: BeforeSendEvent) => ({
        ...event,
        url: withoutLocalUrls(event.url),
      })}
    />
  );
}
