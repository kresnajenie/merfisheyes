import { NextResponse } from "next/server";

import { publishedPresets } from "@/lib/presets/published";

// GET /api/presets — the published preset of every technology that has one.
// Public: the viewer reads it on every open.
export async function GET() {
  return NextResponse.json(
    { presets: await publishedPresets() },
    {
      headers: {
        "Cache-Control": "public, s-maxage=300, stale-while-revalidate=60",
      },
    },
  );
}
