import { randomUUID } from "crypto";

import { NextRequest, NextResponse } from "next/server";

import { auth } from "@/lib/auth";
import { prisma } from "@/lib/prisma";

const corsHeaders = {
  "Access-Control-Allow-Origin": process.env.CORS_ORIGIN || "*",
  "Access-Control-Allow-Methods": "POST, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type",
};

const SID_COOKIE = "mev_sid";

export async function OPTIONS() {
  return NextResponse.json({}, { headers: corsHeaders });
}

/** UTC midnight for the given time — the per-day dedup bucket. */
function utcDay(d = new Date()): Date {
  return new Date(Date.UTC(d.getUTCFullYear(), d.getUTCMonth(), d.getUTCDate()));
}

/**
 * Record a view. Deduped per (dataset, session, UTC day) via a unique
 * constraint — a genuinely new row increments the cached count + daily rollup;
 * repeats are no-ops. No auth required. Anonymous viewers get an httpOnly
 * `mev_sid` cookie as their session key.
 */
export async function POST(
  request: NextRequest,
  { params }: { params: Promise<{ datasetId: string }> },
) {
  const { datasetId } = await params;

  // Session key: signed-in user id, else the anon cookie (mint if absent).
  const session = await auth();
  let sessionKey = session?.user?.id ?? null;
  let setCookie: string | null = null;

  if (!sessionKey) {
    const existing = request.cookies.get(SID_COOKIE)?.value;

    sessionKey = existing ?? randomUUID();
    if (!existing) setCookie = sessionKey;
  }

  const day = utcDay();

  try {
    // The unique (datasetId, sessionKey, day) makes this idempotent: the create
    // throws P2002 on a repeat, which we swallow as "already counted".
    await prisma.$transaction([
      prisma.datasetView.create({
        data: {
          datasetId,
          sessionKey,
          day,
          userId: session?.user?.id ?? null,
        },
      }),
      prisma.dataset.update({
        where: { id: datasetId },
        data: { viewCount: { increment: 1 } },
      }),
      prisma.datasetViewDaily.upsert({
        where: { datasetId_day: { datasetId, day } },
        create: { datasetId, day, count: 1 },
        update: { count: { increment: 1 } },
      }),
    ]);

    const row = await prisma.dataset.findUnique({
      where: { id: datasetId },
      select: { viewCount: true },
    });
    const res = NextResponse.json(
      { counted: true, viewCount: row?.viewCount ?? null },
      { headers: corsHeaders },
    );

    if (setCookie) {
      res.cookies.set(SID_COOKIE, setCookie, {
        httpOnly: true,
        sameSite: "lax",
        path: "/",
        maxAge: 60 * 60 * 24 * 365,
      });
    }

    return res;
  } catch (err: any) {
    // P2002 = duplicate (already counted this session/day). P2003/P2025 =
    // dataset row missing (unregistered) — nothing to count.
    if (err?.code === "P2002") {
      const res = NextResponse.json(
        { counted: false },
        { headers: corsHeaders },
      );

      if (setCookie) {
        res.cookies.set(SID_COOKIE, setCookie, {
          httpOnly: true,
          sameSite: "lax",
          path: "/",
          maxAge: 60 * 60 * 24 * 365,
        });
      }

      return res;
    }

    return NextResponse.json(
      { counted: false, error: "not_registered" },
      { status: 200, headers: corsHeaders },
    );
  }
}
