import { prisma } from "@/lib/prisma";

/**
 * Fixed-window rate limits in shared storage (docs/plans/MCP.md D14): one
 * counter row per key per window, raised by one atomic upsert. Serverless
 * instances share no memory, so a per-instance counter would be bypassed by
 * spreading requests out.
 */

export interface Limit {
  /** What is counted, e.g. `agent:ip:1.2.3.4` — a rule name and a subject. */
  key: string;
  limit: number;
  windowSeconds: number;
}

export type LimitResult =
  | { ok: true }
  | { ok: false; retryAfter: number; key: string };

/** The client address as Vercel reports it (first x-forwarded-for hop). */
export function clientIp(req: Request): string {
  const forwarded = req.headers.get("x-forwarded-for")?.split(",")[0]?.trim();

  return forwarded || req.headers.get("x-real-ip") || "unknown";
}

/** Count one call against every limit; the first one exceeded refuses it.
 * Every limit is counted, refused or not, so hammering does not reset. */
export async function checkRateLimits(
  limits: Limit[],
  now = new Date(),
): Promise<LimitResult> {
  let refused: LimitResult = { ok: true };

  for (const l of limits) {
    const windowMs = l.windowSeconds * 1000;
    const start = new Date(Math.floor(now.getTime() / windowMs) * windowMs);
    const [row] = await prisma.$queryRaw<{ count: number }[]>`
      INSERT INTO rate_limit_counters (key, window_start, count)
      VALUES (${l.key}, ${start}, 1)
      ON CONFLICT (key, window_start)
      DO UPDATE SET count = rate_limit_counters.count + 1
      RETURNING count`;

    if (refused.ok && row.count > l.limit) {
      refused = {
        ok: false,
        key: l.key,
        retryAfter: Math.max(
          1,
          Math.ceil((start.getTime() + windowMs - now.getTime()) / 1000),
        ),
      };
    }
  }
  // Sweep windows older than the longest one (a day) now and then.
  if (Math.random() < 0.01) {
    await prisma.rateLimitCounter.deleteMany({
      where: { windowStart: { lt: new Date(now.getTime() - 2 * 86_400_000) } },
    });
  }

  return refused;
}

const MINUTE = 60;
const HOUR = 3600;
const DAY = 86_400;

/** The limits of D14, by what is being done and who does it. */
export const LIMITS = {
  agentAnonymous: (ip: string): Limit[] => [
    { key: `agent:ip:${ip}`, limit: 60, windowSeconds: MINUTE },
  ],
  agentUser: (who: { grantId: string | null; userId: string }): Limit[] => [
    {
      key: `agent:grant:${who.grantId ?? `session:${who.userId}`}`,
      limit: 300,
      windowSeconds: MINUTE,
    },
    { key: `agent:user:${who.userId}`, limit: 10_000, windowSeconds: DAY },
  ],
  /** Calls that read S3 or ask AWS Batch: a catalog row or dataset with its
   * sidecar files, a processing status. */
  agentExpensive: (subject: string): Limit[] => [
    { key: `agent:expensive:${subject}`, limit: 30, windowSeconds: MINUTE },
  ],
  register: (ip: string): Limit[] => [
    { key: `oauth:register:${ip}`, limit: 10, windowSeconds: HOUR },
  ],
  token: (clientId: string, ip: string): Limit[] => [
    {
      key: `oauth:token:${clientId.slice(0, 200)}:${ip}`,
      limit: 30,
      windowSeconds: MINUTE,
    },
  ],
  consent: (userId: string): Limit[] => [
    { key: `oauth:consent:${userId}`, limit: 20, windowSeconds: HOUR },
  ],
};
