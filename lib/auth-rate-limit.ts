import { prisma } from "@/lib/prisma";

/**
 * Throttling for the email sign-in code flow.
 *
 * Two separate limits, because they defend against different things:
 *
 * - **Issuance** stops someone using the sign-in form to mail-bomb an address
 *   they don't own (and stops us burning SES quota).
 * - **Verification** is the one that actually protects the account. A six-digit
 *   code is only a 10^6 space, and a wrong guess does not consume the stored
 *   token — Auth.js's adapter deletes it only on a *successful* lookup. Without
 *   a cap, an attacker can keep guessing for the code's full lifetime.
 *
 * State lives in Postgres rather than memory because serverless requests are
 * spread across instances; a per-instance counter is bypassed by just sending
 * more requests.
 */

const WINDOW_MS = 15 * 60 * 1000;
const MAX_ISSUES_PER_WINDOW = 5;
const MAX_FAILED_VERIFIES_PER_WINDOW = 8;

export const RATE_LIMIT_WINDOW_MINUTES = WINDOW_MS / 60000;

export type AttemptKind = "issue" | "verify_fail";

export function normalizeIdentifier(email: string): string {
  return email.trim().toLowerCase();
}

async function countRecent(
  identifier: string,
  kind: AttemptKind,
): Promise<number> {
  return prisma.loginAttempt.count({
    where: {
      identifier,
      kind,
      createdAt: { gt: new Date(Date.now() - WINDOW_MS) },
    },
  });
}

export async function recordAttempt(
  identifier: string,
  kind: AttemptKind,
): Promise<void> {
  await prisma.loginAttempt.create({ data: { identifier, kind } });
}

/**
 * Drop rows that have aged out of the window.
 *
 * Opportunistic — called on the write paths so the table cannot grow without
 * bound, which saves standing up a scheduled job for a table that only ever
 * holds a few minutes of history.
 */
export async function pruneAttempts(): Promise<void> {
  try {
    await prisma.loginAttempt.deleteMany({
      where: { createdAt: { lt: new Date(Date.now() - WINDOW_MS) } },
    });
  } catch {
    // Housekeeping only — never fail a sign-in over it.
  }
}

export async function clearAttempts(identifier: string): Promise<void> {
  try {
    await prisma.loginAttempt.deleteMany({ where: { identifier } });
  } catch {
    // Best effort: a stale counter expires on its own within the window.
  }
}

/** Throws when this address has requested too many codes recently. */
export async function assertCanIssue(identifier: string): Promise<void> {
  if ((await countRecent(identifier, "issue")) >= MAX_ISSUES_PER_WINDOW) {
    throw new Error(
      `Too many sign-in codes requested. Try again in ${RATE_LIMIT_WINDOW_MINUTES} minutes.`,
    );
  }
}

/** True when this address still has verification attempts left. */
export async function canAttemptVerify(identifier: string): Promise<boolean> {
  return (
    (await countRecent(identifier, "verify_fail")) <
    MAX_FAILED_VERIFIES_PER_WINDOW
  );
}
