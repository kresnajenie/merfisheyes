import type { NextRequest } from "next/server";

import { handlers, EMAIL_CODE_PROVIDER_ID } from "@/lib/auth";
import {
  RATE_LIMIT_WINDOW_MINUTES,
  canAttemptVerify,
  clearAttempts,
  normalizeIdentifier,
  pruneAttempts,
  recordAttempt,
} from "@/lib/auth-rate-limit";

export const POST = handlers.POST;

const EMAIL_CALLBACK_PATH = `/api/auth/callback/${EMAIL_CODE_PROVIDER_ID}`;

/**
 * Auth.js's GET handler, with an attempt cap in front of the email-code callback.
 *
 * The cap cannot live inside the provider: Auth.js verifies the code internally
 * and the provider is never consulted on the callback leg. It matters because a
 * wrong code leaves the stored token intact (the adapter deletes it only on a
 * successful lookup), so a 10^6 code space would otherwise be guessable for the
 * code's whole lifetime.
 *
 * Success/failure is read off the redirect Auth.js returns — it sends failures
 * to an `?error=` URL. A misread only miscounts an attempt; it cannot let a bad
 * code through, because the verification itself is still Auth.js's.
 */
export async function GET(request: NextRequest) {
  const url = new URL(request.url);

  if (url.pathname !== EMAIL_CALLBACK_PATH) {
    return handlers.GET(request);
  }

  const email = normalizeIdentifier(url.searchParams.get("email") ?? "");

  if (email && !(await canAttemptVerify(email))) {
    return Response.redirect(
      new URL(
        `/auth/signin?error=TooManyAttempts&minutes=${RATE_LIMIT_WINDOW_MINUTES}`,
        url.origin,
      ),
      302,
    );
  }

  const response = await handlers.GET(request);

  if (email) {
    const location = response.headers.get("location") ?? "";

    if (location.includes("error=")) {
      await recordAttempt(email, "verify_fail");
      void pruneAttempts();
    } else {
      // Signed in — drop the counters so a later typo is not judged against
      // attempts from a session that already succeeded.
      void clearAttempts(email);
    }
  }

  return response;
}
