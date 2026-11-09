import { NextResponse } from "next/server";

import { auth } from "@/lib/auth";
import { issuer } from "@/lib/oauth/config";
import { verifySigned } from "@/lib/oauth/crypto";
import { LIMITS, checkRateLimits } from "@/lib/rate-limit";
import {
  type AuthorizeRequest,
  OAuthError,
  approve,
  errorRedirect,
} from "@/lib/oauth/flow";

/**
 * POST /api/oauth/authorize — the consent form's decision. Accepted only from
 * the site's own page (Origin), only for the signed request it rendered, and
 * only for the person it rendered it for.
 */

interface Signed {
  req: AuthorizeRequest;
  userId: string;
  exp: number;
}

const page = (message: string, status: number) =>
  new NextResponse(message, {
    status,
    headers: { "Content-Type": "text/plain; charset=utf-8" },
  });

export async function POST(request: Request) {
  const origin = request.headers.get("origin");

  if (origin !== new URL(issuer()).origin)
    return page("This form must be sent from MERFISHEYES itself.", 403);
  const form = await request.formData();
  const signed = verifySigned<Signed>(String(form.get("request") ?? ""));

  if (!signed || signed.exp < Date.now()) {
    return page("This request expired. Start connecting the agent again.", 400);
  }
  const session = await auth();

  if (session?.user?.id !== signed.userId)
    return page("You are signed in as someone else. Start again.", 403);
  const limited = await checkRateLimits(LIMITS.consent(signed.userId));

  if (!limited.ok) {
    return page(
      `Too many tries. Wait ${limited.retryAfter} s and start again.`,
      429,
    );
  }
  const { req } = signed;

  if (form.get("decision") !== "approve") {
    return NextResponse.redirect(
      errorRedirect(
        req.redirectUri,
        req.state,
        "access_denied",
        "the person declined",
      ),
      303,
    );
  }
  const reach =
    form.get("reach") === "some"
      ? {
          all: false as const,
          projectIds: form.getAll("projectIds").map(String),
        }
      : { all: true as const };

  try {
    return NextResponse.redirect(await approve(signed.userId, req, reach), 303);
  } catch (e) {
    if (e instanceof OAuthError)
      return page(`${e.message}. Go back and choose again.`, 400);
    throw e;
  }
}
