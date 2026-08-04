import { type Page, expect } from "@playwright/test";

/**
 * Sign the current browser context in via the no-verification `dev-email`
 * NextAuth provider (enabled in tests by ALLOW_DEV_EMAIL_LOGIN). Uploads now
 * require a signed-in owner, so the @upload suite calls this before driving the
 * upload modal. Uses the API sign-in flow (CSRF + credentials callback) so the
 * session cookie lands in the shared context — no UI dependency.
 */
export async function signInDev(
  page: Page,
  email = "e2e@example.com",
): Promise<void> {
  const csrf = await (await page.request.get("/api/auth/csrf")).json();

  const res = await page.request.post("/api/auth/callback/dev-email", {
    form: {
      csrfToken: csrf.csrfToken,
      email,
      callbackUrl: "/",
      json: "true",
    },
  });

  expect(res.ok(), "dev-email sign-in request").toBeTruthy();

  // Confirm the session is established before the test proceeds.
  const session = await (await page.request.get("/api/auth/session")).json();

  expect(session?.user?.email, "signed-in session user").toBeTruthy();
}
