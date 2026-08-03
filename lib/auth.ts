import type { Provider } from "next-auth/providers";

import NextAuth from "next-auth";
import Google from "next-auth/providers/google";
import Credentials from "next-auth/providers/credentials";
import { PrismaAdapter } from "@auth/prisma-adapter";

import { prisma } from "@/lib/prisma";
import { sendVerificationCodeEmail } from "@/lib/ses";
import {
  assertCanIssue,
  normalizeIdentifier,
  pruneAttempts,
  recordAttempt,
} from "@/lib/auth-rate-limit";

/** How long a sign-in code stays valid. */
export const EMAIL_CODE_MAX_AGE_SECONDS = 10 * 60;

/** Provider id for the email-code flow; referenced by the sign-in page. */
export const EMAIL_CODE_PROVIDER_ID = "email-code";

/**
 * Dev-only sign-in with no verification at all: type an address, you're in.
 *
 * Exists so the ingestion work can be tested repeatedly without an inbox
 * round-trip per upload. It is a login with no secret, so it must never be
 * enabled in production — set it in Vercel's Preview scope only, where
 * Deployment Protection already keeps the deployment non-public. The sign-in
 * page surfaces it prominently so it cannot be on unnoticed.
 */
export const DEV_EMAIL_LOGIN_ENABLED =
  process.env.ALLOW_DEV_EMAIL_LOGIN === "true" &&
  process.env.VERCEL_ENV !== "production";

function sixDigitCode(): string {
  // crypto.getRandomValues over Math.random: this is an auth secret, and
  // rejection-free modulo bias across 10^6 from a 32-bit draw is negligible
  // but avoided anyway by drawing until in range.
  const max = 1_000_000;
  const limit = Math.floor(0xffffffff / max) * max;
  const buf = new Uint32Array(1);

  let value = 0;

  do {
    crypto.getRandomValues(buf);
    value = buf[0];
  } while (value >= limit);

  return String(value % max).padStart(6, "0");
}

/**
 * Email + 6-digit code.
 *
 * Built as a plain `type: "email"` provider rather than via the `Nodemailer`
 * factory: that factory throws unless given an SMTP `server`, and imports
 * `nodemailer` at module load. We send through the SES client this app already
 * uses, so neither is wanted.
 *
 * Auth.js does the rest — it hashes the code before storing it, expires it, and
 * makes it single-use.
 */
const emailCodeProvider = {
  id: EMAIL_CODE_PROVIDER_ID,
  type: "email" as const,
  name: "Email code",
  from: process.env.SES_FROM_EMAIL || "noreply@merfisheyes.com",
  maxAge: EMAIL_CODE_MAX_AGE_SECONDS,
  generateVerificationToken: sixDigitCode,
  async sendVerificationRequest({
    identifier,
    token,
  }: {
    identifier: string;
    token: string;
  }) {
    const email = normalizeIdentifier(identifier);

    // Throwing here aborts the sign-in, so the code is never mailed. Auth.js
    // still writes the token, but an unsent code is harmless.
    await assertCanIssue(email);
    await recordAttempt(email, "issue");
    void pruneAttempts();

    await sendVerificationCodeEmail({
      email,
      code: token,
      expiresInMinutes: EMAIL_CODE_MAX_AGE_SECONDS / 60,
    });
  },
};

const providers: Provider[] = [
  Google({
    clientId: process.env.GOOGLE_CLIENT_ID!,
    clientSecret: process.env.GOOGLE_CLIENT_SECRET!,
  }),
  emailCodeProvider,
];

if (DEV_EMAIL_LOGIN_ENABLED) {
  providers.push(
    Credentials({
      id: "dev-email",
      name: "Developer email (no verification)",
      credentials: { email: { label: "Email", type: "email" } },
      async authorize(credentials) {
        const email = normalizeIdentifier(String(credentials?.email ?? ""));

        if (!email.includes("@")) return null;

        // Credentials bypasses the adapter, so the User row is created here —
        // Dataset.ownerId is a real foreign key and would not resolve without it.
        const user = await prisma.user.upsert({
          where: { email },
          update: {},
          create: { email },
        });

        return { id: user.id, email: user.email, name: user.name };
      },
    }),
  );
}

export const { handlers, auth, signIn, signOut } = NextAuth({
  adapter: PrismaAdapter(prisma),
  providers,
  session: { strategy: "jwt" },
  callbacks: {
    async jwt({ token, user }) {
      if (user) {
        token.id = user.id;
      }
      // Always refresh role from DB so admin changes take effect immediately
      if (token.id) {
        const dbUser = await prisma.user.findUnique({
          where: { id: token.id as string },
          select: { role: true },
        });

        token.role = dbUser?.role ?? "USER";
      }

      return token;
    },
    async session({ session, token }) {
      session.user.id = token.id as string;
      session.user.role = token.role as "USER" | "ADMIN" | "SUPER_ADMIN";

      return session;
    },
  },
  pages: {
    signIn: "/auth/signin",
  },
});
