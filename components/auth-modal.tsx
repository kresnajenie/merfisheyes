"use client";

import { Button, Input } from "@heroui/react";
import { getSession, signIn, useSession } from "next-auth/react";
import Link from "next/link";
import { useEffect, useState } from "react";
import { createPortal } from "react-dom";

const PROVIDER = "email-code";
const CODE_LENGTH = 6;

interface AuthModalProps {
  isOpen: boolean;
  /** Closed without signing in. */
  onClose: () => void;
  /** Signed in successfully — the caller resumes whatever it was doing. */
  onAuthenticated: () => void;
}

/**
 * Quick email + 6-digit-code sign-in as a popup.
 *
 * Used from the upload cards so a signed-out user can drop a file, verify, and
 * have the upload continue — without navigating away (a full-page redirect
 * would drop the in-memory File). The code is therefore verified in place with
 * a `fetch` to the Auth.js callback (which sets the session cookie); a fresh
 * `getSession()` is the source of truth for success. Google / dev sign-in stay
 * on the full /auth/signin page.
 *
 * Built as a plain overlay (like RawUploadOverlay) rather than the HeroUI Modal,
 * whose backdrop/content chrome doesn't render solid in this app's theme.
 */
export function AuthModal({
  isOpen,
  onClose,
  onAuthenticated,
}: AuthModalProps) {
  const { update } = useSession();
  const [email, setEmail] = useState("");
  const [code, setCode] = useState("");
  const [step, setStep] = useState<"email" | "code">("email");
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState("");

  const reset = () => {
    setEmail("");
    setCode("");
    setStep("email");
    setError("");
    setBusy(false);
  };

  const close = () => {
    if (busy) return;
    reset();
    onClose();
  };

  // Escape closes the modal (unless a request is in flight).
  useEffect(() => {
    if (!isOpen) return;
    const onKey = (e: KeyboardEvent) => {
      if (e.key === "Escape") close();
    };

    window.addEventListener("keydown", onKey);

    return () => window.removeEventListener("keydown", onKey);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [isOpen, busy]);

  const sendCode = async () => {
    setError("");
    setBusy(true);
    try {
      const res = await signIn(PROVIDER, { email, redirect: false });

      if (res?.error) {
        setError(
          "Couldn't send a code to that address. Check it and try again.",
        );

        return;
      }
      setStep("code");
    } catch {
      setError("Something went wrong sending the code. Try again.");
    } finally {
      setBusy(false);
    }
  };

  const verify = async () => {
    setError("");
    setBusy(true);
    try {
      const params = new URLSearchParams({
        token: code.trim(),
        email: email.trim().toLowerCase(),
        callbackUrl: "/",
      });

      // Verify without navigating: the callback sets the session cookie and
      // 302s; we don't follow the redirect (that would reload and lose the
      // dropped file). getSession() then confirms whether it worked.
      await fetch(`/api/auth/callback/${PROVIDER}?${params.toString()}`, {
        method: "GET",
        redirect: "manual",
        credentials: "include",
      });

      const session = await getSession();

      if (!session?.user) {
        setError("That code is wrong or has expired. Request a new one.");

        return;
      }

      await update?.();
      reset();
      onAuthenticated();
    } catch {
      setError("Couldn't verify that code. Try again.");
    } finally {
      setBusy(false);
    }
  };

  if (!isOpen || typeof document === "undefined") return null;

  // Portal to <body>: the upload cards live inside a transformed ancestor
  // (the crossfade layout's translate/scale), which would otherwise anchor
  // this `fixed` overlay to that ancestor instead of the viewport.
  return createPortal(
    <div
      className="fixed inset-0 z-[var(--z-modal-top)] flex items-center justify-center bg-black/70 p-6 backdrop-blur-md"
      role="presentation"
      onClick={(e) => {
        // Backdrop click closes; clicks inside the card don't bubble to here.
        if (e.target === e.currentTarget) close();
      }}
    >
      <div
        aria-label="Sign in to upload"
        aria-modal="true"
        className="glass-panel w-full max-w-sm rounded-2xl p-7"
        role="dialog"
      >
        <h2 className="text-xl font-semibold text-foreground">
          {step === "email" ? "Sign in to upload" : "Enter your code"}
        </h2>
        <p className="mt-1 text-sm text-default-500">
          {step === "email"
            ? "We'll email you a 6-digit code — no password needed."
            : `Sent to ${email}`}
        </p>

        {error ? (
          <p className="mt-3 text-sm text-danger" role="alert">
            {error}
          </p>
        ) : null}

        {step === "email" ? (
          <div className="mt-5 flex flex-col gap-4">
            <Input
              isDisabled={busy}
              label="Email"
              placeholder="you@example.com"
              type="email"
              value={email}
              variant="bordered"
              onKeyDown={(e) => {
                if (e.key === "Enter" && email.includes("@")) sendCode();
              }}
              onValueChange={setEmail}
            />
            <Button
              color="primary"
              isDisabled={!email.includes("@")}
              isLoading={busy}
              size="lg"
              onPress={sendCode}
            >
              Email me a code
            </Button>
            <p className="text-center text-sm text-default-500">
              Prefer Google?{" "}
              <Link className="text-primary underline" href="/auth/signin">
                More sign-in options
              </Link>
            </p>
          </div>
        ) : (
          <div className="mt-5 flex flex-col gap-4">
            <Input
              inputMode="numeric"
              isDisabled={busy}
              label="6-digit code"
              maxLength={CODE_LENGTH}
              placeholder="000000"
              value={code}
              variant="bordered"
              onKeyDown={(e) => {
                if (e.key === "Enter" && code.length === CODE_LENGTH) verify();
              }}
              onValueChange={(v) => setCode(v.replace(/\D/g, ""))}
            />
            <Button
              color="primary"
              isDisabled={code.length !== CODE_LENGTH}
              isLoading={busy}
              size="lg"
              onPress={verify}
            >
              Verify &amp; continue
            </Button>
            <button
              className="text-sm text-default-500 underline hover:text-foreground disabled:opacity-50"
              disabled={busy}
              type="button"
              onClick={() => {
                setStep("email");
                setCode("");
                setError("");
              }}
            >
              Use a different email
            </button>
          </div>
        )}
      </div>
    </div>,
    document.body,
  );
}
