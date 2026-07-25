"use client";

import { Button } from "@heroui/button";
import { Divider } from "@heroui/divider";
import { Input } from "@heroui/input";
import { Spinner } from "@heroui/spinner";
import { getProviders, signIn } from "next-auth/react";
import { useSearchParams } from "next/navigation";
import { Suspense, useEffect, useState } from "react";

const EMAIL_CODE_PROVIDER_ID = "email-code";
const CODE_LENGTH = 6;

/** Google's multicolour "G", so the OAuth button reads as an official one. */
function GoogleIcon() {
  return (
    <svg aria-hidden="true" className="w-4 h-4" viewBox="0 0 48 48">
      <path
        fill="#EA4335"
        d="M24 9.5c3.54 0 6.71 1.22 9.21 3.6l6.85-6.85C35.9 2.38 30.47 0 24 0 14.62 0 6.51 5.38 2.56 13.22l7.98 6.19C12.43 13.72 17.74 9.5 24 9.5z"
      />
      <path
        fill="#4285F4"
        d="M46.98 24.55c0-1.57-.15-3.09-.38-4.55H24v9.02h12.94c-.58 2.96-2.26 5.48-4.78 7.18l7.73 6c4.51-4.18 7.09-10.36 7.09-17.65z"
      />
      <path
        fill="#FBBC05"
        d="M10.53 28.59c-.48-1.45-.76-2.99-.76-4.59s.27-3.14.76-4.59l-7.98-6.19C.92 16.46 0 20.12 0 24c0 3.88.92 7.54 2.56 10.78l7.97-6.19z"
      />
      <path
        fill="#34A853"
        d="M24 48c6.48 0 11.93-2.13 15.89-5.81l-7.73-6c-2.15 1.45-4.92 2.3-8.16 2.3-6.26 0-11.57-4.22-13.47-9.91l-7.98 6.19C6.51 42.62 14.62 48 24 48z"
      />
    </svg>
  );
}

function SignInContent() {
  const searchParams = useSearchParams();
  const callbackUrl = searchParams.get("callbackUrl") || "/";

  const [email, setEmail] = useState("");
  const [code, setCode] = useState("");
  const [step, setStep] = useState<"email" | "code">("email");
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState("");
  // The dev bypass is a server-side env flag, so rather than mirror it into a
  // NEXT_PUBLIC_ twin that can drift, ask Auth.js which providers are actually
  // registered.
  const [devLoginEnabled, setDevLoginEnabled] = useState(false);

  useEffect(() => {
    getProviders()
      .then((providers) =>
        setDevLoginEnabled(Boolean(providers?.["dev-email"])),
      )
      .catch(() => setDevLoginEnabled(false));
  }, []);

  useEffect(() => {
    const param = searchParams.get("error");

    if (param === "TooManyAttempts") {
      setError(
        `Too many incorrect codes. Try again in ${searchParams.get("minutes") || 15} minutes.`,
      );
    } else if (param === "Verification") {
      setError("That code is wrong or has expired. Request a new one.");
    } else if (param) {
      setError("Sign in failed. Please try again.");
    }
  }, [searchParams]);

  const sendCode = async () => {
    setError("");
    setBusy(true);
    try {
      const result = await signIn(EMAIL_CODE_PROVIDER_ID, {
        email,
        callbackUrl,
        redirect: false,
      });

      if (result?.error) {
        setError(
          "Couldn't send a code to that address. Check it and try again.",
        );

        return;
      }
      setStep("code");
    } finally {
      setBusy(false);
    }
  };

  // Auth.js verifies the code on its own callback route, so submitting is just
  // a navigation to that URL with the typed code as the token.
  const submitCode = () => {
    setBusy(true);
    const params = new URLSearchParams({
      callbackUrl,
      token: code.trim(),
      email: email.trim().toLowerCase(),
    });

    window.location.href = `/api/auth/callback/${EMAIL_CODE_PROVIDER_ID}?${params}`;
  };

  return (
    <div className="flex items-center justify-center min-h-[70vh] px-4">
      <div className="glass-panel rounded-3xl w-full max-w-sm p-8 flex flex-col items-center gap-5">
        <div className="flex flex-col items-center gap-1 text-center">
          <span className="text-xs font-semibold tracking-[0.28em] text-white/80">
            MERFISHEYES
          </span>
          <h1 className="text-2xl font-bold">Sign In</h1>
          <p className="text-default-400 text-sm">
            {step === "email"
              ? "Sign in to upload and manage your datasets"
              : `Enter the ${CODE_LENGTH}-digit code sent to ${email}`}
          </p>
        </div>

        {error ? (
          <p className="w-full text-center text-sm text-danger">{error}</p>
        ) : null}

        {step === "email" ? (
          <div className="w-full flex flex-col gap-4">
            <Input
              isDisabled={busy}
              label="Email"
              placeholder="you@example.com"
              type="email"
              value={email}
              onKeyDown={(e) => {
                if (e.key === "Enter" && email.includes("@")) sendCode();
              }}
              onValueChange={setEmail}
            />
            <Button
              className="w-full"
              color="primary"
              isDisabled={!email.includes("@")}
              isLoading={busy}
              size="lg"
              onPress={sendCode}
            >
              Email me a code
            </Button>

            <div className="w-full flex items-center gap-3 py-1">
              <Divider className="flex-1" />
              <span className="text-tiny text-default-400">or</span>
              <Divider className="flex-1" />
            </div>

            <Button
              className="w-full"
              size="lg"
              startContent={<GoogleIcon />}
              variant="flat"
              onPress={() => signIn("google", { callbackUrl })}
            >
              Continue with Google
            </Button>

            {devLoginEnabled ? (
              <div className="w-full mt-1 rounded-medium border border-warning-400/60 bg-warning-50/10 p-3 flex flex-col gap-2">
                <p className="text-tiny text-warning-500">
                  Dev sign-in is enabled: any address signs in with no
                  verification. Never enable this in production.
                </p>
                <Button
                  color="warning"
                  isDisabled={!email.includes("@")}
                  size="sm"
                  variant="flat"
                  onPress={() => signIn("dev-email", { email, callbackUrl })}
                >
                  Sign in without a code
                </Button>
              </div>
            ) : null}
          </div>
        ) : (
          <div className="w-full flex flex-col gap-4">
            <Input
              inputMode="numeric"
              label="Code"
              maxLength={CODE_LENGTH}
              placeholder="000000"
              value={code}
              isDisabled={busy}
              onKeyDown={(e) => {
                if (e.key === "Enter" && code.length === CODE_LENGTH) {
                  submitCode();
                }
              }}
              onValueChange={(v) => setCode(v.replace(/\D/g, ""))}
            />
            <Button
              className="w-full"
              color="primary"
              isDisabled={code.length !== CODE_LENGTH}
              isLoading={busy}
              size="lg"
              onPress={submitCode}
            >
              Sign in
            </Button>
            <Button
              className="w-full"
              isDisabled={busy}
              size="sm"
              variant="light"
              onPress={() => {
                setStep("email");
                setCode("");
                setError("");
              }}
            >
              Use a different email
            </Button>
          </div>
        )}
      </div>
    </div>
  );
}

export default function SignInPage() {
  return (
    <Suspense
      fallback={
        <div className="flex items-center justify-center min-h-[70vh]">
          <Spinner size="lg" />
        </div>
      }
    >
      <SignInContent />
    </Suspense>
  );
}
