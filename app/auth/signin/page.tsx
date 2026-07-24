"use client";

import { Button } from "@heroui/button";
import { Card, CardBody, CardHeader } from "@heroui/card";
import { Divider } from "@heroui/divider";
import { Input } from "@heroui/input";
import { getProviders, signIn } from "next-auth/react";
import { useSearchParams } from "next/navigation";
import { Suspense, useEffect, useState } from "react";

const EMAIL_CODE_PROVIDER_ID = "email-code";
const CODE_LENGTH = 6;

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
    <div className="flex items-center justify-center min-h-[60vh]">
      <Card
        isBlurred
        className="border-none bg-background/60 dark:bg-default-100/50 w-full max-w-md"
        shadow="sm"
      >
        <CardHeader className="flex flex-col gap-1 items-center pt-8">
          <h1 className="text-2xl font-bold">Sign In</h1>
          <p className="text-default-500 text-sm">
            {step === "email"
              ? "Sign in to access MERFISHEYES"
              : `Enter the ${CODE_LENGTH}-digit code sent to ${email}`}
          </p>
        </CardHeader>
        <CardBody className="flex flex-col gap-4 items-center pb-8">
          {error ? (
            <p className="w-full max-w-xs text-center text-sm text-danger">
              {error}
            </p>
          ) : null}

          {step === "email" ? (
            <>
              <Input
                className="max-w-xs"
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
                className="w-full max-w-xs"
                color="primary"
                isDisabled={!email.includes("@")}
                isLoading={busy}
                size="lg"
                onPress={sendCode}
              >
                Email me a code
              </Button>

              <div className="w-full max-w-xs flex items-center gap-3 py-1">
                <Divider className="flex-1" />
                <span className="text-tiny text-default-400">or</span>
                <Divider className="flex-1" />
              </div>

              <Button
                className="w-full max-w-xs"
                size="lg"
                variant="flat"
                onPress={() => signIn("google", { callbackUrl })}
              >
                Continue with Google
              </Button>

              {devLoginEnabled ? (
                <div className="w-full max-w-xs mt-2 rounded-medium border border-warning-400/60 bg-warning-50/40 p-3 flex flex-col gap-2">
                  <p className="text-tiny text-warning-700 dark:text-warning-500">
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
            </>
          ) : (
            <>
              <Input
                className="max-w-xs"
                inputMode="numeric"
                isDisabled={busy}
                label="Code"
                maxLength={CODE_LENGTH}
                placeholder="000000"
                value={code}
                onKeyDown={(e) => {
                  if (e.key === "Enter" && code.length === CODE_LENGTH) {
                    submitCode();
                  }
                }}
                onValueChange={(v) => setCode(v.replace(/\D/g, ""))}
              />
              <Button
                className="w-full max-w-xs"
                color="primary"
                isDisabled={code.length !== CODE_LENGTH}
                isLoading={busy}
                size="lg"
                onPress={submitCode}
              >
                Sign in
              </Button>
              <Button
                className="w-full max-w-xs"
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
            </>
          )}
        </CardBody>
      </Card>
    </div>
  );
}

export default function SignInPage() {
  return (
    <Suspense>
      <SignInContent />
    </Suspense>
  );
}
