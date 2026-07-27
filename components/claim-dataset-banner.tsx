"use client";

import { Button } from "@heroui/button";
import { Input } from "@heroui/input";
import { useSession, signIn } from "next-auth/react";
import { useRouter, useSearchParams } from "next/navigation";
import { useEffect, useRef, useState } from "react";
import { toast } from "react-toastify";

import { useViewerRegistrationStore } from "@/lib/stores/viewerRegistrationStore";
import { glassPanel } from "@/components/primitives";

const EMAIL_CODE_PROVIDER_ID = "email-code";
const CODE_LENGTH = 6;

/**
 * Top banner shown in the from-s3 viewer when the dataset isn't owned yet
 * (unregistered, or registered but ownerless). Loading is never blocked by
 * this — the dataset renders regardless. Signed-in users claim in one click;
 * signed-out users get an inline email-code sign-in, then auto-claim on return.
 */
export function ClaimDatasetBanner() {
  const { data: session, status } = useSession();
  const router = useRouter();
  const searchParams = useSearchParams();

  const { dbId, ownerId, registered, s3Url, set } = useViewerRegistrationStore();

  const [dismissed, setDismissed] = useState(false);
  const [busy, setBusy] = useState(false);
  const [step, setStep] = useState<"idle" | "email" | "code">("idle");
  const [email, setEmail] = useState("");
  const [code, setCode] = useState("");
  const autoClaimedRef = useRef(false);

  const userId = session?.user?.id ?? null;
  // "Not owned yet" = no owner. Owned-by-someone-else is simply hidden.
  const notOwned = !!s3Url && ownerId == null;
  const ownedByOther = !!ownerId && !!userId && ownerId !== userId;

  async function claim() {
    if (!s3Url) return;
    setBusy(true);
    try {
      const res = await fetch("/api/datasets/register-s3", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ url: s3Url }),
      });

      if (!res.ok) {
        toast.error("Couldn't add this dataset.");

        return;
      }
      const { dataset, ownedByOther: taken } = await res.json();

      if (taken) {
        toast.error("This dataset is already owned by another user.");
        set({ dbId: dataset.id, ownerId: dataset.ownerId, registered: true });

        return;
      }
      set({ dbId: dataset.id, ownerId: dataset.ownerId, registered: true });
      toast.success("Dataset added to your account.");
    } catch {
      toast.error("Couldn't add this dataset.");
    } finally {
      setBusy(false);
    }
  }

  // Auto-claim after returning from the inline email-code sign-in.
  useEffect(() => {
    if (
      status === "authenticated" &&
      searchParams.get("claim") === "1" &&
      notOwned &&
      !autoClaimedRef.current
    ) {
      autoClaimedRef.current = true;
      claim().finally(() => {
        // Strip the one-shot flag so a reload doesn't re-trigger.
        const params = new URLSearchParams(window.location.search);

        params.delete("claim");
        const qs = params.toString();

        router.replace(qs ? `?${qs}` : window.location.pathname, {
          scroll: false,
        });
      });
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [status, notOwned, searchParams]);

  if (dismissed || ownedByOther || !notOwned || status === "loading") {
    return null;
  }

  // Where to return after the email-code sign-in completes: back here with the
  // one-shot claim flag so the effect above finishes the claim.
  function claimCallbackUrl() {
    const params = new URLSearchParams(window.location.search);

    params.set("claim", "1");

    return `${window.location.pathname}?${params.toString()}`;
  }

  async function sendCode() {
    setBusy(true);
    try {
      const result = await signIn(EMAIL_CODE_PROVIDER_ID, {
        email,
        callbackUrl: claimCallbackUrl(),
        redirect: false,
      });

      if (result?.error) {
        toast.error("Couldn't send a code to that address.");

        return;
      }
      setStep("code");
    } finally {
      setBusy(false);
    }
  }

  function submitCode() {
    const params = new URLSearchParams({
      callbackUrl: claimCallbackUrl(),
      token: code.trim(),
      email: email.trim().toLowerCase(),
    });

    window.location.href = `/api/auth/callback/${EMAIL_CODE_PROVIDER_ID}?${params}`;
  }

  return (
    <div
      className="absolute top-4 left-1/2 -translate-x-1/2 z-[var(--z-chrome)] w-[min(92vw,560px)]"
      data-ui-overlay
    >
      <div className={`${glassPanel()} px-4 py-3 flex items-center gap-3`}>
        <div className="flex-1 min-w-0">
          <p className="text-sm font-medium text-white">
            This dataset isn&apos;t owned yet
          </p>
          <p className="text-xs text-white/60">
            Add it to your account to track views and save viewer defaults.
          </p>
        </div>

        {userId ? (
          <Button color="primary" isLoading={busy} size="sm" onPress={claim}>
            Add to my account
          </Button>
        ) : step === "idle" ? (
          <Button
            color="primary"
            size="sm"
            variant="flat"
            onPress={() => setStep("email")}
          >
            Add this dataset
          </Button>
        ) : step === "email" ? (
          <div className="flex items-center gap-2">
            <Input
              aria-label="Email"
              className="w-48"
              isDisabled={busy}
              placeholder="you@example.com"
              size="sm"
              type="email"
              value={email}
              onKeyDown={(e) => {
                if (e.key === "Enter" && email.includes("@")) sendCode();
              }}
              onValueChange={setEmail}
            />
            <Button
              color="primary"
              isDisabled={!email.includes("@")}
              isLoading={busy}
              size="sm"
              onPress={sendCode}
            >
              Email code
            </Button>
          </div>
        ) : (
          <div className="flex items-center gap-2">
            <Input
              aria-label="Code"
              className="w-28"
              inputMode="numeric"
              isDisabled={busy}
              maxLength={CODE_LENGTH}
              placeholder="000000"
              size="sm"
              value={code}
              onKeyDown={(e) => {
                if (e.key === "Enter" && code.length === CODE_LENGTH) submitCode();
              }}
              onValueChange={(v) => setCode(v.replace(/\D/g, ""))}
            />
            <Button
              color="primary"
              isDisabled={code.length !== CODE_LENGTH}
              size="sm"
              onPress={submitCode}
            >
              Sign in
            </Button>
          </div>
        )}

        <button
          aria-label="Dismiss"
          className="text-white/50 hover:text-white transition-colors"
          type="button"
          onClick={() => setDismissed(true)}
        >
          <svg className="w-4 h-4" fill="none" stroke="currentColor" strokeWidth={2} viewBox="0 0 24 24">
            <path d="M6 18L18 6M6 6l12 12" strokeLinecap="round" strokeLinejoin="round" />
          </svg>
        </button>
      </div>
    </div>
  );
}
