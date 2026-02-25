"use client";

import { Button } from "@heroui/button";
import { Card, CardBody, CardHeader } from "@heroui/card";
import { signIn } from "next-auth/react";
import { useSearchParams } from "next/navigation";
import { Suspense } from "react";

function SignInContent() {
  const searchParams = useSearchParams();
  const callbackUrl = searchParams.get("callbackUrl") || "/";

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
            Sign in to access MERFISHEYES
          </p>
        </CardHeader>
        <CardBody className="flex flex-col gap-4 items-center pb-8">
          <Button
            className="w-full max-w-xs"
            color="primary"
            size="lg"
            variant="flat"
            onPress={() => signIn("google", { callbackUrl })}
          >
            Continue with Google
          </Button>
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
