"use client";

import { Button } from "@heroui/button";
import { useSession } from "next-auth/react";

export type OwnerValue = "me" | "admin";

/**
 * Admin-only owner picker (Me vs Admin), shown at dataset-creation time so an
 * admin can choose personal vs collective "admin" ownership. Renders nothing
 * for non-admins — their datasets are always personally owned. Controlled.
 */
export function OwnerChoice({
  value,
  onChange,
  className,
}: {
  value: OwnerValue;
  onChange: (v: OwnerValue) => void;
  className?: string;
}) {
  const { data: session } = useSession();
  const isAdmin =
    session?.user?.role === "ADMIN" || session?.user?.role === "SUPER_ADMIN";

  if (!isAdmin) return null;

  return (
    <div className={className}>
      <p className="text-sm font-medium mb-1.5">Owner</p>
      <div className="flex gap-2">
        <Button
          color={value === "me" ? "primary" : "default"}
          size="sm"
          variant={value === "me" ? "solid" : "flat"}
          onPress={() => onChange("me")}
        >
          Me
        </Button>
        <Button
          color={value === "admin" ? "primary" : "default"}
          size="sm"
          variant={value === "admin" ? "solid" : "flat"}
          onPress={() => onChange("admin")}
        >
          Admin
        </Button>
      </div>
      <p className="text-xs text-default-500 mt-1">
        {value === "admin"
          ? "Owned by the admin team (shared across admins)."
          : "Owned by your account."}
      </p>
    </div>
  );
}
