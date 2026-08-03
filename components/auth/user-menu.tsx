"use client";

import {
  Dropdown,
  DropdownTrigger,
  DropdownMenu,
  DropdownItem,
} from "@heroui/dropdown";
import { Avatar } from "@heroui/avatar";
import { Button } from "@heroui/button";
import { useSession, signOut } from "next-auth/react";
import NextLink from "next/link";
import { usePathname } from "next/navigation";

export function UserMenu() {
  const { data: session, status } = useSession();
  const pathname = usePathname();

  if (status === "loading") return null;

  if (!session?.user) {
    // Route to the full sign-in page (Google + email code), preserving where
    // the user was so they land back here after authenticating.
    const callbackUrl = encodeURIComponent(pathname || "/");

    return (
      <Button
        as={NextLink}
        href={`/auth/signin?callbackUrl=${callbackUrl}`}
        size="sm"
        variant="flat"
      >
        Sign In
      </Button>
    );
  }

  const isAdmin = session.user.role === "ADMIN" || session.user.role === "SUPER_ADMIN";

  return (
    <Dropdown placement="bottom-end">
      <DropdownTrigger>
        <Avatar
          isBordered
          as="button"
          className="transition-transform"
          color="primary"
          imgProps={{ style: { opacity: 1 } }}
          name={session.user.name ?? undefined}
          size="sm"
          src={session.user.image ?? undefined}
        />
      </DropdownTrigger>
      <DropdownMenu aria-label="User menu">
        <DropdownItem key="profile" isReadOnly className="h-14 gap-2">
          <p className="font-semibold">{session.user.name}</p>
          <p className="text-xs text-default-500">{session.user.email}</p>
        </DropdownItem>
        <DropdownItem key="account" as={NextLink} href="/account">
          Your datasets
        </DropdownItem>
        {isAdmin ? (
          <DropdownItem key="admin" as={NextLink} href="/admin">
            Admin Panel
          </DropdownItem>
        ) : null}
        <DropdownItem
          key="signout"
          color="danger"
          onPress={() => signOut()}
        >
          Sign Out
        </DropdownItem>
      </DropdownMenu>
    </Dropdown>
  );
}
