"use client";

import {
  Dropdown,
  DropdownTrigger,
  DropdownMenu,
  DropdownItem,
} from "@heroui/dropdown";
import { Avatar } from "@heroui/avatar";
import { Button } from "@heroui/button";
import { useSession, signIn, signOut } from "next-auth/react";
import NextLink from "next/link";

export function UserMenu() {
  const { data: session, status } = useSession();

  if (status === "loading") return null;

  if (!session?.user) {
    return (
      <Button size="sm" variant="flat" onPress={() => signIn("google")}>
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
