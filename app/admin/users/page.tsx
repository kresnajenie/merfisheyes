"use client";

import { useEffect, useState, useCallback } from "react";
import {
  Table,
  TableHeader,
  TableColumn,
  TableBody,
  TableRow,
  TableCell,
} from "@heroui/table";
import { Chip } from "@heroui/chip";
import { Spinner } from "@heroui/spinner";
import { Avatar } from "@heroui/avatar";
import { useSession } from "next-auth/react";
import { useRouter } from "next/navigation";

interface UserItem {
  id: string;
  name: string | null;
  email: string;
  image: string | null;
  role: "USER" | "ADMIN" | "SUPER_ADMIN";
  permissions: string[];
  createdAt: string;
}

export default function AdminUsersPage() {
  const { data: session } = useSession();
  const router = useRouter();
  const [users, setUsers] = useState<UserItem[]>([]);
  const [loading, setLoading] = useState(true);
  const [togglingId, setTogglingId] = useState<string | null>(null);

  // Redirect non-SUPER_ADMIN users
  useEffect(() => {
    if (session && session.user.role !== "SUPER_ADMIN") {
      router.replace("/admin");
    }
  }, [session, router]);

  const fetchUsers = useCallback(async () => {
    setLoading(true);
    const res = await fetch("/api/admin/users");

    if (res.ok) {
      const data = await res.json();

      setUsers(data.users ?? []);
    }
    setLoading(false);
  }, []);

  useEffect(() => {
    fetchUsers();
  }, [fetchUsers]);

  const toggleRole = async (user: UserItem) => {
    if (user.role === "SUPER_ADMIN") return;
    const newRole = user.role === "ADMIN" ? "USER" : "ADMIN";

    setTogglingId(user.id);
    await fetch(`/api/admin/users/${user.id}`, {
      method: "PATCH",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ role: newRole }),
    });
    await fetchUsers();
    setTogglingId(null);
  };

  const togglePresets = async (user: UserItem) => {
    setTogglingId(user.id);
    await fetch(`/api/admin/users/${user.id}`, {
      method: "PATCH",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        permissions: user.permissions.includes("presets:write")
          ? []
          : ["presets:write"],
      }),
    });
    await fetchUsers();
    setTogglingId(null);
  };

  const roleColor = (role: string) => {
    if (role === "SUPER_ADMIN") return "danger" as const;
    if (role === "ADMIN") return "warning" as const;

    return "default" as const;
  };

  return (
    <div>
      <div className="mb-6">
        <h1 className="text-2xl font-bold">User Management</h1>
        <p className="text-default-500 text-sm mt-1">
          Click a role chip to toggle between USER and ADMIN. Click a Presets
          chip to let a user edit technology presets at /presets and nothing
          else; admins always can.
        </p>
      </div>

      {loading ? (
        <div className="flex justify-center py-12">
          <Spinner size="lg" />
        </div>
      ) : users.length === 0 ? (
        <p className="text-default-500 text-center py-12">No users found.</p>
      ) : (
        <Table aria-label="User management">
          <TableHeader>
            <TableColumn>User</TableColumn>
            <TableColumn>Email</TableColumn>
            <TableColumn>Role</TableColumn>
            <TableColumn>Presets</TableColumn>
            <TableColumn>Joined</TableColumn>
          </TableHeader>
          <TableBody>
            {users.map((user) => (
              <TableRow key={user.id}>
                <TableCell>
                  <div className="flex items-center gap-3">
                    <Avatar
                      imgProps={{ style: { opacity: 1 } }}
                      name={user.name ?? undefined}
                      size="sm"
                      src={user.image ?? undefined}
                    />
                    <span className="font-medium">{user.name ?? "—"}</span>
                  </div>
                </TableCell>
                <TableCell>
                  <span className="text-sm text-default-500">{user.email}</span>
                </TableCell>
                <TableCell>
                  <Chip
                    className={
                      user.role !== "SUPER_ADMIN" ? "cursor-pointer" : ""
                    }
                    color={roleColor(user.role)}
                    isDisabled={togglingId === user.id}
                    size="sm"
                    variant="flat"
                    onClick={() => toggleRole(user)}
                  >
                    {user.role === "SUPER_ADMIN"
                      ? "Super Admin"
                      : user.role === "ADMIN"
                        ? "Admin"
                        : "User"}
                  </Chip>
                </TableCell>
                <TableCell>
                  {user.role === "USER" ? (
                    <Chip
                      className="cursor-pointer"
                      color={
                        user.permissions.includes("presets:write")
                          ? "success"
                          : "default"
                      }
                      data-testid={`presets-permission-${user.email}`}
                      isDisabled={togglingId === user.id}
                      size="sm"
                      variant="flat"
                      onClick={() => togglePresets(user)}
                    >
                      {user.permissions.includes("presets:write")
                        ? "Can edit"
                        : "No"}
                    </Chip>
                  ) : (
                    <span className="text-sm text-default-400">Always</span>
                  )}
                </TableCell>
                <TableCell>
                  <span className="text-sm text-default-500">
                    {new Date(user.createdAt).toLocaleDateString()}
                  </span>
                </TableCell>
              </TableRow>
            ))}
          </TableBody>
        </Table>
      )}
    </div>
  );
}
