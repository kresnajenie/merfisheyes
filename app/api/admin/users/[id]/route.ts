import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { requireSuperAdmin } from "@/lib/admin-auth";

// PATCH /api/admin/users/[id] — update user role (SUPER_ADMIN only)
export async function PATCH(
  req: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { error } = await requireSuperAdmin();

  if (error) return error;

  const { id } = await params;
  const body = await req.json();
  const { role } = body;

  // Only allow setting to USER or ADMIN
  if (role !== "USER" && role !== "ADMIN") {
    return NextResponse.json(
      { error: "Can only set role to USER or ADMIN" },
      { status: 400 },
    );
  }

  // Check target user exists and is not SUPER_ADMIN
  const target = await prisma.user.findUnique({
    where: { id },
    select: { role: true },
  });

  if (!target) {
    return NextResponse.json({ error: "User not found" }, { status: 404 });
  }

  if (target.role === "SUPER_ADMIN") {
    return NextResponse.json(
      { error: "Cannot modify a SUPER_ADMIN user" },
      { status: 403 },
    );
  }

  const updated = await prisma.user.update({
    where: { id },
    data: { role },
    select: { id: true, role: true },
  });

  return NextResponse.json(updated);
}
