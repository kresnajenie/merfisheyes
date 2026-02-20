import { redirect } from "next/navigation";
import NextLink from "next/link";

import { auth } from "@/lib/auth";

export default async function AdminLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  const session = await auth();

  if (!session?.user) {
    redirect("/auth/signin?callbackUrl=/admin");
  }

  if (session.user.role !== "ADMIN") {
    redirect("/");
  }

  return (
    <div className="flex gap-6 w-full">
      {/* Sidebar */}
      <aside className="w-48 shrink-0 hidden md:block">
        <nav className="sticky top-24 flex flex-col gap-2">
          <h2 className="text-sm font-semibold text-default-500 uppercase tracking-wider mb-2">
            Admin
          </h2>
          <NextLink
            className="text-sm px-3 py-2 rounded-lg hover:bg-default-100 transition-colors"
            href="/admin/datasets"
          >
            Datasets
          </NextLink>
        </nav>
      </aside>

      {/* Main content */}
      <div className="flex-1 min-w-0">{children}</div>
    </div>
  );
}
