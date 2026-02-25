// Middleware intentionally left minimal.
// Admin auth is handled in app/admin/layout.tsx via server-side auth() check.
export { auth as middleware } from "@/lib/auth";

export const config = {
  matcher: ["/admin/:path*"],
};
