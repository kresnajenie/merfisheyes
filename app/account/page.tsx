import { redirect } from "next/navigation";

import { auth } from "@/lib/auth";
import { AccountWorkspace } from "@/components/account/account-workspace";
import { title, subtitle } from "@/components/primitives";

export const metadata = {
  title: "Your datasets — MERFISHEYES",
};

export default async function AccountPage() {
  const session = await auth();

  if (!session?.user) {
    redirect("/auth/signin?callbackUrl=/account");
  }

  return (
    <div className="w-full max-w-6xl mx-auto py-8">
      <div className="mb-6">
        <h1 className={title({ size: "md" })}>Your datasets</h1>
        <p className={subtitle({ class: "mt-2" })}>
          Manage the datasets you own — edit their details, group them into
          projects, and submit them to Explore.
        </p>
      </div>
      <AccountWorkspace />
    </div>
  );
}
