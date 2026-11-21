import { redirect } from "next/navigation";

import { auth } from "@/lib/auth";
import { canEditPresets } from "@/lib/presets/permission";
import { PresetEditor } from "@/components/presets/preset-editor";
import { title, subtitle } from "@/components/primitives";

export const metadata = {
  title: "Technology presets — MERFISHEYES",
};

export default async function PresetsPage() {
  const session = await auth();

  if (!session?.user) {
    redirect("/auth/signin?callbackUrl=/presets");
  }
  if (!(await canEditPresets(session.user))) {
    redirect("/");
  }

  return (
    <div className="w-full max-w-4xl mx-auto py-8">
      <div className="mb-6">
        <h1 className={title({ size: "md" })}>Technology presets</h1>
        <p className={subtitle({ class: "mt-2" })}>
          What a dataset of each technology opens with in the viewer. Its
          owner&apos;s saved settings and the link still override a preset.
        </p>
      </div>
      <PresetEditor />
    </div>
  );
}
