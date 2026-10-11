import { FormatGuide } from "@/components/format-guide";
import { subtitle, title } from "@/components/primitives";

export default function DocsPage() {
  return (
    <div className="flex flex-col gap-6">
      <div>
        <h1 className={title({ size: "sm" })}>Supported data formats</h1>
        <p className={subtitle({ class: "mt-2" })}>
          What you can upload, and how each format is expected to look.
        </p>
      </div>
      <FormatGuide />
    </div>
  );
}
