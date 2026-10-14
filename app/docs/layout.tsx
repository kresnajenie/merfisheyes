import { DocsNav, DocsPager } from "@/components/docs/docs-nav";

export default function DocsLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <section className="mx-auto flex w-full max-w-6xl flex-col gap-6 px-4 py-8 md:py-10 lg:flex-row lg:gap-10">
      <aside className="lg:w-52 lg:shrink-0">
        <DocsNav />
      </aside>
      <div className="flex min-w-0 flex-1 flex-col gap-10">
        {children}
        <DocsPager />
      </div>
    </section>
  );
}
