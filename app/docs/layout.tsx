export default function DocsLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <section className="mx-auto flex w-full max-w-4xl flex-col gap-4 px-4 py-8 md:py-10">
      {children}
    </section>
  );
}
