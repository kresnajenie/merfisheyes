export default function ExploreLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <section className="flex flex-col items-center justify-center gap-4 py-8 md:py-10">
      <div className="inline-block w-full max-w-7xl justify-center">
        {children}
      </div>
    </section>
  );
}
