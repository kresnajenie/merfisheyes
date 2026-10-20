"use client";

import NextLink from "next/link";
import { usePathname } from "next/navigation";

export const DOC_PAGES = [
  { href: "/docs", label: "Getting started" },
  { href: "/docs/formats", label: "Data formats" },
  { href: "/docs/loading", label: "Loading data" },
  { href: "/docs/prepare", label: "Preparing data" },
  { href: "/docs/viewer", label: "Using the viewer" },
  { href: "/docs/expression", label: "Expression & DEGs" },
  { href: "/docs/sharing", label: "Sharing on Explore" },
];

/**
 * The docs' page list: a sidebar on wide screens, a row of pills that
 * scrolls sideways on narrow ones.
 */
export function DocsNav() {
  const pathname = usePathname();

  return (
    <nav aria-label="Docs" className="lg:sticky lg:top-28">
      <ul className="-mx-4 flex gap-2 overflow-x-auto px-4 pb-1 lg:mx-0 lg:flex-col lg:gap-1 lg:overflow-visible lg:px-0">
        {DOC_PAGES.map((page) => {
          const active = pathname === page.href;

          return (
            <li key={page.href} className="shrink-0">
              <NextLink
                aria-current={active ? "page" : undefined}
                className={`block whitespace-nowrap rounded-full px-3 py-1.5 text-sm transition-colors lg:rounded-lg ${
                  active
                    ? "bg-primary/10 font-medium text-primary"
                    : "text-default-600 hover:bg-default-100 hover:text-foreground"
                }`}
                href={page.href}
              >
                {page.label}
              </NextLink>
            </li>
          );
        })}
      </ul>
    </nav>
  );
}

/** Previous / next links at the foot of a docs page. */
export function DocsPager() {
  const pathname = usePathname();
  const i = DOC_PAGES.findIndex((p) => p.href === pathname);
  const prev = i > 0 ? DOC_PAGES[i - 1] : null;
  const next = i >= 0 && i < DOC_PAGES.length - 1 ? DOC_PAGES[i + 1] : null;

  return (
    <div className="mt-4 flex justify-between gap-4 border-t border-default-200 pt-6 text-sm">
      {prev ? (
        <NextLink
          className="text-default-600 hover:text-primary"
          href={prev.href}
        >
          ← {prev.label}
        </NextLink>
      ) : (
        <span />
      )}
      {next && (
        <NextLink
          className="text-default-600 hover:text-primary"
          href={next.href}
        >
          {next.label} →
        </NextLink>
      )}
    </div>
  );
}
