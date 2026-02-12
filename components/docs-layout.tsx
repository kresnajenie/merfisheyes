"use client";

import { useState } from "react";
import { Link } from "@heroui/react";

import { DocSection } from "@/lib/utils/markdown";

interface DocsLayoutProps {
  mainSections: DocSection[];
  children: React.ReactNode;
}

export function DocsLayout({ mainSections, children }: DocsLayoutProps) {
  const [activeSection, setActiveSection] = useState<string>("");

  return (
    <div className="fixed inset-0 top-16 bg-background">
      <div className="flex h-full w-full">
        {/* Left Sidebar - Navigation */}
        <aside className="w-64 flex-shrink-0 border-r border-divider overflow-y-auto bg-background">
          <nav className="px-4 py-6 space-y-1 text-left">
            <div className="mb-6">
              <h2 className="text-lg font-bold text-foreground mb-4 text-left">
                Documentation
              </h2>
            </div>
            {mainSections && mainSections.length > 0 ? (
              mainSections.map((section) => (
                <div key={section.id}>
                  <Link
                    className={`block px-3 py-2 text-sm rounded-lg hover:bg-default-100 transition-colors font-bold text-left ${
                      activeSection === section.id
                        ? "bg-primary-50 dark:bg-primary-900/20 text-primary"
                        : "text-foreground"
                    }`}
                    href={`#${section.id}`}
                    onClick={() => setActiveSection(section.id)}
                  >
                    {section.title}
                  </Link>
                  {section.children && section.children.length > 0 && (
                    <div className="ml-4 mt-1 space-y-1">
                      {section.children.map((child) => (
                        <Link
                          key={child.id}
                          className={`block px-3 py-1.5 text-xs transition-colors text-left ${
                            activeSection === child.id
                              ? "text-primary font-medium"
                              : "text-default-600 hover:text-primary"
                          }`}
                          href={`#${child.id}`}
                          onClick={() => setActiveSection(child.id)}
                        >
                          {child.title}
                        </Link>
                      ))}
                    </div>
                  )}
                </div>
              ))
            ) : (
              <p className="text-sm text-default-500">No sections found</p>
            )}
          </nav>
        </aside>

        {/* Main Content Area */}
        <main className="flex-1 overflow-y-auto">
          <div className="px-8 py-8">
            <article className="prose prose-neutral dark:prose-invert max-w-none text-left">
              {children}
            </article>
          </div>
        </main>
      </div>
    </div>
  );
}
