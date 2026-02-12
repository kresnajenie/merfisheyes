import {
  getDocumentationContent,
  parseTableOfContents,
} from "@/lib/utils/markdown";
import { DocsLayout } from "@/components/docs-layout";
import { MarkdownContent } from "@/components/markdown-content";

export default async function DocumentationPage() {
  const { content } = getDocumentationContent();
  const toc = parseTableOfContents(content);

  // Extract h2 sections from the hierarchical tree (they may be nested under h1s)
  const mainSections: any[] = [];

  toc.forEach((section) => {
    if (section.level === 2) {
      mainSections.push(section);
    }
    // If this is an h1, check its children for h2 sections
    if (section.level === 1 && section.children) {
      section.children.forEach((child) => {
        if (child.level === 2) {
          mainSections.push(child);
        }
      });
    }
  });

  // Serialize sections for client component (remove circular refs, ensure plain objects)
  const serializedSections = JSON.parse(JSON.stringify(mainSections));

  return (
    <DocsLayout mainSections={serializedSections}>
      <MarkdownContent content={content} />
    </DocsLayout>
  );
}
