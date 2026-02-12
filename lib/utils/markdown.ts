import fs from "fs";
import path from "path";

import matter from "gray-matter";

export interface DocSection {
  id: string;
  title: string;
  level: number;
  children?: DocSection[];
}

export function getDocumentationContent() {
  const docsPath = path.join(process.cwd(), "DOCUMENTATION.md");
  const fileContents = fs.readFileSync(docsPath, "utf8");
  const { data, content } = matter(fileContents);

  return {
    frontmatter: data,
    content,
  };
}

export function parseTableOfContents(content: string): DocSection[] {
  // Remove code blocks first to avoid parsing # inside code
  const contentWithoutCodeBlocks = content.replace(/```[\s\S]*?```/g, "");

  const headingRegex = /^(#{1,3})\s+(.+)$/gm;
  const sections: DocSection[] = [];
  const stack: DocSection[] = [];
  let match;

  while ((match = headingRegex.exec(contentWithoutCodeBlocks)) !== null) {
    const level = match[1].length;
    const title = match[2].trim();
    const id = title
      .toLowerCase()
      .replace(/[^a-z0-9]+/g, "-")
      .replace(/(^-|-$)/g, "");

    const section: DocSection = {
      id,
      title,
      level,
      children: [],
    };

    // Find the appropriate parent based on level
    while (stack.length > 0 && stack[stack.length - 1].level >= level) {
      stack.pop();
    }

    if (stack.length === 0) {
      // Top-level section
      sections.push(section);
    } else {
      // Nested section
      const parent = stack[stack.length - 1];

      if (!parent.children) {
        parent.children = [];
      }
      parent.children.push(section);
    }

    stack.push(section);
  }

  return sections;
}
