"use client";

import { Code, Link, Divider } from "@heroui/react";
import {
  Table,
  TableHeader,
  TableColumn,
  TableBody,
  TableRow,
  TableCell,
} from "@heroui/react";

export const MDXComponents = {
  h1: ({ children, ...props }: any) => (
    <h1
      className="text-4xl font-bold mb-6 mt-8 scroll-mt-32"
      id={children
        ?.toString()
        .toLowerCase()
        .replace(/[^a-z0-9]+/g, "-")
        .replace(/(^-|-$)/g, "")}
      {...props}
    >
      {children}
    </h1>
  ),
  h2: ({ children, ...props }: any) => (
    <h2
      className="text-3xl font-bold mb-4 mt-8 scroll-mt-32"
      id={children
        ?.toString()
        .toLowerCase()
        .replace(/[^a-z0-9]+/g, "-")
        .replace(/(^-|-$)/g, "")}
      {...props}
    >
      {children}
    </h2>
  ),
  h3: ({ children, ...props }: any) => (
    <h3
      className="text-2xl font-bold mb-3 mt-6 scroll-mt-32"
      id={children
        ?.toString()
        .toLowerCase()
        .replace(/[^a-z0-9]+/g, "-")
        .replace(/(^-|-$)/g, "")}
      {...props}
    >
      {children}
    </h3>
  ),
  h4: ({ children, ...props }: any) => (
    <h4
      className="text-xl font-semibold mb-2 mt-4 scroll-mt-32"
      id={children
        ?.toString()
        .toLowerCase()
        .replace(/[^a-z0-9]+/g, "-")
        .replace(/(^-|-$)/g, "")}
      {...props}
    >
      {children}
    </h4>
  ),
  p: ({ children, ...props }: any) => (
    <p className="mb-4 leading-7 text-default-700" {...props}>
      {children}
    </p>
  ),
  a: ({ children, href, ...props }: any) => (
    <Link href={href} isExternal={href?.startsWith("http")} {...props}>
      {children}
    </Link>
  ),
  code: ({ children, className, ...props }: any) => {
    const isInline = !className;

    if (isInline) {
      return <Code size="sm">{children}</Code>;
    }

    return (
      <Code
        className="block p-4 mb-4 whitespace-pre overflow-x-auto"
        size="sm"
        {...props}
      >
        {children}
      </Code>
    );
  },
  pre: ({ children, ...props }: any) => (
    <div className="mb-4 rounded-lg overflow-hidden bg-default-100">
      <pre className="p-4 overflow-x-auto" {...props}>
        {children}
      </pre>
    </div>
  ),
  ul: ({ children, ...props }: any) => (
    <ul className="list-disc list-inside mb-4 space-y-2 ml-4" {...props}>
      {children}
    </ul>
  ),
  ol: ({ children, ...props }: any) => (
    <ol className="list-decimal list-inside mb-4 space-y-2 ml-4" {...props}>
      {children}
    </ol>
  ),
  li: ({ children, ...props }: any) => (
    <li className="leading-7 text-default-700" {...props}>
      {children}
    </li>
  ),
  blockquote: ({ children, ...props }: any) => (
    <blockquote
      className="border-l-4 border-primary pl-4 my-4 italic text-default-600"
      {...props}
    >
      {children}
    </blockquote>
  ),
  hr: () => <Divider className="my-8" />,
  table: ({ children, ...props }: any) => (
    <div className="mb-6 overflow-x-auto">
      <Table aria-label="Documentation table" {...props}>
        {children}
      </Table>
    </div>
  ),
  thead: ({ children, ...props }: any) => (
    <TableHeader {...props}>{children}</TableHeader>
  ),
  tbody: ({ children, ...props }: any) => (
    <TableBody {...props}>{children}</TableBody>
  ),
  tr: ({ children, ...props }: any) => (
    <TableRow {...props}>{children}</TableRow>
  ),
  th: ({ children, ...props }: any) => (
    <TableColumn {...props}>{children}</TableColumn>
  ),
  td: ({ children, ...props }: any) => (
    <TableCell {...props}>{children}</TableCell>
  ),
  strong: ({ children, ...props }: any) => (
    <strong className="font-semibold text-default-900" {...props}>
      {children}
    </strong>
  ),
  em: ({ children, ...props }: any) => (
    <em className="italic" {...props}>
      {children}
    </em>
  ),
};
