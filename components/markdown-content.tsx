"use client";

import ReactMarkdown from "react-markdown";
import remarkGfm from "remark-gfm";
import { Code, Link, Divider } from "@heroui/react";

interface MarkdownContentProps {
  content: string;
}

export function MarkdownContent({ content }: MarkdownContentProps) {
  return (
    <ReactMarkdown
      components={{
        h1: ({ children, ...props }: any) => {
          const id = children
            ?.toString()
            .toLowerCase()
            .replace(/[^a-z0-9]+/g, "-")
            .replace(/(^-|-$)/g, "");

          return (
            <h1
              className="text-4xl font-bold mb-6 mt-8 scroll-mt-20"
              id={id}
              {...props}
            >
              {children}
            </h1>
          );
        },
        h2: ({ children, ...props }: any) => {
          const id = children
            ?.toString()
            .toLowerCase()
            .replace(/[^a-z0-9]+/g, "-")
            .replace(/(^-|-$)/g, "");

          return (
            <h2
              className="text-3xl font-bold mb-4 mt-8 scroll-mt-20 border-b border-divider pb-2"
              id={id}
              {...props}
            >
              {children}
            </h2>
          );
        },
        h3: ({ children, ...props }: any) => {
          const id = children
            ?.toString()
            .toLowerCase()
            .replace(/[^a-z0-9]+/g, "-")
            .replace(/(^-|-$)/g, "");

          return (
            <h3
              className="text-2xl font-bold mb-3 mt-6 scroll-mt-20"
              id={id}
              {...props}
            >
              {children}
            </h3>
          );
        },
        h4: ({ children, ...props }: any) => {
          const id = children
            ?.toString()
            .toLowerCase()
            .replace(/[^a-z0-9]+/g, "-")
            .replace(/(^-|-$)/g, "");

          return (
            <h4
              className="text-xl font-semibold mb-2 mt-4 scroll-mt-20"
              id={id}
              {...props}
            >
              {children}
            </h4>
          );
        },
        p: ({ children, ...props }: any) => (
          <p className="mb-4 leading-7 text-default-700" {...props}>
            {children}
          </p>
        ),
        a: ({ children, href, ...props }: any) => (
          <Link
            className="text-primary hover:underline"
            href={href}
            isExternal={href?.startsWith("http")}
            {...props}
          >
            {children}
          </Link>
        ),
        code: ({ inline, children, className, ...props }: any) => {
          if (inline) {
            return (
              <Code className="px-1.5 py-0.5" size="sm">
                {children}
              </Code>
            );
          }

          return (
            <Code
              className="block p-4 my-4 whitespace-pre overflow-x-auto bg-default-100 dark:bg-default-50 rounded-lg"
              size="sm"
              {...props}
            >
              {children}
            </Code>
          );
        },
        pre: ({ children, ...props }: any) => (
          <pre className="p-4 my-4 overflow-x-auto rounded-lg bg-default-100 dark:bg-default-50" {...props}>
            {children}
          </pre>
        ),
        ul: ({ children, ...props }: any) => (
          <ul className="list-disc list-outside mb-4 space-y-2 ml-6" {...props}>
            {children}
          </ul>
        ),
        ol: ({ children, ...props }: any) => (
          <ol
            className="list-decimal list-outside mb-4 space-y-2 ml-6"
            {...props}
          >
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
            className="border-l-4 border-primary bg-primary-50/50 dark:bg-primary-900/10 rounded-r-lg pl-4 pr-4 py-3 my-4 italic text-default-700"
            {...props}
          >
            {children}
          </blockquote>
        ),
        hr: () => <Divider className="my-8" />,
        table: ({ children, ...props }: any) => (
          <div className="my-6 overflow-x-auto rounded-lg border border-divider">
            <table className="min-w-full border-collapse" {...props}>
              {children}
            </table>
          </div>
        ),
        thead: ({ children, ...props }: any) => (
          <thead className="bg-default-100 dark:bg-default-50" {...props}>
            {children}
          </thead>
        ),
        tbody: ({ children, ...props }: any) => (
          <tbody {...props}>{children}</tbody>
        ),
        tr: ({ children, ...props }: any) => (
          <tr className="border-b border-divider" {...props}>
            {children}
          </tr>
        ),
        th: ({ children, ...props }: any) => (
          <th
            className="px-4 py-3 text-left font-semibold text-sm text-default-900"
            {...props}
          >
            {children}
          </th>
        ),
        td: ({ children, ...props }: any) => (
          <td className="px-4 py-3 text-sm text-default-700" {...props}>
            {children}
          </td>
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
      }}
      remarkPlugins={[remarkGfm]}
    >
      {content}
    </ReactMarkdown>
  );
}
