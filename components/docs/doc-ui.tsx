import type { ReactNode } from "react";

import NextLink from "next/link";
import { Chip } from "@heroui/chip";

import { subtitle, title } from "@/components/primitives";

/** Page heading: a title and a one-line lead. */
export function DocHeader({
  heading,
  lead,
}: {
  heading: string;
  lead: string;
}) {
  return (
    <header className="flex flex-col gap-2">
      <h1 className={title({ size: "sm" })}>{heading}</h1>
      <p
        className={subtitle({ fullWidth: true, class: "text-base lg:text-lg" })}
      >
        {lead}
      </p>
    </header>
  );
}

/** A titled section, linkable by its id. */
export function DocSection({
  id,
  heading,
  children,
}: {
  id: string;
  heading: string;
  children: ReactNode;
}) {
  return (
    <section className="flex scroll-mt-28 flex-col gap-4" id={id}>
      <h2 className="text-xl font-semibold tracking-tight">{heading}</h2>
      {children}
    </section>
  );
}

/** Body text, at a comfortable reading size. */
export function Prose({ children }: { children: ReactNode }) {
  return (
    <div className="flex flex-col gap-3 text-sm leading-relaxed text-default-700 [&_code]:rounded [&_code]:bg-default-100 [&_code]:px-1 [&_code]:py-0.5 [&_code]:text-xs [&_a]:text-primary [&_a:hover]:underline">
      {children}
    </div>
  );
}

/** A link to another place in the docs, styled as a link. */
export function DocLink({
  href,
  children,
}: {
  href: string;
  children: ReactNode;
}) {
  return (
    <NextLink className="text-primary hover:underline" href={href}>
      {children}
    </NextLink>
  );
}

/** Several links in a row, separated by dots. */
export function LinkList({
  links,
}: {
  links: { label: string; href: string }[];
}) {
  return (
    <span className="flex flex-wrap gap-x-2 gap-y-1">
      {links.map((l, i) => (
        <span key={l.href} className="flex gap-2">
          {i > 0 && (
            <span aria-hidden className="text-default-400">
              ·
            </span>
          )}
          <DocLink href={l.href}>{l.label}</DocLink>
        </span>
      ))}
    </span>
  );
}

/** A highlighted aside: a tip, or something to watch out for. */
export function Callout({
  kind = "tip",
  id,
  children,
}: {
  kind?: "tip" | "note";
  /** To link straight to it. */
  id?: string;
  children: ReactNode;
}) {
  return (
    <div
      id={id}
      className={`scroll-mt-28 rounded-xl border px-4 py-3 text-sm leading-relaxed ${
        kind === "tip"
          ? "border-primary/30 bg-primary/5"
          : "border-warning/30 bg-warning/5"
      } [&_code]:rounded [&_code]:bg-default-100 [&_code]:px-1 [&_code]:text-xs [&_a]:text-primary`}
    >
      {children}
    </div>
  );
}

const card = "rounded-xl border border-default-200 bg-default-50/40";

function Arrow({ down = false }: { down?: boolean }) {
  return (
    <span
      aria-hidden
      className={`flex shrink-0 items-center justify-center text-default-400 ${
        down ? "h-6" : "h-6 md:h-auto md:w-6"
      }`}
    >
      <span className={down ? "" : "md:hidden"}>↓</span>
      {!down && <span className="hidden md:inline">→</span>}
    </span>
  );
}

export interface FlowStep {
  title: string;
  detail?: ReactNode;
  href?: string;
}

/**
 * A short pipeline, left to right on wide screens and top to bottom on
 * narrow ones. For up to four steps; longer sequences read better as Steps.
 */
export function Flow({ steps }: { steps: FlowStep[] }) {
  return (
    <ol className="flex flex-col md:flex-row md:items-stretch">
      {steps.map((step, i) => {
        const body = (
          <>
            <span className="text-xs font-medium text-default-500">
              {i + 1}
            </span>
            <span className="font-semibold">{step.title}</span>
            {step.detail && (
              <span className="text-xs leading-relaxed text-default-500">
                {step.detail}
              </span>
            )}
          </>
        );

        return (
          <li key={step.title} className="flex flex-1 flex-col md:flex-row">
            {step.href ? (
              <NextLink
                className={`${card} flex flex-1 flex-col gap-1 p-4 transition-colors hover:border-primary/50`}
                href={step.href}
              >
                {body}
              </NextLink>
            ) : (
              <div className={`${card} flex flex-1 flex-col gap-1 p-4`}>
                {body}
              </div>
            )}
            {i < steps.length - 1 && <Arrow />}
          </li>
        );
      })}
    </ol>
  );
}

/** A numbered vertical timeline, for step-by-step instructions. */
export function Steps({
  steps,
}: {
  steps: { title: string; children: ReactNode }[];
}) {
  return (
    <ol className="flex flex-col">
      {steps.map((step, i) => (
        <li key={step.title} className="relative flex gap-4 pb-6 last:pb-0">
          {i < steps.length - 1 && (
            <span
              aria-hidden
              className="absolute left-4 top-9 bottom-0 w-px bg-default-200"
            />
          )}
          <span className="z-[1] flex h-8 w-8 shrink-0 items-center justify-center rounded-full border border-primary/40 bg-primary/10 text-sm font-semibold text-primary">
            {i + 1}
          </span>
          <div className="flex min-w-0 flex-1 flex-col gap-2 pt-1">
            <h3 className="font-semibold">{step.title}</h3>
            <Prose>{step.children}</Prose>
          </div>
        </li>
      ))}
    </ol>
  );
}

/**
 * A chain of yes/no questions: each "yes" leads to an answer, each "no" to
 * the next question, and the last "no" to a final answer.
 */
export function Decision({
  questions,
  otherwise,
}: {
  questions: { question: string; answer: string; href?: string }[];
  otherwise: { answer: string; href?: string };
}) {
  const answer = (text: string, href?: string) =>
    href ? (
      <NextLink
        className="inline-block rounded-lg border border-success/40 bg-success/10 px-3 py-2 text-sm font-medium text-success-700 hover:border-success dark:text-success"
        href={href}
      >
        {text}
      </NextLink>
    ) : (
      <span className="inline-block rounded-lg border border-success/40 bg-success/10 px-3 py-2 text-sm font-medium text-success-700 dark:text-success">
        {text}
      </span>
    );

  return (
    <div className="flex flex-col">
      {questions.map((q) => (
        <div key={q.question} className="flex flex-col">
          <div className="flex flex-col gap-2 md:flex-row md:items-center">
            <div className={`${card} flex-1 px-4 py-3 text-sm font-medium`}>
              {q.question}
            </div>
            <span className="flex items-center gap-2 text-xs text-default-500 md:w-36 md:shrink-0">
              <span aria-hidden>→</span> yes
            </span>
            <div className="md:w-60 md:shrink-0">
              {answer(q.answer, q.href)}
            </div>
          </div>
          <span
            aria-hidden
            className="flex h-8 items-center gap-2 pl-6 text-xs text-default-500"
          >
            ↓ no
          </span>
        </div>
      ))}
      <div className="md:pl-[calc(100%-15rem)]">
        {answer(otherwise.answer, otherwise.href)}
      </div>
    </div>
  );
}

/** A status as the app shows it, so the guide matches what users see. */
export function Status({
  color,
  children,
}: {
  color: "warning" | "success" | "danger" | "default";
  children: ReactNode;
}) {
  return (
    <Chip color={color} size="sm" variant="flat">
      {children}
    </Chip>
  );
}

/** A comparison table with a label column. */
export function CompareTable({
  headers,
  rows,
  labelHeader,
}: {
  headers: string[];
  rows: ReactNode[][];
  /** Heading of the label column, if it needs one. */
  labelHeader?: string;
}) {
  return (
    <div className={`${card} overflow-x-auto`}>
      <table className="w-full min-w-[36rem] text-sm">
        <thead>
          <tr>
            <th className="w-40 border-b border-default-200 px-3 py-2 text-left text-xs font-medium text-default-500">
              {labelHeader}
            </th>
            {headers.map((h) => (
              <th
                key={h}
                className="border-b border-default-200 px-3 py-2 text-left font-semibold"
              >
                {h}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {rows.map(([label, ...cells], r) => (
            <tr key={r} className="align-top">
              <td className="border-b border-default-100 px-3 py-2 text-xs font-medium text-default-500">
                {label}
              </td>
              {cells.map((cell, c) => (
                <td
                  key={c}
                  className="border-b border-default-100 px-3 py-2 text-default-700"
                >
                  {cell}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

/** A grid of cards linking to other pages. */
export function LinkCards({
  links,
}: {
  links: { title: string; detail: string; href: string }[];
}) {
  return (
    <div className="grid gap-3 sm:grid-cols-2">
      {links.map((l) => (
        <NextLink
          key={l.href}
          className={`${card} flex flex-col gap-1 p-4 transition-colors hover:border-primary/50`}
          href={l.href}
        >
          <span className="font-semibold">{l.title} →</span>
          <span className="text-xs leading-relaxed text-default-500">
            {l.detail}
          </span>
        </NextLink>
      ))}
    </div>
  );
}
