import { redirect } from "next/navigation";

import { prisma } from "@/lib/prisma";
import { auth } from "@/lib/auth";
import { LIFETIME, PUBLIC_SCOPE, SCOPES } from "@/lib/oauth/config";
import { sign } from "@/lib/oauth/crypto";
import { UntrustedRequestError, validateAuthorize } from "@/lib/oauth/flow";

/**
 * /oauth/authorize — where an agent sends the person to allow it
 * (docs/specs/MCP-PROTOCOL.md §6.4): sign in, then consent. The form posts
 * the validated request back signed, so nothing can be changed between this
 * page and the decision (app/api/oauth/authorize/route.ts).
 */

export const dynamic = "force-dynamic";

interface Props {
  searchParams: Promise<Record<string, string | string[] | undefined>>;
}

function Card({ children }: { children: React.ReactNode }) {
  return (
    <div className="flex justify-center px-4 py-12">
      <div className="w-full max-w-lg rounded-2xl border border-default-200 bg-content1 p-6 shadow-sm">
        {children}
      </div>
    </div>
  );
}

export default async function AuthorizePage({ searchParams }: Props) {
  const raw = await searchParams;
  const query = new URLSearchParams(
    Object.entries(raw).flatMap(([k, v]) =>
      typeof v === "string" ? [[k, v]] : [],
    ),
  );

  let validated;

  try {
    validated = await validateAuthorize(query);
  } catch (e) {
    if (!(e instanceof UntrustedRequestError)) throw e;

    return (
      <Card>
        <h1 className="text-xl font-semibold">
          This agent cannot be connected
        </h1>
        <p className="mt-3 text-default-600">{e.message}.</p>
        <p className="mt-3 text-sm text-default-500">
          Nothing was shared. If you were setting up an agent, check its
          MERFISHEYES configuration and try again.
        </p>
      </Card>
    );
  }
  if ("redirect" in validated) redirect(validated.redirect);
  const req = validated.request;

  const session = await auth();

  if (!session?.user?.id) {
    redirect(
      `/auth/signin?callbackUrl=${encodeURIComponent(`/oauth/authorize?${query.toString()}`)}`,
    );
  }
  const projects = await prisma.project.findMany({
    where: { ownerId: session.user.id },
    select: { id: true, title: true },
    orderBy: { updatedAt: "desc" },
    take: 200,
  });
  const signed = sign({
    req,
    userId: session.user.id,
    exp: Date.now() + LIFETIME.consentSeconds * 1000,
  });
  const asked = req.scopes.filter((s) => s !== PUBLIC_SCOPE);
  const needsProjects = asked.some(
    (s) => s !== "viewer:control" && s !== "export",
  );

  return (
    <Card>
      <h1 className="text-xl font-semibold">
        Allow <span className="text-primary">{req.clientName}</span> to use
        MERFISHEYES as you?
      </h1>
      <dl className="mt-4 space-y-1 text-sm text-default-600">
        {req.clientHost && (
          <div>
            <dt className="inline font-medium">Published by: </dt>
            <dd className="inline">{req.clientHost}</dd>
          </div>
        )}
        {req.clientKind === "registered" && (
          <div className="text-warning-600">
            This agent registered itself; its name is not verified.
          </div>
        )}
        <div>
          <dt className="inline font-medium">Returns to: </dt>
          <dd className="inline font-mono">{req.redirectHost}</dd>
        </div>
        {req.loopbackOnly && (
          <div className="text-warning-600">
            This is an app running on your computer. Only continue if you
            started it yourself just now.
          </div>
        )}
      </dl>

      <form
        action="/api/oauth/authorize"
        className="mt-6 space-y-5"
        method="post"
      >
        <input name="request" type="hidden" value={signed} />
        <div>
          <p className="font-medium">It will be able to:</p>
          <ul className="mt-2 list-disc space-y-1 pl-5 text-sm">
            <li>{SCOPES[PUBLIC_SCOPE]}</li>
            {asked.map((s) => (
              <li key={s}>{SCOPES[s]}</li>
            ))}
          </ul>
        </div>

        {needsProjects && (
          <fieldset>
            <legend className="font-medium">Which of your projects?</legend>
            <label className="mt-2 flex items-center gap-2 text-sm">
              <input defaultChecked name="reach" type="radio" value="all" /> All
              my projects and datasets
            </label>
            <label className="mt-1 flex items-center gap-2 text-sm">
              <input
                disabled={projects.length === 0}
                name="reach"
                type="radio"
                value="some"
              />{" "}
              Only these projects
              {projects.length === 0 && (
                <span className="text-default-400">(you have none yet)</span>
              )}
            </label>
            {projects.length > 0 && (
              <div className="ml-6 mt-2 max-h-48 space-y-1 overflow-y-auto rounded-lg border border-default-200 p-2">
                {projects.map((p) => (
                  <label key={p.id} className="flex items-center gap-2 text-sm">
                    <input name="projectIds" type="checkbox" value={p.id} />{" "}
                    {p.title}
                  </label>
                ))}
              </div>
            )}
          </fieldset>
        )}

        <div className="flex justify-end gap-3 pt-2">
          <button
            className="rounded-lg px-4 py-2 text-sm text-default-700 hover:bg-default-100"
            name="decision"
            type="submit"
            value="deny"
          >
            Cancel
          </button>
          <button
            className="rounded-lg bg-primary px-4 py-2 text-sm font-medium text-primary-foreground"
            name="decision"
            type="submit"
            value="approve"
          >
            Allow
          </button>
        </div>
        <p className="text-xs text-default-400">
          You can disconnect this agent at any time from your account page.
        </p>
      </form>
    </Card>
  );
}
