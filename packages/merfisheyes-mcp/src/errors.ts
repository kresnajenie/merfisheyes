/**
 * What an agent API failure means and what to do about it, in words an agent
 * can act on — one text for both servers (MCP.md D23). `remote` is the
 * site's own `/mcp`, where sign-in belongs to the connector, not to a tool.
 */
export function explainStatus(
  status: number,
  body: { error?: string; message?: string },
  headers: Headers,
  opts: { signedIn: boolean; remote?: boolean },
): string {
  const said = body.message ? ` (${body.message})` : "";

  switch (status) {
    case 401:
      if (opts.remote)
        return "The MERFISHEYES connector's sign-in has expired or was disconnected. Ask the person to reconnect MERFISHEYES in their connector settings.";

      return opts.signedIn
        ? "Your MERFISHEYES sign-in has expired or was disconnected. Call meyes_sign_in, open the link it returns, then try again."
        : "This needs you to be signed in to MERFISHEYES. Call meyes_sign_in, open the link it returns, then try again.";
    case 403: {
      const scope = /scope="([^"]+)"/.exec(
        headers.get("www-authenticate") ?? "",
      )?.[1];

      if (opts.remote)
        return `The person did not allow this when they connected MERFISHEYES${scope ? ` (it needs: ${scope})` : ""}. Ask them to reconnect it and allow it.`;

      return `Your sign-in does not allow this${scope ? ` (it needs: ${scope})` : ""}. Call meyes_sign_in with scopes=[${(
        scope ?? ""
      )
        .split(" ")
        .map((s) => `"${s}"`)
        .join(", ")}] to grant it.`;
    }
    case 404:
      return `Not found${said}. Check the id: meyes_search_datasets, meyes_list_my_datasets and meyes_list_projects return valid ones.`;
    case 429:
      return `MERFISHEYES is rate-limiting these calls. Wait ${headers.get("retry-after") ?? "60"} s, then make fewer calls (page through results instead of repeating them).`;
    default:
      return `MERFISHEYES answered HTTP ${status}${said}.`;
  }
}
