/**
 * What a tool needs from wherever it runs (MCP.md D23: one implementation of
 * each tool, two servers). The npm package on a person's computer reaches the
 * site over the agent API and can open their browser; the site's own MCP
 * server (`/mcp`) calls the same API in process and opens nothing.
 */
export interface Host {
  /** The MERFISHEYES site the tools describe and link to, no trailing slash. */
  site: string;
  /** A GET on the agent API (`/api/agent/v1` + path), JSON back; failures
   * throw an Error whose message says what to do. */
  get<T>(
    path: string,
    params: Record<string, string | number | undefined>,
    signal?: AbortSignal,
  ): Promise<T>;
  /** Open a link in the person's browser; only where there is one (the
   * package). Returns whether a browser was likely opened. */
  open?: (url: string) => boolean;
  /** Said after "Not signed in" (where the package keeps sign-ins, how to
   * sign in there). */
  signInHint?: string;
}
