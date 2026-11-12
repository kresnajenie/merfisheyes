import { API_BASE, NAME, VERSION, siteHeaders } from "./config.js";
import { accessToken } from "./auth.js";

/**
 * The one way tools call the MERFISHEYES agent API (/api/agent/v1,
 * MCP-PROTOCOL §4): the access token when signed in, one refresh on a 401,
 * a timeout, and every failure turned into a message that says what to do.
 */

export class ApiError extends Error {
  constructor(
    readonly status: number,
    message: string,
  ) {
    super(message);
  }
}

function explain(
  status: number,
  body: { error?: string; message?: string },
  headers: Headers,
  signedIn: boolean,
): string {
  const said = body.message ? ` (${body.message})` : "";

  switch (status) {
    case 401:
      return signedIn
        ? "Your MERFISHEYES sign-in has expired or was disconnected. Call meyes_sign_in, open the link it returns, then try again."
        : "This needs you to be signed in to MERFISHEYES. Call meyes_sign_in, open the link it returns, then try again.";
    case 403: {
      const scope = /scope="([^"]+)"/.exec(
        headers.get("www-authenticate") ?? "",
      )?.[1];

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

export async function apiGet<T>(
  path: string,
  params: Record<string, string | number | undefined>,
  signal?: AbortSignal,
): Promise<T> {
  const url = new URL(`${API_BASE}${path}`);

  for (const [k, v] of Object.entries(params))
    if (v !== undefined && v !== "") url.searchParams.set(k, String(v));
  const send = async (token: string | null) => {
    try {
      return await fetch(url, {
        headers: {
          accept: "application/json",
          "user-agent": `${NAME}/${VERSION}`,
          ...siteHeaders(),
          ...(token ? { authorization: `Bearer ${token}` } : {}),
        },
        signal: AbortSignal.any([
          AbortSignal.timeout(30_000),
          ...(signal ? [signal] : []),
        ]),
      });
    } catch (e) {
      if (e instanceof DOMException && e.name === "TimeoutError") {
        throw new ApiError(
          504,
          "MERFISHEYES did not answer within 30 s. Try again, or narrow the request.",
        );
      }
      throw new ApiError(
        0,
        `Could not reach MERFISHEYES at ${API_BASE}: ${(e as Error).message}`,
      );
    }
  };

  let token = await accessToken();
  let res = await send(token);

  if (res.status === 401 && token) {
    token = await accessToken(true);
    if (token) res = await send(token);
  }
  if (!res.ok) {
    const body = (await res.json().catch(() => ({}))) as {
      error?: string;
      message?: string;
    };

    throw new ApiError(
      res.status,
      explain(res.status, body, res.headers, !!token),
    );
  }

  return (await res.json()) as T;
}
