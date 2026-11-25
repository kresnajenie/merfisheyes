import { API_BASE, NAME, VERSION, siteHeaders } from "./config.ts";
import { accessToken } from "./auth.ts";
import { explainStatus } from "./errors.ts";

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
      explainStatus(res.status, body, res.headers, { signedIn: !!token }),
    );
  }

  return (await res.json()) as T;
}
