import { revokeToken } from "@/lib/oauth/flow";
import { formParams, oauthJson, preflight } from "@/lib/oauth/http";

// POST /api/oauth/revoke — RFC 7009. Always 200: whether a token existed is
// not revealed.
export async function POST(req: Request) {
  const token = (await formParams(req)).get("token");

  if (token) await revokeToken(token);

  return oauthJson({});
}

export const OPTIONS = preflight;
