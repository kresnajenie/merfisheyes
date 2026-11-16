# ADR: technology presets — where they live, who writes them, where the editor is

Status: accepted, 2026-11-16.

## Context

A technology preset is what a dataset of one technology opens with in the
viewer (colormap, cell opacity, hide zero counts; more later). The first
version kept them in code (`lib/config/technology-presets.ts`), so changing one
needed a deploy. The owner of SpatialData-DB should manage them himself, and
others may be given the same permission later. A separate repository for the
editor, perhaps together with the admin pages, was considered.

What the code has today (read 2026-11-16, develop `df4c2b7`):

- Sign-in is NextAuth with a JWT session; the role (`USER`, `ADMIN`,
  `SUPER_ADMIN`) is re-read from the database on every request
  (`lib/auth.ts`). `requireAdmin()` (`lib/admin-auth.ts`) guards all twelve
  routes under `app/api/admin/` and `app/admin/layout.tsx` guards the pages.
  There is no permission narrower than `ADMIN`.
- Admin is nine pages and twelve routes: catalog entries, community review,
  processing failures and retries, users. Every one reads or writes tables the
  main app also uses, through the same Prisma client.
- The app is an OAuth 2.1 server (`app/api/oauth/*`, `lib/oauth/`), with
  scopes checked by `hasScope()` (`lib/agent/caller.ts`). That is the existing
  way for a program outside the app to act for a user.

## Decision

1. **Presets are data.** Two append-only tables: `technology_preset_versions`
   (a config, its content hash, who saved it, the app version and the store it
   was tried on) and `technology_preset_publishes` (which version a technology
   uses from when, and who set it). The preset in force for a technology is its
   latest publish row; rolling back is publishing an earlier version; turning a
   preset off is a publish row with no version (the technology then opens with
   the viewer's ordinary defaults). Nothing is updated or deleted,
   so the history is the tables themselves.
2. **One contract.** `docs/specs/schemas/preset.schema.json` (JSON Schema
   2020-12, as the other schemas here). The write route refuses a config that
   does not match it (Ajv). The viewer applies the fields it knows and ignores
   the rest, so a newer preset does not break an older viewer.
3. **The built-in table stays as the fallback.** A technology never
   published uses `TECHNOLOGY_PRESETS` in code; so does a viewer that cannot reach
   the API. The viewer never fails to open because presets could not be read.
4. **A permission, not a role.** `users.permissions` (a list of strings) with
   one value so far, `presets:write`. Admins have it implicitly. It is granted
   by a super admin on the Users page. It gives nothing under `/api/admin`.
5. **The editor is a page of this app, `/presets`,** outside `/admin`, guarded
   by the permission.
6. **A draft is previewed by its hash.** `?preset=<hash>` on any viewer link
   makes that version the preset for that open, published or not. Versions are
   readable by hash without sign-in: a config holds nothing private, and a
   preview link has to work for whoever it is sent to.

## The separate repository, costed

| | Page in this app (chosen) | Separate repo, own database access | Separate repo, client of this app's API |
|---|---|---|---|
| Sign-in | the existing session | a second NextAuth setup on the same `users` table, second OAuth client at Google, second cookie domain | OAuth against this app; needs a `presets:write` scope and a consent step |
| Database | one Prisma schema | two copies of `prisma/schema.prisma` that must migrate together; a deploy of one can break the other | none |
| Deploys | one | two, ordered | two, independent |
| Preview | same origin as the viewer | cross-origin link | cross-origin link |
| New code | routes + one page | the same, plus auth and schema sync | the same, plus the OAuth scope and a token-handling client |

A separate repo with its own database access is refused: two programs
migrating one production database is the failure this project's
production-database rules exist to prevent. A separate repo as an API client
is sound, and the routes added here (`/api/presets*`) are the API it would
use; what it would add is the OAuth scope (`lib/oauth/config.ts`) and
`resolveCaller()` in place of the session check. It is not built because it
buys nothing for one form used by one person.

**Splitting all of admin out** is a different question with the same answer
for now. Admin's routes are not a client of anything: they are the write side
of the catalog, community and ingestion tables. Moving the pages without the
routes gives a second deploy and no separation; moving the routes means the
second-Prisma-client problem above, or first turning each into an API the main
app serves. The step that makes a split possible is the one taken here for
presets: a narrow permission and a route per capability. When admin is worth
splitting, do the same per area (catalog, community, failures), then the
pages can move to a client repo one area at a time.

## Consequences

- A publish takes effect on the next open of a dataset; no deploy.
- The viewer makes one more request per open (`GET /api/presets`, cached five
  minutes by the CDN).
- Only technologies the viewer can recognise from the data can have a preset
  that applies: today Visium HD (square bins). The editor offers only those.
  Others need the catalog's platform passed to the viewer.
- A link made while one preset was published can look different after another
  is published, for the fields the link does not carry. That is the same
  behaviour as before, and what a default means.
