# Releasing MERFISHEYES

How we number and cut releases. **Agents: follow this whenever you're asked to
"cut a release", "make a release", "bump the version", or promote `develop` → `main`.**

---

## Versioning

We use [Semantic Versioning](https://semver.org/): **`vMAJOR.MINOR.PATCH`** (e.g. `v0.4.2`).

For this app the bump means:

| Bump | When | Examples |
|------|------|----------|
| **MAJOR** (`1.0.0` → `2.0.0`) | Breaking change to the UX, URLs, data format, or stored dataset layout that isn't backward compatible. | Chunk format change that invalidates existing datasets; removing a viewer route. |
| **MINOR** (`0.1.0` → `0.2.0`) | New user-facing feature, backward compatible. | Server-side ingestion, split screen, a new plot type. |
| **PATCH** (`0.1.0` → `0.1.1`) | Bug fixes and small tweaks only, no new features. | Fixing the column-confirm modal, a rendering bug. |

We are pre-`1.0`: the API/format is still evolving, so breaking changes bump
**MINOR** while MAJOR stays `0`. Bump to `1.0.0` when we decide the platform is
"stable".

The version lives in three places that must stay in sync: the **git tag** on
`main`, `package.json`'s `version`, and the top entry in
[`CHANGELOG.md`](../CHANGELOG.md).

---

## Branch model

- Feature branches → PR into **`develop`** (staging: dev.merfisheyes.com).
- **A release is a promotion of `develop` → `main`** (production:
  www.merfisheyes.com). Every `main` merge that ships user-visible change gets a
  version tag.

---

## Cutting a release

Do this on an up-to-date `develop` that's green in CI.

1. **Pick the version.** Look at what's landed on `develop` since the last tag
   (`git log $(git describe --tags --abbrev=0)..develop --oneline`) and choose the
   bump per the table above.

2. **Update `CHANGELOG.md`.** Rename the `## [Unreleased]` section to
   `## [X.Y.Z] - YYYY-MM-DD` (today's date), leave a fresh empty `## [Unreleased]`
   above it, and add/refresh the `Added` / `Changed` / `Fixed` / `Removed`
   subsections. Update the compare/tag links at the bottom.

3. **Bump `package.json`** `version` to `X.Y.Z` (no `v` prefix there).

   ```bash
   npm version X.Y.Z --no-git-tag-version   # edits package.json only
   ```

4. **Commit** to `develop`:

   ```bash
   git commit -am "chore(release): vX.Y.Z"
   ```

5. **Promote to `main`.** Open a PR `develop` → `main` titled `Release vX.Y.Z`,
   let CI pass, and merge it (merge commit, matching existing history).

6. **Tag `main`** at the merge commit and push the tag:

   ```bash
   git checkout main && git pull
   git tag -a vX.Y.Z -m "vX.Y.Z"
   git push origin vX.Y.Z
   ```

7. **Create the GitHub Release** from the tag, with notes:

   ```bash
   gh release create vX.Y.Z --title "vX.Y.Z" --generate-notes
   ```

   `--generate-notes` builds the notes from merged PRs; paste the curated
   `CHANGELOG.md` highlights on top.

8. **If the worker or `scripts/process_*.py` changed**, rebuild and push the ECR
   image so Fargate runs the new code — see
   [server-side-ingestion §12](server-side-ingestion/README.md#12-deploying-a-worker-code-change).
   This is separate from the app deploy.

---

## Rules for agents

- **Never tag, push a tag, create a GitHub Release, or open a `develop` → `main`
  PR unless the user explicitly asks.** Setting the version in files (steps 2–4)
  is fine to prepare; the outward-facing steps (5–7) need an explicit go-ahead.
  This mirrors the "never commit or push unless asked" rule.
- Keep the git tag, `package.json` version, and `CHANGELOG.md` top entry
  identical.
- One tag per `main` release; don't retag. A mistake gets a new PATCH.
- The tag is `vX.Y.Z` (with `v`); `package.json` is `X.Y.Z` (without).
