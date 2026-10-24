import { readFileSync } from "node:fs";

import { describe, expect, it } from "vitest";

import { releaseNotes } from "../../scripts/release-notes.mjs";

const CHANGELOG = `# Changelog

## [Unreleased]

## [0.2.0] - 2026-01-02

The summary.

### Added

- A thing.

## [0.1.0] - 2026-01-01

- The first.

[Unreleased]: https://example.com/compare/v0.2.0...HEAD
[0.2.0]: https://example.com/compare/v0.1.0...v0.2.0
`;

describe("release notes from the changelog", () => {
  it("takes one version's section, up to the next version", () => {
    expect(releaseNotes(CHANGELOG, "0.2.0")).toBe(
      "The summary.\n\n### Added\n\n- A thing.",
    );
  });

  it("stops the last section at the link references", () => {
    expect(releaseNotes(CHANGELOG, "0.1.0")).toBe("- The first.");
  });

  it("has nothing for a version that is not there", () => {
    expect(releaseNotes(CHANGELOG, "0.3.0")).toBeNull();
    // A prefix is not a version
    expect(releaseNotes(CHANGELOG, "0.2")).toBeNull();
  });

  // What the release workflow relies on when main moves: the version in
  // package.json always has its notes, so its tag can be made
  it("package.json's version has a section in CHANGELOG.md", () => {
    const { version } = JSON.parse(readFileSync("package.json", "utf8"));
    const notes = releaseNotes(readFileSync("CHANGELOG.md", "utf8"), version);

    expect(notes, `CHANGELOG.md has no "## [${version}]"`).toBeTruthy();
  });
});
