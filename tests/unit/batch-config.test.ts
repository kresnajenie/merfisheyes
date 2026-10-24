import { afterEach, describe, expect, it, vi } from "vitest";

import { isBatchConfigured } from "@/lib/batch";

describe("isBatchConfigured", () => {
  afterEach(() => vi.unstubAllEnvs());

  const configured = () => {
    vi.stubEnv("CALLBACK_SECRET", "secret");
    vi.stubEnv("INGEST_CALLBACK_BASE_URL", "https://dev.example.com");
    vi.stubEnv("AWS_S3_ENDPOINT", "");
    vi.stubEnv("S3_ENDPOINT", "");
  };

  it("submits jobs with a callback secret and base URL", () => {
    configured();
    expect(isBatchConfigured()).toBe(true);
  });

  it("never submits against a custom S3 endpoint (the local test stack)", () => {
    configured();
    vi.stubEnv("AWS_S3_ENDPOINT", "http://localhost:9000");
    expect(isBatchConfigured()).toBe(false);
  });

  it("needs the callback secret", () => {
    configured();
    vi.stubEnv("CALLBACK_SECRET", "");
    expect(isBatchConfigured()).toBe(false);
  });
});
