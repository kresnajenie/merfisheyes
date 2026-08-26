import { describe, it, expect } from "vitest";

import {
  computeTierMemoryGb,
  nextComputeTier,
  normalizeComputeTier,
} from "@/lib/batch";
import { planOomEscalation } from "@/lib/ingest/escalation";
import { OUT_OF_MEMORY_LABEL } from "@/lib/ingest/error-classification";

describe("compute tiers", () => {
  it("walks the escalation ladder standard → large → xlarge → none", () => {
    expect(nextComputeTier("standard")).toBe("large");
    expect(nextComputeTier("large")).toBe("xlarge");
    expect(nextComputeTier("xlarge")).toBeNull();
  });

  it("treats unknown/legacy values as standard", () => {
    expect(normalizeComputeTier(null)).toBe("standard");
    expect(normalizeComputeTier("huge")).toBe("standard");
    expect(nextComputeTier(undefined)).toBe("large");
  });

  it("reports tier memory in GB for messages", () => {
    expect(computeTierMemoryGb("standard")).toBe(16);
    expect(computeTierMemoryGb("large")).toBe(32);
    expect(computeTierMemoryGb("xlarge")).toBe(64);
  });
});

describe("planOomEscalation", () => {
  it("escalates an OOM one tier up", () => {
    expect(
      planOomEscalation({
        label: OUT_OF_MEMORY_LABEL,
        currentTier: "standard",
        attempts: 1,
      }),
    ).toEqual({ nextTier: "large" });
  });

  it("does not escalate non-OOM failures", () => {
    expect(
      planOomEscalation({
        label: "Needs review",
        currentTier: "standard",
        attempts: 1,
      }),
    ).toBeNull();
  });

  it("stops at the top tier", () => {
    expect(
      planOomEscalation({
        label: OUT_OF_MEMORY_LABEL,
        currentTier: "xlarge",
        attempts: 2,
      }),
    ).toBeNull();
  });

  it("caps automatic attempts", () => {
    expect(
      planOomEscalation({
        label: OUT_OF_MEMORY_LABEL,
        currentTier: "standard",
        attempts: 3,
      }),
    ).toBeNull();
  });
});
