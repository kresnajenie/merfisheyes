import { describe, expect, it, vi } from "vitest";

import { keepErrorsOpen } from "@/lib/utils/keep-errors-open";

describe("errors stay until they are closed", () => {
  it("whatever a toast.error call asked for", () => {
    const original = vi.fn();
    const toast = { error: original, update: vi.fn() };

    keepErrorsOpen(toast as never);
    toast.error("Upload failed", { autoClose: 3000 });

    expect(original).toHaveBeenCalledWith("Upload failed", {
      autoClose: false,
    });
  });

  it("when a loading toast is updated to an error", () => {
    const original = vi.fn();
    const toast = { error: vi.fn(), update: original };

    keepErrorsOpen(toast as never);
    toast.update("load-gene", {
      render: "Failed to load gene",
      type: "error",
      isLoading: false,
      autoClose: 3000,
    });

    expect(original.mock.calls[0][1]).toMatchObject({
      type: "error",
      autoClose: false,
    });
  });

  it("but leaves other updates as they asked", () => {
    const original = vi.fn();
    const toast = { error: vi.fn(), update: original };

    keepErrorsOpen(toast as never);
    toast.update("load-gene", { type: "success", autoClose: 1500 });

    expect(original.mock.calls[0][1]).toEqual({
      type: "success",
      autoClose: 1500,
    });
  });
});
