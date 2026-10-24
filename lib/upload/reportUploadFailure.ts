/**
 * Tell the server an upload failed in the browser, so the signed-in user gets
 * the reason by email and its dataset, when one was made, is marked FAILED
 * (`/api/uploads/failed`). Never throws: the error is already on screen, and
 * a report that does not get through must not hide it.
 */
export async function reportUploadFailure(failure: {
  datasetId?: string;
  title?: string;
  error: unknown;
}): Promise<void> {
  try {
    await fetch("/api/uploads/failed", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      // A failure as the page is left still gets through
      keepalive: true,
      body: JSON.stringify({
        datasetId: failure.datasetId,
        title: failure.title,
        error:
          failure.error instanceof Error
            ? failure.error.message
            : String(failure.error ?? "Upload failed"),
      }),
    });
  } catch {
    // Nothing more to do: the error is shown
  }
}
