/**
 * Tell the server a save uploaded as it is has finished, so the signed-in
 * user gets one "ready" email for it (`/api/uploads/ready`). `datasetIds` are
 * every dataset the save made, the one to open first first. Never throws: the
 * save is done, and an email that does not get through must not undo that.
 */
export async function reportUploadReady(datasetIds: string[]): Promise<void> {
  if (!datasetIds.length) return;
  try {
    await fetch("/api/uploads/ready", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      // A save finished as the page is left still gets its email
      keepalive: true,
      body: JSON.stringify({ datasetIds }),
    });
  } catch {
    // Nothing more to do: the save itself succeeded
  }
}
