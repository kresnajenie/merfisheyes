/**
 * IndexedDB service for persisting local dataset metadata across browser sessions.
 * Stores only metadata (not actual data) to allow re-upload prompts when revisiting URLs.
 */

const DB_NAME = "merfish-eyes-local";
const STORE_NAME = "dataset-metadata";
const DB_VERSION = 1;
const MAX_ENTRIES = 3;

export interface LocalDatasetMetadata {
  id: string;
  name: string;
  type: string; // "h5ad" | "xenium" | "merscope" | "parquet" | "csv" | "processed_chunked"
  datasetCategory: "cell" | "molecule";
  createdAt: number;
  pointCount?: number;
  geneCount?: number;
  clusterColumns?: string[];
  spatialDimensions?: number;
  moleculeCount?: number;
  uniqueGeneCount?: number;
}

function openDB(): Promise<IDBDatabase> {
  return new Promise((resolve, reject) => {
    const request = indexedDB.open(DB_NAME, DB_VERSION);

    request.onupgradeneeded = () => {
      const db = request.result;

      if (!db.objectStoreNames.contains(STORE_NAME)) {
        db.createObjectStore(STORE_NAME, { keyPath: "id" });
      }
    };
    request.onsuccess = () => resolve(request.result);
    request.onerror = () => reject(request.error);
  });
}

/**
 * Save or update local dataset metadata. Evicts oldest entries if >MAX_ENTRIES total.
 */
export async function saveLocalDatasetMeta(
  meta: LocalDatasetMetadata,
): Promise<void> {
  try {
    const db = await openDB();

    // Upsert the entry
    await new Promise<void>((resolve, reject) => {
      const tx = db.transaction(STORE_NAME, "readwrite");
      const store = tx.objectStore(STORE_NAME);

      store.put(meta);
      tx.oncomplete = () => resolve();
      tx.onerror = () => reject(tx.error);
    });

    // Evict oldest if over limit
    const all = await getAllEntries(db);

    if (all.length > MAX_ENTRIES) {
      all.sort((a, b) => a.createdAt - b.createdAt);
      const toEvict = all.slice(0, all.length - MAX_ENTRIES);

      const tx = db.transaction(STORE_NAME, "readwrite");
      const store = tx.objectStore(STORE_NAME);

      for (const entry of toEvict) {
        store.delete(entry.id);
      }
      await new Promise<void>((resolve, reject) => {
        tx.oncomplete = () => resolve();
        tx.onerror = () => reject(tx.error);
      });
    }

    db.close();
  } catch {
    // Silent fail — IndexedDB may be unavailable (private browsing)
  }
}

/**
 * Get local dataset metadata by ID. Returns null if not found.
 */
export async function getLocalDatasetMeta(
  id: string,
): Promise<LocalDatasetMetadata | null> {
  try {
    const db = await openDB();
    const result = await new Promise<LocalDatasetMetadata | undefined>(
      (resolve, reject) => {
        const tx = db.transaction(STORE_NAME, "readonly");
        const store = tx.objectStore(STORE_NAME);
        const req = store.get(id);

        req.onsuccess = () => resolve(req.result);
        req.onerror = () => reject(req.error);
      },
    );

    db.close();

    return result ?? null;
  } catch {
    return null;
  }
}

/**
 * Delete local dataset metadata by ID.
 */
export async function deleteLocalDatasetMeta(id: string): Promise<void> {
  try {
    const db = await openDB();

    await new Promise<void>((resolve, reject) => {
      const tx = db.transaction(STORE_NAME, "readwrite");
      const store = tx.objectStore(STORE_NAME);

      store.delete(id);
      tx.oncomplete = () => resolve();
      tx.onerror = () => reject(tx.error);
    });
    db.close();
  } catch {
    // Silent fail
  }
}

/**
 * Returns true if ID does NOT start with S3 prefixes (ds_ or sm_).
 * S3-uploaded datasets use server-generated IDs: ds_* (cell), sm_* (molecule).
 * Local datasets use format-based IDs: h5ad_*, xenium_*, merscope_*, etc.
 */
export function isLocalDatasetId(id: string): boolean {
  return !id.startsWith("ds_") && !id.startsWith("sm_");
}

// Internal helper
function getAllEntries(db: IDBDatabase): Promise<LocalDatasetMetadata[]> {
  return new Promise((resolve, reject) => {
    const tx = db.transaction(STORE_NAME, "readonly");
    const store = tx.objectStore(STORE_NAME);
    const req = store.getAll();

    req.onsuccess = () => resolve(req.result);
    req.onerror = () => reject(req.error);
  });
}
