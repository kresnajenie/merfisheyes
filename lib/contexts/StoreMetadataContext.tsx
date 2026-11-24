"use client";

import type { StoreMetadata } from "@/lib/spatialdata/metadata";

import { createContext, useContext } from "react";

/**
 * The metadata of the store a viewer (or a split panel) shows, as
 * `readStoreMetadata` read it from the store's root attributes. Provided by
 * the page or panel that opened the store; null where nothing was opened
 * from a store, or the store carries none — the info button then hides.
 */
export const StoreMetadataContext = createContext<StoreMetadata | null>(null);

export const useStoreMetadata = () => useContext(StoreMetadataContext);
