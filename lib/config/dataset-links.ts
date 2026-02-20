export interface DatasetLinkConfig {
  linkColumn: string; // Cluster column to read (e.g., "batch")
  links: Record<string, string>; // Column value → SM S3 base URL
}

// Keyed by the dataset's customS3BaseUrl (from dataset.metadata.customS3BaseUrl)
export const DATASET_LINK_REGISTRY: Record<string, DatasetLinkConfig> = {
  "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/12_16_2025_all_sets_meyes":
    {
      linkColumn: "batch",
      links: {
        set1: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/12_16_2025_sm_mereyes_set1",
        set2: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/12_16_2025_sm_mereyes_set2",
        set3: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/12_16_2025_sm_mereyes_set3",
        set4: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/12_16_2025_sm_mereyes_set4",
        set5: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/12_16_2025_sm_mereyes_set5",
        set6: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/12_16_2025_sm_mereyes_set6",
        set7: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/12_16_2025_sm_mereyes_set7",
        set8: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/12_16_2025_sm_mereyes_set8",
        set9: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/12_16_2025_sm_mereyes_set9",
        set10:
          "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/12_16_2025_sm_mereyes_set10",
      },
    },
};

export function getDatasetLinkConfig(dataset: {
  metadata?: Record<string, any>;
}): DatasetLinkConfig | null {
  const s3Url = dataset.metadata?.customS3BaseUrl;

  return s3Url ? (DATASET_LINK_REGISTRY[s3Url] ?? null) : null;
}
