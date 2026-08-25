# Benchmark — linux-x64

**Machine:** AMD Ryzen Threadripper PRO 5975WX 32-Cores · 64 cores · 503.52 GB RAM · Chromium 149.0.7827.55

`process` = parse + build the chunked layout. Byte transfer and Fargate
provisioning are excluded from it and listed separately.

| dataset | format | size | shape | local | local peak | server | server peak |
|---|---|---:|---|---|---:|---|---:|
| MouseSpinalCord_Run8_Slide2_Mus_M7_F7_ace-pie-cot | xenium | 0.01 GB | 182,698 cells × 541 | 6.3s | 2.26 GB | 1.2m | 0.48 GB |
| MouseSpinalCord_Run3_Slide2_Mus_M2_F2_ace-pie-bud | xenium | 0.01 GB | 203,848 cells × 541 | 6.5s | 2.41 GB | 1.3m | 0.53 GB |
| MouseSpinalCord_Run5_Slide2_Mus_M4_F4_ace-pie-cab | xenium | 0.01 GB | 200,677 cells × 541 | 6.8s | 2.35 GB | 1.4m | 0.57 GB |
| MouseSpinalCord_Run6_Slide1_Mus_M5_F5_ace-pie-can | xenium | 0.01 GB | 195,045 cells × 541 | 6.6s | 2.31 GB | 1.2m | 0.52 GB |
| MouseSpinalCord_Run6_Slide2_Mus_M5_F5_ace-pie-cap | xenium | 0.01 GB | 202,387 cells × 541 | 6.5s | 2.31 GB | 1.2m | 0.51 GB |
| MouseSpinalCord_Run7_Slide1_Mus_M6_F6_ace-pie-car | xenium | 0.01 GB | 201,328 cells × 541 | 6.1s | 2.22 GB | 54.0s | 0.43 GB |
| MouseSpinalCord_Run7_Slide2_Mus_M6_F6_ace-pie-cat | xenium | 0.01 GB | 197,012 cells × 541 | 6.3s | 2.27 GB | 55.4s | 0.44 GB |
| MouseSpinalCord_Run8_Slide1_Mus_M7_F7_ace-pie-cop | xenium | 0.01 GB | 193,138 cells × 541 | 6.6s | 2.29 GB | 54.2s | 0.51 GB |
| MouseSpinalCord_Run1_Slide1_Mus_M1_F1_ace-pie-bin | xenium | 0.01 GB | 208,799 cells × 541 | 7.0s | 2.37 GB | 1.3m | 0.54 GB |
| MouseSpinalCord_Run3_Slide1_Mus_M2_F2_ace-pie-box | xenium | 0.01 GB | 214,194 cells × 541 | 6.8s | 2.39 GB | 1.3m | 0.52 GB |
| MouseSpinalCord_Run4_Slide1_Mus_M3_F3_ace-pie-bug | xenium | 0.01 GB | 214,349 cells × 541 | 6.9s | 2.39 GB | 1.4m | 0.55 GB |
| MouseSpinalCord_Run5_Slide1_Mus_M4_F4_ace-pie-bus | xenium | 0.01 GB | 211,359 cells × 541 | 7.2s | 2.49 GB | 1.1m | 0.61 GB |
| MouseSpinalCord_Run1_Slide2_Mus_M1_F1_ace-pie-bit | xenium | 0.01 GB | 226,204 cells × 541 | 7.0s | 2.52 GB | 1.5m | 0.58 GB |
| MouseSpinalCord_Run2_Slide2_Mus_M2_F2_ace-pie-boo | xenium | 0.01 GB | 239,887 cells × 541 | 7.1s | 2.57 GB | 1.3m | 0.55 GB |
| MouseSpinalCord_Run4_Slide2_Mus_M3_F3_ace-pie-bun | xenium | 0.01 GB | 238,228 cells × 541 | 7.1s | 2.63 GB | 1.6m | 0.61 GB |
| MouseSpinalCord_Run2_Slide1_Mus_M2_F2_ace-pie-bog | xenium | 0.01 GB | 258,147 cells × 541 | 7.5s | 2.68 GB | 1.4m | 0.6 GB |
| merscope-800m-ace-dry-dog-detected_transcripts3-csv | single-molecule | 0.04 GB | 499,999 molecules | 3.1s | 1.01 GB | 2.2s | 0.25 GB |
| h5ad__140MB_cbm2_combined_neuropath | h5ad | 0.13 GB | 57,489 cells × 223 | 2.1s | 1.28 GB | 9.0s | 0.27 GB |
| Human_Ovarian_Cancer_150mb | merscope | 0.16 GB | 71,381 cells × 550 | 14.2s | 3.1 GB | 42.1s | 0.83 GB |
| Mouse_receptors_rep3_1_190mb | merscope | 0.19 GB | 70,844 cells × 649 | 14.3s | 3.26 GB | 35.6s | 0.97 GB |
| Mouse_receptors_rep1_1_210mb | merscope | 0.20 GB | 78,329 cells × 649 | 16.0s | 3.2 GB | 40.7s | 1.05 GB |
| Mouse_receptors_rep2_1_220mb | merscope | 0.22 GB | 83,461 cells × 649 | 16.7s | 3.64 GB | 45.0s | 1.1 GB |
| Mouse_receptors_rep3_2_220mb | merscope | 0.22 GB | 83,461 cells × 649 | 17.1s | 3.61 GB | 44.6s | 1.11 GB |
| Mouse_receptors_rep2_2_220mb | merscope | 0.22 GB | 84,172 cells × 649 | 17.1s | 3.61 GB | 42.0s | 1.12 GB |
| Mouse_receptors_rep1_3_220mb | merscope | 0.22 GB | 84,636 cells × 649 | 17.3s | 3.68 GB | 43.1s | 1.12 GB |
| Mouse_receptors_rep2_3_230mb | merscope | 0.23 GB | 85,958 cells × 649 | 16.6s | 3.72 GB | 44.9s | 1.14 GB |
| h5ad__250MB_CBM7_labeled | h5ad | 0.23 GB | 139,997 cells × 223 | 3.4s | 1.59 GB | 9.6s | 0.37 GB |
| Mouse_receptors_rep1_2_240mb | merscope | 0.23 GB | 88,884 cells × 649 | 17.7s | 3.81 GB | 45.1s | 1.23 GB |
| hum_brain2_280mb | merscope | 0.30 GB | 137,627 cells × 1000 | **oom** | 4.66 GB | 1.7m | 2.69 GB |
| xenium-Human_Lung_Cancer_FFPE-parquet | single-molecule | 0.33 GB | 32,073,729 molecules | 55.5s | 5.37 GB | 46.8s | 2.99 GB |
| xenium-Human_Melanoma_10gb-parquet | single-molecule | 0.39 GB | 23,651,907 molecules | 39.3s | 4.54 GB | 34.6s | 2.32 GB |
| ed_lein_ace-irk-sag | merscope | 0.45 GB | 878,735 cells × 315 | **failed** | 5.58 GB | **error** | – |
| h5ad__cbm5_combined_neuropath | h5ad | 0.51 GB | 152,568 cells × 223 | 4.3s | 1.98 GB | 16.4s | 0.39 GB |
| merscope-800m-ace-dry-dog-detected_transcripts2-csv | single-molecule | 0.63 GB | 8,044,719 molecules | **timeout** | 1.6 GB | 49.6s | 1.77 GB |
| h5ad__scdata_2000gnMERFISH_DCBB_5_9_2024__set3 | h5ad | 0.75 GB | 163,800 cells × 275 | 6.8s | 3.01 GB | 26.8s | 0.57 GB |
| merscope-800m-ace-dry-dog-csv | single-molecule | 0.79 GB | 10,107,749 molecules | **timeout** | 1.66 GB | 57.2s | 2.14 GB |
| xenium-Human_Renal_Carcinoma_50gb-parquet | single-molecule | 0.95 GB | 77,056,936 molecules | **oom** | 7.35 GB | 2.0m | 7.16 GB |
| Human_Uterine_Cancer_800mb | merscope | 1.05 GB | 739,735 cells × 550 | **oom** | 4.81 GB | 9.6m | 9.54 GB |
| WTA_Preview_FFPE_Breast_Cancer | xenium | 1.40 GB | 170,057 cells × 27104 | 6.5s ⚠️ | 5.85 GB | 31.1m | 10.61 GB |
| Human_Colon_Cancer_1-2gb | merscope | 1.43 GB | 706,468 cells × 900 | **oom** | 4.69 GB | 4.8m | 14.25 GB |
| Xenium_Prime_Mouse_Pup_FFPE | xenium | 1.51 GB | 1,298,870 cells × 13780 | 9.7s ⚠️ | 2.94 GB | 14.3m | 7.84 GB |
| merscope-ed_lein_ace-irk-sag-H22-30-001-CX-13-01-0-csv | single-molecule | 1.53 GB | 18,644,913 molecules | **timeout** | 3.01 GB | 1.5m | 3.96 GB |
| h5ad__year1s_labelled_edu | h5ad | 1.81 GB | 381,811 cells × 275 | 30.0s | 5.77 GB | 1.4m | 0.88 GB |
| xenium-Human_Melanoma_10gb-csv | single-molecule | 1.93 GB | 23,390,260 molecules | **timeout** | 4.28 GB | 32.9s | 4.39 GB |
| merscope-ed_lein_ace-irk-sag-csv | single-molecule | 2.13 GB | 26,196,203 molecules | **timeout** | 0.68 GB | 2.4m | 5.36 GB |
| merscope-2g-ace-ear-nap-202204011403_2022040113532-csv | single-molecule | 2.44 GB | 32,146,915 molecules | **timeout** | 0.67 GB | 3.0m | 5.91 GB |
| h5ad__retroB30_all_BRBB_Axon_Virus_AnnotBRBBTransfer | h5ad | 2.49 GB | 3,009,652 cells × 992 | **failed** | 0.63 GB | 65.9m | 13.49 GB |
| xenium-Human_Lymph_Node-parquet | single-molecule | 2.78 GB | 232,650,139 molecules | **failed** | 0.81 GB | **oom** → ok @32 GB 5.8m (peak 20.8 GB) | – |
| WTA_Preview_FFPE_Cervical_Cancer | xenium | 2.92 GB | 717,576 cells × 27104 | 9.0s ⚠️ | 4.71 GB | **oom** → ok @32 GB 43.9m (peak 16.26 GB) | – |
| h5ad__scdata_2000gnMERFISH_DCBB_5_9_2024__set1_FULL_corrected__wit | h5ad | 2.94 GB | 93,022 cells × 2088 | **failed** | 0.64 GB | 1.6m | 1.1 GB |
| h5ad__4GB_2kset2_labelled_clusters | h5ad | 3.97 GB | 125,665 cells × 2088 | **failed** | 0.63 GB | 2.1m | 1.38 GB |
| merscope-5g-ace-dry-dog-csv | single-molecule | 5.43 GB | 69,965,035 molecules | **timeout** | 1.8 GB | 6.5m | 13.72 GB |
| xenium-Whole_mouse_pup_65gb-parquet | single-molecule | 8.01 GB | 669,251,333 molecules | **failed** | 0.88 GB | **oom** → ok @64 GB 13.5m (peak 59.44 GB) | – |
| xenium-WTA_Preview_FFPE_Breast_Cancer-parquet | single-molecule | 9.41 GB | 740,442,119 molecules | **failed** | 1.93 GB | **oom** → oom @64 GB | – |
| xenium-WTA_Preview_FFPE_Cervical_Cancer-parquet | single-molecule | 15.62 GB | 1,235,087,993 molecules | **failed** | 0.8 GB | **oom** → oom @64 GB | – |

⚠️ = finished but loaded **no real gene expression**. Neither path reads
`cell_feature_matrix.h5`, which is how these Xenium exports ship it — the
browser loads zero genes and the server publishes a 1-gene placeholder.
Xenium rows are cells-only on both paths, not full comparisons.

## Ceiling

- Largest processed **in the browser**: `h5ad__year1s_labelled_edu` — 1.81 GB, 381,811 cells × 275, 30.0s, peak 5.77 GB
- Smallest **browser crash**: `merscope__hum_brain2_280mb` — 0.30 GB, 137,627 cells × 1000
- Browser crashes peaked at 4.66–7.35 GB on a machine with 503.52 GB of RAM — the cap is per-renderer, not the hardware.
- Largest processed **on the server**: `single-molecule__merscope-5g-ace-dry-dog-csv` — 5.43 GB, 69,965,035 molecules, 6.5m, peak 13.72 GB

## Where the server wins (same datasets, both paths succeeded)

| dataset | local process | server process | local peak | server peak |
|---|---:|---:|---:|---:|
| MouseSpinalCord_Run8_Slide2_Mus_M7_F7_ace-pie-cot | 6.3s | 1.2m | 2.26 GB | 0.48 GB |
| MouseSpinalCord_Run3_Slide2_Mus_M2_F2_ace-pie-bud | 6.5s | 1.3m | 2.41 GB | 0.53 GB |
| MouseSpinalCord_Run5_Slide2_Mus_M4_F4_ace-pie-cab | 6.8s | 1.4m | 2.35 GB | 0.57 GB |
| MouseSpinalCord_Run6_Slide1_Mus_M5_F5_ace-pie-can | 6.6s | 1.2m | 2.31 GB | 0.52 GB |
| MouseSpinalCord_Run6_Slide2_Mus_M5_F5_ace-pie-cap | 6.5s | 1.2m | 2.31 GB | 0.51 GB |
| MouseSpinalCord_Run7_Slide1_Mus_M6_F6_ace-pie-car | 6.1s | 54.0s | 2.22 GB | 0.43 GB |
| MouseSpinalCord_Run7_Slide2_Mus_M6_F6_ace-pie-cat | 6.3s | 55.4s | 2.27 GB | 0.44 GB |
| MouseSpinalCord_Run8_Slide1_Mus_M7_F7_ace-pie-cop | 6.6s | 54.2s | 2.29 GB | 0.51 GB |
| MouseSpinalCord_Run1_Slide1_Mus_M1_F1_ace-pie-bin | 7.0s | 1.3m | 2.37 GB | 0.54 GB |
| MouseSpinalCord_Run3_Slide1_Mus_M2_F2_ace-pie-box | 6.8s | 1.3m | 2.39 GB | 0.52 GB |
| MouseSpinalCord_Run4_Slide1_Mus_M3_F3_ace-pie-bug | 6.9s | 1.4m | 2.39 GB | 0.55 GB |
| MouseSpinalCord_Run5_Slide1_Mus_M4_F4_ace-pie-bus | 7.2s | 1.1m | 2.49 GB | 0.61 GB |
| MouseSpinalCord_Run1_Slide2_Mus_M1_F1_ace-pie-bit | 7.0s | 1.5m | 2.52 GB | 0.58 GB |
| MouseSpinalCord_Run2_Slide2_Mus_M2_F2_ace-pie-boo | 7.1s | 1.3m | 2.57 GB | 0.55 GB |
| MouseSpinalCord_Run4_Slide2_Mus_M3_F3_ace-pie-bun | 7.1s | 1.6m | 2.63 GB | 0.61 GB |
| MouseSpinalCord_Run2_Slide1_Mus_M2_F2_ace-pie-bog | 7.5s | 1.4m | 2.68 GB | 0.6 GB |
| merscope-800m-ace-dry-dog-detected_transcripts3-csv | 3.1s | 2.2s | 1.01 GB | 0.25 GB |
| h5ad__140MB_cbm2_combined_neuropath | 2.1s | 9.0s | 1.28 GB | 0.27 GB |
| Human_Ovarian_Cancer_150mb | 14.2s | 42.1s | 3.1 GB | 0.83 GB |
| Mouse_receptors_rep3_1_190mb | 14.3s | 35.6s | 3.26 GB | 0.97 GB |
| Mouse_receptors_rep1_1_210mb | 16.0s | 40.7s | 3.2 GB | 1.05 GB |
| Mouse_receptors_rep2_1_220mb | 16.7s | 45.0s | 3.64 GB | 1.1 GB |
| Mouse_receptors_rep3_2_220mb | 17.1s | 44.6s | 3.61 GB | 1.11 GB |
| Mouse_receptors_rep2_2_220mb | 17.1s | 42.0s | 3.61 GB | 1.12 GB |
| Mouse_receptors_rep1_3_220mb | 17.3s | 43.1s | 3.68 GB | 1.12 GB |
| Mouse_receptors_rep2_3_230mb | 16.6s | 44.9s | 3.72 GB | 1.14 GB |
| h5ad__250MB_CBM7_labeled | 3.4s | 9.6s | 1.59 GB | 0.37 GB |
| Mouse_receptors_rep1_2_240mb | 17.7s | 45.1s | 3.81 GB | 1.23 GB |
| xenium-Human_Lung_Cancer_FFPE-parquet | 55.5s | 46.8s | 5.37 GB | 2.99 GB |
| xenium-Human_Melanoma_10gb-parquet | 39.3s | 34.6s | 4.54 GB | 2.32 GB |
| h5ad__cbm5_combined_neuropath | 4.3s | 16.4s | 1.98 GB | 0.39 GB |
| h5ad__scdata_2000gnMERFISH_DCBB_5_9_2024__set3 | 6.8s | 26.8s | 3.01 GB | 0.57 GB |
| WTA_Preview_FFPE_Breast_Cancer | 6.5s | 31.1m | 5.85 GB | 10.61 GB |
| Xenium_Prime_Mouse_Pup_FFPE | 9.7s | 14.3m | 2.94 GB | 7.84 GB |
| h5ad__year1s_labelled_edu | 30.0s | 1.4m | 5.77 GB | 0.88 GB |

## Excluded from `process` (reported, not hidden)

| dataset | upload | MB/s | queue (Fargate) |
|---|---:|---:|---:|
| MouseSpinalCord_Run8_Slide2_Mus_M7_F7_ace-pie-cot | 2.0s | 13 | 41.6s |
| MouseSpinalCord_Run3_Slide2_Mus_M2_F2_ace-pie-bud | 2.6s | 11 | 37.9s |
| MouseSpinalCord_Run5_Slide2_Mus_M4_F4_ace-pie-cab | 3.3s | 9.3 | 34.1s |
| MouseSpinalCord_Run6_Slide1_Mus_M5_F5_ace-pie-can | 2.1s | 12.9 | 33.8s |
| MouseSpinalCord_Run6_Slide2_Mus_M5_F5_ace-pie-cap | 2.1s | 12.9 | 39.4s |
| MouseSpinalCord_Run7_Slide1_Mus_M6_F6_ace-pie-car | 2.1s | 11.5 | 57.7s |
| MouseSpinalCord_Run7_Slide2_Mus_M6_F6_ace-pie-cat | 2.2s | 10.7 | 31.5s |
| MouseSpinalCord_Run8_Slide1_Mus_M7_F7_ace-pie-cop | 2.4s | 11.6 | 29.8s |
| MouseSpinalCord_Run1_Slide1_Mus_M1_F1_ace-pie-bin | 2.4s | 12.4 | 32.1s |
| MouseSpinalCord_Run3_Slide1_Mus_M2_F2_ace-pie-box | 3.7s | 7.8 | 39.5s |
| MouseSpinalCord_Run4_Slide1_Mus_M3_F3_ace-pie-bug | 2.3s | 13.1 | 33.0s |
| MouseSpinalCord_Run5_Slide1_Mus_M4_F4_ace-pie-bus | 2.7s | 12.2 | 35.8s |
| MouseSpinalCord_Run1_Slide2_Mus_M1_F1_ace-pie-bit | 2.1s | 15.3 | 33.5s |
| MouseSpinalCord_Run2_Slide2_Mus_M2_F2_ace-pie-boo | 3.0s | 10.4 | 38.6s |
| MouseSpinalCord_Run4_Slide2_Mus_M3_F3_ace-pie-bun | 2.6s | 13.1 | 32.1s |
| MouseSpinalCord_Run2_Slide1_Mus_M2_F2_ace-pie-bog | 2.5s | 13.6 | 44.9s |
| merscope-800m-ace-dry-dog-detected_transcripts3-csv | 2.0s | 20.5 | 36.8s |
| h5ad__140MB_cbm2_combined_neuropath | 2.8s | 50 | 34.9s |
| Human_Ovarian_Cancer_150mb | 3.8s | 44 | 32.7s |
| Mouse_receptors_rep3_1_190mb | 6.8s | 29.5 | 40.1s |
| Mouse_receptors_rep1_1_210mb | 4.0s | 55.4 | 36.3s |
| Mouse_receptors_rep2_1_220mb | 6.4s | 36.6 | 37.2s |
| Mouse_receptors_rep3_2_220mb | 8.3s | 28.2 | 34.4s |
| Mouse_receptors_rep2_2_220mb | 4.0s | 59.3 | 28.6s |
| Mouse_receptors_rep1_3_220mb | 3.4s | 70.3 | 36.4s |
| Mouse_receptors_rep2_3_230mb | 4.7s | 51.6 | 30.8s |
| h5ad__250MB_CBM7_labeled | 2.9s | 86.3 | 37.7s |
| Mouse_receptors_rep1_2_240mb | 3.5s | 71 | 40.7s |
| hum_brain2_280mb | 7.8s | 41.4 | 32.6s |
| xenium-Human_Lung_Cancer_FFPE-parquet | 3.4s | 103.7 | 40.4s |
| xenium-Human_Melanoma_10gb-parquet | 3.7s | 113.5 | 37.9s |
| ed_lein_ace-irk-sag | 11.1s | 43.6 | 1.7m |
| h5ad__cbm5_combined_neuropath | 5.5s | 100 | 28.7s |
| merscope-800m-ace-dry-dog-detected_transcripts2-csv | 6.0s | 111.5 | 42.9s |
| h5ad__scdata_2000gnMERFISH_DCBB_5_9_2024__set3 | 6.0s | 135.1 | 32.3s |
| merscope-800m-ace-dry-dog-csv | 6.2s | 135.4 | 42.3s |
| xenium-Human_Renal_Carcinoma_50gb-parquet | 9.2s | 111.4 | 37.8s |
| Human_Uterine_Cancer_800mb | 10.1s | 111.4 | 31.7s |
| WTA_Preview_FFPE_Breast_Cancer | 9.6s | 42.5 | 33.3s |
| Human_Colon_Cancer_1-2gb | 15.1s | 101.3 | 29.4s |
| Xenium_Prime_Mouse_Pup_FFPE | 6.0s | 77.2 | 42.8s |
| merscope-ed_lein_ace-irk-sag-H22-30-001-CX-13-01-0-csv | 12.9s | 127.4 | 34.8s |
| h5ad__year1s_labelled_edu | 18.8s | 103.2 | 35.7s |
| xenium-Human_Melanoma_10gb-csv | 14.8s | 140.5 | 36.8s |
| merscope-ed_lein_ace-irk-sag-csv | 19.7s | 116.1 | 36.5s |
| merscope-2g-ace-ear-nap-202204011403_2022040113532-csv | 20.0s | 131.1 | 33.9s |
| h5ad__retroB30_all_BRBB_Axon_Virus_AnnotBRBBTransfer | 23.8s | 112.6 | 29.5s |
| xenium-Human_Lymph_Node-parquet | 35.6s | 83.8 | 2.8m |
| WTA_Preview_FFPE_Cervical_Cancer | 11.3s | 79 | 2.6m |
| h5ad__scdata_2000gnMERFISH_DCBB_5_9_2024__set1_FULL_corrected__wit | 26.1s | 121 | 44.0s |
| h5ad__4GB_2kset2_labelled_clusters | 38.5s | 110.8 | 31.5s |
| merscope-5g-ace-dry-dog-csv | 47.4s | 123.2 | 40.3s |
| xenium-Whole_mouse_pup_65gb-parquet | 1.1m | 127.2 | 3.5m |
| xenium-WTA_Preview_FFPE_Breast_Cancer-parquet | 1.1m | 151.4 | 3.5m |
| xenium-WTA_Preview_FFPE_Cervical_Cancer-parquet | 3.0m | 92.7 | 5.0m |

Upload time is network-bound and varies by more than 10× between machines,
which is why it is never part of the headline number.

