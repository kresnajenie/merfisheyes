/**
 * Seed script for BIL (Brain Image Library) datasets.
 * Run with: npx tsx prisma/seed-bil-data.ts
 */
import { PrismaClient } from "@prisma/client";

const prisma = new PrismaClient();

interface SeedEntry {
  label: string;
  datasetType: "single_cell" | "single_molecule";
  s3BaseUrl?: string;
  datasetId?: string;
}

interface SeedDataset {
  title: string;
  description: string;
  species: string;
  tissue: string;
  platform: string;
  externalLink: string;
  entries: SeedEntry[];
}

const bilDatasets: SeedDataset[] = [
  {
    title: "MERSCOPE Imaging of WT, Trem2-R47H, 5xFAD and Trem2-R47H; 5xFAD samples",
    description:
      "This dataset contains MERSCOPE Imaging data of 19 coronal hemispheres in WT, Trem2-R47H, 5xFAD, and Trem2-R47H; 5xFAD mice. Sections were imaged with a 300 gene panel, combined with DAPI imaging for localization of cellular nuclei.",
    species: "Mouse",
    tissue: "Brain",
    platform: "MERFISH",
    externalLink: "https://bil.merfisheyes.com/dataset/ace-ear-nap",
    entries: [
      { label: "All Single Cell", datasetType: "single_cell", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/MERFISH_Data_corrected_orientation_meyes2" },
      { label: "TREM2_5xFAD4", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/xiangmin_ace-ear-nap-processed/ace-ear-nap/repository/Data_Repository/202204011403_20220401135320220401TREM5xADF12mo4hemip_VMSC00101/Region_0_single_molecule" },
      { label: "TREM2_4", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/xiangmin_ace-ear-nap-processed/ace-ear-nap/repository/Data_Repository/202204011403_20220401135320220401TREM5xADF12mo4hemip_VMSC00101/Region_2_single_molecule" },
      { label: "TREM2_5xFAD3", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/xiangmin_ace-ear-nap-processed/ace-ear-nap/repository/Data_Repository/202209041036_20220904-5xFAD-TREM2-AnteHip-300x05_VMSC00101/Region_0_single_molecule" },
      { label: "WT3", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/xiangmin_ace-ear-nap-processed/ace-ear-nap/repository/Data_Repository/202209041036_20220904-5xFAD-TREM2-AnteHip-300x05_VMSC00101/Region_1_single_molecule" },
      { label: "TREM2_3", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/xiangmin_ace-ear-nap-processed/ace-ear-nap/repository/Data_Repository/202209041036_20220904-5xFAD-TREM2-AnteHip-300x05_VMSC00101/Region_2_single_molecule" },
      { label: "5xFAD3", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/xiangmin_ace-ear-nap-processed/ace-ear-nap/repository/Data_Repository/202209041036_20220904-5xFAD-TREM2-AnteHip-300x05_VMSC00101/Region_3_single_molecule" },
      { label: "TREM2_5xFAD5", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/xiangmin_ace-ear-nap-processed/ace-ear-nap/repository/Data_Repository/202302031055_20230203TREM5xFAD-VZG171x04_VMSC05201/Region_0_single_molecule" },
      { label: "5xFAD5", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/xiangmin_ace-ear-nap-processed/ace-ear-nap/repository/Data_Repository/202302031055_20230203TREM5xFAD-VZG171x04_VMSC05201/Region_1_single_molecule" },
      { label: "TREM2_5", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/xiangmin_ace-ear-nap-processed/ace-ear-nap/repository/Data_Repository/202302031055_20230203TREM5xFAD-VZG171x04_VMSC05201/Region_2_single_molecule" },
      { label: "WT5", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/xiangmin_ace-ear-nap-processed/ace-ear-nap/repository/Data_Repository/202302031055_20230203TREM5xFAD-VZG171x04_VMSC05201/Region_3_single_molecule" },
      { label: "WT1", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/xiangmin_ace-ear-nap-processed/ace-ear-nap/repository/Data_Repository/Kim_202201111602_20220111-WT-5xFAD10518pHip300GP_VMSC00101/region_0_single_molecule" },
      { label: "5xFAD1", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/xiangmin_ace-ear-nap-processed/ace-ear-nap/repository/Data_Repository/Kim_202201111602_20220111-WT-5xFAD10518pHip300GP_VMSC00101/region_1_single_molecule" },
      { label: "WT2", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/xiangmin_ace-ear-nap-processed/ace-ear-nap/repository/Data_Repository/Kim_202201111602_20220111-WT-5xFAD10518pHip300GP_VMSC00101/region_2_single_molecule" },
      { label: "TREM2_2", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/xiangmin_ace-ear-nap-processed/ace-ear-nap/repository/Data_Repository/Kim2_202112171955_12172021TREM2-5x12Mo300GP_VMSC00101/Region_2_single_molecule" },
    ],
  },
  {
    title: "Spatial Transcriptomic Atlas of Adult Human Basal Ganglia: H22.30.001.CX.13",
    description:
      "This dataset is part of a whole-brain atlas effort for human and non-human primate. The data collected for this submission represent spatial transcriptomics data from core human basal ganglia regions collected using a gene panel designed to identify all known cell types identified in the snRNAseq dataset from human generated from similar regions.",
    species: "Human",
    tissue: "Brain",
    platform: "MERFISH",
    externalLink: "https://bil.merfisheyes.com/dataset/ace-irk-sag",
    entries: [
      { label: "Single Cell", datasetType: "single_cell", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/other_merfisheyes" },
      { label: "Single Molecule", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/other_smeyes" },
    ],
  },
  {
    title: "VA00265 300-gene MERFISH results of human visual cortex (GW24 & GW34)",
    description:
      "This dataset is part of a set of MERFISH experiments to identify the cell types and their locations in the developing human visual cortex. Includes samples from both GW24 and GW34 gestational stages.",
    species: "Human",
    tissue: "Brain",
    platform: "MERFISH",
    externalLink: "https://bil.merfisheyes.com/dataset/ace-dry",
    entries: [
      { label: "GW24 SC", datasetType: "single_cell", datasetId: "ds_hXFp5BVbKe" },
      { label: "GW24 SM", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/ace-dip-dog-sm" },
      { label: "GW34 SC", datasetType: "single_cell", datasetId: "ds_0uyEF2TwZ2" },
      { label: "GW34 SM", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/ace-dry-dry-sm" },
    ],
  },
  {
    title: "MERFISH Imaging of mesocorticolimbic regions in the mouse brain",
    description:
      "This dataset is generated by MERFISH experiments with a 300 gene panel on whole brain coronal sections covering major mesocorticolimbic areas (e.g. PFC, ventral STR, TH). After classifying transcriptomic cell types from the MERFISH expression of cell type marker genes, cell types showing differential expressions of cocaine-related genes in brain regions of interest were identified.",
    species: "Mouse",
    tissue: "Brain",
    platform: "MERFISH",
    externalLink: "https://bil.merfisheyes.com/dataset/ace-dip-tub",
    entries: [
      { label: "1231122257 SC", datasetType: "single_cell", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/hongkui_ace-dip-tub-processed/ace-dip-tub/1231122257/merfish_output/cellpose_cyto2_nuclei_single_cell" },
      { label: "1231122257 SM", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/hongkui_ace-dip-tub-processed/ace-dip-tub/1231122257/merfish_output/cellpose_cyto2_nuclei_single_molecule" },
      { label: "1231122260 SC", datasetType: "single_cell", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/hongkui_ace-dip-tub-processed/ace-dip-tub/1231122260/merfish_output/cellpose_cyto2_nuclei_single_cell" },
      { label: "1231122263 SC", datasetType: "single_cell", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/hongkui_ace-dip-tub-processed/ace-dip-tub/1231122263/merfish_output/cellpose_cyto2_nuclei_single_cell" },
      { label: "1231122263 SM", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/hongkui_ace-dip-tub-processed/ace-dip-tub/1231122263/merfish_output/cellpose_cyto2_nuclei_single_molecule" },
      { label: "1231122266 SC", datasetType: "single_cell", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/hongkui_ace-dip-tub-processed/ace-dip-tub/1231122266/merfish_output/cellpose_cyto2_nuclei_single_cell" },
      { label: "1231122266 SM", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/hongkui_ace-dip-tub-processed/ace-dip-tub/1231122266/merfish_output/cellpose_cyto2_nuclei_single_molecule" },
    ],
  },
  {
    title: "Subcortical Spatial Transcriptomic Atlas of the Adult Rhesus Macaque",
    description:
      "MERSCOPE imaging of 10\u00b5m sections of core basal ganglia regions (caudate nucleus, putamen, substantia nigra, and nucleus accumbens) of fresh-frozen macaca mulatta brain slabs.",
    species: "Rhesus Macaque",
    tissue: "Brain",
    platform: "MERFISH",
    externalLink: "https://bil.merfisheyes.com/dataset/ace-dud-was",
    entries: [
      { label: "ace-dud-wag SC", datasetType: "single_cell", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/monkeys/ace-dud-wag/ace-dud-wag-sc" },
      { label: "ace-dud-wag SM", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/monkeys/ace-dud-wag/ace-dud-wag-sm" },
      { label: "ace-dud-vex SC", datasetType: "single_cell", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/monkeys/processed-ace-dud-vex" },
      { label: "ace-dud-vex SM", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/monkeys/ace-dud-vex-sm" },
      { label: "ace-dud-was SC", datasetType: "single_cell", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/monkeys/ace-dud-was/ace-dud-was-processed" },
      { label: "ace-dud-was SM", datasetType: "single_molecule", s3BaseUrl: "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/monkeys/ace-dud-was/ace-dud-was-sm" },
    ],
  },
];

async function main() {
  console.log("Seeding BIL datasets...");

  for (const ds of bilDatasets) {
    // Check for existing by title (idempotent)
    const existing = await prisma.catalogDataset.findFirst({
      where: { title: ds.title, isBil: true },
    });

    if (existing) {
      console.log(`  Skipping (already exists): ${ds.title.slice(0, 60)}...`);
      continue;
    }

    const created = await prisma.catalogDataset.create({
      data: {
        title: ds.title,
        description: ds.description,
        species: ds.species,
        tissue: ds.tissue,
        platform: ds.platform,
        externalLink: ds.externalLink,
        isBil: true,
        isPublished: true,
        entries: {
          create: ds.entries.map((e, i) => ({
            label: e.label,
            datasetType: e.datasetType,
            s3BaseUrl: e.s3BaseUrl ?? null,
            datasetId: e.datasetId ?? null,
            sortOrder: i,
          })),
        },
      },
      include: { entries: true },
    });

    console.log(`  Created: ${created.title.slice(0, 60)}... (${created.entries.length} entries)`);
  }

  const total = await prisma.catalogDataset.count({ where: { isBil: true } });
  console.log(`\nDone. Total BIL datasets: ${total}`);
}

main()
  .catch((e) => {
    console.error(e);
    process.exit(1);
  })
  .finally(() => prisma.$disconnect());
