import type { Dataset } from "@/components/dataset-card";

export const exampleDatasets: Dataset[] = [
  {
    id: "mouse-brain-merscope",
    title: "Mouse Brain",
    description:
      "High-resolution spatial transcriptomics of mouse brain with 2000+ genes",
    image:
      "https://images.unsplash.com/photo-1559757175-0eb30cd8c063?w=400&h=300&fit=crop",
    type: "MERFISHPLUS",
    link: "https://roymaimonel.merfisheyes.com",
  },
  {
    id: "human-heart-xenium",
    title: "Fetal Human Heart",
    description: "Spatial profiling of fetal human heart with 2000+ genes",
    image:
      "https://images.unsplash.com/photo-1530026405186-ed1f139313f8?w=400&h=300&fit=crop",
    type: "MERFISHPLUS",
    link: "https://zhu.merfisheyes.com",
  },
  {
    id: "mouse-liver-organoid",
    title: "Single Cell Gastrulation to Organogenesis Zebrafish",
    description:
      "Single Cell Spatial transcriptomics profiling of zebrafish (50pe, 75pe, 6-somite)",
    image:
      "https://images.unsplash.com/photo-1576086213369-97a306d36557?w=400&h=300&fit=crop",
    type: "weMERFISH & scMultiome",
    link: "https://schier.merfisheyes.com",
  },
  {
    id: "mouse-liver-organodid",
    title: "Single Molecule Zebrafish",
    description:
      "Single Molecule Spatial transcriptomics profiling of zebrafish (50pe, 75pe, 6-somite)",
    image:
      "https://images.unsplash.com/photo-1576086213369-97a306d36557?w=400&h=300&fit=crop",
    type: "weMERFISH",
    link: "https://sm-schier.merfisheyes.com",
  },
];
