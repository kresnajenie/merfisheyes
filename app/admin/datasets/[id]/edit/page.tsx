"use client";

import { useEffect, useState } from "react";
import { useParams } from "next/navigation";
import { Spinner } from "@heroui/spinner";

import { CatalogForm, CatalogFormData } from "@/components/admin/catalog-form";

export default function EditCatalogDatasetPage() {
  const { id } = useParams<{ id: string }>();
  const [initialData, setInitialData] = useState<Partial<CatalogFormData> | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    (async () => {
      const res = await fetch(`/api/admin/catalog/${id}`);
      if (!res.ok) {
        setLoading(false);
        return;
      }
      const item = await res.json();
      setInitialData({
        title: item.title ?? "",
        description: item.description ?? "",
        datasetType: item.datasetType ?? "single_cell",
        s3BaseUrl: item.s3BaseUrl ?? "",
        datasetId: item.datasetId ?? "",
        externalLink: item.externalLink ?? "",
        species: item.species ?? "",
        disease: item.disease ?? "",
        institute: item.institute ?? "",
        tissue: item.tissue ?? "",
        platform: item.platform ?? "",
        tags: (item.tags ?? []).join(", "),
        thumbnailUrl: item.thumbnailUrl ?? "",
        numCells: item.numCells != null ? String(item.numCells) : "",
        numGenes: item.numGenes != null ? String(item.numGenes) : "",
        isPublished: item.isPublished ?? false,
        isFeatured: item.isFeatured ?? false,
        sortOrder: String(item.sortOrder ?? 0),
      });
      setLoading(false);
    })();
  }, [id]);

  if (loading) {
    return (
      <div className="flex justify-center py-12">
        <Spinner size="lg" />
      </div>
    );
  }

  if (!initialData) {
    return <p className="text-danger">Dataset not found</p>;
  }

  return (
    <div>
      <h1 className="text-2xl font-bold mb-6">Edit Catalog Dataset</h1>
      <CatalogForm
        initialData={initialData}
        method="PATCH"
        submitUrl={`/api/admin/catalog/${id}`}
      />
    </div>
  );
}
