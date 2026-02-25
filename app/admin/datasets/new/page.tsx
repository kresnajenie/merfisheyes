"use client";

import { CatalogForm } from "@/components/admin/catalog-form";

export default function NewCatalogDatasetPage() {
  return (
    <div>
      <h1 className="text-2xl font-bold mb-6">Add Catalog Dataset</h1>
      <CatalogForm method="POST" submitUrl="/api/admin/catalog" />
    </div>
  );
}
