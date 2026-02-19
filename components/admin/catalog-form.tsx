"use client";

import { useState } from "react";
import { Input, Textarea } from "@heroui/input";
import { Button } from "@heroui/button";
import { Select, SelectItem } from "@heroui/select";
import { Switch } from "@heroui/switch";
import { useRouter } from "next/navigation";

export interface CatalogFormData {
  title: string;
  description: string;
  datasetType: string;
  s3BaseUrl: string;
  datasetId: string;
  externalLink: string;
  species: string;
  disease: string;
  institute: string;
  tissue: string;
  platform: string;
  tags: string;
  thumbnailUrl: string;
  numCells: string;
  numGenes: string;
  isPublished: boolean;
  isFeatured: boolean;
  sortOrder: string;
}

const emptyForm: CatalogFormData = {
  title: "",
  description: "",
  datasetType: "single_cell",
  s3BaseUrl: "",
  datasetId: "",
  externalLink: "",
  species: "",
  disease: "",
  institute: "",
  tissue: "",
  platform: "",
  tags: "",
  thumbnailUrl: "",
  numCells: "",
  numGenes: "",
  isPublished: false,
  isFeatured: false,
  sortOrder: "0",
};

interface CatalogFormProps {
  initialData?: Partial<CatalogFormData>;
  /** PATCH for edit, POST for create */
  submitUrl: string;
  method: "POST" | "PATCH";
}

export function CatalogForm({ initialData, submitUrl, method }: CatalogFormProps) {
  const router = useRouter();
  const [form, setForm] = useState<CatalogFormData>({ ...emptyForm, ...initialData });
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const set = (key: keyof CatalogFormData) => (val: string | boolean) => {
    setForm((prev) => ({ ...prev, [key]: val }));
  };

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!form.title.trim()) {
      setError("Title is required");
      return;
    }
    setSaving(true);
    setError(null);

    const payload = {
      title: form.title.trim(),
      description: form.description.trim() || null,
      datasetType: form.datasetType,
      s3BaseUrl: form.s3BaseUrl.trim() || null,
      datasetId: form.datasetId.trim() || null,
      externalLink: form.externalLink.trim() || null,
      species: form.species.trim() || null,
      disease: form.disease.trim() || null,
      institute: form.institute.trim() || null,
      tissue: form.tissue.trim() || null,
      platform: form.platform.trim() || null,
      tags: form.tags
        .split(",")
        .map((t) => t.trim())
        .filter(Boolean),
      thumbnailUrl: form.thumbnailUrl.trim() || null,
      numCells: form.numCells ? Number(form.numCells) : null,
      numGenes: form.numGenes ? Number(form.numGenes) : null,
      isPublished: form.isPublished,
      isFeatured: form.isFeatured,
      sortOrder: Number(form.sortOrder) || 0,
    };

    const res = await fetch(submitUrl, {
      method,
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });

    if (!res.ok) {
      setError("Failed to save");
      setSaving(false);
      return;
    }

    router.push("/admin/datasets");
    router.refresh();
  };

  return (
    <form className="flex flex-col gap-5 max-w-2xl" onSubmit={handleSubmit}>
      {error && (
        <div className="text-danger text-sm bg-danger-50 dark:bg-danger-100/10 px-4 py-2 rounded-lg">
          {error}
        </div>
      )}

      <Input
        isRequired
        label="Title"
        value={form.title}
        onValueChange={set("title")}
      />

      <Textarea
        label="Description"
        minRows={3}
        value={form.description}
        onValueChange={set("description")}
      />

      <Select
        label="Dataset Type"
        selectedKeys={[form.datasetType]}
        onSelectionChange={(keys) => {
          const val = Array.from(keys)[0] as string;
          if (val) set("datasetType")(val);
        }}
      >
        <SelectItem key="single_cell">Single Cell</SelectItem>
        <SelectItem key="single_molecule">Single Molecule</SelectItem>
      </Select>

      <div className="border border-default-200 rounded-xl p-4 flex flex-col gap-4">
        <p className="text-sm font-semibold text-default-500">Data Source (set one)</p>
        <Input
          description="Public S3 URL for from-s3 loading"
          label="S3 Base URL"
          value={form.s3BaseUrl}
          onValueChange={set("s3BaseUrl")}
        />
        <Input
          description="ID of dataset uploaded through the app"
          label="Dataset ID"
          value={form.datasetId}
          onValueChange={set("datasetId")}
        />
        <Input
          description="External link (opens in new tab)"
          label="External Link"
          value={form.externalLink}
          onValueChange={set("externalLink")}
        />
      </div>

      <div className="grid grid-cols-2 gap-4">
        <Input label="Species" value={form.species} onValueChange={set("species")} />
        <Input label="Disease" value={form.disease} onValueChange={set("disease")} />
        <Input label="Institute" value={form.institute} onValueChange={set("institute")} />
        <Input label="Tissue" value={form.tissue} onValueChange={set("tissue")} />
        <Input label="Platform" value={form.platform} onValueChange={set("platform")} />
        <Input
          description="Comma-separated"
          label="Tags"
          value={form.tags}
          onValueChange={set("tags")}
        />
      </div>

      <Input
        description="Direct URL to thumbnail image"
        label="Thumbnail URL"
        value={form.thumbnailUrl}
        onValueChange={set("thumbnailUrl")}
      />

      <div className="grid grid-cols-3 gap-4">
        <Input
          label="Num Cells / Molecules"
          type="number"
          value={form.numCells}
          onValueChange={set("numCells")}
        />
        <Input
          label="Num Genes"
          type="number"
          value={form.numGenes}
          onValueChange={set("numGenes")}
        />
        <Input
          label="Sort Order"
          type="number"
          value={form.sortOrder}
          onValueChange={set("sortOrder")}
        />
      </div>

      <div className="flex gap-8">
        <Switch isSelected={form.isPublished} onValueChange={set("isPublished")}>
          Published
        </Switch>
        <Switch isSelected={form.isFeatured} onValueChange={set("isFeatured")}>
          Featured
        </Switch>
      </div>

      <div className="flex gap-3">
        <Button color="primary" isLoading={saving} type="submit">
          {method === "POST" ? "Create" : "Save"}
        </Button>
        <Button variant="flat" onPress={() => router.back()}>
          Cancel
        </Button>
      </div>
    </form>
  );
}
