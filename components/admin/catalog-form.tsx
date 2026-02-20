"use client";

import { useState } from "react";
import { Input, Textarea } from "@heroui/input";
import { Button } from "@heroui/button";
import { Select, SelectItem } from "@heroui/select";
import { Switch } from "@heroui/switch";
import { useRouter } from "next/navigation";

export interface CatalogEntryFormData {
  label: string;
  datasetType: string;
  s3BaseUrl: string;
  datasetId: string;
  sortOrder: number;
}

export interface CatalogFormData {
  title: string;
  description: string;
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
  isBil: boolean;
  sortOrder: string;
  entries: CatalogEntryFormData[];
}

const emptyEntry: CatalogEntryFormData = {
  label: "",
  datasetType: "single_cell",
  s3BaseUrl: "",
  datasetId: "",
  sortOrder: 0,
};

const emptyForm: CatalogFormData = {
  title: "",
  description: "",
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
  isBil: false,
  sortOrder: "0",
  entries: [],
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

  const set = (key: keyof Omit<CatalogFormData, "entries">) => (val: string | boolean) => {
    setForm((prev) => ({ ...prev, [key]: val }));
  };

  const addEntry = () => {
    setForm((prev) => ({
      ...prev,
      entries: [...prev.entries, { ...emptyEntry, sortOrder: prev.entries.length }],
    }));
  };

  const removeEntry = (index: number) => {
    setForm((prev) => ({
      ...prev,
      entries: prev.entries.filter((_, i) => i !== index),
    }));
  };

  const updateEntry = (index: number, key: keyof CatalogEntryFormData, value: string | number) => {
    setForm((prev) => ({
      ...prev,
      entries: prev.entries.map((e, i) => (i === index ? { ...e, [key]: value } : e)),
    }));
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
      isBil: form.isBil,
      sortOrder: Number(form.sortOrder) || 0,
      entries: form.entries.map((entry, i) => ({
        label: entry.label.trim() || (entry.datasetType === "single_cell" ? "Single Cell" : "Single Molecule"),
        datasetType: entry.datasetType,
        s3BaseUrl: entry.s3BaseUrl.trim() || null,
        datasetId: entry.datasetId.trim() || null,
        sortOrder: entry.sortOrder ?? i,
      })),
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

      <Input
        description="Original source page (BIL, publication, data repository)"
        label="Source Link"
        value={form.externalLink}
        onValueChange={set("externalLink")}
      />

      {/* Data Entries */}
      <div className="border border-default-200 rounded-xl p-4 flex flex-col gap-4">
        <div className="flex items-center justify-between">
          <p className="text-sm font-semibold text-default-500">
            Data Entries ({form.entries.length})
          </p>
          <Button size="sm" variant="flat" onPress={addEntry}>
            Add Entry
          </Button>
        </div>

        {form.entries.length === 0 && (
          <p className="text-xs text-default-400 text-center py-2">
            No entries yet. Add entries to link this catalog item to viewer datasets.
          </p>
        )}

        {form.entries.map((entry, i) => (
          <div
            key={i}
            className="border border-default-100 rounded-lg p-3 flex flex-col gap-3 bg-default-50/50"
          >
            <div className="flex items-center justify-between">
              <span className="text-xs font-medium text-default-500">Entry {i + 1}</span>
              <Button
                color="danger"
                size="sm"
                variant="light"
                onPress={() => removeEntry(i)}
              >
                Remove
              </Button>
            </div>
            <div className="grid grid-cols-2 gap-3">
              <Input
                label="Label"
                placeholder="e.g. Brain Section 1"
                size="sm"
                value={entry.label}
                onValueChange={(v) => updateEntry(i, "label", v)}
              />
              <Select
                label="Type"
                selectedKeys={[entry.datasetType]}
                size="sm"
                onSelectionChange={(keys) => {
                  const val = Array.from(keys)[0] as string;
                  if (val) updateEntry(i, "datasetType", val);
                }}
              >
                <SelectItem key="single_cell">Single Cell</SelectItem>
                <SelectItem key="single_molecule">Single Molecule</SelectItem>
              </Select>
            </div>
            <Input
              description="Public S3 URL for from-s3 loading"
              label="S3 Base URL"
              size="sm"
              value={entry.s3BaseUrl}
              onValueChange={(v) => updateEntry(i, "s3BaseUrl", v)}
            />
            <Input
              description="ID of dataset uploaded through the app"
              label="Dataset ID"
              size="sm"
              value={entry.datasetId}
              onValueChange={(v) => updateEntry(i, "datasetId", v)}
            />
          </div>
        ))}
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
        <Switch isSelected={form.isBil} onValueChange={set("isBil")}>
          BIL
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
