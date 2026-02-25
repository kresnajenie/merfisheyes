"use client";

import { useEffect, useState, useCallback } from "react";
import {
  Table,
  TableHeader,
  TableColumn,
  TableBody,
  TableRow,
  TableCell,
} from "@heroui/table";
import { Button } from "@heroui/button";
import { Chip } from "@heroui/chip";
import { Spinner } from "@heroui/spinner";
import NextLink from "next/link";

interface CatalogEntry {
  id: string;
  label: string;
  datasetType: string;
}

interface CatalogItem {
  id: string;
  title: string;
  entries: CatalogEntry[];
  species: string | null;
  tissue: string | null;
  platform: string | null;
  isPublished: boolean;
  isFeatured: boolean;
  isInternal: boolean;
  numCells: number | null;
  numGenes: number | null;
}

export default function AdminDatasetsPage() {
  const [items, setItems] = useState<CatalogItem[]>([]);
  const [loading, setLoading] = useState(true);

  const fetchItems = useCallback(async () => {
    setLoading(true);
    const res = await fetch("/api/admin/catalog");
    const data = await res.json();
    setItems(data.items ?? []);
    setLoading(false);
  }, []);

  useEffect(() => {
    fetchItems();
  }, [fetchItems]);

  const toggleField = async (id: string, field: "isPublished" | "isFeatured" | "isInternal", current: boolean) => {
    await fetch(`/api/admin/catalog/${id}`, {
      method: "PATCH",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ [field]: !current }),
    });
    fetchItems();
  };

  const deleteItem = async (id: string) => {
    if (!confirm("Delete this catalog entry?")) return;
    await fetch(`/api/admin/catalog/${id}`, { method: "DELETE" });
    fetchItems();
  };

  return (
    <div>
      <div className="flex items-center justify-between mb-6">
        <h1 className="text-2xl font-bold">Catalog Datasets</h1>
        <Button as={NextLink} color="primary" href="/admin/datasets/new">
          Add Dataset
        </Button>
      </div>

      {loading ? (
        <div className="flex justify-center py-12">
          <Spinner size="lg" />
        </div>
      ) : items.length === 0 ? (
        <p className="text-default-500 text-center py-12">
          No catalog datasets yet. Click &quot;Add Dataset&quot; to create one.
        </p>
      ) : (
        <Table aria-label="Catalog datasets">
          <TableHeader>
            <TableColumn>Title</TableColumn>
            <TableColumn>Entries</TableColumn>
            <TableColumn>Species</TableColumn>
            <TableColumn>Published</TableColumn>
            <TableColumn>Featured</TableColumn>
            <TableColumn>Internal</TableColumn>
            <TableColumn>Actions</TableColumn>
          </TableHeader>
          <TableBody>
            {items.map((item) => (
              <TableRow key={item.id}>
                <TableCell>
                  <span className="font-medium">{item.title}</span>
                </TableCell>
                <TableCell>
                  <div className="flex flex-wrap gap-1">
                    {item.entries.length === 0 ? (
                      <span className="text-default-400 text-xs">None</span>
                    ) : (
                      item.entries.map((entry) => (
                        <Chip
                          key={entry.id}
                          color={entry.datasetType === "single_cell" ? "primary" : "secondary"}
                          size="sm"
                          variant="flat"
                        >
                          {entry.label || (entry.datasetType === "single_cell" ? "SC" : "SM")}
                        </Chip>
                      ))
                    )}
                  </div>
                </TableCell>
                <TableCell>{item.species ?? "—"}</TableCell>
                <TableCell>
                  <Chip
                    className="cursor-pointer"
                    color={item.isPublished ? "success" : "default"}
                    size="sm"
                    variant="flat"
                    onClick={() => toggleField(item.id, "isPublished", item.isPublished)}
                  >
                    {item.isPublished ? "Published" : "Draft"}
                  </Chip>
                </TableCell>
                <TableCell>
                  <Chip
                    className="cursor-pointer"
                    color={item.isFeatured ? "warning" : "default"}
                    size="sm"
                    variant="flat"
                    onClick={() => toggleField(item.id, "isFeatured", item.isFeatured)}
                  >
                    {item.isFeatured ? "Featured" : "—"}
                  </Chip>
                </TableCell>
                <TableCell>
                  <Chip
                    className="cursor-pointer"
                    color={item.isInternal ? "secondary" : "default"}
                    size="sm"
                    variant="flat"
                    onClick={() => toggleField(item.id, "isInternal", item.isInternal)}
                  >
                    {item.isInternal ? "Internal" : "—"}
                  </Chip>
                </TableCell>
                <TableCell>
                  <div className="flex gap-2">
                    <Button
                      as={NextLink}
                      href={`/admin/datasets/${item.id}/edit`}
                      size="sm"
                      variant="flat"
                    >
                      Edit
                    </Button>
                    <Button
                      color="danger"
                      size="sm"
                      variant="flat"
                      onPress={() => deleteItem(item.id)}
                    >
                      Delete
                    </Button>
                  </div>
                </TableCell>
              </TableRow>
            ))}
          </TableBody>
        </Table>
      )}
    </div>
  );
}
