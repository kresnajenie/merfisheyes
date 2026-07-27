"use client";

import {
  Table,
  TableHeader,
  TableColumn,
  TableBody,
  TableRow,
  TableCell,
  Button,
  Input,
  Chip,
  Select,
  SelectItem,
  Spinner,
} from "@heroui/react";
import NextLink from "next/link";
import { useCallback, useEffect, useState } from "react";
import { toast } from "react-toastify";

interface Row {
  id: string;
  title: string | null;
  status: string;
  datasetType: string | null;
  numCells: number;
  numGenes: number;
  viewCount: number;
  ingestSource: string | null;
  createdAt: string;
  viewerUrl: string | null;
}

type SortKey = "date" | "molecules" | "name";

const SORT_OPTIONS: { key: SortKey; label: string }[] = [
  { key: "date", label: "Date" },
  { key: "molecules", label: "Molecules" },
  { key: "name", label: "Name" },
];

const statusColor: Record<string, "success" | "warning" | "danger" | "default"> = {
  COMPLETE: "success",
  PROCESSING: "warning",
  QUEUED: "warning",
  FAILED: "danger",
};

function fmtDate(iso: string): string {
  try {
    return new Date(iso).toLocaleDateString(undefined, {
      year: "numeric",
      month: "short",
      day: "numeric",
    });
  } catch {
    return "—";
  }
}

export function AccountDatasets() {
  const [rows, setRows] = useState<Row[]>([]);
  const [loading, setLoading] = useState(true);
  const [sort, setSort] = useState<SortKey>("date");
  const [dir, setDir] = useState<"asc" | "desc">("desc");
  const [editingId, setEditingId] = useState<string | null>(null);
  const [draftTitle, setDraftTitle] = useState("");

  const load = useCallback(async () => {
    setLoading(true);
    try {
      const res = await fetch(`/api/ingest/mine?sort=${sort}&dir=${dir}`, {
        cache: "no-store",
      });

      if (res.ok) {
        const j = await res.json();

        setRows(j.datasets ?? []);
      }
    } finally {
      setLoading(false);
    }
  }, [sort, dir]);

  useEffect(() => {
    load();
  }, [load]);

  const saveTitle = async (id: string) => {
    const next = draftTitle.trim();

    setEditingId(null);
    const current = rows.find((r) => r.id === id);

    if (!current || next === (current.title ?? "")) return;

    // Optimistic update.
    setRows((prev) => prev.map((r) => (r.id === id ? { ...r, title: next } : r)));
    const res = await fetch(`/api/ingest/mine/${id}`, {
      method: "PATCH",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ title: next }),
    });

    if (!res.ok) {
      toast.error("Couldn't rename dataset.");
      load();
    }
  };

  return (
    <div className="flex flex-col gap-4">
      <div className="flex items-center gap-2">
        <Select
          aria-label="Sort by"
          className="w-40"
          selectedKeys={[sort]}
          size="sm"
          label="Sort by"
          onChange={(e) => setSort(e.target.value as SortKey)}
        >
          {SORT_OPTIONS.map((o) => (
            <SelectItem key={o.key}>{o.label}</SelectItem>
          ))}
        </Select>
        <Button
          size="sm"
          variant="flat"
          onPress={() => setDir((d) => (d === "asc" ? "desc" : "asc"))}
        >
          {dir === "asc" ? "Ascending ↑" : "Descending ↓"}
        </Button>
      </div>

      <Table aria-label="Your datasets">
        <TableHeader>
          <TableColumn>NAME</TableColumn>
          <TableColumn>TYPE</TableColumn>
          <TableColumn>MOLECULES</TableColumn>
          <TableColumn>GENES</TableColumn>
          <TableColumn>VIEWS</TableColumn>
          <TableColumn>CREATED</TableColumn>
          <TableColumn>STATUS</TableColumn>
          <TableColumn> </TableColumn>
        </TableHeader>
        <TableBody
          emptyContent={
            loading ? " " : "You don't own any datasets yet."
          }
          isLoading={loading}
          items={rows}
          loadingContent={<Spinner label="Loading…" />}
        >
          {(item) => (
            <TableRow key={item.id}>
              <TableCell>
                {editingId === item.id ? (
                  <Input
                    autoFocus
                    className="max-w-xs"
                    size="sm"
                    value={draftTitle}
                    onBlur={() => saveTitle(item.id)}
                    onKeyDown={(e) => {
                      if (e.key === "Enter") saveTitle(item.id);
                      if (e.key === "Escape") setEditingId(null);
                    }}
                    onValueChange={setDraftTitle}
                  />
                ) : (
                  <button
                    className="text-left hover:underline"
                    type="button"
                    onClick={() => {
                      setEditingId(item.id);
                      setDraftTitle(item.title ?? "");
                    }}
                  >
                    {item.title || "Untitled dataset"}
                  </button>
                )}
              </TableCell>
              <TableCell className="text-default-500">
                {item.datasetType === "single_molecule"
                  ? "Molecule"
                  : "Cell"}
              </TableCell>
              <TableCell>{item.numCells.toLocaleString()}</TableCell>
              <TableCell>{item.numGenes.toLocaleString()}</TableCell>
              <TableCell>{item.viewCount.toLocaleString()}</TableCell>
              <TableCell className="text-default-500">
                {fmtDate(item.createdAt)}
              </TableCell>
              <TableCell>
                <Chip
                  color={statusColor[item.status] ?? "default"}
                  size="sm"
                  variant="flat"
                >
                  {item.status.toLowerCase()}
                </Chip>
              </TableCell>
              <TableCell>
                {item.viewerUrl ? (
                  <Button
                    as={NextLink}
                    color="primary"
                    href={item.viewerUrl}
                    size="sm"
                    variant="flat"
                  >
                    Open
                  </Button>
                ) : (
                  <span className="text-default-400 text-xs">—</span>
                )}
              </TableCell>
            </TableRow>
          )}
        </TableBody>
      </Table>
    </div>
  );
}
