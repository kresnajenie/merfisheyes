"use client";

import type { SortDescriptor } from "@heroui/table";

import { useCallback, useEffect, useState } from "react";
import { Chip } from "@heroui/chip";
import { Spinner } from "@heroui/spinner";
import { Pagination } from "@heroui/pagination";
import {
  Table,
  TableHeader,
  TableColumn,
  TableBody,
  TableRow,
  TableCell,
} from "@heroui/table";
import { toast } from "react-toastify";

type Fault = "USER" | "PLATFORM" | "UNKNOWN";

interface Failure {
  id: string;
  title: string | null;
  datasetType: string | null;
  errorMessage: string | null;
  faultCategory: Fault | null;
  autoLabel: string;
  createdAt: string;
  completedAt: string | null;
  owner: { name: string | null; email: string | null } | null;
}

const PAGE_LIMIT = 20;

const FAULTS: {
  key: Fault;
  label: string;
  color: "warning" | "danger" | "default";
}[] = [
  { key: "USER", label: "User", color: "warning" },
  { key: "PLATFORM", label: "Us", color: "danger" },
  { key: "UNKNOWN", label: "?", color: "default" },
];

function fmtDate(iso: string | null): string {
  if (!iso) return "—";

  return new Date(iso).toLocaleString(undefined, {
    month: "short",
    day: "numeric",
    hour: "2-digit",
    minute: "2-digit",
  });
}

export default function AdminFailuresPage() {
  const [failures, setFailures] = useState<Failure[]>([]);
  const [total, setTotal] = useState(0);
  const [page, setPage] = useState(1);
  const [loading, setLoading] = useState(true);
  // Default: most recent uploads first.
  const [sort, setSort] = useState<SortDescriptor>({
    column: "created",
    direction: "descending",
  });
  // Error cells render clamped; clicking one expands the full message.
  const [expanded, setExpanded] = useState<Set<string>>(new Set());

  const load = useCallback(async () => {
    setLoading(true);
    const params = new URLSearchParams({
      page: String(page),
      limit: String(PAGE_LIMIT),
      sort: String(sort.column),
      dir: sort.direction === "ascending" ? "asc" : "desc",
    });
    const res = await fetch(`/api/admin/failures?${params}`, {
      cache: "no-store",
    });

    if (res.ok) {
      const j = await res.json();

      setFailures(j.failures ?? []);
      setTotal(j.total ?? 0);
    } else {
      toast.error("Couldn't load failures.");
    }
    setLoading(false);
  }, [page, sort]);

  useEffect(() => {
    load();
  }, [load]);

  const onSortChange = (next: SortDescriptor) => {
    setSort(next);
    setPage(1);
  };

  const toggleExpanded = (id: string) => {
    setExpanded((prev) => {
      const next = new Set(prev);

      if (next.has(id)) next.delete(id);
      else next.add(id);

      return next;
    });
  };

  const setFault = async (id: string, next: Fault | null) => {
    // Optimistic — the tag flips instantly.
    setFailures((rows) =>
      rows.map((f) => (f.id === id ? { ...f, faultCategory: next } : f)),
    );
    const res = await fetch(`/api/admin/processing/${id}`, {
      method: "PATCH",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ faultCategory: next }),
    });

    if (!res.ok) {
      toast.error("Couldn't update the fault tag.");
      load();
    }
  };

  const pages = Math.max(1, Math.ceil(total / PAGE_LIMIT));

  return (
    <div className="space-y-6">
      <div>
        <h1 className="text-2xl font-bold">Failures ({total})</h1>
        <p className="mt-1 text-default-500">
          Every failed upload, most recent first. Click a column header to sort;
          click an error to expand the full message.
        </p>
      </div>

      {loading ? (
        <div className="flex justify-center py-16">
          <Spinner label="Loading…" size="lg" />
        </div>
      ) : failures.length === 0 ? (
        <p className="py-12 text-center text-default-500">No failures. 🎉</p>
      ) : (
        <>
          <Table
            removeWrapper
            aria-label="Failed uploads"
            sortDescriptor={sort}
            onSortChange={onSortChange}
          >
            <TableHeader>
              <TableColumn key="title" allowsSorting>
                DATASET
              </TableColumn>
              <TableColumn key="owner" allowsSorting>
                OWNER
              </TableColumn>
              <TableColumn key="why">WHY</TableColumn>
              <TableColumn key="error">ERROR</TableColumn>
              <TableColumn key="fault" allowsSorting>
                FAULT
              </TableColumn>
              <TableColumn key="created" allowsSorting>
                UPLOADED
              </TableColumn>
              <TableColumn key="failedAt" allowsSorting>
                FAILED
              </TableColumn>
            </TableHeader>
            <TableBody>
              {failures.map((f) => (
                <TableRow key={f.id}>
                  <TableCell>
                    <div className="flex max-w-[220px] flex-col">
                      <span
                        className="truncate text-sm font-medium"
                        title={f.title ?? undefined}
                      >
                        {f.title || "Untitled"}
                      </span>
                      <span className="text-[11px] text-default-400">
                        {f.datasetType ?? "—"} · {f.id}
                      </span>
                    </div>
                  </TableCell>
                  <TableCell>
                    {f.owner?.email ? (
                      <a
                        className="text-sm text-primary hover:underline"
                        href={`mailto:${f.owner.email}?subject=${encodeURIComponent(
                          `Your MERFISHEYES upload "${f.title || f.id}"`,
                        )}`}
                      >
                        {f.owner.name || f.owner.email}
                      </a>
                    ) : (
                      <span className="text-sm text-default-400">—</span>
                    )}
                  </TableCell>
                  <TableCell>
                    <Chip size="sm" variant="flat">
                      {f.autoLabel}
                    </Chip>
                  </TableCell>
                  <TableCell>
                    {f.errorMessage ? (
                      <button
                        className="block max-w-[420px] cursor-pointer text-left"
                        title={
                          expanded.has(f.id)
                            ? "Click to collapse"
                            : "Click to expand"
                        }
                        type="button"
                        onClick={() => toggleExpanded(f.id)}
                      >
                        <span
                          className={`block whitespace-pre-wrap [overflow-wrap:anywhere] text-xs text-default-500 ${
                            expanded.has(f.id) ? "" : "line-clamp-2"
                          }`}
                        >
                          {f.errorMessage}
                        </span>
                      </button>
                    ) : (
                      <span className="text-xs text-default-500">—</span>
                    )}
                  </TableCell>
                  <TableCell>
                    <div className="flex gap-1">
                      {FAULTS.map((opt) => {
                        const active = f.faultCategory === opt.key;

                        return (
                          <button
                            key={opt.key}
                            type="button"
                            onClick={() =>
                              setFault(f.id, active ? null : opt.key)
                            }
                          >
                            <Chip
                              className="cursor-pointer"
                              color={active ? opt.color : "default"}
                              size="sm"
                              variant={active ? "solid" : "bordered"}
                            >
                              {opt.label}
                            </Chip>
                          </button>
                        );
                      })}
                    </div>
                  </TableCell>
                  <TableCell>
                    <span className="whitespace-nowrap text-xs text-default-400">
                      {fmtDate(f.createdAt)}
                    </span>
                  </TableCell>
                  <TableCell>
                    <span className="whitespace-nowrap text-xs text-default-400">
                      {fmtDate(f.completedAt)}
                    </span>
                  </TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>

          {pages > 1 && (
            <div className="flex justify-center">
              <Pagination
                showControls
                page={page}
                size="sm"
                total={pages}
                onChange={setPage}
              />
            </div>
          )}
        </>
      )}
    </div>
  );
}
