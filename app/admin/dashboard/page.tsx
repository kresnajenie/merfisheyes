"use client";

import { useCallback, useEffect, useState } from "react";
import { Card, CardBody } from "@heroui/card";
import { Chip } from "@heroui/chip";
import { Spinner } from "@heroui/spinner";
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

interface PerUser {
  userId: string;
  name: string | null;
  email: string | null;
  total: number;
  complete: number;
  failed: number;
  inProgress: number;
}

interface Stats {
  statusCounts: Record<string, number>;
  total: number;
  successRate: number | null;
  failedCount: number;
  uploaders: number;
  recent: {
    created7: number;
    created30: number;
    active7: number;
    active30: number;
  };
  failures: Failure[];
  perUser: PerUser[];
}

const STATUS_ORDER = [
  "UPLOADING",
  "QUEUED",
  "PROCESSING",
  "COMPLETE",
  "FAILED",
] as const;

const STATUS_COLOR: Record<
  string,
  "default" | "primary" | "secondary" | "success" | "warning" | "danger"
> = {
  UPLOADING: "default",
  QUEUED: "secondary",
  PROCESSING: "primary",
  COMPLETE: "success",
  FAILED: "danger",
};

const FAULTS: {
  key: Fault;
  label: string;
  color: "warning" | "danger" | "default";
}[] = [
  { key: "USER", label: "User", color: "warning" },
  { key: "PLATFORM", label: "Us", color: "danger" },
  { key: "UNKNOWN", label: "?", color: "default" },
];

function StatTile({
  label,
  value,
  hint,
}: {
  label: string;
  value: string | number;
  hint?: string;
}) {
  return (
    <Card shadow="sm">
      <CardBody className="gap-1">
        <span className="text-xs uppercase tracking-wider text-default-400">
          {label}
        </span>
        <span className="text-2xl font-semibold">{value}</span>
        {hint && <span className="text-xs text-default-400">{hint}</span>}
      </CardBody>
    </Card>
  );
}

function fmtDate(iso: string | null): string {
  if (!iso) return "—";

  return new Date(iso).toLocaleString(undefined, {
    month: "short",
    day: "numeric",
    hour: "2-digit",
    minute: "2-digit",
  });
}

export default function AdminDashboardPage() {
  const [data, setData] = useState<Stats | null>(null);
  const [loading, setLoading] = useState(true);

  const load = useCallback(async () => {
    const res = await fetch("/api/admin/processing", { cache: "no-store" });

    if (res.ok) setData(await res.json());
    setLoading(false);
  }, []);

  useEffect(() => {
    load();
  }, [load]);

  const setFault = async (id: string, next: Fault | null) => {
    // Optimistic — the tag flips instantly.
    setData((d) =>
      d
        ? {
            ...d,
            failures: d.failures.map((f) =>
              f.id === id ? { ...f, faultCategory: next } : f,
            ),
          }
        : d,
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

  if (loading) {
    return (
      <div className="flex justify-center py-16">
        <Spinner label="Loading…" size="lg" />
      </div>
    );
  }

  if (!data) {
    return (
      <p className="py-12 text-center text-default-500">
        Couldn&apos;t load stats.
      </p>
    );
  }

  const successPct =
    data.successRate == null ? "—" : `${Math.round(data.successRate * 100)}%`;

  return (
    <div className="space-y-8">
      <div>
        <h1 className="text-2xl font-bold">Processing dashboard</h1>
        <p className="mt-1 text-default-500">
          Upload &amp; processing health across all datasets.
        </p>
      </div>

      {/* Stat tiles */}
      <div className="grid grid-cols-2 gap-4 sm:grid-cols-3 lg:grid-cols-5">
        <StatTile label="Total datasets" value={data.total} />
        <StatTile
          hint={`${data.statusCounts.COMPLETE} ok · ${data.failedCount} failed`}
          label="Success rate"
          value={successPct}
        />
        <StatTile label="Failed" value={data.failedCount} />
        <StatTile label="Uploaders" value={data.uploaders} />
        <StatTile
          hint={`${data.recent.created30} uploads`}
          label="Active (30d)"
          value={data.recent.active30}
        />
      </div>

      {/* Status breakdown */}
      <div className="flex flex-wrap gap-2">
        {STATUS_ORDER.map((s) => (
          <Chip key={s} color={STATUS_COLOR[s]} variant="flat">
            {s.toLowerCase()}: {data.statusCounts[s] ?? 0}
          </Chip>
        ))}
      </div>

      {/* Recent activity */}
      <div className="grid grid-cols-2 gap-4 sm:grid-cols-4">
        <StatTile label="Uploads (7d)" value={data.recent.created7} />
        <StatTile label="Uploads (30d)" value={data.recent.created30} />
        <StatTile label="Active users (7d)" value={data.recent.active7} />
        <StatTile label="Active users (30d)" value={data.recent.active30} />
      </div>

      {/* Failures */}
      <section className="space-y-3">
        <h2 className="text-lg font-semibold">
          Recent failures ({data.failures.length})
        </h2>
        {data.failures.length === 0 ? (
          <p className="py-6 text-center text-default-500">No failures. 🎉</p>
        ) : (
          <Table removeWrapper aria-label="Recent failures">
            <TableHeader>
              <TableColumn>DATASET</TableColumn>
              <TableColumn>OWNER</TableColumn>
              <TableColumn>WHY</TableColumn>
              <TableColumn>ERROR</TableColumn>
              <TableColumn>FAULT</TableColumn>
              <TableColumn>WHEN</TableColumn>
            </TableHeader>
            <TableBody>
              {data.failures.map((f) => (
                <TableRow key={f.id}>
                  <TableCell>
                    <div className="flex flex-col">
                      <span className="max-w-[220px] truncate text-sm font-medium">
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
                    <span
                      className="block max-w-[280px] truncate text-xs text-default-500"
                      title={f.errorMessage ?? ""}
                    >
                      {f.errorMessage || "—"}
                    </span>
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
                      {fmtDate(f.completedAt ?? f.createdAt)}
                    </span>
                  </TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        )}
      </section>

      {/* Per-user activity */}
      <section className="space-y-3">
        <h2 className="text-lg font-semibold">
          Uploaders ({data.perUser.length})
        </h2>
        {data.perUser.length === 0 ? (
          <p className="py-6 text-center text-default-500">No uploaders yet.</p>
        ) : (
          <Table removeWrapper aria-label="Per-user activity">
            <TableHeader>
              <TableColumn>USER</TableColumn>
              <TableColumn>TOTAL</TableColumn>
              <TableColumn>COMPLETE</TableColumn>
              <TableColumn>FAILED</TableColumn>
              <TableColumn>IN PROGRESS</TableColumn>
            </TableHeader>
            <TableBody>
              {data.perUser.map((u) => (
                <TableRow key={u.userId}>
                  <TableCell>
                    <div className="flex flex-col">
                      <span className="text-sm font-medium">
                        {u.name || "—"}
                      </span>
                      {u.email && (
                        <a
                          className="text-[11px] text-default-400 hover:text-primary"
                          href={`mailto:${u.email}`}
                        >
                          {u.email}
                        </a>
                      )}
                    </div>
                  </TableCell>
                  <TableCell>{u.total}</TableCell>
                  <TableCell>
                    <span className="text-success">{u.complete}</span>
                  </TableCell>
                  <TableCell>
                    <span className={u.failed ? "text-danger" : ""}>
                      {u.failed}
                    </span>
                  </TableCell>
                  <TableCell>{u.inProgress}</TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        )}
      </section>
    </div>
  );
}
