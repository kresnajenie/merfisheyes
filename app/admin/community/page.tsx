"use client";

import { useCallback, useEffect, useState } from "react";
import { Card, CardBody } from "@heroui/card";
import { Chip } from "@heroui/chip";
import { Button } from "@heroui/button";
import { Textarea } from "@heroui/input";
import { Tabs, Tab } from "@heroui/tabs";
import { Spinner } from "@heroui/spinner";
import NextLink from "next/link";
import { toast } from "react-toastify";

type ReviewStatus = "PENDING" | "APPROVED" | "REJECTED";

interface Entry {
  id: string;
  label: string;
  datasetType: string;
  s3BaseUrl: string | null;
  datasetId: string | null;
}

interface CommunityItem {
  id: string;
  title: string;
  description: string | null;
  species: string | null;
  tissue: string | null;
  platform: string | null;
  thumbnailUrl: string | null;
  reviewStatus: ReviewStatus;
  reviewNote: string | null;
  sourceDatasetId: string | null;
  sourceProjectId: string | null;
  updatedAt: string;
  entries: Entry[];
  submitter?: { name: string | null; email: string | null } | null;
}

export default function AdminCommunityPage() {
  const [status, setStatus] = useState<ReviewStatus>("PENDING");
  const [items, setItems] = useState<CommunityItem[]>([]);
  const [loading, setLoading] = useState(true);
  const [notes, setNotes] = useState<Record<string, string>>({});
  const [busyId, setBusyId] = useState<string | null>(null);

  const load = useCallback(async () => {
    setLoading(true);
    try {
      const res = await fetch(`/api/admin/community?status=${status}`, {
        cache: "no-store",
      });

      if (res.ok) {
        const j = await res.json();

        setItems(j.items ?? []);
      }
    } finally {
      setLoading(false);
    }
  }, [status]);

  useEffect(() => {
    load();
  }, [load]);

  const act = async (id: string, action: "approve" | "reject") => {
    setBusyId(id);
    try {
      const res = await fetch(`/api/admin/community/${id}`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ action, note: notes[id] ?? "" }),
      });

      if (res.ok) {
        toast.success(
          action === "approve" ? "Approved & published." : "Rejected.",
        );
        await load();
      } else {
        toast.error("Action failed.");
      }
    } finally {
      setBusyId(null);
    }
  };

  return (
    <div>
      <h1 className="text-2xl font-bold mb-6">Community review</h1>

      <Tabs
        aria-label="Review status"
        className="mb-4"
        selectedKey={status}
        variant="underlined"
        onSelectionChange={(k) => setStatus(k as ReviewStatus)}
      >
        <Tab key="PENDING" title="Pending" />
        <Tab key="APPROVED" title="Approved" />
        <Tab key="REJECTED" title="Rejected" />
      </Tabs>

      {loading ? (
        <div className="flex justify-center py-12">
          <Spinner size="lg" />
        </div>
      ) : items.length === 0 ? (
        <p className="text-default-500 text-center py-12">
          Nothing {status.toLowerCase()}.
        </p>
      ) : (
        <div className="flex flex-col gap-4">
          {items.map((item) => (
            <Card key={item.id} shadow="sm">
              <CardBody className="flex flex-col md:flex-row gap-4">
                {item.thumbnailUrl && (
                  // eslint-disable-next-line @next/next/no-img-element
                  <img
                    alt={item.title}
                    className="w-full md:w-40 h-28 object-cover rounded-lg"
                    src={item.thumbnailUrl}
                  />
                )}
                <div className="flex-1 flex flex-col gap-2 min-w-0">
                  <div className="flex items-center gap-2 flex-wrap">
                    <h3 className="font-semibold">{item.title}</h3>
                    <Chip size="sm" variant="flat">
                      {item.sourceProjectId ? "Project" : "Dataset"}
                    </Chip>
                    <Chip size="sm" variant="flat">
                      {item.entries.length} entr
                      {item.entries.length === 1 ? "y" : "ies"}
                    </Chip>
                  </div>

                  {item.submitter && (
                    <p className="text-xs text-default-400">
                      Submitted by{" "}
                      {item.submitter.name || item.submitter.email || "unknown"}
                    </p>
                  )}

                  {item.description && (
                    <p className="text-sm text-default-500 line-clamp-2">
                      {item.description}
                    </p>
                  )}

                  <div className="flex flex-wrap gap-1">
                    {item.species && (
                      <Chip size="sm" variant="flat">
                        {item.species}
                      </Chip>
                    )}
                    {item.tissue && (
                      <Chip size="sm" variant="flat">
                        {item.tissue}
                      </Chip>
                    )}
                    {item.platform && (
                      <Chip size="sm" variant="flat">
                        {item.platform}
                      </Chip>
                    )}
                  </div>

                  <div className="flex flex-wrap gap-1">
                    {item.entries.map((e) => {
                      const base =
                        e.datasetType === "single_molecule"
                          ? "/sm-viewer"
                          : "/viewer";
                      const href = e.s3BaseUrl
                        ? `${base}/from-s3?url=${encodeURIComponent(e.s3BaseUrl)}`
                        : e.datasetId
                          ? `${base}/${e.datasetId}`
                          : null;

                      return href ? (
                        <Button
                          key={e.id}
                          as={NextLink}
                          href={href}
                          size="sm"
                          target="_blank"
                          variant="flat"
                        >
                          {e.label}
                        </Button>
                      ) : (
                        <Chip key={e.id} size="sm" variant="flat">
                          {e.label}
                        </Chip>
                      );
                    })}
                  </div>

                  {status === "PENDING" && (
                    <Textarea
                      className="mt-1"
                      label="Rejection note (optional)"
                      minRows={1}
                      value={notes[item.id] ?? ""}
                      onValueChange={(v) =>
                        setNotes((n) => ({ ...n, [item.id]: v }))
                      }
                    />
                  )}

                  {item.reviewNote && status !== "PENDING" && (
                    <p className="text-xs text-danger">
                      Note: {item.reviewNote}
                    </p>
                  )}

                  <div className="flex gap-2 mt-1">
                    <Button
                      as={NextLink}
                      href={`/admin/datasets/${item.id}/edit`}
                      size="sm"
                      variant="flat"
                    >
                      Edit
                    </Button>
                    {status !== "APPROVED" && (
                      <Button
                        color="success"
                        isLoading={busyId === item.id}
                        size="sm"
                        onPress={() => act(item.id, "approve")}
                      >
                        Approve & publish
                      </Button>
                    )}
                    {status !== "REJECTED" && (
                      <Button
                        color="danger"
                        isLoading={busyId === item.id}
                        size="sm"
                        variant="flat"
                        onPress={() => act(item.id, "reject")}
                      >
                        Reject
                      </Button>
                    )}
                  </div>
                </div>
              </CardBody>
            </Card>
          ))}
        </div>
      )}
    </div>
  );
}
