"use client";

import type { MetadataDraft } from "./metadata-edit-modal";
import type { DatasetRow, ProjectRow } from "./types";

import { useCallback, useEffect, useMemo, useState } from "react";
import { Tabs, Tab } from "@heroui/tabs";
import { Button } from "@heroui/button";
import { Spinner } from "@heroui/spinner";
import { toast } from "react-toastify";

import { AccountDatasetCard } from "./account-dataset-card";
import { AccountProjectCard } from "./account-project-card";
import { MetadataEditModal } from "./metadata-edit-modal";
import { AddToProjectModal } from "./add-to-project-modal";
import {
  AccountFilterBar,
  DEFAULT_FILTERS,
  applyDatasetFilters,
  type DatasetFilters,
} from "./account-filter-bar";

type EditTarget =
  | { kind: "dataset"; data: DatasetRow }
  | { kind: "project"; data: ProjectRow }
  | { kind: "new-project" }
  | null;

export function AccountWorkspace() {
  const [datasets, setDatasets] = useState<DatasetRow[]>([]);
  const [projects, setProjects] = useState<ProjectRow[]>([]);
  const [loading, setLoading] = useState(true);
  const [editTarget, setEditTarget] = useState<EditTarget>(null);
  const [addToProjectFor, setAddToProjectFor] = useState<DatasetRow | null>(
    null,
  );
  const [filters, setFilters] = useState<DatasetFilters>(DEFAULT_FILTERS);

  const loadDatasets = useCallback(async () => {
    const res = await fetch("/api/ingest/mine?sort=date&dir=desc", {
      cache: "no-store",
    });

    if (res.ok) {
      const j = await res.json();

      setDatasets(j.datasets ?? []);
    }
  }, []);

  const loadProjects = useCallback(async () => {
    const res = await fetch("/api/projects", { cache: "no-store" });

    if (res.ok) {
      const j = await res.json();

      setProjects(j.projects ?? []);
    }
  }, []);

  const loadAll = useCallback(async () => {
    setLoading(true);
    try {
      await Promise.all([loadDatasets(), loadProjects()]);
    } finally {
      setLoading(false);
    }
  }, [loadDatasets, loadProjects]);

  useEffect(() => {
    loadAll();
  }, [loadAll]);

  // Map dataset id → project titles it belongs to (for the card subtitle).
  const projectNamesByDataset = useMemo(() => {
    const map = new Map<string, string[]>();

    for (const p of projects) {
      for (const d of p.datasets) {
        const list = map.get(d.id) ?? [];

        list.push(p.title);
        map.set(d.id, list);
      }
    }

    return map;
  }, [projects]);

  const filteredDatasets = useMemo(
    () => applyDatasetFilters(datasets, filters),
    [datasets, filters],
  );

  // ---- Metadata save (edit dataset / edit project / create project) ----
  const handleSave = async (draft: MetadataDraft): Promise<boolean> => {
    if (!editTarget) return false;

    try {
      if (editTarget.kind === "dataset") {
        const res = await fetch(`/api/ingest/mine/${editTarget.data.id}`, {
          method: "PATCH",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            title: draft.title,
            description: draft.description,
            species: draft.species,
            disease: draft.disease,
            tissue: draft.tissue,
            platform: draft.platform,
            institute: draft.institute,
            tags: draft.tags,
            externalLink: draft.externalLink,
            publicationLink: draft.publicationLink,
            metadata: draft.metadata,
          }),
        });

        if (!res.ok) return false;
        await loadDatasets();

        return true;
      }

      const isNew = editTarget.kind === "new-project";
      const url = isNew
        ? "/api/projects"
        : `/api/projects/${(editTarget as { data: ProjectRow }).data.id}`;
      const res = await fetch(url, {
        method: isNew ? "POST" : "PATCH",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(draft),
      });

      if (!res.ok) return false;
      await loadProjects();

      return true;
    } catch {
      return false;
    }
  };

  // ---- Submit / withdraw ----
  const submit = async (type: "dataset" | "project", id: string) => {
    const res = await fetch("/api/community/submit", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ type, id }),
    });

    if (res.ok) {
      const j = await res.json();

      toast.success(
        j.reviewStatus === "APPROVED"
          ? "Submission updated — live on Explore."
          : "Submitted for review. An admin will approve it soon.",
      );
      // Only the affected list needs to refresh its submission status.
      await (type === "project" ? loadProjects() : loadDatasets());
    } else {
      const j = await res.json().catch(() => ({}));

      toast.error(j.error ?? "Couldn't submit.");
    }
  };

  const withdraw = async (type: "dataset" | "project", id: string) => {
    if (!confirm("Remove this from Explore / cancel the submission?")) return;
    const res = await fetch("/api/community/submit", {
      method: "DELETE",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ type, id }),
    });

    if (res.ok) {
      toast.success("Withdrawn.");
      await (type === "project" ? loadProjects() : loadDatasets());
    } else {
      toast.error("Couldn't withdraw.");
    }
  };

  // ---- Project membership ----
  // A lightweight member summary for optimistic updates (matches the shape the
  // API returns for a project's `datasets`).
  const summaryFor = (datasetId: string) => {
    const d = datasets.find((x) => x.id === datasetId);

    return d
      ? {
          id: d.id,
          title: d.title,
          datasetType: d.datasetType,
          thumbnailUrl: d.thumbnailUrl,
          status: d.status,
        }
      : null;
  };

  const toggleMembership = async (
    projectId: string,
    datasetId: string,
    member: boolean,
  ): Promise<boolean> => {
    const summary = summaryFor(datasetId);

    // Optimistic: reflect the change locally right away — the card badge and
    // the modal checkbox update instantly, no full refetch on the happy path.
    setProjects((prev) =>
      prev.map((p) => {
        if (p.id !== projectId) return p;
        if (member) {
          if (p.datasets.some((d) => d.id === datasetId)) return p;

          return {
            ...p,
            datasets: summary ? [...p.datasets, summary] : p.datasets,
            datasetCount: p.datasetCount + 1,
          };
        }

        return {
          ...p,
          datasets: p.datasets.filter((d) => d.id !== datasetId),
          datasetCount: Math.max(0, p.datasetCount - 1),
        };
      }),
    );

    const res = await fetch(
      member
        ? `/api/projects/${projectId}/datasets`
        : `/api/projects/${projectId}/datasets/${datasetId}`,
      {
        method: member ? "POST" : "DELETE",
        headers: { "Content-Type": "application/json" },
        body: member ? JSON.stringify({ datasetId }) : undefined,
      },
    );

    if (!res.ok) {
      toast.error("Couldn't update project.");
      // Reconcile with server truth on failure.
      await loadProjects();

      return false;
    }

    return true;
  };

  const createAndAdd = async (
    title: string,
    datasetId: string,
  ): Promise<boolean> => {
    const res = await fetch("/api/projects", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ title }),
    });

    if (!res.ok) {
      toast.error("Couldn't create project.");

      return false;
    }
    const project = await res.json();
    const summary = summaryFor(datasetId);

    // Optimistically insert the new project (with this dataset) — no refetch.
    setProjects((prev) => [
      {
        ...project,
        tags: project.tags ?? [],
        datasetCount: summary ? 1 : 0,
        datasets: summary ? [summary] : [],
        submission: null,
      },
      ...prev,
    ]);

    const addRes = await fetch(`/api/projects/${project.id}/datasets`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ datasetId }),
    });

    if (!addRes.ok) {
      // Project was created but the dataset couldn't be attached — reconcile.
      await loadProjects();
    }

    return true;
  };

  const deleteProject = async (project: ProjectRow) => {
    if (
      !confirm(
        `Delete project "${project.title}"? The datasets themselves are not deleted.`,
      )
    )
      return;
    const res = await fetch(`/api/projects/${project.id}`, {
      method: "DELETE",
    });

    if (res.ok) {
      toast.success("Project deleted.");
      await loadProjects();
    } else {
      toast.error("Couldn't delete project.");
    }
  };

  const editInitial: Partial<MetadataDraft> =
    editTarget && editTarget.kind !== "new-project"
      ? {
          ...editTarget.data,
          title: editTarget.data.title ?? "",
          metadata: editTarget.data.metadata ?? {},
        }
      : {};

  return (
    <>
      <Tabs aria-label="Account tabs" variant="underlined">
        <Tab key="datasets" title={`Datasets (${datasets.length})`}>
          {loading ? (
            <div className="flex justify-center py-16">
              <Spinner label="Loading…" />
            </div>
          ) : datasets.length === 0 ? (
            <p className="text-default-500 py-12 text-center">
              You don&apos;t own any datasets yet.
            </p>
          ) : (
            <div className="mt-2">
              <AccountFilterBar
                datasets={datasets}
                filters={filters}
                shown={filteredDatasets.length}
                onChange={setFilters}
              />
              {filteredDatasets.length === 0 ? (
                <p className="text-default-500 py-12 text-center">
                  No datasets match your filters.
                </p>
              ) : (
                <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
                  {filteredDatasets.map((d) => (
                    <AccountDatasetCard
                      key={d.id}
                      dataset={d}
                      projectNames={projectNamesByDataset.get(d.id) ?? []}
                      onAddToProject={setAddToProjectFor}
                      onEdit={(data) =>
                        setEditTarget({ kind: "dataset", data })
                      }
                      onSubmit={(data) => submit("dataset", data.id)}
                      onWithdraw={(data) => withdraw("dataset", data.id)}
                    />
                  ))}
                </div>
              )}
            </div>
          )}
        </Tab>

        <Tab key="projects" title={`Projects (${projects.length})`}>
          <div className="flex justify-end mt-2 mb-4">
            <Button
              color="primary"
              size="sm"
              onPress={() => setEditTarget({ kind: "new-project" })}
            >
              New Project
            </Button>
          </div>
          {loading ? (
            <div className="flex justify-center py-16">
              <Spinner label="Loading…" />
            </div>
          ) : projects.length === 0 ? (
            <p className="text-default-500 py-12 text-center">
              No projects yet. Create one to group several datasets into a
              single Explore card.
            </p>
          ) : (
            <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
              {projects.map((p) => (
                <AccountProjectCard
                  key={p.id}
                  project={p}
                  onDelete={deleteProject}
                  onEdit={(data) => setEditTarget({ kind: "project", data })}
                  onSubmit={(data) => submit("project", data.id)}
                  onWithdraw={(data) => withdraw("project", data.id)}
                />
              ))}
            </div>
          )}
        </Tab>
      </Tabs>

      <MetadataEditModal
        heading={
          editTarget?.kind === "dataset"
            ? "Edit dataset"
            : editTarget?.kind === "project"
              ? "Edit project"
              : "New project"
        }
        initial={editInitial}
        isOpen={editTarget !== null}
        showThumbnail={editTarget?.kind !== "dataset"}
        onClose={() => setEditTarget(null)}
        onSave={handleSave}
      />

      <AddToProjectModal
        dataset={addToProjectFor}
        isOpen={addToProjectFor !== null}
        projects={projects}
        onClose={() => setAddToProjectFor(null)}
        onCreateAndAdd={createAndAdd}
        onToggle={toggleMembership}
      />
    </>
  );
}
