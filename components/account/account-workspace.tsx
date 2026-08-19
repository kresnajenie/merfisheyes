"use client";

import type { MetadataDraft } from "./metadata-edit-modal";
import type { DatasetRow, ProjectRow } from "./types";
import type { AccountFacets } from "./account-filter-bar";

import { useCallback, useEffect, useRef, useState } from "react";
import { Tabs, Tab } from "@heroui/tabs";
import { Button } from "@heroui/button";
import { Input } from "@heroui/input";
import { Spinner } from "@heroui/spinner";
import { Pagination } from "@heroui/pagination";
import { useRouter, useSearchParams } from "next/navigation";
import { toast } from "react-toastify";

import { AccountDatasetCard } from "./account-dataset-card";
import { AccountProjectCard } from "./account-project-card";
import { MetadataEditModal } from "./metadata-edit-modal";
import { AddToProjectModal } from "./add-to-project-modal";
import { AccountFilterBar } from "./account-filter-bar";

type EditTarget =
  | { kind: "dataset"; data: DatasetRow }
  | { kind: "project"; data: ProjectRow }
  | { kind: "new-project" }
  | null;

const PAGE_LIMIT = 20;
const EMPTY_FACETS: AccountFacets = { species: [], tissues: [], platforms: [] };

export function AccountWorkspace() {
  const router = useRouter();
  const searchParams = useSearchParams();

  // ---- URL-backed view / filter / page state ----
  const view =
    searchParams.get("view") === "projects" ? "projects" : "datasets";
  const [search, setSearch] = useState(searchParams.get("q") ?? "");
  const [debouncedSearch, setDebouncedSearch] = useState(
    searchParams.get("q") ?? "",
  );

  useEffect(() => {
    const t = setTimeout(() => setDebouncedSearch(search), 250);

    return () => clearTimeout(t);
  }, [search]);

  // Gene search — debounced like the text search, same contract as Explore.
  const [geneSearch, setGeneSearch] = useState(searchParams.get("gene") ?? "");
  const [debouncedGeneSearch, setDebouncedGeneSearch] = useState(
    searchParams.get("gene") ?? "",
  );

  useEffect(() => {
    const t = setTimeout(() => setDebouncedGeneSearch(geneSearch), 250);

    return () => clearTimeout(t);
  }, [geneSearch]);

  const [type, setType] = useState(searchParams.get("type") ?? "");
  const [status, setStatus] = useState(searchParams.get("status") ?? "");
  const [species, setSpecies] = useState(searchParams.get("species") ?? "");
  const [tissue, setTissue] = useState(searchParams.get("tissue") ?? "");
  const [platform, setPlatform] = useState(searchParams.get("platform") ?? "");
  const [page, setPage] = useState(Number(searchParams.get("page")) || 1);

  // ---- Data ----
  const [datasets, setDatasets] = useState<DatasetRow[]>([]);
  const [datasetsTotal, setDatasetsTotal] = useState(0);
  const [facets, setFacets] = useState<AccountFacets>(EMPTY_FACETS);
  const [projects, setProjects] = useState<ProjectRow[]>([]);
  const [projectsTotal, setProjectsTotal] = useState(0);
  // Full (unpaginated) project list for the "Add to project" picker.
  const [allProjects, setAllProjects] = useState<ProjectRow[]>([]);
  const [loading, setLoading] = useState(true);

  const [editTarget, setEditTarget] = useState<EditTarget>(null);
  const [addToProjectFor, setAddToProjectFor] = useState<DatasetRow | null>(
    null,
  );

  // Reset page to 1 on a real filter/view change (snapshot compare, so mount
  // and Strict Mode remounts don't clobber a deep-linked ?page=N).
  const prevKey = useRef<string | null>(null);
  const filterKey = JSON.stringify([
    view,
    debouncedSearch,
    debouncedGeneSearch,
    type,
    status,
    species,
    tissue,
    platform,
  ]);

  useEffect(() => {
    if (prevKey.current === null) {
      prevKey.current = filterKey;

      return;
    }
    if (prevKey.current !== filterKey) {
      prevKey.current = filterKey;
      setPage(1);
    }
  }, [filterKey]);

  // ---- URL sync (replace, not push) ----
  const isInitialMount = useRef(true);

  useEffect(() => {
    if (isInitialMount.current) {
      isInitialMount.current = false;

      return;
    }
    const params = new URLSearchParams();

    if (view !== "datasets") params.set("view", view);
    if (debouncedSearch) params.set("q", debouncedSearch);
    if (debouncedGeneSearch) params.set("gene", debouncedGeneSearch);
    if (type) params.set("type", type);
    if (status) params.set("status", status);
    if (species) params.set("species", species);
    if (tissue) params.set("tissue", tissue);
    if (platform) params.set("platform", platform);
    if (page > 1) params.set("page", String(page));

    const qs = params.toString();

    router.replace(`/account${qs ? `?${qs}` : ""}`, { scroll: false });
  }, [
    view,
    debouncedSearch,
    debouncedGeneSearch,
    type,
    status,
    species,
    tissue,
    platform,
    page,
    router,
  ]);

  // ---- Fetchers ----
  const loadDatasets = useCallback(async () => {
    const params = new URLSearchParams();

    params.set("page", String(page));
    params.set("limit", String(PAGE_LIMIT));
    if (debouncedSearch) params.set("search", debouncedSearch);
    if (debouncedGeneSearch) params.set("genes", debouncedGeneSearch);
    if (type) params.set("type", type);
    if (status) params.set("status", status);
    if (species) params.set("species", species);
    if (tissue) params.set("tissue", tissue);
    if (platform) params.set("platform", platform);

    const res = await fetch(`/api/ingest/mine?${params}`, {
      cache: "no-store",
    });

    if (res.ok) {
      const j = await res.json();

      setDatasets(j.datasets ?? []);
      setDatasetsTotal(j.total ?? 0);
      setFacets(j.filters ?? EMPTY_FACETS);
    }
  }, [
    page,
    debouncedSearch,
    debouncedGeneSearch,
    type,
    status,
    species,
    tissue,
    platform,
  ]);

  const loadProjects = useCallback(async () => {
    const params = new URLSearchParams();

    params.set("page", String(page));
    params.set("limit", String(PAGE_LIMIT));
    if (debouncedSearch) params.set("search", debouncedSearch);

    const res = await fetch(`/api/projects?${params}`, { cache: "no-store" });

    if (res.ok) {
      const j = await res.json();

      setProjects(j.projects ?? []);
      setProjectsTotal(j.total ?? 0);
    }
  }, [page, debouncedSearch]);

  // The full project list (no pagination) feeds the Add-to-project picker.
  const loadAllProjects = useCallback(async () => {
    const res = await fetch("/api/projects", { cache: "no-store" });

    if (res.ok) {
      const j = await res.json();

      setAllProjects(j.projects ?? []);
    }
  }, []);

  useEffect(() => {
    loadAllProjects();
  }, [loadAllProjects]);

  // Reload the active view whenever its inputs change.
  useEffect(() => {
    let active = true;

    setLoading(true);
    (view === "projects" ? loadProjects() : loadDatasets()).finally(() => {
      if (active) setLoading(false);
    });

    return () => {
      active = false;
    };
  }, [view, loadDatasets, loadProjects]);

  const setView = (next: "datasets" | "projects") => {
    const params = new URLSearchParams(searchParams.toString());

    if (next === "datasets") params.delete("view");
    else params.set("view", next);
    params.delete("page");
    router.replace(`/account${params.toString() ? `?${params}` : ""}`, {
      scroll: false,
    });
  };

  const datasetsPages = Math.max(1, Math.ceil(datasetsTotal / PAGE_LIMIT));
  const projectsPages = Math.max(1, Math.ceil(projectsTotal / PAGE_LIMIT));

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
      await Promise.all([loadProjects(), loadAllProjects()]);

      return true;
    } catch {
      return false;
    }
  };

  // ---- Submit / withdraw ----
  const submit = async (kind: "dataset" | "project", id: string) => {
    const res = await fetch("/api/community/submit", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ type: kind, id }),
    });

    if (res.ok) {
      const j = await res.json();

      toast.success(
        j.reviewStatus === "APPROVED"
          ? "Submission updated — live on Explore."
          : "Submitted for review. An admin will approve it soon.",
      );
      await (kind === "project" ? loadProjects() : loadDatasets());
    } else {
      const j = await res.json().catch(() => ({}));

      toast.error(j.error ?? "Couldn't submit.");
    }
  };

  const withdraw = async (kind: "dataset" | "project", id: string) => {
    if (!confirm("Remove this from Explore / cancel the submission?")) return;
    const res = await fetch("/api/community/submit", {
      method: "DELETE",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ type: kind, id }),
    });

    if (res.ok) {
      toast.success("Withdrawn.");
      await (kind === "project" ? loadProjects() : loadDatasets());
    } else {
      toast.error("Couldn't withdraw.");
    }
  };

  // ---- Project membership (optimistic) ----
  // Patch the local caches (both project lists + the dataset's badge names)
  // so the UI flips instantly; the API call settles in the background and a
  // failure applies the inverse patch. No refetches.
  const applyMembership = useCallback(
    (projectId: string, datasetId: string, member: boolean) => {
      const ds = datasets.find((d) => d.id === datasetId);
      const patchProject = (p: ProjectRow): ProjectRow => {
        if (p.id !== projectId) return p;
        const next = member
          ? p.datasets.some((d) => d.id === datasetId)
            ? p.datasets
            : [
                ...p.datasets,
                {
                  id: datasetId,
                  title: ds?.title ?? null,
                  datasetType: ds?.datasetType ?? null,
                  thumbnailUrl: ds?.thumbnailUrl ?? null,
                  status: ds?.status ?? "COMPLETE",
                },
              ]
          : p.datasets.filter((d) => d.id !== datasetId);

        return { ...p, datasets: next, datasetCount: next.length };
      };

      setProjects((ps) => ps.map(patchProject));
      setAllProjects((ps) => ps.map(patchProject));

      const projectTitle = allProjects.find((p) => p.id === projectId)?.title;

      if (projectTitle) {
        setDatasets((rows) =>
          rows.map((d) => {
            if (d.id !== datasetId) return d;
            const names = member
              ? d.projectNames.includes(projectTitle)
                ? d.projectNames
                : [...d.projectNames, projectTitle]
              : d.projectNames.filter((n) => n !== projectTitle);

            return { ...d, projectNames: names };
          }),
        );
      }
    },
    [datasets, allProjects],
  );

  const toggleMembership = async (
    projectId: string,
    datasetId: string,
    member: boolean,
  ): Promise<boolean> => {
    applyMembership(projectId, datasetId, member);

    const res = await fetch(
      member
        ? `/api/projects/${projectId}/datasets`
        : `/api/projects/${projectId}/datasets/${datasetId}`,
      {
        method: member ? "POST" : "DELETE",
        headers: { "Content-Type": "application/json" },
        body: member ? JSON.stringify({ datasetId }) : undefined,
      },
    ).catch(() => null);

    if (!res?.ok) {
      applyMembership(projectId, datasetId, !member);
      toast.error("Couldn't update project.");

      return false;
    }

    // After ADDING a dataset, jump straight to the project page (the API has
    // settled, so the server-rendered page includes the new member).
    if (member) {
      setAddToProjectFor(null);
      router.push(`/account/projects/${projectId}`);
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

    await fetch(`/api/projects/${project.id}/datasets`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ datasetId }),
    });
    // Go straight to the new project's page — no list refetches needed, we're
    // leaving this view.
    setAddToProjectFor(null);
    router.push(`/account/projects/${project.id}`);

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
      await Promise.all([loadProjects(), loadAllProjects()]);
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

  const pager = (pages: number) =>
    pages > 1 ? (
      <div className="flex justify-center py-6">
        <Pagination
          showControls
          page={page}
          size="sm"
          total={pages}
          onChange={setPage}
        />
      </div>
    ) : null;

  return (
    <>
      <Tabs
        aria-label="Account tabs"
        selectedKey={view}
        variant="underlined"
        onSelectionChange={(k) => setView(k as "datasets" | "projects")}
      >
        <Tab key="datasets" title="Datasets">
          <div className="mt-2">
            <AccountFilterBar
              facets={facets}
              geneSearch={geneSearch}
              platform={platform}
              search={search}
              species={species}
              status={status}
              tissue={tissue}
              type={type}
              onGeneSearchChange={setGeneSearch}
              onPlatformChange={setPlatform}
              onSearchChange={setSearch}
              onSpeciesChange={setSpecies}
              onStatusChange={setStatus}
              onTissueChange={setTissue}
              onTypeChange={setType}
            />

            <p className="mb-4 text-sm text-default-500">
              {datasetsTotal} dataset{datasetsTotal === 1 ? "" : "s"}
            </p>

            {loading ? (
              <div className="flex justify-center py-16">
                <Spinner label="Loading…" />
              </div>
            ) : datasets.length === 0 ? (
              <p className="py-12 text-center text-default-500">
                No datasets match your filters.
              </p>
            ) : (
              <>
                <div className="grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4">
                  {datasets.map((d) => (
                    <AccountDatasetCard
                      key={d.id}
                      dataset={d}
                      geneHighlight={debouncedGeneSearch || undefined}
                      projectNames={d.projectNames ?? []}
                      onAddToProject={setAddToProjectFor}
                      onEdit={(data) =>
                        setEditTarget({ kind: "dataset", data })
                      }
                      onSubmit={(data) => submit("dataset", data.id)}
                      onWithdraw={(data) => withdraw("dataset", data.id)}
                    />
                  ))}
                </div>
                {pager(datasetsPages)}
              </>
            )}
          </div>
        </Tab>

        <Tab key="projects" title="Projects">
          <div className="mt-2">
            <div className="mb-4 flex items-center gap-2">
              <Input
                classNames={{ inputWrapper: "bg-default-100" }}
                placeholder="Search your projects…"
                value={search}
                onValueChange={setSearch}
              />
              <Button
                color="primary"
                onPress={() => setEditTarget({ kind: "new-project" })}
              >
                New Project
              </Button>
            </div>

            <p className="mb-4 text-sm text-default-500">
              {projectsTotal} project{projectsTotal === 1 ? "" : "s"}
            </p>

            {loading ? (
              <div className="flex justify-center py-16">
                <Spinner label="Loading…" />
              </div>
            ) : projects.length === 0 ? (
              <p className="py-12 text-center text-default-500">
                No projects yet. Create one to group several datasets into a
                single Explore card.
              </p>
            ) : (
              <>
                <div className="grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4">
                  {projects.map((p) => (
                    <AccountProjectCard
                      key={p.id}
                      project={p}
                      onDelete={deleteProject}
                      onEdit={(data) =>
                        setEditTarget({ kind: "project", data })
                      }
                      onSubmit={(data) => submit("project", data.id)}
                      onWithdraw={(data) => withdraw("project", data.id)}
                    />
                  ))}
                </div>
                {pager(projectsPages)}
              </>
            )}
          </div>
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
        projects={allProjects}
        onClose={() => setAddToProjectFor(null)}
        onCreateAndAdd={createAndAdd}
        onToggle={toggleMembership}
      />
    </>
  );
}
