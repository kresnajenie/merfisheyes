"use client";

import { useCallback, useEffect, useMemo, useState } from "react";
import { Button } from "@heroui/button";
import { Chip } from "@heroui/chip";

import {
  RECOGNISED_TECHNOLOGIES,
  TECHNOLOGY_PRESETS,
  readPreset,
  type TechnologyPreset,
} from "@/lib/config/technology-presets";
import { COLORMAPS } from "@/lib/utils/colormaps";
import { decodeCellVizState } from "@/lib/utils/url-visualization-state";

interface Person {
  name: string | null;
  email: string;
}
interface Version {
  hash: string;
  technology: string;
  config: TechnologyPreset;
  note: string | null;
  testedOn: string | null;
  appVersion: string;
  createdAt: string;
  createdBy: Person | null;
}
interface Publish {
  id: string;
  technology: string;
  versionHash: string | null;
  createdAt: string;
  publishedBy: Person | null;
}

const field =
  "w-full rounded-lg border border-default-200 bg-transparent px-3 py-2 text-sm";

function describe(config: TechnologyPreset): string {
  const parts = [
    config.colormap && `colormap ${config.colormap}`,
    config.hideZeroExpression !== undefined &&
      (config.hideZeroExpression ? "zero counts hidden" : "zero counts shown"),
    config.cellOpacity !== undefined &&
      `opacity ${Math.round(config.cellOpacity * 100)}%`,
  ].filter(Boolean);

  return parts.length ? parts.join(", ") : "nothing set";
}

const by = (p: Person | null) => p?.name || p?.email || "someone";
const when = (iso: string) => new Date(iso).toLocaleString();

/** A viewer link with ?preset=<hash> added (or replaced). */
function previewLink(viewerLink: string, hash: string): string | null {
  try {
    const url = new URL(viewerLink, window.location.origin);

    url.searchParams.set("preset", hash);

    return url.toString();
  } catch {
    return null;
  }
}

export function PresetEditor() {
  const [technology, setTechnology] = useState(RECOGNISED_TECHNOLOGIES[0]);
  const [versions, setVersions] = useState<Version[]>([]);
  const [publishes, setPublishes] = useState<Publish[]>([]);

  // The form. "" = not set: the field is left out of the config.
  const [colormap, setColormap] = useState("");
  const [hideZero, setHideZero] = useState("");
  const [opacity, setOpacity] = useState("");
  const [note, setNote] = useState("");
  const [testedOn, setTestedOn] = useState("");

  const [message, setMessage] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);

  const refresh = useCallback(async () => {
    const r = await fetch("/api/presets/versions");

    if (!r.ok) return;
    const j = await r.json();

    setVersions(j.versions);
    setPublishes(j.publishes);
  }, []);

  useEffect(() => {
    refresh();
  }, [refresh]);

  const fill = (config: TechnologyPreset) => {
    setColormap(config.colormap ?? "");
    setHideZero(
      config.hideZeroExpression === undefined
        ? ""
        : String(config.hideZeroExpression),
    );
    setOpacity(
      config.cellOpacity === undefined ? "" : String(config.cellOpacity),
    );
  };

  const config = useMemo(() => {
    const c: TechnologyPreset = {};

    if (colormap) c.colormap = colormap;
    if (hideZero) c.hideZeroExpression = hideZero === "true";
    if (opacity !== "") c.cellOpacity = Number(opacity);

    return c;
  }, [colormap, hideZero, opacity]);

  const mine = versions.filter((v) => v.technology === technology);
  const log = publishes.filter((p) => p.technology === technology);
  // The latest publish row decides; none at all = the built-in preset.
  const inForce = log[0];
  const builtIn = TECHNOLOGY_PRESETS[technology];

  const post = async (path: string, body: object) => {
    setBusy(true);
    setMessage(null);
    const r = await fetch(path, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
    const j = await r.json().catch(() => ({}));

    setBusy(false);
    if (!r.ok) {
      setMessage(j.error ?? `Failed (${r.status})`);

      return null;
    }
    await refresh();

    return j;
  };

  const save = async () => {
    const j = await post("/api/presets/versions", {
      technology,
      config,
      note,
      testedOn,
    });

    if (j) setMessage(`Saved as draft ${j.version.hash}. Preview it below.`);
  };

  const publish = async (hash: string | null) => {
    const j = await post("/api/presets/publish", { technology, hash });

    if (j)
      setMessage(
        hash
          ? `${hash} is now what ${technology} opens with.`
          : `${technology} now opens with no preset.`,
      );
  };

  const fillFromLink = (link: string) => {
    try {
      const v = new URL(link, window.location.origin).searchParams.get("v");
      const state = v ? decodeCellVizState(v) : null;

      if (!state) {
        setMessage("That link carries no viewer settings.");

        return;
      }
      fill(
        readPreset({
          colormap: state.cm,
          hideZeroExpression: state.hz,
          cellOpacity: state.op,
        }),
      );
      setTestedOn(link);
      setMessage(
        "Filled from the link. A setting the link does not carry is left not set.",
      );
    } catch {
      setMessage("That is not a link.");
    }
  };

  const fillFromFile = async (file: File) => {
    try {
      fill(readPreset(JSON.parse(await file.text())));
      setMessage(`Filled from ${file.name}.`);
    } catch {
      setMessage(`${file.name} is not JSON.`);
    }
  };

  return (
    <div className="flex flex-col gap-8">
      <section className="flex flex-col gap-2">
        <label className="text-sm font-medium" htmlFor="preset-technology">
          Technology
        </label>
        <select
          className={field}
          data-testid="preset-technology"
          id="preset-technology"
          value={technology}
          onChange={(e) => setTechnology(e.target.value)}
        >
          {RECOGNISED_TECHNOLOGIES.map((t) => (
            <option key={t} value={t}>
              {t}
            </option>
          ))}
        </select>
        <p className="text-sm text-default-500" data-testid="preset-in-force">
          In force now:{" "}
          {!inForce
            ? builtIn
              ? `the built-in preset (${describe(builtIn)})`
              : "none"
            : inForce.versionHash
              ? `${inForce.versionHash} (${describe(
                  mine.find((v) => v.hash === inForce.versionHash)?.config ??
                    {},
                )})`
              : "turned off"}
        </p>
      </section>

      <section className="flex flex-col gap-4">
        <h2 className="text-lg font-semibold">New version</h2>
        <div className="grid grid-cols-1 gap-4 sm:grid-cols-3">
          <label className="flex flex-col gap-1 text-sm">
            Colormap
            <select
              className={field}
              data-testid="preset-colormap"
              value={colormap}
              onChange={(e) => setColormap(e.target.value)}
            >
              <option value="">Not set</option>
              {Object.values(COLORMAPS).map((c) => (
                <option key={c.name} value={c.name}>
                  {c.label}
                </option>
              ))}
            </select>
          </label>
          <label className="flex flex-col gap-1 text-sm">
            Zero counts
            <select
              className={field}
              data-testid="preset-hide-zero"
              value={hideZero}
              onChange={(e) => setHideZero(e.target.value)}
            >
              <option value="">Not set</option>
              <option value="true">Hidden</option>
              <option value="false">Shown</option>
            </select>
          </label>
          <label className="flex flex-col gap-1 text-sm">
            Cell opacity (0.05 to 1)
            <input
              className={field}
              data-testid="preset-opacity"
              max={1}
              min={0.05}
              placeholder="Not set"
              step={0.05}
              type="number"
              value={opacity}
              onChange={(e) => setOpacity(e.target.value)}
            />
          </label>
        </div>
        <label className="flex flex-col gap-1 text-sm">
          Viewer link of a store to try it on (kept with the version)
          <input
            className={field}
            data-testid="preset-tested-on"
            placeholder="https://…/spatialdata-viewer/…"
            value={testedOn}
            onChange={(e) => setTestedOn(e.target.value)}
          />
        </label>
        <label className="flex flex-col gap-1 text-sm">
          Note
          <input
            className={field}
            data-testid="preset-note"
            value={note}
            onChange={(e) => setNote(e.target.value)}
          />
        </label>
        <div className="flex flex-wrap items-center gap-3">
          <Button
            color="primary"
            data-testid="preset-save"
            isDisabled={busy}
            onPress={save}
          >
            Save draft
          </Button>
          <Button
            isDisabled={!testedOn}
            variant="flat"
            onPress={() => fillFromLink(testedOn)}
          >
            Fill from that link&apos;s settings
          </Button>
          <label className="cursor-pointer text-sm text-primary">
            Fill from a JSON file
            <input
              accept="application/json,.json"
              className="hidden"
              data-testid="preset-file"
              type="file"
              onChange={(e) => {
                const f = e.target.files?.[0];

                if (f) fillFromFile(f);
                e.target.value = "";
              }}
            />
          </label>
        </div>
        {message && (
          <p className="text-sm" data-testid="preset-message" role="status">
            {message}
          </p>
        )}
      </section>

      <section className="flex flex-col gap-3">
        <div className="flex items-center justify-between">
          <h2 className="text-lg font-semibold">Versions</h2>
          <Button
            data-testid="preset-turn-off"
            isDisabled={busy || (!!inForce && !inForce.versionHash)}
            size="sm"
            variant="flat"
            onPress={() => publish(null)}
          >
            Turn the preset off
          </Button>
        </div>
        {mine.length === 0 && (
          <p className="text-sm text-default-500">No versions saved yet.</p>
        )}
        {mine.map((v) => {
          const link = v.testedOn && previewLink(v.testedOn, v.hash);
          const current = inForce?.versionHash === v.hash;

          return (
            <div
              key={v.hash}
              className="flex flex-col gap-2 rounded-xl border border-default-200 p-4"
              data-testid={`preset-version-${v.hash}`}
            >
              <div className="flex flex-wrap items-center gap-2">
                <code className="text-sm">{v.hash}</code>
                {current && (
                  <Chip color="success" size="sm" variant="flat">
                    In force
                  </Chip>
                )}
                <span className="text-sm">{describe(v.config)}</span>
              </div>
              <p className="text-xs text-default-500">
                Saved by {by(v.createdBy)} on {when(v.createdAt)}, app version{" "}
                {v.appVersion}
                {v.testedOn ? `, tried on ${v.testedOn}` : ""}
                {v.note ? `. ${v.note}` : ""}
              </p>
              <div className="flex flex-wrap gap-2">
                {link ? (
                  <Button
                    as="a"
                    data-testid={`preset-preview-${v.hash}`}
                    href={link}
                    rel="noreferrer"
                    size="sm"
                    target="_blank"
                    variant="flat"
                  >
                    Preview
                  </Button>
                ) : (
                  <span className="self-center text-xs text-default-500">
                    To preview, add ?preset={v.hash} to any viewer link.
                  </span>
                )}
                <Button size="sm" variant="flat" onPress={() => fill(v.config)}>
                  Edit a copy
                </Button>
                <Button
                  color="primary"
                  data-testid={`preset-publish-${v.hash}`}
                  isDisabled={busy || current}
                  size="sm"
                  onPress={() => publish(v.hash)}
                >
                  Publish
                </Button>
              </div>
            </div>
          );
        })}
      </section>

      <section className="flex flex-col gap-2">
        <h2 className="text-lg font-semibold">History</h2>
        {log.length === 0 && (
          <p className="text-sm text-default-500">Nothing published yet.</p>
        )}
        <ul
          className="flex flex-col gap-1 text-sm"
          data-testid="preset-history"
        >
          {log.map((p) => (
            <li key={p.id}>
              {when(p.createdAt)}: {by(p.publishedBy)}{" "}
              {p.versionHash ? `published ${p.versionHash}` : "turned it off"}
            </li>
          ))}
        </ul>
      </section>
    </div>
  );
}
