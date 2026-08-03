"use client";

import { useEffect, useMemo, useState } from "react";
import {
  Modal,
  ModalContent,
  ModalHeader,
  ModalBody,
  ModalFooter,
} from "@heroui/modal";
import { Button } from "@heroui/button";
import { Chip } from "@heroui/chip";
import { Spinner } from "@heroui/spinner";
import { useSession } from "next-auth/react";

import { OwnerChoice, type OwnerValue } from "@/components/owner-choice";

type Mode = "add" | "overwrite" | "replace_all";

interface ImportResult {
  mode: Mode;
  created: { id: string; title: string }[];
  updated: { id: string; title: string }[];
  skipped: { title: string; reason: string }[];
  errors: { title?: string; error: string }[];
  deleted?: { id: string; title: string }[];
}

interface CatalogImportModalProps {
  isOpen: boolean;
  onClose: () => void;
  /** Called after a successful import so the parent table can refresh. */
  onImported?: () => void;
}

/**
 * Small bank of low-ambiguity words used for type-to-confirm prompts. Kept
 * short so the typing isn't punishing, but distinct enough that the user
 * actually has to *read* each prompt instead of muscle-memory-clicking.
 */
const CONFIRM_WORDS = [
  "orange",
  "falcon",
  "river",
  "mango",
  "violet",
  "ember",
  "harbor",
  "lantern",
  "pebble",
  "willow",
  "saffron",
  "marble",
  "thunder",
  "amber",
  "raven",
  "comet",
  "compass",
  "summit",
  "garnet",
  "hazel",
];

/** Per-mode number of type-to-confirm prompts before Import is enabled. */
const CONFIRM_COUNT: Record<Exclude<Mode, "add">, number> = {
  overwrite: 3,
  replace_all: 5,
};

function pickWords(n: number): string[] {
  const pool = [...CONFIRM_WORDS];

  // In-place shuffle, then slice — avoids picking the same word twice.
  for (let i = pool.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));

    [pool[i], pool[j]] = [pool[j], pool[i]];
  }

  return pool.slice(0, n);
}

/**
 * Import-from-JSON modal for the admin catalog page.
 *
 *   - add (default)   skips rows that already exist
 *   - overwrite       updates matched rows in place; needs 3 typed confirms
 *   - replace_all     same as overwrite + deletes rows NOT in the JSON;
 *                     SUPER_ADMIN only; needs 5 typed confirms
 *
 * Each confirmation step asks the user to type a distinct random word.
 * Picking a different file or mode wipes the in-progress confirmations.
 */
export function CatalogImportModal({
  isOpen,
  onClose,
  onImported,
}: CatalogImportModalProps) {
  const { data: session } = useSession();
  const isSuperAdmin = session?.user?.role === "SUPER_ADMIN";

  const [file, setFile] = useState<File | null>(null);
  const [parsedItemCount, setParsedItemCount] = useState<number | null>(null);
  const [parseError, setParseError] = useState<string | null>(null);
  const [dragOver, setDragOver] = useState(false);
  const [mode, setMode] = useState<Mode>("add");
  const [owner, setOwner] = useState<OwnerValue>("me");

  // Confirmation state. `words` is the (frozen-per-mode) list of words the
  // user must type, and `step` is how many they've already cleared (0..N).
  const [words, setWords] = useState<string[]>([]);
  const [step, setStep] = useState(0);
  const [typed, setTyped] = useState("");

  const [submitting, setSubmitting] = useState(false);
  const [result, setResult] = useState<ImportResult | null>(null);
  const [submitError, setSubmitError] = useState<string | null>(null);

  // Whenever the user enters a destructive mode (or relands on it), generate
  // a fresh set of random words. Resets in-flight typing too.
  useEffect(() => {
    if (mode === "add") {
      setWords([]);
      setStep(0);
      setTyped("");

      return;
    }
    setWords(pickWords(CONFIRM_COUNT[mode]));
    setStep(0);
    setTyped("");
  }, [mode]);

  const requiredConfirms = mode === "add" ? 0 : CONFIRM_COUNT[mode];
  const allConfirmed = step >= requiredConfirms;

  const resetAll = () => {
    setFile(null);
    setParsedItemCount(null);
    setParseError(null);
    setDragOver(false);
    setMode("add");
    setOwner("me");
    setWords([]);
    setStep(0);
    setTyped("");
    setSubmitting(false);
    setResult(null);
    setSubmitError(null);
  };

  const handleClose = () => {
    resetAll();
    onClose();
  };

  const handleFile = async (f: File | null) => {
    setFile(f);
    setParsedItemCount(null);
    setParseError(null);
    setResult(null);
    setSubmitError(null);
    // Picking a new file restarts the confirmation chain too.
    setStep(0);
    setTyped("");
    if (mode !== "add") setWords(pickWords(CONFIRM_COUNT[mode]));
    if (!f) return;
    try {
      const text = await f.text();
      const parsed = JSON.parse(text);
      const count = Array.isArray(parsed) ? parsed.length : 1;

      setParsedItemCount(count);
    } catch (e) {
      setParseError(`Could not parse JSON: ${(e as Error).message}`);
    }
  };

  const handleAdvanceConfirm = () => {
    if (mode === "add") return;
    const expected = words[step];

    if (typed.trim().toLowerCase() === expected?.toLowerCase()) {
      setStep((p) => p + 1);
      setTyped("");
    }
  };

  const canImport =
    !!file &&
    parseError === null &&
    parsedItemCount !== null &&
    parsedItemCount > 0 &&
    !submitting &&
    (mode === "add" || allConfirmed);

  const submit = async () => {
    if (!file || !canImport) return;
    setSubmitting(true);
    setSubmitError(null);
    try {
      const text = await file.text();
      const parsed = JSON.parse(text);
      const res = await fetch("/api/admin/catalog/import", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ mode, owner, items: parsed }),
      });

      if (!res.ok) {
        const body = await res.json().catch(() => ({}));

        throw new Error(body.error ?? `HTTP ${res.status}`);
      }
      const data: ImportResult = await res.json();

      setResult(data);
      onImported?.();
    } catch (e) {
      setSubmitError((e as Error).message);
    } finally {
      setSubmitting(false);
    }
  };

  return (
    <Modal isOpen={isOpen} size="2xl" onClose={handleClose}>
      <ModalContent>
        <ModalHeader>Import Catalog JSON</ModalHeader>
        <ModalBody className="flex flex-col gap-5">
          {/* File picker — click to browse or drag-and-drop a .json file */}
          <div className="flex flex-col gap-2">
            <span className="text-sm font-medium">JSON file</span>
            <label
              className={`flex flex-col items-center justify-center gap-1 rounded-xl border-2 border-dashed px-4 py-6 text-center cursor-pointer transition-colors ${
                dragOver
                  ? "border-primary bg-primary-50 dark:bg-primary-100/10"
                  : "border-default-300 hover:border-default-400"
              }`}
              htmlFor="catalog-json-file"
              onDragLeave={() => setDragOver(false)}
              onDragOver={(e) => {
                e.preventDefault();
                setDragOver(true);
              }}
              onDrop={(e) => {
                e.preventDefault();
                setDragOver(false);
                handleFile(e.dataTransfer.files?.[0] ?? null);
              }}
            >
              <svg
                className="w-6 h-6 text-default-400"
                fill="none"
                stroke="currentColor"
                strokeWidth={1.5}
                viewBox="0 0 24 24"
              >
                <path
                  d="M12 16.5V6m0 0-3.5 3.5M12 6l3.5 3.5M4.5 15.75v1.5A2.25 2.25 0 0 0 6.75 19.5h10.5a2.25 2.25 0 0 0 2.25-2.25v-1.5"
                  strokeLinecap="round"
                  strokeLinejoin="round"
                />
              </svg>
              <span className="text-sm font-medium">
                {file ? file.name : "Drop a JSON file here, or click to browse"}
              </span>
              <span className="text-[11px] text-default-400">
                .json — a single object or an array
              </span>
            </label>
            <input
              accept="application/json"
              className="hidden"
              id="catalog-json-file"
              type="file"
              onChange={(e) => handleFile(e.target.files?.[0] ?? null)}
            />
            {file && parseError === null && parsedItemCount !== null && (
              <p className="text-xs text-default-500">
                {file.name} — {parsedItemCount} item
                {parsedItemCount === 1 ? "" : "s"} parsed.
              </p>
            )}
            {parseError && <p className="text-xs text-danger">{parseError}</p>}
            <p className="text-[11px] text-default-400">
              Schema: see <code>scripts/sample-catalog.json</code>. Accepts a
              single object or an array. Paste a viewer/sm-viewer URL into{" "}
              <code>s3BaseUrl</code> and the importer auto-fills{" "}
              <code>datasetType</code> + <code>datasetId</code>.
            </p>
          </div>

          {/* Mode toggle */}
          <div className="flex flex-col gap-2">
            <p className="text-sm font-medium">Mode</p>
            <div className="flex flex-wrap gap-2">
              <Button
                color={mode === "add" ? "primary" : "default"}
                size="sm"
                variant={mode === "add" ? "solid" : "bordered"}
                onPress={() => setMode("add")}
              >
                Add (skip duplicates)
              </Button>
              <Button
                color={mode === "overwrite" ? "danger" : "default"}
                size="sm"
                variant={mode === "overwrite" ? "solid" : "bordered"}
                onPress={() => setMode("overwrite")}
              >
                Overwrite existing
              </Button>
              {isSuperAdmin && (
                <Button
                  color={mode === "replace_all" ? "danger" : "default"}
                  size="sm"
                  variant={mode === "replace_all" ? "solid" : "bordered"}
                  onPress={() => setMode("replace_all")}
                >
                  Replace All (SUPER_ADMIN)
                </Button>
              )}
            </div>

            {mode !== "add" && (
              <TypeWordConfirm
                allConfirmed={allConfirmed}
                mode={mode}
                step={step}
                totalSteps={requiredConfirms}
                typed={typed}
                wordList={words}
                onAdvance={handleAdvanceConfirm}
                onTypedChange={setTyped}
              />
            )}
          </div>

          {/* Owner choice (admins only) — owns the datasets this import registers */}
          <OwnerChoice value={owner} onChange={setOwner} />

          {/* Submit error */}
          {submitError && (
            <p className="text-sm text-danger bg-danger-50 dark:bg-danger-100/10 px-3 py-2 rounded-md">
              {submitError}
            </p>
          )}

          {/* Result summary */}
          {result && (
            <div className="flex flex-col gap-3 border-t border-default-200 pt-4">
              <p className="text-sm font-medium">
                Import complete ({result.mode})
              </p>
              <div className="flex flex-wrap gap-2">
                <Chip color="success" size="sm" variant="flat">
                  Created: {result.created.length}
                </Chip>
                <Chip color="primary" size="sm" variant="flat">
                  Updated: {result.updated.length}
                </Chip>
                <Chip color="default" size="sm" variant="flat">
                  Skipped: {result.skipped.length}
                </Chip>
                <Chip color="danger" size="sm" variant="flat">
                  Errors: {result.errors.length}
                </Chip>
                {result.deleted && result.deleted.length > 0 && (
                  <Chip color="danger" size="sm" variant="solid">
                    Deleted: {result.deleted.length}
                  </Chip>
                )}
              </div>
              {result.skipped.length > 0 && (
                <ResultList
                  items={result.skipped.map((s) => `${s.title} — ${s.reason}`)}
                  title="Skipped"
                />
              )}
              {result.errors.length > 0 && (
                <ResultList
                  isError
                  items={result.errors.map(
                    (e) => `${e.title ?? "(no title)"} — ${e.error}`,
                  )}
                  title="Errors"
                />
              )}
              {result.deleted && result.deleted.length > 0 && (
                <ResultList
                  isError
                  items={result.deleted.map((d) => d.title)}
                  title="Deleted (not in JSON)"
                />
              )}
            </div>
          )}
        </ModalBody>
        <ModalFooter>
          <Button variant="light" onPress={handleClose}>
            {result ? "Close" : "Cancel"}
          </Button>
          {!result && (
            <Button
              color={mode === "add" ? "primary" : "danger"}
              isDisabled={!canImport}
              onPress={submit}
            >
              {submitting ? (
                <span className="flex items-center gap-2">
                  <Spinner size="sm" /> Importing…
                </span>
              ) : mode === "replace_all" ? (
                "Replace All & Import"
              ) : mode === "overwrite" ? (
                "Overwrite & Import"
              ) : (
                "Import"
              )}
            </Button>
          )}
        </ModalFooter>
      </ModalContent>
    </Modal>
  );
}

const MODE_MESSAGES: Record<
  Exclude<Mode, "add">,
  { headline: string; subhead: string }
> = {
  overwrite: {
    headline: "⚠ Overwrite",
    subhead:
      "Existing rows matching the JSON will be REPLACED. Entries are deleted and recreated.",
  },
  replace_all: {
    headline: "☢ Replace All",
    subhead:
      "Same as Overwrite — AND every catalog row NOT in the JSON will be DELETED. This wipes anything you forgot to include.",
  },
};

function TypeWordConfirm({
  mode,
  step,
  totalSteps,
  wordList,
  typed,
  onTypedChange,
  onAdvance,
  allConfirmed,
}: {
  mode: Exclude<Mode, "add">;
  step: number;
  totalSteps: number;
  wordList: string[];
  typed: string;
  onTypedChange: (v: string) => void;
  onAdvance: () => void;
  allConfirmed: boolean;
}) {
  const { headline, subhead } = MODE_MESSAGES[mode];
  const expected = wordList[step] ?? "";

  const matches = useMemo(
    () => typed.trim().toLowerCase() === expected.toLowerCase(),
    [typed, expected],
  );

  return (
    <div className="flex flex-col gap-3 bg-danger-50 dark:bg-danger-100/10 border border-danger-200 dark:border-danger-200/30 px-3 py-3 rounded-md text-xs">
      <div className="flex items-start gap-2">
        <span className="text-danger font-semibold whitespace-nowrap">
          {headline}
        </span>
        <span className="text-default-700 flex-1">{subhead}</span>
      </div>

      {allConfirmed ? (
        <p className="text-danger font-medium">
          ✓ All {totalSteps} confirmations cleared. Import is enabled.
        </p>
      ) : (
        <>
          <p className="text-default-600">
            Step <span className="font-medium">{step + 1}</span> of{" "}
            <span className="font-medium">{totalSteps}</span>. Type the word
            below to continue. Each step uses a different word.
          </p>
          <div className="flex items-center gap-2">
            <code className="px-2 py-1 rounded bg-default-100 dark:bg-default-200/20 text-default-800 select-all font-mono text-sm">
              {expected}
            </code>
            <input
              className="flex-1 bg-default-100 dark:bg-default-200/20 rounded-md px-2 py-1.5 text-sm outline-none focus:ring-1 focus:ring-danger/50"
              placeholder="Type the word above…"
              type="text"
              value={typed}
              onChange={(e) => onTypedChange(e.target.value)}
              onKeyDown={(e) => {
                if (e.key === "Enter" && matches) {
                  e.preventDefault();
                  onAdvance();
                }
              }}
            />
            <Button
              color="danger"
              isDisabled={!matches}
              size="sm"
              variant="flat"
              onPress={onAdvance}
            >
              Confirm
            </Button>
          </div>
        </>
      )}
    </div>
  );
}

function ResultList({
  title,
  items,
  isError,
}: {
  title: string;
  items: string[];
  isError?: boolean;
}) {
  return (
    <div className="flex flex-col gap-1">
      <p
        className={`text-xs font-medium ${
          isError ? "text-danger" : "text-default-600"
        }`}
      >
        {title}
      </p>
      <ul className="text-xs text-default-500 max-h-32 overflow-y-auto pl-4 list-disc">
        {items.slice(0, 50).map((line) => (
          <li key={line}>{line}</li>
        ))}
        {items.length > 50 && (
          <li className="text-default-400">
            …{items.length - 50} more not shown
          </li>
        )}
      </ul>
    </div>
  );
}
