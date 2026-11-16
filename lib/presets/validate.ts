/**
 * The preset contract, enforced: a config is accepted only if it matches
 * docs/specs/schemas/preset.schema.json. Server-only (Ajv, node:crypto).
 */
import type { TechnologyPreset } from "@/lib/config/technology-presets";

import { createHash } from "node:crypto";

import Ajv2020 from "ajv/dist/2020";

import schema from "@/docs/specs/schemas/preset.schema.json";

const validate = new Ajv2020({ allErrors: true }).compile(schema);

/** Null when the config matches the schema, else what is wrong with it. */
export function presetConfigError(config: unknown): string | null {
  if (validate(config)) return null;

  return (validate.errors ?? [])
    .map((e) => `${e.instancePath || "config"} ${e.message}`)
    .join("; ");
}

/**
 * The name of a version: the same technology and config always give the same
 * hash, whatever order the fields were written in.
 */
export function presetHash(
  technology: string,
  config: TechnologyPreset,
): string {
  const sorted = Object.fromEntries(
    Object.entries(config).sort(([a], [b]) => (a < b ? -1 : 1)),
  );

  return createHash("sha256")
    .update(JSON.stringify([technology, sorted]))
    .digest("hex")
    .slice(0, 16);
}
