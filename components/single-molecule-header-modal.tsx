"use client";

import type { MoleculeColumnMapping } from "@/lib/config/moleculeColumnMappings";

import { useEffect, useMemo, useState } from "react";
import {
  Modal,
  ModalContent,
  ModalHeader,
  ModalBody,
  ModalFooter,
} from "@heroui/modal";
import { Button } from "@heroui/button";

const NONE = "__none__";

type Role = "gene" | "x" | "y" | "z" | "cellId";

const FIELDS: {
  role: Role;
  label: string;
  required: boolean;
  help?: string;
}[] = [
  { role: "gene", label: "Gene / feature", required: true },
  { role: "x", label: "X coordinate", required: true },
  { role: "y", label: "Y coordinate", required: true },
  { role: "z", label: "Z coordinate (optional)", required: false },
  {
    role: "cellId",
    label: "Cell ID (optional)",
    required: false,
    help: "Used to split assigned vs. unassigned molecules (e.g. Xenium's UNASSIGNED).",
  },
];

// Static (purge-safe) highlight classes per mapped role.
const HL: Record<Role, { th: string; td: string }> = {
  gene: { th: "text-primary bg-primary/10", td: "bg-primary/5 font-medium" },
  x: { th: "text-success bg-success/10", td: "bg-success/5 font-medium" },
  y: { th: "text-success bg-success/10", td: "bg-success/5 font-medium" },
  z: { th: "text-success bg-success/10", td: "bg-success/5 font-medium" },
  cellId: { th: "text-warning bg-warning/10", td: "bg-warning/5 font-medium" },
};

interface Props {
  isOpen: boolean;
  fileName: string;
  columns: string[];
  rows: Record<string, unknown>[];
  /** Auto-detected mapping used to seed the dropdowns. */
  autoMapping: MoleculeColumnMapping;
  onConfirm: (mapping: MoleculeColumnMapping) => void;
  onCancel: () => void;
}

function cell(value: unknown): string {
  if (value === null || value === undefined || value === "") return "—";

  return String(value);
}

export function SingleMoleculeHeaderModal({
  isOpen,
  fileName,
  columns,
  rows,
  autoMapping,
  onConfirm,
  onCancel,
}: Props) {
  const seed = (col: string | null) =>
    col && columns.includes(col) ? col : NONE;

  const [sel, setSel] = useState<Record<Role, string>>({
    gene: NONE,
    x: NONE,
    y: NONE,
    z: NONE,
    cellId: NONE,
  });

  // Reseed from the auto-detected mapping whenever the modal (re)opens.
  useEffect(() => {
    if (!isOpen) return;
    setSel({
      gene: seed(autoMapping.gene),
      x: seed(autoMapping.x),
      y: seed(autoMapping.y),
      z: seed(autoMapping.z),
      cellId: seed(autoMapping.cellId),
    });
  }, [isOpen]);

  const roleByColumn = useMemo(() => {
    const m = new Map<string, Role>();

    (Object.keys(sel) as Role[]).forEach((role) => {
      if (sel[role] !== NONE) m.set(sel[role], role);
    });

    return m;
  }, [sel]);

  const valid = sel.gene !== NONE && sel.x !== NONE && sel.y !== NONE;

  const confirm = () => {
    if (!valid) return;
    onConfirm({
      gene: sel.gene,
      x: sel.x,
      y: sel.y,
      z: sel.z === NONE ? "" : sel.z,
      cellId: sel.cellId === NONE ? null : sel.cellId,
    });
  };

  return (
    <Modal
      isOpen={isOpen}
      scrollBehavior="inside"
      size="4xl"
      onClose={onCancel}
    >
      <ModalContent>
        <ModalHeader className="flex flex-col gap-1">
          <span>Confirm columns</span>
          <span className="text-xs font-normal text-default-400">
            {fileName} — {columns.length} columns · showing {rows.length} sample
            rows
          </span>
        </ModalHeader>
        <ModalBody className="gap-5">
          {/* Mapping dropdowns */}
          <div className="grid grid-cols-2 gap-3 sm:grid-cols-3">
            {FIELDS.map((f) => (
              <label key={f.role} className="flex flex-col gap-1">
                <span className="text-xs font-medium text-default-600">
                  {f.label}
                  {f.required && <span className="text-danger"> *</span>}
                </span>
                <select
                  className="rounded-lg border border-default-300 bg-default-100 px-2 py-1.5 text-sm outline-none focus:border-primary"
                  value={sel[f.role]}
                  onChange={(e) =>
                    setSel((s) => ({ ...s, [f.role]: e.target.value }))
                  }
                >
                  {!f.required && <option value={NONE}>(none)</option>}
                  {f.required && sel[f.role] === NONE && (
                    <option value={NONE}>Select a column…</option>
                  )}
                  {columns.map((c) => (
                    <option key={c} value={c}>
                      {c}
                    </option>
                  ))}
                </select>
                {f.help && (
                  <span className="text-[11px] text-default-400">{f.help}</span>
                )}
              </label>
            ))}
          </div>

          {/* Sample rows — mapped columns highlighted so you can see values. */}
          <div className="overflow-x-auto rounded-lg border border-default-200">
            <table className="min-w-full text-xs">
              <thead>
                <tr className="border-b border-default-200 bg-default-100/60">
                  {columns.map((c) => {
                    const role = roleByColumn.get(c);

                    return (
                      <th
                        key={c}
                        className={`whitespace-nowrap px-2 py-1.5 text-left font-medium ${
                          role ? HL[role].th : "text-default-500"
                        }`}
                      >
                        {c}
                        {role && (
                          <span className="ml-1 text-[9px] uppercase opacity-70">
                            {role}
                          </span>
                        )}
                      </th>
                    );
                  })}
                </tr>
              </thead>
              <tbody>
                {rows.map((row, i) => (
                  <tr key={i} className="border-b border-default-100">
                    {columns.map((c) => {
                      const role = roleByColumn.get(c);

                      return (
                        <td
                          key={c}
                          className={`whitespace-nowrap px-2 py-1 ${
                            role ? HL[role].td : ""
                          }`}
                        >
                          {cell(row[c])}
                        </td>
                      );
                    })}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </ModalBody>
        <ModalFooter>
          <Button variant="light" onPress={onCancel}>
            Cancel
          </Button>
          <Button color="primary" isDisabled={!valid} onPress={confirm}>
            Looks good — process
          </Button>
        </ModalFooter>
      </ModalContent>
    </Modal>
  );
}
