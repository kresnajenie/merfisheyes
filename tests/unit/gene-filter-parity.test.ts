import { describe, it, expect } from "vitest";
import { readFileSync } from "fs";
import path from "path";

import { shouldFilterGene } from "@/lib/utils/gene-filters";

/**
 * The control-probe filter must be the same rule in both pipelines. The
 * Python scripts keep a literal list of regex sources; this test reads them
 * out of the source files and checks the JS function agrees on a panel of
 * gene names — so an edit to one side without the other fails CI.
 */
const PANEL = [
  "Gad1", "Slc17a7", "MT-CO1", "Blank-12", "BLANK_0001", "NegControlProbe_3",
  "NegControlCodeword_0502", "Negative Control Probe 1", "neg-ctrl-7",
  "Unassigned Codeword 12", "DeprecatedCodeword_3", "blankish_gene", "Codeword_9",
  "Vip", "Pvalb", "UNASSIGNED", "",
];

function pythonPatterns(file: string): RegExp[] {
  const src = readFileSync(path.join(__dirname, "..", "..", "scripts", file), "utf-8");
  const out: RegExp[] = [];
  // r"pattern" literals inside the filter list(s)
  const re = /r"((?:[^"\\]|\\.)*)"/g;
  let m: RegExpExecArray | null;

  while ((m = re.exec(src))) {
    const p = m[1];

    if (/negative|neg|unassigned|deprecated|codeword|blank|negcontrol/i.test(p)) {
      out.push(new RegExp(p.replace(/\\s/g, "\\s"), "i"));
    }
  }

  return out;
}

function pythonFilters(patterns: RegExp[], gene: string): boolean {
  return patterns.some((p) => p.test(gene));
}

describe("control-probe gene filter parity (JS vs Python)", () => {
  for (const script of ["process_spatial_data.py", "process_single_molecule.py"]) {
    it(`agrees with ${script} on every panel gene`, () => {
      const patterns = pythonPatterns(script);

      expect(patterns.length).toBeGreaterThanOrEqual(7);
      for (const gene of PANEL) {
        expect(
          { gene, python: pythonFilters(patterns, gene) },
        ).toEqual({ gene, python: shouldFilterGene(gene) });
      }
    });
  }
});
