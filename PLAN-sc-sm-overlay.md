# SC + SM Same-Space Overlay Plan

## Status: WIP — Proof of concept committed

The basic rendering works: when `mapping.json` has `linkColumn: "__all__"`, SM molecules load and render in the same Three.js scene as single cell points.

## What's Done

- Detect `__all__` mapping in `three-scene.tsx`
- Load SM dataset from `links["__all__"]` URL in background
- Auto-select 3 default genes
- Render SM point clouds in same scene (`renderOrder: -1`, below cell points)
- Loading indicator while SM dataset loads
- Cleanup on gene deselection

## What's Left

### 1. SM Controls in Left Sidebar
- Add a collapsible "Molecules" section below existing Celltype/Gene/Size controls
- Gene search + checkbox list (same pattern as `single-molecule-controls.tsx`)
- Per-gene color chips showing current color
- Master toggle to show/hide all SM molecules at once

### 2. Per-Gene Controls
- Color picker (per gene)
- Opacity slider (per gene or global)
- Size slider (global SM size multiplier)
- Visibility toggle (per gene)

### 3. Render Order
- SC points should render on top — currently using `renderOrder` but may need
  to add a small Z offset to SC points to ensure they're visually on top in 2D view

### 4. Store Integration
- Either reuse `singleMoleculeVisualizationStore` for the overlay, or create
  a lightweight overlay-specific store
- Need: selectedGenes, gene colors, gene visibility, global scale, master toggle
- Consider: should overlay genes persist in URL state?

### 5. URL State
- Encode SM overlay gene selections in the cell viz URL state
- New field in `CellVizUrlState` for overlay genes

### 6. Interaction
- Tooltips: when hovering an SM molecule point, show gene name + coordinates
- Decide: should raycaster work on SM points too, or only cell points?

### 7. Performance
- Lazy-load SM genes on-demand (already supported by `fromCustomS3`)
- Cache loaded genes (already built into `SingleMoleculeDataset`)
- Consider: limit max simultaneous SM genes to avoid memory pressure

## Architecture Notes

- SM dataset is loaded via `SingleMoleculeDataset.fromCustomS3(url)`
- Gene coordinates come from `getCoordinatesByGene(gene)` → `Float32Array`
- Coordinates are raw microns — same space as SC raw coordinates
- Each gene gets its own `THREE.Points` mesh with `PointsMaterial`
- SC uses custom `ShaderMaterial`, SM uses `PointsMaterial` — both coexist in same scene

## Dependencies

- This branch should be merged AFTER the aesthetics/visualization revamp
  to avoid merge conflicts in `three-scene.tsx` and control components
