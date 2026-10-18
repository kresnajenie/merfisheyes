# How MERFISHEYES is built

These are the rules a change is judged by. They exist so that the next
format, the next viewer and the next contributor fit without a rewrite.

## 1. One implementation per concern

Every concern has exactly one home, and every surface uses it. Molecules are
drawn by one layer whether they sit on their own, over cells, in a split
panel or from a folder on disk. A dataset is opened by one opener whether it
comes by id, by URL or from disk. If a second copy of a concern appears, the
copy is the bug.

## 2. An element is a kind, and a kind is a contract

Data arrives as elements: a table of cells, a points element of molecules, a
shapes element of polygons, a multiscale image, a label image. Each kind is
one entry in a registry that names its opener, its layer, its controls, its
URL state and its tests. Adding a kind is adding an entry; it never edits a
viewer. Every layer implements the same lifecycle: attach, update, tick,
report, abort, dispose, within a stated byte and memory budget
(`docs/specs/LAYER-CONTRACT.md`).

## 3. The format is the product

What we write to disk is a specification first and code second: versioned,
with a JSON Schema for its metadata, a validator anyone can run, and a
reference reader and writer that agree byte for byte. A format is finished
when a library we do not maintain (the `spatialdata` Python package, an
OME-Zarr reader) opens it, and when someone outside the team has reviewed
the spec. Formats build on published conventions (OME-NGFF, SpatialData,
zarr-conventions, GeoParquet) and add only what those do not say.

## 4. Read the view, not the dataset

Every read is bounded by the screen, never by the size of the data: a level
whose molecules in view fit the budget, the tiles under the camera, the one
gene asked for. Levels are uniform samples of each other so density is
honest at every zoom. Zooming or panning back fetches nothing. What is on
screen stays until what replaces it has arrived. Nothing is drawn twice.

## 5. Human time is the budget

A frame takes 16 ms while the camera moves. An interaction shows its result
or a visible start within 100 ms. A new view shows content within 300 ms and
settles within 1 s. Anything past 10 s has a percent-done indicator and can
be cancelled; every download can be cancelled. These are acceptance
criteria measured on real data in the perf gate (`docs/specs/BUDGETS.md`),
not aspirations. Sources: Nielsen's 0.1 / 1 / 10 s limits (Miller 1968,
Card et al. 1991), the RAIL frame budget, Liu & Heer 2014 on 500 ms.

## 6. Correctness is shown from outside

A change is right when something other than its author says so: the
validator and the reader agree, the upstream library round-trips the store,
the conformance suite passes for every layer, the perf budgets hold on the
largest real store, and a fresh implementer can add a kind from the written
contract alone. Tests run in CI on every push; a browser test exercises every
surface a feature touches.

## 7. Start from the raw material

Design from the code, the on-disk metadata of real stores, the upstream
specifications and the network traces, not from summaries. Write the
measurement before the optimisation. Record decisions, with their reasons,
where the next person will look (the spec, the ADR, CLAUDE.md).

## 8. The largest dataset we have is the design target

Today that is the 11.24 billion molecule, ten-section mouse brain store.
A design that explores it seamlessly is the bar for the next five years; a
design that only works on small data is not finished.
