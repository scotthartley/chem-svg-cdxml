---
name: chem-svg-cdxml
description: Draw chemical structures, Lewis structures, reaction schemes, and arrow-pushing mechanisms in the ChemDraw ACS Document 1996 style as editable SVG or CDXML, with a PNG preview Claude inspects before delivering. Use when the user asks to draw a molecule, a Lewis structure, a reaction scheme, a mechanism with curved arrows, resonance structures, polymers (repeat-unit brackets, copolymers), or a catalytic cycle.
---

# chem-svg-cdxml

`<repo>` is this skill's directory: the target of the project's `.claude/skills/chem-svg-cdxml` symlink.
Run every command as `uv run --project <repo> python <repo>/scripts/<script> ...`.
If the Read tool refuses `<repo>` files (outside the working directory, symlink or not), read them with `cat` / `sed`.

## Workflow

1. **Write a spec** (YAML). The core is below; for each feature the request needs, read its
   schema section first (table under "Features"). Full schema: `references/scheme-schema.md`.
   ```yaml
   preset: acs1996          # or acs1996_x0.8, acs1996_x0.6, or any in presets/
   species:
     ester:   {smiles: "COC(=O)c1ccc(O)cc1", abbreviate: [OMe], orientation: {rotate: 120}}
     bnbr:    {smiles: "BrCc1ccccc1"}
     ether:   {smiles: "COC(=O)c1ccc(OCc2ccccc2)cc1", abbreviate: [OMe, OBn]}
     alcohol: {smiles: "OCc1ccc(OCc2ccccc2)cc1", abbreviate: [OBn]}
   steps:
     - {from: [ester, bnbr], to: ether, above: K2CO3, below: "acetone, reflux"}
     - {to: alcohol, above: LiAlH4, below: "THF, 0 °C"}
   numbers: {ester: 1, bnbr: 2, ether: 3, alcohol: 4}
   ```
   - Use only the conditions the user gave. A reagent named in the request (BnBr) goes over
     the arrow unless the user wants it drawn and numbered.
   - Conditions text: digits after a letter or `)` become subscripts (K2CO3, Pd(PPh3)4);
     a trailing `+`/`−` is a superscript (H3O+; `−` is U+2212, `Br-` is an error); `*n*-BuLi` is italic; `\` escapes.
     Several lines: `above: [line1, line2]`; numbered stages: `stages: ["*n*-BuLi, THF, −78 °C", "H3O+"]`.
   - Products are aligned to the species they come from. Orient only the first species of each
     sequence (and separate reactants) with `orientation: {rotate: …, flip_h: …}` (degrees counter-clockwise
     on screen); `align: false` opts out.
   - Without `steps`, the species are drawn side by side (single structures).
2. **Validate:** `scripts/validate.py spec.yaml`. Fix the spec until it prints OK (schema,
   RDKit chemistry, abbreviations, electron bookkeeping, overlaps, curved arrows crossing
   bonds or labels, dots on bonds, conditions wider than their arrow).
3. **Render:** `scripts/render.py spec.yaml --out DIR [--crops] [--maps]` writes `spec.svg`,
   `spec.png` (4 px/pt; `--px-per-pt 10` for close inspection) and `spec.cdxml`; `--crops` adds
   close-ups of arrows, charged labels and bracket corners.
4. **Inspect** the PNG (and crops): orientation of the first species, label sides, double-bond
   sides, wedges, charges, conditions centered on their arrows, numbers on one line; for
   mechanisms, that each arrow leaves its lone pair or bond and points at its target.
5. **Adjust** the spec only (orientation, abbreviations, line breaks, `align`, `overrides`),
   never the SVG. At most 3 tries per curved arrow and 5 render rounds overall; then deliver
   and report what is left.
6. **Deliver** `.svg` (editable text, named groups; `--affinity` writes for Affinity: separate
   arrowhead shapes, bold/italic face names), `.yaml`, `.png`, and `.cdxml` for ChemDraw (every spec
   kind; `--no-cdxml` skips it). Report stderr warnings (fallback font; `warning: no CDXML:` means
   that spec has no `.cdxml`).

## Features

Read the named section of `references/scheme-schema.md` before using a feature. Which features
combine (and which raise an error): `references/feature-matrix.md`.

| Need | Key (minimal form) | Schema section |
|---|---|---|
| Abbreviated groups (OTf, Boc, *t*-Bu, TBS, …) | `abbreviate: [OTf, t-Bu]`; own groups under `abbreviations:` | Schemes (`abbreviate`) |
| Generic atoms, ligand text (X, R, L, Lₙ, dppf) | stand-in atom `[Cl:9]` / `[*:1]` + top-level `labels: {9: X}` | Schemes (`labels`) |
| Reversible step | `arrow_type: equilibrium` (or `equilibrium_products` / `_reactants`) | Schemes |
| CH₃ labels, Kekulé form, molfile input | `show_carbons: terminal`, `kekule:`, `molfile:` | Single species |
| Charges: circled, custom text | `charges: {circled: true}`; `charge_labels: {1: "+3"}` | Single species |
| Lewis structure (every atom labeled) | `lewis: true` \| `{bonds: lines}` \| `{electrons: formal_charge}` | Lewis structures |
| Chair conformer | `chair: true` \| `{flip: true}` \| `{axial: [maps]}` | Chairs |
| Polymer, copolymer | `polymer: {units: ["*CC(*)c1ccccc1"], start, end, kind, ratio}` | Polymers |
| Monomer coefficient (n styrene →) | `coefficients: {styrene: n}` on the step | Polymers |
| Catalytic cycle | `cycle:` in place of `steps:` | Catalytic cycles |
| Curved arrows | `arrows:` on a step | Mechanisms (and below) |
| Resonance, more rows | `resonance:`, `new_row: true` | Resonance, sequences, rows |

## Mechanisms (curved arrows)

Worked example and common mistakes: `references/workflow.md`.
1. Write the steps in words (who attacks what, what leaves).
2. Assign map numbers once, in a table, to every atom an arrow touches (`[O-:1]`); keep them
   for the whole scheme. Map each H that moves as an atom (`[H:4][CH2:3]…`). Writing one of a
   carbon's H's as an atom draws all of them as bonds, so write only the H that moves.
3. Write the reactants and each step's arrows (drawn on the step's reactants).
4. Write rough products, run `validate.py --suggest`, and paste the `suggested:` lines. Every
   step's arrows must turn its reactants into its products exactly (charges, radicals, bond
   orders, H moved); the validator's message names the arrow and the atom.
5. **Account for every product:** an acid/base that only moves a proton goes `over:` the arrow
   (its conjugate and the one regenerated last are `omit:`ted, without text); omit other products
   the next step doesn't use, noting a loss that matters on the arrow (`above: "−H2O"`; never a
   loss below ⇌). Why: `references/workflow.md` ("A product that disappears").
6. Render with `--maps` (map numbers beside the atoms, preview only): it shows which end of
   each species faces its partner.
7. Orient, render, inspect. `render.py` prints each arrow's default `bow`; if an arrow crosses
   something, try `overrides` (the other `bow`, then `tail_side`, `curvature`, `shape: arc|tangent`), or turn a species
   (`orientation`) or swap the reactants (`from: [mebr, hydroxide]`) so the atoms face each other.

```yaml
species:
  hydroxide: {smiles: "[H][O-:1]"}
  mebr:      {smiles: "[CH3:2][Br:3]", show_carbons: terminal}
  methanol:  {smiles: "[OH:1][CH3:2]", show_carbons: terminal}
  bromide:   {smiles: "[Br-:3]"}
steps:
  - from: [hydroxide, mebr]
    to: [methanol, bromide]
    arrows:
      - {from: "lp(1)", to: "between(1,2)"}            # type: two_electron (default) | fishhook
      - {from: "bond(2-3)", to: "atom(3)"}
overrides:
  arrows: [{step: 1, index: 1, bow: right}]           # index counts the step's arrows from 0
```
- Endpoints: from `lp(n)`, `bond(a-b)` / `bond(a=b)` / `bond(a#b)`, `electron(n)` (a radical);
  to `atom(n)` or `between(a,b)` (a new or stronger bond). Full vocabulary and drawing rules:
  `references/scheme-schema.md` (Mechanisms).
- Lone pairs the arrows don't involve: `show_lone_pairs: [n]` per species.
- `[C]`/`[CH2]`/`[N]` are singlet carbenes/nitrenes (a drawn pair, `lp(n)`); `triplet: [id]` for two unpaired dots.

## Resonance, sequences, rows

```yaml
layout: {width: 504}                                  # optional; pt
steps:
  - {from: [benzene, nitronium], to: sigma_a, arrows: [...]}
  - resonance: [sigma_a, sigma_b, sigma_c]            # first = the previous step's first product
    arrows: [[{from: "bond(13=14)", to: "between(12,13)"}], [...]]   # list k: contributor k -> k+1
  - {from: [sigma_c, bisulfate], to: [nitrobenzene, sulfuric], omit: [sulfuric], arrows: [...]}
  - {new_row: true, from: [a, b], to: c}             # another sequence, on a new row
```
- Contributors keep their atom maps and connectivity; they align to the first (orient only it).
- Reuse a species by name when it appears again: each later occurrence is its own copy.
- `omit:` products stay in the bookkeeping, out of the drawing; a later step that names one
  adds it. A mechanism step's added reactants (a base) and its `over:` sit above its arrow.
- Rows break only between steps; if a step plus a bracket can't fit, raise `layout.width`.

## Known limits (raise errors rather than guess)

- Explicit CH₂/CH labels on internal carbons; condensed labels that aren't abbreviations (CH₃Cl as text); aromatic circles.
- Lewis structures: rings, more than 4 neighbors (PCl₅), Lewis species in schemes.
- Chairs: fused or spiro rings (decalins, steroids), sp² ring atoms, boats.

## Rules

- Never edit `src/`, `scripts/`, `presets/`, or `references/` from a use session.
  A project's own style, shared groups, molfiles and scripts go in the project: `references/extending.md`.
- Keep the user's preset (default acs1996), even when rows wrap; raise `layout.width` only on request.
- Problems and workarounds go to `<repo>/issues/inbox.md` (create it if missing), appended, with:
  - the date,
  - a spec excerpt,
  - the symptom,
  - the workaround.
- Validation proves internal consistency, not mechanistic correctness. Say so for
  research-level mechanisms.
