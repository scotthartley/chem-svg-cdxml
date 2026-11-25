# Feature combinations

Which species features work in which contexts. A feature that isn't allowed raises a schema
error that names it: it never draws something wrong without telling you. Keys and details:
`references/scheme-schema.md`.

**✓** = covered by an example or gallery case (named). **✓ p** = tried on 2026-11-25 (validates
and renders correctly), but no example or gallery case covers it yet. **✗** = refused with an
error. **?** = not tried: render it and inspect closely, and log anything odd to the inbox.

| Feature | Single species | Scheme step | Mechanism (curved arrows) | Resonance | Cycle |
|---|---|---|---|---|---|
| `abbreviate` | ✓ m1_abbrev, m1_abbrev2 | ✓ m2_sequence | ✓ p (MeOTf + phenoxide) | ✓ m4_resonance (NO₂) | ✓ m6_suzuki_cycle (PMe₃) |
| `labels` (X, R, L, dppf) | ✓ m1_abbrev2 | ? | ✓ tests/specs/phenol_x_abstraction | ? | ✓ sonogashira_cycle |
| `charges.circled` | ✓ m1_charges | ? | ? | ✓ p (acetate) | ✓ p |
| `charge_labels` | ✓ m1_charges | ? | ? (display only: the bookkeeping uses the SMILES charge) | ? | ✓ p |
| `lewis` | ✓ m1_lewis, m1_lewis_lines | ✗ | ✗ | ✗ | ✗ |
| `chair` | ✓ m1_chair | ✓ p (ketone reduction) | ✓ m3_e2_chair | ? | ? |
| `polymer` | ✓ m7_polymers, m7_raft | ✓ m7_scheme; ✓ p with `new_row` | ? (`chain_ends` on a SMILES species instead: ✓ m7_propagation) | ? | ? |
| `kekule` | ? | ? | ✓ tests/specs/phenol_x_abstraction | ? | ? |
| `molfile` | ✓ | ✓ m2_scheme_basic, m4_wrap | ? | ? | ? |
| `triplet` / carbenes | ✓ m1_carbenes | ? | ? | ? | ? |
| `arrow_type: equilibrium` | n/a | ✓ m2_equilibrium | ✓ m3_proton_transfer, m4_perkin | ✗ (↔ only) | ✓ m6_enamine_cycle |

## Exclusions on one species

- `lewis`: SMILES only; not with `kekule`, `abbreviate`, `show_carbons` or `show_lone_pairs`
  (the last two are implied); `orientation.rotate` a multiple of 90.
- `chair`: SMILES only; not with `lewis`, `molfile` or `orientation`.
- `polymer`: not with `kekule`, `orientation`, `lewis`, `chair` or `align_to`; it aligns to nothing.
- `cycle`: at least 3 steps.

## Known rough edges

- `abbreviate: [OTf]` needs a bonded O–Tf. On the free triflate anion it reports "no such
  terminal group". `abbreviate: [Tf]` works, but draws the O with its charge and lone pairs
  (⁻:Ö–Tf), not as TfO⁻. (2026-11-25)

## Adding a feature

Fill in its row: try it in each context, either cover the combination with an example or
gallery case (✓), refuse it with a schema error (✗), or leave it as ? and say so in the commit.
