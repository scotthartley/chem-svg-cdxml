# Extending from a project folder

A project that has the skill linked in (`.claude/skills/chem-svg-cdxml`) can add its own style,
groups, structures and scripts without touching the skill repo. Everything below lives in the
project; the repo's `src/`, `scripts/`, `presets/` and `references/` stay as they are.

## Layout

```
my-project/
  styles/
    mystyle.yaml               # a preset (copied from presets/ or written by new_preset.py)
    cdxml/
      mystyle.xml              # its CDXML header (root attributes, fonts, colors, page)
      mystyle.nicknames.xml    # its nickname templates (may be empty: generated fragments)
  groups.yaml                  # shared abbreviations
  structures/
    taxol_core.mol             # structures drawn in ChemDraw, saved as molfiles
  figures/
    fig1.yaml                  # specs
  scripts/
    make_figures.py            # optional: batch rendering
```

Paths in a spec are relative to the spec's own folder (`../styles/...` from `figures/`).

## Your own style: `preset:` as a path

```yaml
preset: ../styles/mystyle.yaml
```

A value ending in `.yaml` (or `.yml`) is a file; anything else is a name in the repo's `presets/`.
The CDXML export reads the header and nicknames from `cdxml/<stem>.xml` and
`cdxml/<stem>.nicknames.xml` next to the preset file. If a file is missing, the error gives the path
it looked for. The SVG and PNG don't need either file.

Two ways to make one:

- **From ChemDraw** (recommended): the style probe (style-spec.md §16), written into the project:

  ```bash
  uv run --project ~/Coding/productivity/chem-svg-cdxml python ~/Coding/productivity/chem-svg-cdxml/scripts/new_preset.py probe --out probe
  # open probe/style_probe.cdxml in ChemDraw, apply your style, redo signs and dots, save as CDXML + SVG
  uv run --project ~/Coding/productivity/chem-svg-cdxml python ~/Coding/productivity/chem-svg-cdxml/scripts/new_preset.py \
      build --name mystyle --probe probe/style_probe.cdxml --out-dir styles
  ```

  This writes `styles/mystyle.yaml`, `styles/cdxml/mystyle.xml` and an empty
  `styles/cdxml/mystyle.nicknames.xml`.
- **By hand**: copy `presets/acs1996.yaml` and `presets/cdxml/acs1996.xml` / `.nicknames.xml`
  (renamed to the new stem) and edit them. Keep the header's `BondLength`, `LineWidth` and label
  settings in step with the preset, or ChemDraw will show other settings than the drawing uses.
  Changing one length does not rescale the others: for a smaller or larger style, use the probe.

## Shared groups: `abbreviations:` from files

```yaml
abbreviations: ../groups.yaml                     # one file
abbreviations: [../groups.yaml, {Xy: {kind: text}}]  # files and inline entries, merged in order
```

Each file has the same shape as the inline block (scheme-schema.md, `abbreviations`):

```yaml
OCy:  {smiles: "*OC1CCCCC1"}
TBDMS: {smiles: "*[Si](C)(C)C(C)(C)C", kind: group, text: "TBDMS"}
dtbpy: {kind: text}
```

Later entries replace earlier ones with the same key, and every entry adds to or replaces the
built-in registry (`src/chem_svg_cdxml/data/abbreviations.yaml`).

## Reusable structures: `molfile:`

```yaml
species:
  core: {molfile: ../structures/taxol_core.mol}
```

Molfiles resolve against the spec's folder first, then the repo. The layout is kept and scaled to
the preset's bond length; wedges and atom aliases naming abbreviations are kept
(scheme-schema.md, `molfile`).

## Scripts

A project script runs in the skill's environment:

```bash
uv run --project ~/Coding/productivity/chem-svg-cdxml python scripts/make_figures.py
```

```python
from pathlib import Path

from chem_svg_cdxml.figure import build_figure
from chem_svg_cdxml.preview import rasterize
from chem_svg_cdxml.svg_writer import to_svg

for spec in sorted(Path("figures").glob("*.yaml")):
    fig = build_figure(spec)                 # validates and lays out the spec
    svg, offset = to_svg(fig.root)
    spec.with_suffix(".svg").write_text(svg)
    spec.with_suffix(".png").write_bytes(rasterize(svg, px_per_pt=4))
    spec.with_suffix(".cdxml").write_text(fig.cdxml(offset))  # CdxmlUnsupported if it can't
```

The public calls are `figure.build_figure`, `Figure.root` / `Figure.cdxml`, `svg_writer.to_svg`,
`preview.rasterize` and `spec.load_spec`. Everything else in `chem_svg_cdxml` is internal and may
change. For the electron and geometry checks, run `scripts/validate.py` on the spec first.
