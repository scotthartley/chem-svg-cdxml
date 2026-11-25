"""scripts/render.py writes CDXML by default, and a spec it can't export still gets SVG + PNG."""

from pathlib import Path

import render
from render import main

from chem_svg_cdxml.cdxml_writer import CdxmlUnsupported
from chem_svg_cdxml.figure import Figure

ROOT = Path(__file__).resolve().parents[1]
SCHEME = ROOT / "examples" / "m3" / "m3_carbonyl_addition_x1.0.yaml"
SPECIES = ROOT / "examples" / "m1" / "m1_charges_x1.0.yaml"


def test_cdxml_is_written_by_default(tmp_path):
    for spec in (SCHEME, SPECIES):
        paths = render.render(spec, tmp_path)
        assert [p.suffix for p in paths] == [".svg", ".png", ".cdxml"]
        assert paths[2].read_text().lstrip().startswith("<?xml")


def test_no_cdxml_writes_none_and_the_same_svg(tmp_path):
    for spec in (SCHEME, SPECIES):
        on = render.render(spec, tmp_path / "on")
        off = render.render(spec, tmp_path / "off", with_cdxml=False)
        assert [p.suffix for p in off] == [".svg", ".png"]
        assert not list((tmp_path / "off").glob("*.cdxml"))
        assert on[0].read_text() == off[0].read_text()


def test_flags(tmp_path, capsys):
    assert main([str(SPECIES), "--out", str(tmp_path / "a"), "--no-cdxml"]) == 0
    assert not list((tmp_path / "a").glob("*.cdxml"))
    assert main([str(SPECIES), "--out", str(tmp_path / "b"), "--cdxml"]) == 0  # old commands still work
    assert list((tmp_path / "b").glob("*.cdxml"))


def test_unsupported_cdxml_warns_and_still_renders(tmp_path, monkeypatch, capsys):
    def refuse(*args, **kwargs):
        raise CdxmlUnsupported("test: not exportable")
    monkeypatch.setattr(Figure, "cdxml", refuse)
    assert main([str(SCHEME), "--out", str(tmp_path)]) == 0
    assert sorted(p.suffix for p in tmp_path.iterdir()) == [".png", ".svg"]
    assert capsys.readouterr().err.count("warning: no CDXML: test: not exportable") == 1


def test_crops_cover_single_species(tmp_path):
    """--crops on a spec without steps: charged / dotted Lewis atoms and polymer bracket corners
    (the r2 trials got no crops for these)."""
    lewis = render.render(ROOT / "examples" / "benchmarks" / "ref_lewis.yaml", tmp_path, with_crops=True, with_cdxml=False)
    assert any(p.name == "ref_lewis_crop_nitrate-atom2.png" for p in lewis)
    poly = render.render(ROOT / "examples" / "m7" / "m7_polymers_x1.0.yaml", tmp_path, with_crops=True, with_cdxml=False)
    names = {p.name for p in poly}
    assert {f"m7_polymers_x1.0_crop_pmma-bracket-1-{s}-{c}.png" for s in ("left", "right") for c in ("top", "bottom")} <= names


def test_render_prints_each_chairs_axial_and_equatorial_groups(tmp_path, capsys):
    """r2 chair trial: a cis/trans slip in the SMILES validated and drew without a word; render.py
    now says which substituents went axial."""
    render.render(ROOT / "examples" / "m1" / "m1_chair_x1.0.yaml", tmp_path, with_cdxml=False)
    err = capsys.readouterr().err
    assert "chair cis_4_tert_butylcyclohexanol: axial: O (atom 9); equatorial: C (atom 2)" in err
    assert "chair menthol: axial: none;" in err


def test_a_plain_coefficient_is_italic_like_the_subscript(tmp_path):
    """User choice 2026-11-25: `coefficients: {styrene: n}` sets n italic, as the bracket's n."""
    (svg, *_) = render.render(ROOT / "examples" / "m7" / "m7_scheme_x1.0.yaml", tmp_path, with_cdxml=False)
    text = svg.read_text()
    group = text[text.index('id="coefficient-styrene"'):]
    assert 'font-style="italic"' in group[:group.index("</g>")]
