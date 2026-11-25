"""A project's own style, groups and structures beside its specs (references/extending.md): a preset
by path with its CDXML header beside it, `abbreviations:` from files, molfiles relative to the spec."""

import math
import shutil
from pathlib import Path

import pytest
import yaml
from lxml import etree

from chem_svg_cdxml import abbrev, cdxml_writer
from chem_svg_cdxml.figure import build_figure, spec_cdxml
from chem_svg_cdxml.preset import PRESET_DIR, load_preset, resolve
from chem_svg_cdxml.spec import load_spec

REPO = Path(__file__).resolve().parents[1]
BOND = 20.0


def _project(tmp: Path, header: bool = True) -> Path:
    """tmp/ as a user's project: styles/mystyle.yaml (acs1996, bond length 20) with styles/cdxml/,
    groups.yaml, structures/alanine.mol and figures/fig.yaml using all three."""
    styles = tmp / "styles"
    (styles / "cdxml").mkdir(parents=True)
    raw = yaml.safe_load((PRESET_DIR / "acs1996.yaml").read_text())
    raw["name"] = "mystyle"
    raw["lengths"]["bond_length"] = BOND
    (styles / "mystyle.yaml").write_text(yaml.safe_dump(raw))
    if header:
        text = (PRESET_DIR / "cdxml" / "acs1996.xml").read_text().replace('BondLength="14.40"', f'BondLength="{BOND:.2f}"')
        (styles / "cdxml" / "mystyle.xml").write_text(text)
        shutil.copy(PRESET_DIR / "cdxml" / "acs1996.nicknames.xml", styles / "cdxml" / "mystyle.nicknames.xml")
    (tmp / "groups.yaml").write_text('OCy: {smiles: "*OC1CCCCC1"}\nXy: {kind: text}\n')
    (tmp / "structures").mkdir()
    shutil.copy(REPO / "templates" / "m1" / "m1_set_x1.0__3.mol", tmp / "structures" / "alanine.mol")
    (tmp / "figures").mkdir()
    spec = tmp / "figures" / "fig.yaml"
    spec.write_text(
        "preset: ../styles/mystyle.yaml\n"
        "abbreviations: [../groups.yaml, {Xy: {kind: text, text: '*Xy*'}}]\n"
        "species:\n"
        "  ether: {smiles: 'COc1ccc(OC2CCCCC2)cc1', abbreviate: [OCy]}\n"
        "  ala: {molfile: ../structures/alanine.mol}\n")
    return spec


def _bond_lengths(mol):
    return [math.dist((mol.atoms[b.a].x, mol.atoms[b.a].y), (mol.atoms[b.b].x, mol.atoms[b.b].y)) for b in mol.bonds]


def test_user_preset_groups_and_structures(tmp_path):
    spec_path = _project(tmp_path)
    spec = load_spec(spec_path)
    assert Path(spec.preset) == (tmp_path / "styles" / "mystyle.yaml").resolve()
    assert "OCy" in spec.registry and spec.registry.get("Xy").text == "*Xy*"  # inline after the file
    assert "OCy" not in abbrev.default()
    fig = build_figure(spec_path)
    assert {p.name for p in fig.placed} == {"ether", "ala"}
    for p in fig.placed:  # the molfile is rescaled to the style's bond length too
        assert _bond_lengths(p.mol) == pytest.approx([BOND] * len(p.mol.bonds), abs=0.05), p.name
    root = etree.fromstring(spec_cdxml(spec_path, warn=lambda m: None).encode())
    assert root.get("BondLength") == f"{BOND:.2f}"  # styles/cdxml/mystyle.xml, not presets/cdxml/
    assert "OCy" in ["".join(t.itertext()).strip() for t in root.iter("t")]  # groups.yaml's label, a generated fragment
    assert any(n.get("NodeType") == "Fragment" for n in root.iter("n"))


def test_resolve_names_and_paths(tmp_path):
    assert resolve("acs1996") == (PRESET_DIR / "acs1996.yaml", PRESET_DIR / "cdxml")
    p = tmp_path / "styles" / "x.yaml"
    assert resolve(str(p)) == (p, tmp_path / "styles" / "cdxml")


def test_missing_header_names_the_expected_file(tmp_path):
    spec_path = _project(tmp_path, header=False)
    build_figure(spec_path)  # the SVG needs no header
    with pytest.raises(FileNotFoundError, match=r"styles/cdxml/mystyle\.xml"):
        spec_cdxml(spec_path, warn=lambda m: None)


def test_preset_path_relative_to_the_spec_not_the_cwd(tmp_path, monkeypatch):
    spec_path = _project(tmp_path)
    monkeypatch.chdir(REPO)
    assert load_preset(load_spec(spec_path).preset).lengths.bond_length == BOND


@pytest.fixture(autouse=True)
def _clear_caches():
    yield
    cdxml_writer.load_nicknames.cache_clear()
