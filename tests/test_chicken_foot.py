"""Acyclic four-coordinate centers with three terminal heavy atoms (tBu, SiMe3, CF3) drawn as a
"chicken foot" (chem._chicken_feet), and a two-letter symbol centered when two of its bonds are
opposite at any orientation (draw_mol.centered_label). The user, 2026-11-25; unmeasured.
inbox 2026-11-17."""

import math
from pathlib import Path

import pytest

from chem_svg_cdxml import chem
from chem_svg_cdxml.draw_mol import Drawer
from chem_svg_cdxml.preset import load_preset
from chem_svg_cdxml.spec import load_spec, scheme_molecules

P = load_preset("acs1996")
BL = P.lengths.bond_length
SPECS = Path(__file__).parent / "specs"
TMS_ARYL = "Cc1ccc(cc1)C#C[Si](C)(C)C"


def _deg(m, i: int, j: int) -> float:
    a, b = m.atoms[i], m.atoms[j]
    return math.degrees(math.atan2(-(b.y - a.y), b.x - a.x))


def _feet(m, i: int) -> list[float]:
    """Directions of the terminal bonds round center i, relative to its one non-terminal bond."""
    nbs = [j for j, _ in m.neighbors(i)]
    (inc,) = [j for j in nbs if m.degree(j) > 1]
    t = _deg(m, i, inc)
    return sorted((_deg(m, i, j) - t) % 360 for j in nbs if j != inc)


@pytest.mark.parametrize("smiles, center", [
    (TMS_ARYL, "Si"), ("C[Si](C)(C)C#C[Cu]", "Si"), ("CC(C)(C)c1ccccc1", "C"),
    ("FC(F)(F)c1ccccc1", "C"), ("C[Si](C)(C)OC1CCCCC1", "Si"),
])
def test_chicken_foot(smiles, center):
    """With the bond in at 0°, the three terminal bonds at 120°, 180° and 240°, all 14.4 pt."""
    m = chem.from_smiles(smiles, "m", BL)
    (i,) = [a.idx for a in m.atoms if a.element == center and m.degree(a.idx) == 4]
    assert _feet(m, i) == pytest.approx([120.0, 180.0, 240.0], abs=0.5)
    assert [math.dist(m.atoms[i].p, m.atoms[j].p) for j, _ in m.neighbors(i)] == pytest.approx([BL] * 4, abs=0.01)


def test_tms_label_centered():
    m = chem.from_smiles(TMS_ARYL, "m", BL)
    d = Drawer(m, P)
    si = next(a for a in m.atoms if a.element == "Si")
    lab = d.labels[si.idx]
    assert lab.h_side == "center"
    assert lab.elem_origin[0] + d.face.advance("Si", lab.size) / 2 == pytest.approx(si.x, abs=1e-9)


def test_tilted_centered_label_cdxml_on_atom():
    """A centered label on a tilted bond pair writes Center justification with p.x on the atom,
    as the cross does (cdxml_writer)."""
    from lxml import etree
    from chem_svg_cdxml.cdxml_writer import to_cdxml
    m = chem.from_smiles(TMS_ARYL, "m", BL)
    d = Drawer(m, P)
    si = next(a for a in m.atoms if a.element == "Si")
    root = etree.fromstring(to_cdxml([(m, d, (0.0, 0.0))], "acs1996", warn=lambda msg: None).encode())
    t = next(t for t in root.iter("t") if "".join(x.text for x in t.iter("s")) == "Si")
    assert t.get("LabelJustification") == t.get("Justification") == "Center"
    assert float(t.get("p").split()[0]) == pytest.approx(si.x, abs=0.006)


def test_ring_quaternary_carbon_unchanged():
    """1,1-Dimethylcyclohexane: the center is in a ring, so CoordGen's layout stands."""
    from rdkit.Chem import rdDepictor
    smi = "CC1(C)CCCCC1"
    mol = chem.parse_smiles(smi)
    chem.layout_2d(mol)
    ref = chem.parse_smiles(smi)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(ref)
    assert chem._coords(mol) == pytest.approx(chem._coords(ref), abs=1e-9)


def test_sonogashira_cycle_si_centered():
    """Every Si in the Sonogashira cycle (inbox 2026-11-17), transm pinned to oxadd included, is a
    chicken foot with its symbol centered."""
    path = SPECS / "sonogashira_cycle_x1.0.yaml"
    mols = scheme_molecules(load_spec(path), path.parent)
    seen = 0
    for name, m in mols.items():
        d = Drawer(m, P)
        for a in m.atoms:
            if a.element == "Si":
                seen += 1
                assert _feet(m, a.idx) == pytest.approx([120.0, 180.0, 240.0], abs=0.5), name
                assert d.labels[a.idx].h_side == "center", name
    assert seen == 4


def test_species_override_keeps_coordgen():
    """`chicken_foot: false` on a species keeps CoordGen's layout (m4_nitration's H2SO4, drawn
    like bisulfate; the user, 2026-11-25)."""
    from chem_svg_cdxml.spec import Spec, molecules
    smi = "[O:5]=[S:6](=[O:7])([OH:8])[O:9][H:10]"

    def s_dirs(flag):
        spec = Spec.model_validate({"species": {"a": {"smiles": smi, "chicken_foot": flag}}})
        (m,) = molecules(spec, Path("."))
        (s,) = [a.idx for a in m.atoms if a.element == "S"]
        return _feet(m, s)

    assert s_dirs(True) == pytest.approx([120.0, 180.0, 240.0], abs=0.5)
    assert s_dirs(False) != pytest.approx([120.0, 180.0, 240.0], abs=0.5)
