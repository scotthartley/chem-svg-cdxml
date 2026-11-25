"""Polymers (style-spec.md §17; polymer.py): assembly, the backbone layout, stretched crossing bonds,
brackets, subscripts, the random copolymer's slash and ratio, block junctions, wavy chain ends,
per-bond lengths and reactant coefficients. Against ChemDraw: test_polymer_vs_chemdraw.py."""

import math
from pathlib import Path

import pytest
import yaml
from layouts import example

from chem_svg_cdxml import polymer
from chem_svg_cdxml.figure import build_figure
from chem_svg_cdxml.layout import build_scheme
from chem_svg_cdxml.preset import load_preset
from chem_svg_cdxml.spec import Polymer, Spec, molecules

P = load_preset("acs1996")
L = P.lengths.bond_length


def _spec(species: dict, **extra) -> Spec:
    return Spec.model_validate({"species": species, **extra})


@pytest.fixture(scope="module")
def mols():
    spec = Spec.model_validate(yaml.safe_load(example("m7_polymers_x1.0").read_text()))
    return {m.name: m for m in molecules(spec, Path("."))}


def _crossings(m):
    """(length / L, angle deg (y up) from the inner atom) of each crossing bond, left to right."""
    out = []
    for pr in m.polymer.pairs:
        for bi, inner in pr.crossing:
            b = m.bonds[bi]
            o = b.b if b.a == inner else b.a
            d = (m.atoms[o].x - m.atoms[inner].x, m.atoms[o].y - m.atoms[inner].y)
            out.append((math.hypot(*d) / L, math.degrees(math.atan2(-d[1], d[0]))))
    return out


def test_assemble_joins_head_to_tail_with_bare_ends():
    mol = polymer.assemble(Polymer(units=["*CC(*)c1ccccc1"]))
    ends = [a for a in mol.GetAtoms() if a.HasProp(polymer.POLY_END)]
    assert [a.GetAtomicNum() for a in ends] == [0, 0]
    assert sorted(a.GetIntProp(polymer.POLY_SEG) for a in ends) == [0, 2]
    assert sum(a.GetAtomicNum() == 6 for a in mol.GetAtoms()) == 8


def test_assemble_end_groups_replace_the_bare_ends():
    mol = polymer.assemble(Polymer(start="N#CC(C)(C)*", units=["*CC(*)c1ccccc1"], end="Br*"))
    assert not any(a.GetAtomicNum() == 0 for a in mol.GetAtoms())
    assert {a.GetSymbol() for a in mol.GetAtoms() if a.GetIntProp(polymer.POLY_SEG) == 2} == {"Br"}


@pytest.mark.parametrize("p, msg", [
    (dict(units=["*CCc1ccccc1"]), "needs two '*'"),
    (dict(units=["*CC(*)C"], start="CC"), "needs one '*'"),
    (dict(units=["*[CH2:1]C(*)C"], end="[CH3:1]*"), "unique across"),
])
def test_assemble_errors(p, msg):
    with pytest.raises(ValueError, match=msg):
        polymer.assemble(Polymer(**p))


def test_backbone_horizontal_pendants_below(mols):
    m = mols["ps_x"]
    ring = [a for a in m.atoms if a.seg == 1 and any(a.idx in r for r in m.rings)]
    backbone = [a for a in m.atoms if a.seg == 1 and not any(a.idx in r for r in m.rings)]
    assert min(a.y for a in ring) > max(a.y for a in backbone)  # y down: the phenyl hangs below


def test_crossing_bonds_kept_or_stretched_as_polymer_x1_0(mols):
    """PS: the left crossing bond stays 1 x at the zigzag's 150 deg, the bracket at its midpoint; the
    right one passes the phenyl, so it is 2 x with its rise kept (polymer_x1.0: 28.8 pt at 15 deg)."""
    (l_len, l_ang), (r_len, r_ang) = _crossings(mols["ps_x"])
    assert l_len == pytest.approx(1.0, abs=1e-3) and l_ang == pytest.approx(-150.0, abs=0.5) or \
        l_ang == pytest.approx(150.0, abs=0.5)
    assert r_len == pytest.approx(P.polymer.crossing_stretch, abs=1e-3)
    assert abs(r_ang) == pytest.approx(math.degrees(math.asin(0.25)), abs=0.5)  # rise L/2 over 2L: 14.5 deg
    m = mols["ps_x"]
    pr = m.polymer.pairs[0]
    bi, inner = pr.crossing[0]
    b = m.bonds[bi]
    mid = (m.atoms[b.a].x + m.atoms[b.b].x) / 2
    assert pr.x[0] == pytest.approx(mid, abs=1e-6)


def test_ring_backbone_unstretched(mols):
    assert [round(x, 3) for x, _ in _crossings(mols["ppp"])] == [1.0, 1.0]


def test_brackets_share_height_and_lip(mols):
    for name in ("ps_ends", "pmma", "random", "block"):
        m = mols[name]
        tops = {round(p.top, 6) for p in m.polymer.pairs}
        bottoms = {round(p.bottom, 6) for p in m.polymer.pairs}
        assert len(tops) == len(bottoms) == 1
        enclosed = [m.atoms[i] for p in m.polymer.pairs for i in p.atoms]
        assert m.polymer.pairs[0].top < min(a.y for a in enclosed)
        assert m.polymer.pairs[0].bottom > max(a.y for a in enclosed)


def test_subscript_at_the_right_bracket(mols):
    pr = mols["ps_n1"].polymer.pairs[0]
    assert pr.label_origin == pytest.approx((pr.x[1] + P.polymer.label_dx, pr.bottom + P.polymer.label_drop))
    runs = polymer.label_runs("n-1")
    assert [(r.text, r.style) for r in runs] == [("n", "italic"), ("−1", "regular")]
    assert [(r.text, r.style) for r in polymer.label_runs("10")] == [("10", "regular")]


def test_random_slash_and_ratio(mols):
    m = mols["random"]
    (bi, a, b), = m.polymer.slashes
    bond = m.bonds[bi]
    p, q = m.atoms[bond.a].p, m.atoms[bond.b].p
    assert math.dist(p, q) == pytest.approx(P.polymer.crossing_stretch * L, abs=0.05)
    assert ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2) == pytest.approx(((p[0] + q[0]) / 2, (p[1] + q[1]) / 2))
    assert math.dist(a, b) == pytest.approx(P.polymer.slash_length)
    left, right = sorted((p, q))
    turn = math.degrees(math.atan2(b[1] - a[1], b[0] - a[0]) - math.atan2(right[1] - left[1], right[0] - left[0])) % 180
    assert turn == pytest.approx(P.polymer.slash_angle, abs=1e-6)
    assert len(m.polymer.pairs) == 1
    ink = polymer.ratio_line(m.polymer, P).ink()
    assert (ink[0] + ink[2]) / 2 == pytest.approx(sum(m.polymer.pairs[0].x) / 2, abs=1e-6)


def test_block_brackets_apart_on_the_junction(mols):
    m = mols["block"]
    a, b = m.polymer.pairs
    assert [p.label for p in (a, b)] == ["m", "n"]
    assert b.x[0] - a.x[1] >= P.polymer.block_gap - 1e-6
    assert a.crossing[1][0] == b.crossing[0][0]  # one bond crosses both brackets


def test_marks_drawn_in_the_species_group():
    fig = build_figure(example("m7_polymers_x1.0"))
    ids = [g.id for sp in fig.root.children for g in sp.children]
    assert "species-random-polymer" in ids
    poly = next(g for sp in fig.root.children for g in sp.children if g.id == "species-random-polymer")
    srcs = {leaf.source for leaf in poly.leaves()}
    assert any(s.startswith("polymer-slash:") for s in srcs) and "polymer-ratio" in srcs


def test_wavy_chain_end():
    spec = Spec.model_validate(yaml.safe_load(example("m7_chain_end_x1.0").read_text()))
    m, = molecules(spec, Path("."))
    end = next(a for a in m.atoms if a.end == "wavy")
    assert end.element == "C" and not end.num_h
    path, = polymer.wavy_ends(m, P)
    box = path.ink_bbox()
    nb = m.atoms[m.neighbors(end.idx)[0][0]]
    bond = (end.x - nb.x, end.y - nb.y)  # the squiggle runs across the bond, centered on its end
    chord = (path.segments[-1][-1][0] - path.segments[0][1][0], path.segments[-1][-1][1] - path.segments[0][1][1])
    assert abs(bond[0] * chord[0] + bond[1] * chord[1]) / math.hypot(*bond) / math.hypot(*chord) < 0.05  # ends on a crest, as ChemDraw's
    span = math.hypot(box[2] - box[0], box[3] - box[1])
    assert span == pytest.approx(math.hypot(P.polymer.wave_length, 2 * P.polymer.wave_amplitude), abs=1.5)
    assert ((box[0] + box[2]) / 2, (box[1] + box[3]) / 2) == pytest.approx(end.p, abs=0.6)


def test_bond_lengths_override():
    spec = _spec({"x": {"smiles": "[CH3:1][CH2:2][CH2:3]c1ccccc1", "bond_lengths": {"2-3": 2.0}}})
    m, = molecules(spec, Path("."))
    by = {a.map: a for a in m.atoms if a.map}
    assert math.dist(by[2].p, by[3].p) == pytest.approx(2 * L, abs=1e-6)
    assert math.dist(by[1].p, by[2].p) == pytest.approx(L, abs=0.01)


def test_bond_lengths_ring_bond_refused():
    spec = _spec({"x": {"smiles": "[CH2:1]1[CH2:2]CCCC1", "bond_lengths": {"1-2": 2.0}}})
    with pytest.raises(ValueError, match="acyclic"):
        molecules(spec, Path("."))


def test_coefficient_before_its_reactant():
    lay = build_scheme(Spec.model_validate(yaml.safe_load(example("m7_scheme_x1.0").read_text())), Path("."))
    ln, = lay.coefficients
    ink = ln.ink()
    sty = lay.species["styrene"].group.ink_bbox()
    assert sty[0] - ink[2] == pytest.approx(P.polymer.coefficient_gap)
    assert (ink[1] + ink[3]) / 2 == pytest.approx(lay.steps[0].tail[1], abs=1e-6)  # centered on the arrow axis
    assert [r.style for r in ln.runs] == ["normal"]  # roman, as polymer_x1.0's n


@pytest.mark.parametrize("p, msg", [
    (dict(units=["*CC(*)C", "*CC(*)C"]), "kind homo takes one unit"),
    (dict(kind="random", units=["*CC(*)C"]), "two or more units"),
    (dict(units=["*CC(*)C"], ratio="1:1"), "random copolymer only"),
    (dict(kind="block", units=["*CC(*)C", "*CC(*)C"], label="n"), "go on its units"),
])
def test_polymer_schema(p, msg):
    with pytest.raises(ValueError, match=msg):
        Polymer(**p)


def _angles(m, i):
    """Bond directions (deg, y up, 0-360) round atom i."""
    a = m.atoms[i]
    return sorted(math.degrees(math.atan2(-(m.atoms[j].y - a.y), m.atoms[j].x - a.x)) % 360 for j, _ in m.neighbors(i))


def _on_grid(m, ids):
    """Every bond round the atoms `ids` on the 30 deg grid, and no two bonds at one atom under 60 deg."""
    for i in ids:
        ang = _angles(m, i)
        assert all(min(x % 30, 30 - x % 30) < 0.5 for x in ang), (m.atoms[i].element, ang)
        gaps = [(b - a) for a, b in zip(ang, ang[1:] + [ang[0] + 360])]
        assert min(gaps) > 59.5, (m.atoms[i].element, ang)


def test_ester_pendants_at_120(mols):
    """PMMA's ester as polymer_x1.0 draws it: C=O down-left, O down-right, its methyl straight down,
    every unstretched bond on the 30 deg grid; the quaternary carbon's methyl up, between the
    backbone bonds at 60 deg."""
    for name in ("pmma", "random", "block"):
        m = mols[name]
        stretched = {bi for pr in m.polymer.pairs for bi, _ in pr.crossing} | {bi for bi, *_ in m.polymer.slashes}
        ids = [a.idx for a in m.atoms if a.end is None and not any(b.idx in stretched for b in m.bonds
                                                                   if a.idx in (b.a, b.b))]
        _on_grid(m, ids)
    m = mols["pmma"]
    o2 = next(a for a in m.atoms if a.element == "O" and any(b.order == 2 for b in m.bonds if a.idx in (b.a, b.b)))
    c = m.atoms[m.neighbors(o2.idx)[0][0]]
    assert sorted(_angles(m, c.idx)) == pytest.approx([90, 210, 330], abs=0.5)


def test_start_group_continues_the_zigzag(mols):
    """AIBN's quaternary carbon: the nitrile continues the backbone zigzag (straight through its
    sp carbon), one methyl between the backbone bonds at 60 deg, the other 120 deg from both."""
    m = mols["ps_ends"]
    cq = next(a for a in m.atoms if a.seg == 0 and len(m.neighbors(a.idx)) == 4)
    ang = _angles(m, cq.idx)
    assert ang == pytest.approx([30, 90, 150, 270], abs=0.5)
    cn = next(a for a in m.atoms if a.seg == 0 and any(b.order == 3 for b in m.bonds if a.idx in (b.a, b.b))
              and a.element == "C")
    assert _angles(m, cn.idx) == pytest.approx([150, 330], abs=0.5)


def test_chain_end_drawn_as_a_growing_chain():
    """polymer_x1.0's ~~CH2-CH(Ph)•: the chain end on the left, the CH2 up, Ph straight down."""
    spec = Spec.model_validate(yaml.safe_load(example("m7_chain_end_x1.0").read_text()))
    m, = molecules(spec, Path("."))
    end = next(a for a in m.atoms if a.end == "wavy")
    ch2 = m.atoms[m.neighbors(end.idx)[0][0]]
    ch = next(m.atoms[j] for j, _ in m.neighbors(ch2.idx) if j != end.idx)
    ph = next(m.atoms[j] for j, _ in m.neighbors(ch.idx) if j != ch2.idx)
    assert end.x < ch2.x < ch.x and ch2.y < end.y
    assert ph.x == pytest.approx(ch.x, abs=1e-6) and ph.y > ch.y


def test_chain_end_pendants_kept_apart_by_a_stretched_bond():
    """~~CH2CH(Ph)CH2CH(Ph)•: both Ph straight down, the bond leaving the first CH(Ph) CHAIN_STRETCH x
    long (rise kept) so the rings don't meet."""
    spec = _spec({"x": {"smiles": "[*:1][CH2:2][CH:3](c1ccccc1)[CH2:4][CH:5]c1ccccc1", "chain_ends": [1]}})
    m, = molecules(spec, Path("."))
    by = {a.map: a for a in m.atoms if a.map}
    for k in (3, 5):
        ph = next(m.atoms[j] for j, _ in m.neighbors(by[k].idx) if m.atoms[j].map is None or m.atoms[j].map == 0)
        assert ph.x == pytest.approx(by[k].x, abs=1e-6) and ph.y > by[k].y
    assert math.dist(by[3].p, by[4].p) == pytest.approx(polymer.CHAIN_STRETCH * L, abs=0.01)
    assert math.dist(by[4].p, by[5].p) == pytest.approx(L, abs=0.01)


def test_chain_stretch_is_the_presets():
    for name in ("acs1996", "acs1996_x0.8", "acs1996_x0.6", "presentation_x1.0"):
        assert load_preset(name).polymer.crossing_stretch == polymer.CHAIN_STRETCH
