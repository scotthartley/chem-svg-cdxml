"""CDXML export, phase 7 (style-spec.md §17): polymer brackets as ChemDraw writes them in polymer_x1.0
(a ']' <graphic> with BracketUsage MultipleGroup and its two objecttags, a '[', a <bracketedgroup> with
the enclosed atoms and crossing bonds), bare chain ends as carbons, wavy ends as ExternalConnectionPoint
nodes, the random copolymer's slash as a <curve> in its fragment, ratio and coefficient captions."""

import pytest
from lxml import etree

import measure_refs as mr
from cdxml_helpers import PARSER, exported, needs_refs

POLYMERS = "examples/m7/m7_polymers_x1.0.yaml"


def _root(spec):
    return etree.fromstring(exported(spec).encode(), PARSER)


def test_bracketed_groups():
    root = _root(POLYMERS)
    page = root.find("page")
    els = {el.get("id"): el for el in page.iter() if el.get("id") and el.tag != "objecttag"}
    groups = list(page.iter("bracketedgroup"))
    assert len(groups) == 9  # one per bracket pair: 7 homopolymers and random, 2 blocks
    counts = sorted(g.get("RepeatCount") for g in groups)
    assert counts == ["1"] * 8 + ["10"]
    for g in groups:
        assert g.get("BracketUsage") == "MultipleGroup"
        inner = set(g.get("BracketedObjectIDs").split())
        assert all(els[i].tag == "n" for i in inner)
        atts = g.findall("bracketattachment")
        assert len(atts) == 2
        left, right = (els[a.get("GraphicID")] for a in atts)
        assert left.get("BracketUsage") is None and right.get("BracketUsage") == "MultipleGroup"
        lx, ly0, _, ly1 = (float(v) for v in left.get("BoundingBox").split())
        rx, ry0, _, ry1 = (float(v) for v in right.get("BoundingBox").split())
        assert lx < rx and ly0 > ly1 and ry0 < ry1  # '[' bottom first, ']' top first (polymer_x1.0)
        for a in atts:
            cb = a.find("crossingbond")
            b = els[cb.get("BondID")]
            assert b.tag == "b" and cb.get("InnerAtomID") in inner
            assert len({b.get("B"), b.get("E")} & inner) == 1
        tags = {t.get("Name"): t for t in right.findall("objecttag")}
        assert set(tags) == {"bracketusage", "parameterizedBracketLabel"}
        label = tags["parameterizedBracketLabel"]
        assert label.get("PositioningType") == "offset" and label.get("Visible") == "no"
        dx, dy = (float(v) for v in label.get("PositioningOffset").split())
        assert dx == pytest.approx(2.0, abs=0.02) and dy == pytest.approx(ry1 - ry0 + 2.0, abs=0.02)
    crossing = [{a.find("crossingbond").get("BondID") for a in g.findall("bracketattachment")} for g in groups]
    shared = [(i, j) for i in range(len(groups)) for j in range(i + 1, len(groups)) if crossing[i] & crossing[j]]
    assert shared == [(4, 5)]  # the block copolymer's junction bond crosses both blocks' brackets


def test_subscript_runs():
    root = _root(POLYMERS)
    fonts = {f.get("id"): f.get("name") for f in root.iter("font")}
    labels = [t for t in root.iter("objecttag") if t.get("Name") == "parameterizedBracketLabel"]
    runs = [[(s.text, fonts[s.get("font")], s.get("face"), s.get("size")) for s in t.iter("s")] for t in labels]
    assert [("n", "Helvetica Oblique", "2", "7.5"), ("−1", "Helvetica", None, "7.5")] in runs
    assert [("10", "Helvetica", None, "7.5")] in runs


def test_chain_ends_slash_and_captions():
    root = _root(POLYMERS)
    slash = [c for f in root.iter("fragment") for c in f.findall("curve")]
    assert len(slash) == 1 and len(slash[0].get("CurvePoints").split()) == 12
    texts = ["".join(s.text for s in t.iter("s")) for t in root.find("page").findall("t")]
    assert "70:30" in texts
    wavy = _root("examples/m7/m7_chain_end_x1.0.yaml")
    ecp, = [n for n in wavy.iter("n") if n.get("ExternalConnectionType")]  # not Ph's own ECP
    assert ecp.get("ExternalConnectionType") == "Wavy" and ecp.get("Element") is None
    scheme = _root("examples/m7/m7_scheme_x1.0.yaml")
    step, = scheme.iter("step")
    page_t = {t.get("id"): "".join(s.text for s in t.iter("s")) for t in scheme.find("page").findall("t")}
    assert "n" in page_t.values()
    cited = " ".join(step.get(k) or "" for k in step.keys() if k.startswith("ReactionStep")).split()
    assert not set(cited) & {i for i, v in page_t.items() if v == "n"}  # a free caption, as polymer_x1.0's


@needs_refs
def test_reference_writes_what_we_write():
    """polymer_x1.0's own bracket objects carry the attributes we write."""
    ref = etree.parse(str(mr.REF_DIR / "polymer_x1.0.cdxml"), PARSER).getroot()
    ours = _root(POLYMERS)
    for tag in ("bracketedgroup", "bracketattachment", "crossingbond"):
        assert {k for e in ours.iter(tag) for k in e.keys()} == {k for e in ref.iter(tag) for k in e.keys()}
    for name in ("bracketusage", "parameterizedBracketLabel"):
        r = {k for t in ref.iter("objecttag") if t.get("Name") == name for k in t.keys()}
        o = {k for t in ours.iter("objecttag") if t.get("Name") == name for k in t.keys()}
        assert o == r
    r = {k for g in ref.iter("graphic") if g.get("GraphicType") == "Bracket" for k in g.keys()}
    o = {k for g in ours.iter("graphic") if g.get("GraphicType") == "Bracket" for k in g.keys()}
    assert o == r
    assert {(n.get("NodeType"), n.get("ExternalConnectionType")) for n in ref.iter("n") if n.get("ExternalConnectionType")} == \
        {("ExternalConnectionPoint", "Wavy")}
