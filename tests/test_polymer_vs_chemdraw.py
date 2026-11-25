"""Polymers against polymer_x1.0 (style-spec.md §17): the measured values in preset.polymer (the wave,
the lip rule, the stretch, the subscript size) and the hand-placed defaults within 0.5 pt of what the
reference shows where it is consistent (measure_refs.py --polymer)."""

import pytest

from chemdraw_refs.polymer import polymer_summary
from chem_svg_cdxml.preset import load_preset
from cdxml_helpers import needs_refs

P = load_preset("acs1996").polymer


@pytest.fixture(scope="module")
def ref():
    return polymer_summary("polymer_x1.0")


@needs_refs
def test_chemdraw_generated_values(ref):
    s = ref["summary"]
    assert s["lip_per_height"]["min"] == pytest.approx(P.lip_per_height, abs=0.001)
    assert s["lip_per_height"]["max"] == pytest.approx(P.lip_per_height, abs=0.001)
    assert s["stretched_bonds"]["median"] == pytest.approx(P.crossing_stretch, abs=0.002)
    w, = ref["wavy"]
    assert w["length"] == pytest.approx(P.wave_length, abs=0.05)
    assert w["half_period"] == pytest.approx(P.wave_half_period, abs=0.002)
    assert w["amplitude"] == pytest.approx(P.wave_amplitude, abs=0.005)
    subs = [b["subscript"] for g in ref["groups"] for b in g["brackets"] if "subscript" in b]
    assert {x["size"] for x in subs} == {P.label_size}


@needs_refs
def test_placement_defaults(ref):
    s = ref["summary"]
    assert s["subscript_drop"]["median"] == pytest.approx(P.label_drop, abs=0.5)
    assert s["subscript_dx"]["median"] == pytest.approx(P.label_dx, abs=0.5)
    assert s["crossing_fraction"]["median"] == pytest.approx(0.5, abs=0.05)  # unstretched: the midpoint
    sl, = ref["slashes"]
    assert sl["length"] == pytest.approx(P.slash_length, abs=0.5)
    assert sl["angle_from_bond_clockwise"] == pytest.approx(P.slash_angle, abs=0.5)
    assert sl["on_bond_fraction"] == pytest.approx(0.5, abs=0.01)
    ratio = next(c for c in ref["captions"] if c["text"] == "70:30")
    assert ratio["baseline_below_brackets"] == pytest.approx(P.ratio_drop, abs=0.5)
    assert abs(ratio["center_vs_brackets"]) < 1.0  # centered, by eye
    assert ref["coefficient"]["gap_to_species"] == pytest.approx(P.coefficient_gap, abs=0.5)
    assert abs(ref["coefficient"]["center_above_axis"]) < 0.05
