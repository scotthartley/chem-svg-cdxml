"""Polymers (style-spec.md §17): polymer_x1.0's repeat-unit brackets, subscripts, the random copolymer's
slash and ratio, the wavy chain end and the scheme's coefficient. The brackets there were placed by eye
(round defaults in preset.polymer); the wave and the CDXML objects are ChemDraw's own.

All positions are in the SVG frame (pt, y down); the windows scale with the file's BondLength and
LabelSize, so the same code reads a probe at another style."""

from __future__ import annotations

import math
from pathlib import Path

from lxml import etree

from chem_svg_cdxml.geometry import length, sub

from .core import REF_DIR, assign, label_glyphs, load, stats, union
from .measurements import _brackets

POLYMER_REFS = ["polymer_x1.0"]


def _ink(doc, asg, ids: set[str]) -> tuple:
    """The drawn ink of the atoms `ids` and the bonds between them (assigned shapes, label glyphs)."""
    boxes = []
    for b in doc.bonds:
        if b.a in ids and b.b in ids and b.id in asg:
            r = asg[b.id]
            boxes += [ln.shape.bbox for ln in r.get("lines", [])] + [sh.bbox for sh, _ in r.get("shapes", [])]
    for a in ids:
        if a in doc.atoms:
            boxes += [g.box for g in label_glyphs(doc, doc.atoms[a])[1]]
    return union(boxes)


def polymer_summary(name: str, ref_dir: Path = REF_DIR) -> dict:
    """Per <bracketedgroup> (left, then right bracket): the bar's x, top, bottom and height, the lip
    and lip / height, the gap from the enclosed ink (x) and above / below it, the crossing bond's
    length (in bond lengths), its angle and where the bracket crosses it (fraction from the inner
    atom), the subscript (text, size, italic runs) and its origin from the ']' bar and bottom.
    Then the slash (a <curve> in a fragment: length, angle to its bond, where on the bond), the
    ratio caption (center vs the brackets', baseline below them), the wavy chain end (length,
    half period, amplitude, center vs the node) and the coefficient (gap to its species, center
    vs the arrow axis)."""
    doc = load(name, ref_dir)
    bl = float(doc.settings["BondLength"])
    off = doc.offset
    root = etree.parse(str(ref_dir / f"{name}.cdxml")).getroot()
    graphics = {g.get("id"): g for g in root.iter("graphic")}
    asg = assign(doc)
    drawn = _brackets(doc)
    bonds = {b.id: b for b in doc.bonds}
    groups = []
    for bg in root.iter("bracketedgroup"):
        ids = set(bg.get("BracketedObjectIDs").split())
        ink = _ink(doc, asg, ids)
        out = {"repeat_count": bg.get("RepeatCount"), "usage": bg.get("BracketUsage"), "ink": [round(v, 2) for v in ink],
               "brackets": []}
        for att in bg.iter("bracketattachment"):
            g = graphics[att.get("GraphicID")]
            bb = [float(v) for v in g.get("BoundingBox").split()]
            x, top, bottom = bb[0] + off[0], min(bb[1], bb[3]) + off[1], max(bb[1], bb[3]) + off[1]
            seg = next((d for d in drawn if abs(d["x"] - x) < 0.05 and abs(d["y_top"] - top) < 0.05), None)
            cb = att.find("crossingbond")
            bd = bonds[cb.get("BondID")]
            inner = cb.get("InnerAtomID")
            outer = bd.b if bd.a == inner else bd.a
            pi, po = doc.atoms[inner].p, doc.atoms[outer].p
            side = "left" if po[0] < pi[0] else "right"
            item = {"side": side, "x": round(x, 3), "top": round(top, 3), "bottom": round(bottom, 3),
                    "height": round(bottom - top, 3),
                    "lip": seg["lip_length"]["median"] if seg else None,
                    "lip_per_height": round(seg["lip_length"]["median"] / (bottom - top), 4) if seg else None,
                    "stroke_width": seg["stroke_width"] if seg else None,
                    "gap_x": round(ink[0] - x if side == "left" else x - ink[2], 3),
                    "gap_above": round(ink[1] - top, 3), "gap_below": round(bottom - ink[3], 3),
                    "crossing_bond_per_bond_length": round(length(sub(po, pi)) / bl, 3),
                    "crossing_angle": round(math.degrees(math.atan2(-(po[1] - pi[1]), po[0] - pi[0])), 1),
                    "crossing_fraction": round((x - pi[0]) / (po[0] - pi[0]), 3) if abs(po[0] - pi[0]) > 1e-6 else None,
                    "outer_label": doc.atoms[outer].label}
            tag = g.find(".//objecttag[@Name='parameterizedBracketLabel']")
            if tag is not None:
                t = tag.find("t")
                p = [float(v) for v in t.get("p").split()]
                runs = [(s.text, float(s.get("size")), s.get("face") == "2") for s in t.iter("s")]
                item["subscript"] = {"text": "".join(r[0] for r in runs), "size": runs[0][1],
                                     "italic_runs": [r[0] for r in runs if r[2]],
                                     "origin_dx": round(p[0] + off[0] - x, 3), "baseline_drop": round(p[1] + off[1] - bottom, 3),
                                     "positioning_offset": tag.get("PositioningOffset"),
                                     "bracketusage_text": "".join(s.text for s in g.find(".//objecttag[@Name='bracketusage']").iter("s"))}
            out["brackets"].append(item)
        groups.append(out)
    res = {"settings": {k: doc.settings[k] for k in ("BondLength", "LineWidth", "LabelSize")}, "groups": groups}
    br = [b for g in groups for b in g["brackets"]]
    res["summary"] = {k: stats([b[k] for b in br]) for k in
                      ("lip_per_height", "gap_x", "gap_above", "gap_below", "crossing_fraction")}
    res["summary"]["stretched_bonds"] = stats([b["crossing_bond_per_bond_length"] for b in br
                                               if 1.5 < b["crossing_bond_per_bond_length"] < 2.5])
    subs = [b["subscript"] for b in br if "subscript" in b]
    res["summary"]["subscript_dx"] = stats([s["origin_dx"] for s in subs])
    res["summary"]["subscript_drop"] = stats([s["baseline_drop"] for s in subs])
    res["slashes"] = _slashes(doc, root, off)
    res["wavy"] = _wavy(doc, root, off)
    res["captions"] = []
    for el in root.iter("t"):
        if el.getparent().tag != "page":
            continue
        bb = [float(v) for v in el.get("BoundingBox").split()]
        bb = [bb[0] + off[0], bb[1] + off[1], bb[2] + off[0], bb[3] + off[1]]
        p = [float(v) for v in el.get("p").split()]
        c = {"text": "".join(s.text or "" for s in el.iter("s")), "bbox": [round(v, 2) for v in bb],
             "p": [round(p[0] + off[0], 2), round(p[1] + off[1], 2)]}
        for g in groups:  # a ratio below a bracket pair: center vs theirs, baseline below their bottom
            b = g["brackets"]
            if b[0]["x"] < (bb[0] + bb[2]) / 2 < b[-1]["x"] and 0 < bb[1] - b[0]["bottom"] < 20:
                c["center_vs_brackets"] = round((bb[0] + bb[2]) / 2 - (b[0]["x"] + b[-1]["x"]) / 2, 3)
                c["baseline_below_brackets"] = round(c["p"][1] - b[0]["bottom"], 3)
        res["captions"].append(c)
    res["coefficient"] = _coefficient(doc, groups)
    return res


def _slashes(doc, root, off) -> list[dict]:
    out = []
    for frag in root.iter("fragment"):
        for c in frag.findall("curve"):
            v = [float(t) for t in c.get("CurvePoints").split()]
            p0, p1 = (v[0] + off[0], v[1] + off[1]), (v[-2] + off[0], v[-1] + off[1])
            mid = ((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2)
            best = None
            for b in doc.bonds:
                a, e = doc.atoms[b.a].p, doc.atoms[b.b].p
                if a[0] > e[0]:
                    a, e = e, a
                L = length(sub(e, a))
                t = ((mid[0] - a[0]) * (e[0] - a[0]) + (mid[1] - a[1]) * (e[1] - a[1])) / L ** 2
                d = abs((e[0] - a[0]) * (mid[1] - a[1]) - (e[1] - a[1]) * (mid[0] - a[0])) / L
                if 0 <= t <= 1 and (best is None or d < best[0]):
                    best = (d, t, a, e, L)
            d, t, a, e, L = best
            bond_ang = math.degrees(math.atan2(e[1] - a[1], e[0] - a[0]))
            sl_ang = math.degrees(math.atan2(p1[1] - p0[1], p1[0] - p0[0]))
            out.append({"length": round(length(sub(p1, p0)), 3), "angle_from_horizontal": round(sl_ang, 1),
                        "angle_from_bond_clockwise": round((sl_ang - bond_ang) % 180, 1),
                        "on_bond_fraction": round(t, 3), "off_bond": round(d, 3),
                        "bond_per_bond_length": round(L / float(doc.settings["BondLength"]), 3),
                        "points": c.get("CurvePoints")})
    return out


def _wavy(doc, root, off) -> list[dict]:
    out = []
    for n in root.iter("n"):
        if n.get("ExternalConnectionType") != "Wavy":
            continue
        c = doc.atoms[n.get("id")].p
        nb = next(b.a if b.b == n.get("id") else b.b for b in doc.bonds if n.get("id") in (b.a, b.b))
        d = sub(c, doc.atoms[nb].p)
        u = (-d[1] / length(d), d[0] / length(d))  # across the bond
        sq = [s for s in doc.shapes if s.kind == "stroke" and s.curved
              and min(length(sub(p, c)) for p in s.pts) < 2.0]
        if not sq:
            continue
        w = (d[0] / length(d), d[1] / length(d))  # along the bond
        loc = [((p[0] - c[0]) * u[0] + (p[1] - c[1]) * u[1], (p[0] - c[0]) * w[0] + (p[1] - c[1]) * w[1])
               for p in sq[0].pts]
        anchors = loc[0::3]
        zeros = [a for a, v in anchors if abs(v) < 0.05]
        peaks = [abs(v) for _, v in anchors if abs(v) > 0.05]
        s = [a for a, _ in loc]
        out.append({"length": round(max(s) - min(s), 3), "center_off_node": round((max(s) + min(s)) / 2, 3),
                    "half_period": round((zeros[-1] - zeros[0]) / (len(zeros) - 1), 4) if len(zeros) > 1 else None,
                    "amplitude": round(sum(peaks) / len(peaks), 3) if peaks else None,
                    "half_waves": len(loc) // 3 / 2, "stroke_width": round(sq[0].width, 3),
                    "node": n.get("NodeType")})
    return out


def _coefficient(doc, groups) -> dict | None:
    """A one-letter caption left of a scheme's first reactant (polymer_x1.0's n): its ink to the
    reactant's fragment box, and its ink center against the arrow axis."""
    if not doc.arrows:
        return None
    axis = doc.arrows[0]["tail"][1]
    for c in doc.captions:
        if c["formula"] or len(c["text"]) > 3:
            continue
        right = [fb for fid, fb in doc.fragments.items() if fb[0] > c["bbox"][2] and fb[1] < axis < fb[3]]
        if not right:
            continue
        fb = min(right, key=lambda b: b[0])
        glyphs = [g for g in doc.glyphs if g.box and abs(g.base - c["p"][1]) < 0.5 and g.ch == c["text"][0]
                  and abs(g.x - c["p"][0]) < 1.0]
        ink = glyphs[0].box if glyphs else c["bbox"]
        return {"text": c["text"], "gap_to_species": round(fb[0] - ink[2], 3),
                "center_above_axis": round(axis - (ink[1] + ink[3]) / 2, 3)}
    return None
