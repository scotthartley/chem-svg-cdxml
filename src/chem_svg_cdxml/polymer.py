"""Polymers (style-spec.md §17; polymer_x1.0): repeat units in brackets, copolymers, chain ends.

A spec's `polymer:` (spec.Polymer) becomes one RDKit mol (assemble): the start group, the units
head to tail and the end group, each unit's first `*` its head (left) and its second its tail. A
missing end group leaves the unit's `*` as a bare chain end (POLY_END "dangling": drawn as a bond
to nowhere, written to CDXML as a plain carbon, as ChemDraw's references do). Each atom carries
its segment (POLY_SEG): 0 the start side, k unit k, n + 1 the end side; the bonds between
segments are the crossing bonds (and a random copolymer's slashed bonds).

orient() lays the backbone out horizontally, head on the left, pendants below (polymer_x1.0); the
backbone runs on through the start and end groups along their longest chain (_backbone).

finish() works on the drawn Molecule (pt), where label ink is known. Rules (§17), from
polymer_x1.0, whose brackets were placed by eye; round values in preset.polymer:
- A bracket crosses its crossing bond at the bond's midpoint (measured 0.44-0.52 of 14.4 pt)
  when that leaves min_gap_x to the enclosed ink. Otherwise the bond is stretched to
  crossing_stretch x bond_length keeping its rise (28.8 pt at 15 deg in the reference) and the
  bracket goes as far beyond the ink as the pair's other bracket at its midpoint (within
  min_gap_x..gap_x; else gap_x), the bond reaching at least half its old run past it.
- A random copolymer's bonds between units are crossing_stretch x long too, with a slash across
  their middle (slash_length, slash_angle clockwise from the bond); one bracket pair round all
  units, `ratio` centered below.
- A block copolymer has a bracket pair per unit; on a junction bond the two brackets are
  block_gap apart (more if the left subscript needs it), each placed from its own block. A
  `bracket: false` unit between blocks has none: the brackets beside it take the end-group rule.
- Both brackets of the polymer reach from the highest enclosed ink (or crossing point) less
  gap_above to the lowest plus gap_below; the lip is lip_per_height x the bar (ChemDraw's).
- The subscript (n, m, 10, n-1: letters italic) sits label_dx right of each ']', its baseline
  label_drop below the bracket's bottom.

stretch_bond() is the general per-bond length override (spec bond_lengths).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

from rdkit import Chem

from . import chem, richtext
from .chem import POLY_END, POLY_SEG  # atom props: segment; "dangling" / "wavy" chain end
from .geometry import Point, add, length, mul, sub, union_box, unit
from .preset import Preset
from .scene import Group, Path, Polyline, TextRun

KAPPA = 0.5523  # cubic handle of a quarter ellipse
CHAIN_STRETCH = 2.0  # chain_layout's stretched last bond: polymer.crossing_stretch (polymer_x1.0, 2.000; every preset)


# --------------------------------------------------------------------------- assembly


def _piece(smiles: str, n_stars: int, what: str) -> tuple[Chem.Mol, list[int]]:
    """A start group, unit or end group, and its connection `*` atoms (unmapped) in SMILES order."""
    try:
        mol = chem.parse_smiles(smiles)
    except ValueError as e:
        raise ValueError(f"polymer {what}: {e}") from None
    stars = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and not a.GetAtomMapNum()]
    if len(stars) != n_stars:
        need = "two '*' (head, then tail)" if n_stars == 2 else "one '*' (the bond to the chain)"
        raise ValueError(f"polymer {what}: {smiles!r} has {len(stars)} '*'; it needs {need}")
    for s in stars:
        a = mol.GetAtomWithIdx(s)
        if a.GetDegree() != 1 or a.GetBonds()[0].GetBondType() != Chem.BondType.SINGLE:
            raise ValueError(f"polymer {what}: {smiles!r}: each '*' needs exactly one single bond")
    return mol, stars


def assemble(p) -> Chem.Mol:
    """The polymer (spec.Polymer) as one Kekulé mol with POLY_SEG on every atom and POLY_END on
    bare chain ends. Atom maps must be unique across the pieces."""
    pieces = []  # (mol, head star, tail star, seg)
    if p.start:
        m, (s,) = _piece(p.start, 1, "start")
        pieces.append((m, None, s, 0))
    for k, u in enumerate(p.units, 1):
        m, (h, t) = _piece(u.smiles, 2, f"unit {k}")
        pieces.append((m, h, t, k))
    if p.end:
        m, (s,) = _piece(p.end, 1, "end")
        pieces.append((m, s, None, len(p.units) + 1))
    seen: dict[int, str] = {}
    for m, _, _, seg in pieces:
        for a in m.GetAtoms():
            if a.GetAtomMapNum():
                if a.GetAtomMapNum() in seen:
                    raise ValueError(f"polymer: map number {a.GetAtomMapNum()} is used twice; "
                                     "atom maps must be unique across the start, units and end")
                seen[a.GetAtomMapNum()] = str(seg)
    combo, offs = Chem.Mol(), []
    for m, _, _, seg in pieces:
        for a in m.GetAtoms():
            a.SetIntProp(POLY_SEG, seg)
        offs.append(combo.GetNumAtoms())
        combo = Chem.CombineMols(combo, m) if combo.GetNumAtoms() else Chem.Mol(m)
    rw = Chem.RWMol(combo)
    drop = []
    for (m1, _, t1, _), (m2, h2, _, _), o1, o2 in zip(pieces, pieces[1:], offs, offs[1:]):
        a = rw.GetAtomWithIdx(o1 + t1).GetNeighbors()[0].GetIdx()
        b = rw.GetAtomWithIdx(o2 + h2).GetNeighbors()[0].GetIdx()
        rw.AddBond(a, b, Chem.BondType.SINGLE)
        drop += [o1 + t1, o2 + h2]
    first, last = pieces[0], pieces[-1]
    for (m, h, t, seg), o, star, side in ((first, offs[0], first[1], 0), (last, offs[-1], last[2], len(p.units) + 1)):
        if star is not None:  # no end group on this side: the unit's '*' is a bare chain end
            a = rw.GetAtomWithIdx(o + star)
            a.SetProp(POLY_END, "dangling")
            a.SetIntProp(POLY_SEG, side)
    for i in sorted(drop, reverse=True):
        rw.RemoveAtom(i)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol, chem._SANITIZE)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol


def segs(mol: Chem.Mol) -> list[int]:
    return [a.GetIntProp(POLY_SEG) if a.HasProp(POLY_SEG) else -1 for a in mol.GetAtoms()]


def is_polymer(mol: Chem.Mol) -> bool:
    return mol.GetNumAtoms() > 0 and mol.GetAtomWithIdx(0).HasProp(POLY_SEG)


def _junction(mol: Chem.Mol, seg: list[int], s: int) -> tuple[int, int]:
    """(atom in segment s, atom in segment s + 1) of the bond joining them."""
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if {seg[i], seg[j]} == {s, s + 1}:
            return (i, j) if seg[i] == s else (j, i)
    raise ValueError(f"polymer: no bond joins segment {s} to {s + 1}")


def _chainy(mol: Chem.Mol, i: int) -> bool:
    """An atom a drawn chain may run through: not H, not in a ring, not an abbreviation label."""
    a = mol.GetAtomWithIdx(i)
    return a.GetAtomicNum() != 1 and not mol.GetRingInfo().NumAtomRings(i) and not a.HasProp(chem.ABBREV_PROP)


def _longest(mol: Chem.Mol, u: int, prev: int, keep=lambda i: True) -> list[int]:
    """The longest path from u (reached from prev) through chainy atoms that `keep` allows (a
    tree): ties to the larger branch, then to one not leaving through a carbonyl carbon (an ester
    stays a pendant, as polymer_x1.0's PMMA)."""
    def rank(path):
        if len(path) == 1:
            return (1, -1, False)
        w = path[1]
        carbonyl = any(b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(mol.GetAtomWithIdx(w)).GetAtomicNum() == 8
                       for b in mol.GetAtomWithIdx(w).GetBonds())
        return (len(path), len(_branch(mol, w, u)), not carbonyl)

    best = [u]
    for nb in mol.GetAtomWithIdx(u).GetNeighbors():
        w = nb.GetIdx()
        if w != prev and _chainy(mol, w) and keep(w):
            cand = [u] + _longest(mol, w, u, keep)
            if rank(cand) > rank(best):
                best = cand
    return best


def _backbone(mol: Chem.Mol) -> list[int]:
    """Atoms from the head-side chain end to the tail-side one: junction to junction, then on into
    the start and end groups along their longest chain (_longest), up to an sp atom (style-spec.md §17)."""
    seg = segs(mol)
    n = max(seg)
    left, _ = _junction(mol, seg, 0)
    _, right = _junction(mol, seg, n - 1)
    path = list(Chem.GetShortestPath(mol, left, right))
    on = set(path)

    def into(end: int, side: int) -> list[int]:  # path[0] or path[-1] onward, excluding it
        if not _chainy(mol, end):
            return []
        ext = _longest(mol, end, -1, lambda w: seg[w] == side and w not in on)[1:]
        cut = next((k for k, i in enumerate(ext) if _linear(mol.GetAtomWithIdx(i))), len(ext))
        return ext[:cut]

    return into(left, 0)[::-1] + path + into(right, n)


def layout(mol: Chem.Mol) -> None:
    """2D coordinates: an acyclic backbone pinned to a horizontal zigzag at 120 deg (the phase that
    puts more pendant-bearing atoms on its lower vertices), the rest laid out by CoordGen round
    it; a backbone through a ring (poly(p-phenylene)) laid out freely. Then orient()."""
    from rdkit import Geometry
    from rdkit.Chem import rdDepictor

    path = _backbone(mol)
    ring = mol.GetRingInfo()
    if any(ring.NumAtomRings(i) for i in path):
        chem.layout_2d(mol)
        orient(mol)
        return
    on = set(path)
    pend = [k for k, i in enumerate(path) if any(nb.GetIdx() not in on and nb.GetAtomicNum() != 1
                                                 for nb in mol.GetAtomWithIdx(i).GetNeighbors())]
    phase = 0 if sum(k % 2 == 0 for k in pend) >= sum(k % 2 == 1 for k in pend) else 1
    b = chem.STANDARD_BOND
    cmap = {i: Geometry.Point2D(k * b * math.cos(math.radians(30)), -b / 2 if (k + phase) % 2 == 0 else 0.0)
            for k, i in enumerate(path)}  # y up: the phase's vertices low
    m = chem._dative_as_single(mol)
    prefer = rdDepictor.GetPreferCoordGen()
    rdDepictor.SetPreferCoordGen(False)  # RDKit's depictor holds a coordMap exactly, at 1.5 A bonds
    try:
        rdDepictor.Compute2DCoords(m, coordMap=cmap)
    finally:
        rdDepictor.SetPreferCoordGen(prefer)
    chem._take_conformer(mol, m)
    chem._chicken_feet(mol, fixed=frozenset(cmap))
    if not _tree_layout(mol, path, phase):
        _upright_pairs(mol, path, phase)
    orient(mol, level=True)


def _linear(atom) -> bool:
    """An sp atom (a triple bond, or two double bonds): its bonds run straight through it."""
    kinds = [b.GetBondType() for b in atom.GetBonds()]
    return Chem.BondType.TRIPLE in kinds or kinds.count(Chem.BondType.DOUBLE) >= 2


def _ring_systems(mol: Chem.Mol) -> list[set[int]]:
    systems: list[set[int]] = []
    for r in mol.GetRingInfo().AtomRings():
        s = set(r)
        for o in [o for o in systems if o & s]:
            s |= o
            systems.remove(o)
        systems.append(s)
    return systems


def _tree_layout(mol: Chem.Mol, path: list[int], phase: int, end_pendant: bool = False) -> bool:
    """Everything off the pinned backbone at 120 deg (style-spec.md §17): a backbone atom's one
    pendant points straight away from the chain; of two, the larger points away and the smaller
    straight the other way, between the backbone bonds at 60 deg to each (polymer_x1.0's PMMA).
    Off the chain, each branch zigzags (the larger branch continues the zigzag; at a carbonyl the
    chain heads toward the tail, as polymer_x1.0's esters), a fourth bond sits between the parent
    bond and the main one, sp atoms run straight, rings keep the depictor's shape, turned so the
    bond to them bisects their outside angle. A start or end group continues the backbone's
    zigzag (its two other groups an upright pair); with `end_pendant`, a last atom whose one other
    group is a ring or a label has it straight away from the chain instead, as a pendant (polymer_x1.0's
    ~~CH2-CH(Ph)•). Returns False (coordinates untouched) when the
    backbone has an sp atom inside, an atom has more than four bonds, or two atoms would clash."""
    conf = mol.GetConformer()
    old = {i: (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y) for i in range(mol.GetNumAtoms())}
    b = chem.STANDARD_BOND
    if any(_linear(mol.GetAtomWithIdx(i)) for i in path[1:-1]):
        return False
    pos = {i: old[i] for i in path}  # y up, the pinned zigzag
    systems = {i: s for s in _ring_systems(mol) for i in s}

    def ang(p, q):
        return math.degrees(math.atan2(q[1] - p[1], q[0] - p[0]))

    def at(p, a):
        return (p[0] + b * math.cos(math.radians(a)), p[1] + b * math.sin(math.radians(a)))

    def size(start, block):
        return len(_branch(mol, start, block))

    def out_of_ring(r, sys):
        nbs = [n.GetIdx() for n in mol.GetAtomWithIdx(r).GetNeighbors() if n.GetIdx() in sys]
        pts = [pos[n] for n in nbs] if len(nbs) == 2 else [pos[n] for n in sys]
        cx, cy = sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)
        return math.degrees(math.atan2(pos[r][1] - cy, pos[r][0] - cx))

    def place_ring(c, target, a_in):
        sys = systems[c]
        nbs = [n.GetIdx() for n in mol.GetAtomWithIdx(c).GetNeighbors() if n.GetIdx() in sys]
        pts = [old[n] for n in nbs] if len(nbs) == 2 else [old[n] for n in sys]
        cx, cy = sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)
        have = math.atan2(old[c][1] - cy, old[c][0] - cx)
        r = math.radians(a_in + 180) - have
        cr, sr = math.cos(r), math.sin(r)
        for x in sys:
            dx, dy = old[x][0] - old[c][0], old[x][1] - old[c][1]
            pos[x] = (target[0] + dx * cr - dy * sr, target[1] + dx * sr + dy * cr)
        for x in sys:
            for nb in mol.GetAtomWithIdx(x).GetNeighbors():
                if nb.GetIdx() not in pos:
                    grow(x, [nb.GetIdx()], out_of_ring(x, sys), 0)

    def grow(u, kids, a_in, turn):
        """Place kids of u (u placed, reached along a_in), main branch first."""
        if not kids:
            return
        if _linear(mol.GetAtomWithIdx(u)) and len(kids) == 1:
            rel = [0]
        else:
            if len(kids) > 3:
                raise _Clash
            kids = sorted(kids, key=lambda k: -size(k, u))
            t = -turn
            if not t:  # toward the tail (+x) if that is a side; else counterclockwise
                t = 1 if math.cos(math.radians(a_in + 60)) >= math.cos(math.radians(a_in - 60)) - 1e-9 else -1
            rel = [60 * t, -60 * t, 120 * t][:len(kids)]
        for k, r in zip(kids, rel):
            a = a_in + r
            p = at(pos[u], a)
            if k in systems:
                place_ring(k, p, a)
                continue
            pos[k] = p
            nxt = [n.GetIdx() for n in mol.GetAtomWithIdx(k).GetNeighbors() if n.GetIdx() not in pos]
            grow(k, nxt, a, (r > 0) - (r < 0))

    def free(i):
        return [n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() not in pos]

    try:
        last = len(path) - 1
        for k, i in enumerate(path):
            nbs = free(i)
            if not nbs:
                continue
            if 0 < k < last:
                # the zigzag is horizontal: a low vertex's pendants point straight down (even beside a
                # stretched bond, whose rise is kept), a high one's straight up
                low = pos[i][1] < (pos[path[k - 1]][1] + pos[path[k + 1]][1]) / 2
                inward = 90.0 if low else -90.0
                if len(nbs) > 2:
                    raise _Clash
                nbs = sorted(nbs, key=lambda j: -size(j, i))
                for j, a in zip(nbs, (inward + 180, inward)):
                    p = at(pos[i], a)
                    if j in systems:
                        place_ring(j, p, a)
                    else:
                        pos[j] = p
                        grow(j, free(j), a, 0)
            else:  # a start or end group: continue the zigzag
                nb = path[1] if k == 0 else path[-2]
                if end_pendant and len(nbs) == 1 and _pendant_end(mol, nbs[0]):
                    a = -90.0 if pos[i][1] < pos[nb][1] else 90.0  # a low vertex's pendants point down
                    p = at(pos[i], a)
                    if nbs[0] in systems:
                        place_ring(nbs[0], p, a)
                    else:
                        pos[nbs[0]] = p
                        grow(nbs[0], free(nbs[0]), a, 0)
                    continue
                a_in = ang(pos[nb], pos[i])
                turn = 0
                if len(path) > 2:
                    nb2 = path[2] if k == 0 else path[-3]
                    d = (a_in - ang(pos[nb2], pos[nb]) + 180) % 360 - 180
                    turn = (d > 0) - (d < 0)
                grow(i, nbs, a_in, turn)
    except _Clash:
        return False
    if len(pos) != mol.GetNumAtoms():
        return False
    bonded = {frozenset((x.GetBeginAtomIdx(), x.GetEndAtomIdx())) for x in mol.GetBonds()}
    seg = segs(mol)
    n = max(seg)

    def apart(x, y):  # two units of a copolymer: finish() stretches or shifts the bond between them
        return seg[x] != seg[y] and 1 <= seg[x] < n and 1 <= seg[y] < n

    ids = list(pos)
    for x in range(len(ids)):
        for y in range(x + 1, len(ids)):
            if frozenset((ids[x], ids[y])) not in bonded and not apart(ids[x], ids[y]) \
                    and math.dist(pos[ids[x]], pos[ids[y]]) < 0.6 * b:
                return False
    for i, (x, y) in pos.items():
        conf.SetAtomPosition(i, (x, y, 0.0))
    return True


class _Clash(Exception):
    pass


def _pendant_end(mol: Chem.Mol, i: int) -> bool:
    a = mol.GetAtomWithIdx(i)
    return mol.GetRingInfo().NumAtomRings(i) > 0 or a.HasProp(chem.ABBREV_PROP)


def chain_layout(mol: Chem.Mol) -> bool:
    """A species with one wavy chain end (spec chain_ends; style-spec.md §17) drawn as polymer_x1.0
    draws a growing chain: the chain end on the left, the longest acyclic chain from it (not
    through rings or abbreviation labels) on a horizontal 120 deg zigzag, pendants below, and a
    ring or label at its far end hanging straight down as a pendant; if it would meet the pendant
    before it, the bond leaving that pendant's atom is CHAIN_STRETCH x long with its rise kept
    (failing that, the zigzag continues); the rest as _tree_layout.
    Returns False (layout untouched) for any other species or when the drawing would clash."""
    from rdkit import Geometry
    from rdkit.Chem import rdDepictor

    ends = [a.GetIdx() for a in mol.GetAtoms() if a.HasProp(POLY_END) and a.GetProp(POLY_END) == "wavy"]
    if len(ends) != 1 or is_polymer(mol):
        return False
    path = _longest(mol, ends[0], -1)
    if len(path) < 3 or any(_linear(mol.GetAtomWithIdx(i)) for i in path[1:-1]):
        return False
    on = set(path)
    pend = [k for k, i in enumerate(path) if any(nb.GetIdx() not in on and nb.GetAtomicNum() != 1
                                                 for nb in mol.GetAtomWithIdx(i).GetNeighbors())]
    phase = 0 if sum(k % 2 == 0 for k in pend) >= sum(k % 2 == 1 for k in pend) else 1
    b = chem.STANDARD_BOND
    keep = Chem.Mol(mol)

    last = len(path) - 1
    before = [k for k in pend if k < last]
    cut = before[-1] + 1 if before else last  # the bond leaving the pendant before the far end's

    def pinned(stretch: float) -> None:
        """The zigzag, the bond into path[cut] `stretch` x long with its rise kept (as a stretched
        crossing bond). The depictor lays out the unstretched zigzag (it bends rings round a long
        pinned bond); the stretch moves the pinned atoms only, so the rings keep their shape."""
        run = b * math.cos(math.radians(30))
        cmap = {i: Geometry.Point2D(k * run, -b / 2 if (k + phase) % 2 == 0 else 0.0) for k, i in enumerate(path)}
        m = chem._dative_as_single(mol)
        prefer = rdDepictor.GetPreferCoordGen()
        rdDepictor.SetPreferCoordGen(False)
        try:
            rdDepictor.Compute2DCoords(m, coordMap=cmap)
        finally:
            rdDepictor.SetPreferCoordGen(prefer)
        chem._take_conformer(mol, m)
        chem._chicken_feet(mol, fixed=frozenset(cmap))
        dx = math.sqrt((stretch * b) ** 2 - (b / 2) ** 2) - run
        conf = mol.GetConformer()
        for k in range(cut, len(path)):
            q = conf.GetAtomPosition(path[k])
            conf.SetAtomPosition(path[k], (q.x + dx, q.y, 0.0))

    pinned(1.0)
    if _tree_layout(mol, path, phase, end_pendant=True):
        return True
    pinned(CHAIN_STRETCH)  # the far end's pendant would meet the one before it: stretch the bond leaving that one
    if _tree_layout(mol, path, phase, end_pendant=True):
        return True
    pinned(1.0)
    if _tree_layout(mol, path, phase):  # still clashing: continue the zigzag instead
        return True
    if keep.GetNumConformers():
        chem._take_conformer(mol, keep)
    return False


def _branch(mol: Chem.Mol, start: int, block: int) -> list[int]:
    seen, stack = {start}, [start]
    while stack:
        u = stack.pop()
        for nb in mol.GetAtomWithIdx(u).GetNeighbors():
            w = nb.GetIdx()
            if w != block and w not in seen:
                seen.add(w)
                stack.append(w)
    return sorted(seen)


def _upright_pairs(mol: Chem.Mol, path: list[int], phase: int) -> None:
    """A backbone atom with two pendant branches (PMMA's C(CH3)(CO2Me)) has them upright, as
    polymer_x1.0 draws them: the larger away from the chain (below a low vertex), the smaller
    straight the other way, between the backbone bonds. Each branch turns rigidly."""
    pts = chem._coords(mol)  # y down
    on = set(path)
    for k, i in enumerate(path):
        nbs = [nb.GetIdx() for nb in mol.GetAtomWithIdx(i).GetNeighbors()
               if nb.GetIdx() not in on and nb.GetAtomicNum() != 1]
        if len(nbs) != 2:
            continue
        low = (k + phase) % 2 == 0
        branches = sorted((_branch(mol, j, i) for j in nbs), key=len)
        c = pts[i]
        for br, down in ((branches[1], low), (branches[0], not low)):
            want = math.atan2(1.0 if down else -1.0, 0.0)
            j = next(x for x in br if x in nbs)
            have = math.atan2(pts[j][1] - c[1], pts[j][0] - c[0])
            r = want - have
            cr, sr = math.cos(r), math.sin(r)
            for x in br:
                dx, dy = pts[x][0] - c[0], pts[x][1] - c[1]
                pts[x] = (c[0] + dx * cr - dy * sr, c[1] + dx * sr + dy * cr)
    chem._set_coords(mol, pts)


def orient(mol: Chem.Mol, level: bool = False) -> None:
    """Turn the laid-out conformer so the backbone (head chain end to tail chain end) runs left to
    right along x (its principal axis; `level`: it already does, a pinned zigzag), and mirror it
    top to bottom if most pendant atoms sit above it."""
    seg = segs(mol)
    n = max(seg)
    path = _backbone(mol)
    pts = chem._coords(mol)
    pp = [pts[i] for i in path]
    cx, cy = sum(p[0] for p in pp) / len(pp), sum(p[1] for p in pp) / len(pp)
    sxx = sum((p[0] - cx) ** 2 for p in pp)
    syy = sum((p[1] - cy) ** 2 for p in pp)
    sxy = sum((p[0] - cx) * (p[1] - cy) for p in pp)
    ang = 0.0 if level else 0.5 * math.atan2(2 * sxy, sxx - syy)  # principal axis, y down
    d = (math.cos(ang), math.sin(ang))
    if (pp[-1][0] - pp[0][0]) * d[0] + (pp[-1][1] - pp[0][1]) * d[1] < 0:
        d = (-d[0], -d[1])
    c, s = d
    out = [((x - cx) * c + (y - cy) * s, -(x - cx) * s + (y - cy) * c) for x, y in pts]  # d -> +x
    on_path = set(path)
    pend = [out[i][1] for i in range(len(out)) if i not in on_path and 1 <= seg[i] <= n - 1]
    if pend and sum(pend) / len(pend) < 0:  # pendants above the backbone: mirror (wedges come later)
        out = [(x, -y) for x, y in out]
    chem._set_coords(mol, out)


# --------------------------------------------------------------------------- drawn geometry


@dataclass
class BracketPair:
    """One pair of repeat-unit brackets, in molecule coordinates: the bars' stroke centers x, the
    lips' top and bottom, the subscript, the enclosed atoms and the crossing bonds (bond idx,
    inner atom idx), left then right."""
    x: tuple[float, float]
    top: float
    bottom: float
    label: str
    atoms: list[int]
    crossing: list[tuple[int, int]]
    label_origin: Point = (0.0, 0.0)

    def shifted(self, d: Point) -> "BracketPair":
        return replace(self, x=(self.x[0] + d[0], self.x[1] + d[0]), top=self.top + d[1], bottom=self.bottom + d[1],
                       label_origin=add(self.label_origin, d))


@dataclass
class PolymerMarks:
    """What a polymer draws besides its atoms and bonds (Molecule.polymer), in molecule coordinates."""
    kind: str
    pairs: list[BracketPair]
    slashes: list[tuple[int, Point, Point]] = field(default_factory=list)  # (bond idx, end, end)
    ratio: str | None = None
    ratio_origin: Point = (0.0, 0.0)  # its first glyph origin and baseline


def label_runs(text: str) -> list[richtext.Run]:
    """A subscript as ChemDraw sets it (polymer_x1.0): letters italic, digits and signs regular;
    an ASCII hyphen becomes a minus (n-1 -> n−1)."""
    text = text.replace("-", "−")
    runs: list[richtext.Run] = []
    for ch in text:
        style = "italic" if ch.isalpha() else "regular"
        if runs and runs[-1].style == style:
            runs[-1] = richtext.Run(runs[-1].text + ch, style)
        else:
            runs.append(richtext.Run(ch, style))
    return runs


def _shift_from(m: chem.Molecule, seg: dict[int, int], s: int, dx: float) -> None:
    for a in m.atoms:
        if seg[a.idx] >= s:
            a.x += dx


def _stretch_level(m: chem.Molecule, seg: dict[int, int], i: int, j: int, new_len: float) -> None:
    """Make bond i-j (j right of i, in the next segment) new_len long, keeping its rise: the
    segments from j's on move along x."""
    a, b = m.atoms[i], m.atoms[j]
    rise = b.y - a.y
    run = math.sqrt(max(new_len ** 2 - rise ** 2, 0.0))
    _shift_from(m, seg, seg[j], a.x + run - b.x)


def _ink(m: chem.Molecule, preset: Preset, atoms: set[int]) -> tuple:
    """Ink of the atoms' labels and of the bonds between them (their drawn lines and polygons)."""
    from .draw_mol import Drawer

    d = Drawer(m, preset)
    hw = preset.lengths.line_width / 2
    boxes = []
    for b in m.bonds:
        if b.a not in atoms or b.b not in atoms:
            continue
        g = d.bond_geometry(b)
        for i in (b.a, b.b):  # the main line (drawn as chains), to the atom or its label
            if i not in d.labels:
                x, y = m.atoms[i].p
                boxes.append((x - hw, y - hw, x + hw, y + hw))
        for v, u0, u1 in g.lines:
            p, q = g.frame.at(u0, v), g.frame.at(u1, v)
            boxes.append((min(p[0], q[0]) - hw, min(p[1], q[1]) - hw, max(p[0], q[0]) + hw, max(p[1], q[1]) + hw))
        if g.polygon:
            xs, ys = [p[0] for p in g.polygon], [p[1] for p in g.polygon]
            boxes.append((min(xs), min(ys), max(xs), max(ys)))
    for i in atoms:
        if i in d.labels:
            boxes += [r.bbox() for r in d.labels[i].runs if r.bbox()]
    if not boxes:  # a one-atom unit drawn bare (*C*): its point
        boxes = [(m.atoms[i].x - hw, m.atoms[i].y - hw, m.atoms[i].x + hw, m.atoms[i].y + hw) for i in atoms]
    return union_box(boxes)


def _bond_idx(m: chem.Molecule, i: int, j: int) -> int:
    return next(b.idx for b in m.bonds if {b.a, b.b} == {i, j})


def label_advance(text: str, preset: Preset) -> float:
    return richtext.layout_line(label_runs(text), 0.0, 0.0, preset.polymer.label_size, preset).advance


def finish(m: chem.Molecule, p, preset: Preset) -> chem.Molecule:
    """Stretch the crossing bonds and place the brackets, subscripts, slashes and ratio (module
    docstring) on a drawn polymer Molecule; returns it with Molecule.polymer set."""
    pp = preset.polymer
    L = preset.lengths.bond_length
    seg = {a.idx: a.seg for a in m.atoms}
    n = len(p.units)

    def junction(s: int) -> tuple[int, int]:
        for b in m.bonds:
            if {seg[b.a], seg[b.b]} == {s, s + 1}:
                return (b.a, b.b) if seg[b.a] == s else (b.b, b.a)
        raise ValueError(f"polymer: no bond joins segment {s} to {s + 1}")

    def along(i: int, j: int, dist: float) -> float:
        """x of the point dist from atom i toward atom j."""
        a, b = m.atoms[i].p, m.atoms[j].p
        return add(a, mul(unit(sub(b, a)), dist))[0]

    if p.kind != "block":
        groups = [list(range(1, n + 1))]
    else:  # a pair per bracketed unit; a bracket: false unit sits between two pairs
        groups = [[k] for k in range(1, n + 1) if p.units[k - 1].bracket]
    bracketed = {s for S in groups for s in S}
    labels = p.labels()
    if p.kind == "random":
        for k in range(1, n):
            _stretch_level(m, seg, *junction(k), pp.crossing_stretch * L)
    xs: list[tuple[float, float]] = []
    for g, S in enumerate(groups):
        lo, hi = S[0], S[-1]
        inside = {i for i, s in seg.items() if lo <= s <= hi}
        box = _ink(m, preset, inside)

        def mid_gap(i: int, j: int, own_len: bool) -> float:
            """Ink to a bracket at the midpoint of bond i (inside) - j, unstretched."""
            half = (length(sub(m.atoms[j].p, m.atoms[i].p)) if own_len else L) / 2
            x = along(i, j, half)
            return box[0] - x if x < box[0] else x - box[2]

        def beside_stretch(other: float) -> float:
            """A bracket beside a stretched bond: as far from the ink as the pair's other bracket at
            its bond's midpoint, within [min_gap_x, gap_x]; else gap_x."""
            return other if pp.min_gap_x <= other <= pp.gap_x else pp.gap_x

        left_mid = mid_gap(*junction(lo - 1)[::-1], lo - 1 not in bracketed)
        right_mid = mid_gap(*junction(hi), hi + 1 not in bracketed)
        # left bracket: on the bond from segment lo - 1
        outer, inner = junction(lo - 1)
        if lo - 1 not in bracketed:  # the start side or a bracketless segment
            run0 = abs(m.atoms[inner].x - m.atoms[outer].x)
            xl = along(inner, outer, length(sub(m.atoms[outer].p, m.atoms[inner].p)) / 2)
            if box[0] - xl < pp.min_gap_x:
                before = m.atoms[inner].x
                _stretch_level(m, seg, outer, inner, pp.crossing_stretch * L)
                dx = m.atoms[inner].x - before
                xl = box[0] + dx - beside_stretch(right_mid)
                short = m.atoms[outer].x - (xl - run0 / 2)
                if short > 0:  # the bond must still reach half its old run past the bracket
                    _shift_from(m, seg, lo, short)
                    xl += short
            if lo - 1 > 0 or p.start:  # the drawn atoms outside keep min_gap_x from the bracket too
                ob = _ink(m, preset, {i for i, s in seg.items() if s == lo - 1})
                need = ob[2] + pp.min_gap_x - xl
                if need > 0:
                    _shift_from(m, seg, lo, need)
                    xl += need
        else:  # a block junction: this block's own place for '[', then the block moved block_gap from ']'
            own = along(inner, outer, L / 2)
            if box[0] - own < pp.min_gap_x:
                own = box[0] - beside_stretch(right_mid)
            gap = max(pp.block_gap, pp.label_dx + label_advance(labels[g - 1], preset) + pp.gap_x)
            dx = xs[-1][1] + gap - own
            _shift_from(m, seg, lo, dx)
            xl = own + dx
        # right bracket: on the bond to segment hi + 1
        inner, outer = junction(hi)
        box = _ink(m, preset, inside)
        run0 = abs(m.atoms[outer].x - m.atoms[inner].x)
        loose = hi + 1 not in bracketed  # the end side or a bracketless segment
        xr = along(inner, outer, (length(sub(m.atoms[outer].p, m.atoms[inner].p)) if loose else L) / 2)
        if xr - box[2] < pp.min_gap_x:
            xr = box[2] + beside_stretch(left_mid)
            if loose:
                _stretch_level(m, seg, inner, outer, pp.crossing_stretch * L)
                short = (xr + run0 / 2) - m.atoms[outer].x
                if short > 0:
                    _shift_from(m, seg, hi + 1, short)
        if loose and (hi + 1 <= n or p.end):  # the drawn atoms outside keep min_gap_x from it
            ob = _ink(m, preset, {i for i, s in seg.items() if s == hi + 1})
            need = xr + pp.min_gap_x - ob[0]
            if need > 0:
                _shift_from(m, seg, hi + 1, need)
        xs.append((xl, xr))
    # heights: one span for every bracket of the polymer
    ys = []
    pairs = []
    for g, S in enumerate(groups):
        inside = {i for i, s in seg.items() if S[0] <= s <= S[-1]}
        box = _ink(m, preset, inside)
        ys += [box[1], box[3]]
        crossing = []
        for (i, j), x in ((junction(S[0] - 1)[::-1], xs[g][0]), (junction(S[-1]), xs[g][1])):
            a, b = m.atoms[i].p, m.atoms[j].p  # i inside
            t = (x - a[0]) / (b[0] - a[0]) if abs(b[0] - a[0]) > 1e-9 else 0.5
            ys.append(a[1] + t * (b[1] - a[1]))
            crossing.append((_bond_idx(m, i, j), i))
        pairs.append(BracketPair(xs[g], 0.0, 0.0, labels[g], sorted(inside), crossing))
    top, bottom = min(ys) - pp.gap_above, max(ys) + pp.gap_below
    for pr in pairs:
        pr.top, pr.bottom = top, bottom
        pr.label_origin = (pr.x[1] + pp.label_dx, bottom + pp.label_drop)
    marks = PolymerMarks(p.kind, pairs)
    if p.kind == "random":
        for k in range(1, n):
            i, j = junction(k)
            a, b = m.atoms[i].p, m.atoms[j].p
            mid = mul(add(a, b), 0.5)
            u = unit(sub(b, a))
            t = math.radians(pp.slash_angle)
            d = (u[0] * math.cos(t) - u[1] * math.sin(t), u[0] * math.sin(t) + u[1] * math.cos(t))
            h = mul(d, pp.slash_length / 2)
            marks.slashes.append((_bond_idx(m, i, j), sub(mid, h), add(mid, h)))
        if p.ratio:
            ln = richtext.layout_line(str(p.ratio), 0.0, 0.0, preset.text.caption_size, preset)
            ink = ln.ink()
            cx = (pairs[0].x[0] + pairs[0].x[1]) / 2
            marks.ratio = str(p.ratio)
            marks.ratio_origin = (cx - (ink[0] + ink[2]) / 2, bottom + pp.ratio_drop)
    m.polymer = marks
    return m


# --------------------------------------------------------------------------- drawing


def bracket_lip(top: float, bottom: float, preset: Preset) -> float:
    return preset.polymer.lip_per_height * (bottom - top)


def subscript_line(pr: BracketPair, preset: Preset, source: str = "") -> richtext.Line:
    return richtext.layout_line(label_runs(pr.label), pr.label_origin[0], pr.label_origin[1],
                                preset.polymer.label_size, preset, source)


def ratio_line(marks: PolymerMarks, preset: Preset, source: str = "") -> richtext.Line:
    return richtext.layout_line(marks.ratio, marks.ratio_origin[0], marks.ratio_origin[1],
                                preset.text.caption_size, preset, source)


def draw_marks(m: chem.Molecule, preset: Preset, name: str) -> Group:
    """`species-<name>-polymer`: the brackets (sources polymer-bracket:<crossing bond ids>), the
    subscripts, the slashes (polymer-slash:<bond id>) and the ratio."""
    from .brackets import square_bracket

    marks = m.polymer
    out = Group(f"species-{name}-polymer")
    lw = preset.lengths.line_width
    for k, pr in enumerate(marks.pairs, 1):
        src = "polymer-bracket:" + ",".join(m.bonds[bi].id for bi, _ in pr.crossing)
        lip = bracket_lip(pr.top, pr.bottom, preset)
        for x, side in zip(pr.x, ("left", "right")):
            br = square_bracket(x, pr.top, pr.bottom, side, preset, f"species-{name}-bracket-{k}-{side}", lip=lip)
            for leaf in br.children:
                leaf.source = src
            out.children.append(br)
        out.children.append(Group(f"species-{name}-subscript-{k}", subscript_line(pr, preset, f"polymer-label:{k}").runs))
    for bi, a, b in marks.slashes:
        out.children.append(Polyline([a, b], lw, source=f"polymer-slash:{m.bonds[bi].id}"))
    if marks.ratio:
        out.children.append(Group(f"species-{name}-ratio", ratio_line(marks, preset, "polymer-ratio").runs))
    return out


def wave(center: Point, along_dir: Point, preset: Preset, source: str = "") -> Path:
    """A wavy chain end (polymer_x1.0, ChemDraw's ExternalConnectionType Wavy): half-ellipses
    wave_half_period long and wave_amplitude high, alternating sides, across the bond end,
    wave_length long and centered on it. `along_dir`: the squiggle's direction."""
    pp = preset.polymer
    u = unit(along_dir)
    nv = (-u[1], u[0])
    q = max(2, round(2 * pp.wave_length / pp.wave_half_period))  # quarter arcs
    hp, amp = pp.wave_half_period, pp.wave_amplitude
    total = q * hp / 2

    def at(s: float, v: float) -> Point:
        return add(center, add(mul(u, s - total / 2), mul(nv, v)))
    segs: list[tuple] = [("M", at(0.0, 0.0))]
    s, side = 0.0, 1
    for k in range(q):
        if k % 2 == 0:  # axis -> peak
            p0, p1 = (s, 0.0), (s + hp / 2, side * amp)
            c1, c2 = (s, side * amp * KAPPA), (s + hp / 2 * (1 - KAPPA), side * amp)
        else:  # peak -> axis
            p0, p1 = (s, side * amp), (s + hp / 2, 0.0)
            c1, c2 = (s + hp / 2 * KAPPA, side * amp), (s + hp / 2, side * amp * KAPPA)
            side = -side
        segs.append(("C", at(*c1), at(*c2), at(*p1)))
        s = p1[0]
    return Path(segs, fill=False, stroke=preset.lengths.line_width, closed=False, source=source)


def wavy_ends(m: chem.Molecule, preset: Preset) -> list[Path]:
    out = []
    for a in m.atoms:
        if a.end != "wavy":
            continue
        nbs = m.neighbors(a.idx)
        if len(nbs) != 1:
            raise ValueError(f"chain_ends: atom {a.id} must have exactly one bond")
        d = unit(sub(a.p, m.atoms[nbs[0][0]].p))
        out.append(wave(a.p, (-d[1], d[0]), preset, f"atom:{a.id}"))
    return out


# --------------------------------------------------------------------------- per-bond lengths


def stretch_bond(m: chem.Molecule, i: int, j: int, factor: float, bond_length: float) -> None:
    """Bond i-j made factor x bond_length long along its own direction: the side with fewer
    atoms moves. A ring bond can't be stretched."""
    b = next((bd for bd in m.bonds if {bd.a, bd.b} == {i, j}), None)
    if b is None:
        raise ValueError(f"atoms {m.atoms[i].id} and {m.atoms[j].id} are not bonded")
    if b.ring:
        raise ValueError(f"bond {m.atoms[i].id}-{m.atoms[j].id} is in a ring; only acyclic bonds can be stretched")
    side_j = _side(m, j, i)
    side_i = _side(m, i, j)
    move, fixed, mover = (side_j, i, j) if len(side_j) <= len(side_i) else (side_i, j, i)
    d = sub(m.atoms[mover].p, m.atoms[fixed].p)
    delta = mul(unit(d), factor * bond_length - length(d))
    for k in move:
        m.atoms[k].x += delta[0]
        m.atoms[k].y += delta[1]


def _side(m: chem.Molecule, start: int, block: int) -> set[int]:
    seen, stack = {start}, [start]
    while stack:
        u = stack.pop()
        for w, _ in m.neighbors(u):
            if w != block and w not in seen:
                seen.add(w)
                stack.append(w)
    return seen
