"""Scheme layouts and validations shared by the test files, each built once per session.

The arrow search makes some layouts take seconds (m4_perkin, the nitrations,
tests/specs/fischer_over_x1.0.yaml, tests/specs/phenol_x_abstraction.yaml), and several test
files ask for the same ones. What these return is shared: a test that changes a layout works
on copy.deepcopy(...) of it.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import yaml
from validate import checked_layout

from chem_svg_cdxml.layout import SchemeLayout, build_scheme
from chem_svg_cdxml.spec import Spec

REPO = Path(__file__).resolve().parents[1]
EX = REPO / "examples"


def example(stem: str) -> Path:
    """examples/m<k>/<stem>.yaml."""
    return EX / stem[:2] / f"{stem}.yaml"


@lru_cache(maxsize=None)
def _validated(path: Path) -> tuple[tuple[str, ...], SchemeLayout | None]:
    problems, _, lay = checked_layout(path)
    return tuple(problems), lay


def problems(path: Path) -> list[str]:
    """validate(path): its layout is kept for scheme(path)."""
    return list(_validated(Path(path))[0])


def scheme(path: Path) -> SchemeLayout:
    """build_scheme of the spec at path, as render.py and validate.py build it (base dir: the
    spec's folder; name: its stem): validate's own layout, so asking for both builds it once."""
    lay = _validated(Path(path))[1]
    if lay is None:
        raise ValueError(f"no scheme layout for {path}: {problems(path)}")
    return lay


@lru_cache(maxsize=None)
def overridden(path: Path, arrows: tuple[tuple[tuple[str, object], ...], ...]) -> SchemeLayout:
    """build_scheme of the spec at path with overrides.arrows set to `arrows` (each a sorted
    tuple of an override's items, see arrow_overrides)."""
    path = Path(path)
    raw = yaml.safe_load(path.read_text())
    raw.setdefault("overrides", {})["arrows"] = [dict(a) for a in arrows]
    return build_scheme(Spec.model_validate(raw), path.parent, path.stem)


def arrow_overrides(*overrides: dict) -> tuple:
    """The hashable form of overrides.arrows entries, for overridden()."""
    return tuple(tuple(sorted(o.items())) for o in overrides)
