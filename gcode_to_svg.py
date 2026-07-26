#!/usr/bin/env python3
"""Render a G-code toolpath's spray (and optionally travel) moves as an SVG.

Existing to *validate* :mod:`svg_to_gcode`: converting an SVG's paths to
G-code (``svg_to_gcode.py``) and this file's G-code back to SVG should
reproduce the original artwork, since :func:`prepare_subpaths`'s default
Y-flip is the exact inverse of ``svg_to_gcode.convert``'s (both reflect about
the *spray* geometry's own bounding box, so applying it twice is the
identity). Round-trip test: parse an SVG -> ``svg_to_gcode.convert(...,
normalize=False)`` -> :func:`spraysim.gcode.load_moves` on the emitted G-code
-> :func:`convert` here -> compare the returned points against the original
SVG's parsed path points. Curves are already linearised by the time they
reach G-code, so the round trip is exact for polylines and within
``svg_to_gcode``'s ``--tolerance-mm`` of the true curve otherwise.

Each :class:`spraysim.gcode.Move` is a straight line, so no curve-fitting is
needed here — moves are simply grouped into polylines: consecutive moves
whose endpoints connect become one SVG subpath (``M`` + ``L``\\ s); any gap
(a travel move, or any coordinate discontinuity) starts a new one. This is
computed from move *positions*, not from G-code structure, so it groups
correctly regardless of which tool produced the file.

Usage
-----
    python gcode_to_svg.py path.gcode                     # -> path.svg
    python gcode_to_svg.py path.gcode -o out.svg --show-travel
    python gcode_to_svg.py path.gcode --no-flip-y          # keep raw Y

(run from the repo root; travel moves, if shown, are drawn dashed and grey,
spray moves solid and black.)
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Callable

from spraysim import gcode

Point = tuple[float, float]
Subpath = list[Point]

DEFAULT_MARGIN_MM = 5.0
DEFAULT_STROKE_WIDTH_MM = 0.5
_CONNECT_EPS_MM = 1.0e-6


def moves_to_subpaths(moves: list[gcode.Move],
                       predicate: Callable[[gcode.Move], bool]) -> list[Subpath]:
    """Group the moves matching ``predicate`` into connected polylines (mm).

    A new subpath starts whenever a selected move's start doesn't match the
    previous selected move's end (a travel move in between, or any other
    discontinuity) — geometry-based, not G-code-structure-based.
    """
    segments = [
        ((m.start[0] * 1.0e3, m.start[1] * 1.0e3), (m.end[0] * 1.0e3, m.end[1] * 1.0e3))
        for m in moves if predicate(m)
    ]
    subpaths: list[Subpath] = []
    cur: Subpath | None = None
    prev_end: Point | None = None
    for start, end in segments:
        connected = prev_end is not None and math.hypot(
            start[0] - prev_end[0], start[1] - prev_end[1]
        ) <= _CONNECT_EPS_MM
        if not connected:
            cur = [start]
            subpaths.append(cur)
        cur.append(end)
        prev_end = end
    return subpaths


def prepare_subpaths(
    moves: list[gcode.Move], *, show_travel: bool = False, flip_y: bool = True
) -> tuple[list[Subpath], list[Subpath]]:
    """Return (spray_subpaths, travel_subpaths) in mm, flipped if requested.

    The flip reflects about the *spray* geometry's bounding box only (falling
    back to the travel bbox if there's no spray) — this exactly undoes
    ``svg_to_gcode.convert``'s own flip, which is computed the same way
    before any travel moves exist.
    """
    spray = moves_to_subpaths(moves, lambda m: m.spray_on)
    travel = moves_to_subpaths(moves, lambda m: not m.spray_on) if show_travel else []

    if flip_y:
        basis = [p for sp in spray for p in sp] or [p for sp in travel for p in sp]
        if basis:
            ys = [y for _, y in basis]
            ymin, ymax = min(ys), max(ys)

            def fy(y: float) -> float:
                return ymin + ymax - y

            spray = [[(x, fy(y)) for x, y in sp] for sp in spray]
            travel = [[(x, fy(y)) for x, y in sp] for sp in travel]

    return spray, travel


def _bbox(subpaths: list[Subpath]) -> tuple[float, float, float, float]:
    xs = [x for sp in subpaths for x, _ in sp]
    ys = [y for sp in subpaths for _, y in sp]
    return min(xs), max(xs), min(ys), max(ys)


def render_svg(
    spray: list[Subpath],
    travel: list[Subpath],
    out_path: str | Path,
    *,
    margin_mm: float = DEFAULT_MARGIN_MM,
    stroke_width_mm: float = DEFAULT_STROKE_WIDTH_MM,
    precision: int = 3,
) -> None:
    """Write ``spray``/``travel`` (already flipped/scaled, in mm) as an SVG."""
    all_subpaths = spray + travel
    if not all_subpaths:
        raise ValueError("nothing to draw: no spray or travel segments")

    xmin, xmax, ymin, ymax = _bbox(all_subpaths)
    vb_xmin, vb_ymin = xmin - margin_mm, ymin - margin_mm
    width, height = (xmax - xmin) + 2 * margin_mm, (ymax - ymin) + 2 * margin_mm

    def fmt(v: float) -> str:
        return f"{v:.{precision}f}"

    def path_d(subpaths: list[Subpath]) -> str:
        parts = []
        for sp in subpaths:
            if len(sp) < 2:
                continue
            x0, y0 = sp[0]
            parts.append(f"M{fmt(x0)},{fmt(y0)}")
            parts.extend(f"L{fmt(x)},{fmt(y)}" for x, y in sp[1:])
        return " ".join(parts)

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{fmt(width)}mm" '
        f'height="{fmt(height)}mm" '
        f'viewBox="{fmt(vb_xmin)} {fmt(vb_ymin)} {fmt(width)} {fmt(height)}">'
    ]
    travel_d = path_d(travel)
    if travel_d:
        lines.append(
            f'  <path d="{travel_d}" fill="none" stroke="#999999" '
            f'stroke-width="{fmt(stroke_width_mm / 2)}" '
            f'stroke-dasharray="{fmt(stroke_width_mm * 2)},{fmt(stroke_width_mm * 2)}"/>'
        )
    spray_d = path_d(spray)
    if spray_d:
        lines.append(
            f'  <path d="{spray_d}" fill="none" stroke="#000000" '
            f'stroke-width="{fmt(stroke_width_mm)}" '
            'stroke-linecap="round" stroke-linejoin="round"/>'
        )
    lines.append("</svg>")

    Path(out_path).write_text("\n".join(lines) + "\n")


def convert(
    source: str,
    out_path: str | Path,
    *,
    show_travel: bool = False,
    flip_y: bool = True,
    margin_mm: float = DEFAULT_MARGIN_MM,
    stroke_width_mm: float = DEFAULT_STROKE_WIDTH_MM,
    precision: int = 3,
) -> tuple[list[Subpath], list[Subpath]]:
    """Load G-code ``source`` (file path or inline text), write ``out_path`` as
    an SVG, and return the (spray, travel) subpaths actually written (mm) —
    handy for round-trip comparisons against the original artwork."""
    moves = gcode.load_moves(source)
    spray, travel = prepare_subpaths(moves, show_travel=show_travel, flip_y=flip_y)
    render_svg(spray, travel, out_path, margin_mm=margin_mm,
               stroke_width_mm=stroke_width_mm, precision=precision)
    return spray, travel


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def main() -> None:
    p = argparse.ArgumentParser(
        description="Render a G-code toolpath's spray (and optionally travel) moves as an SVG."
    )
    p.add_argument("gcode", help="G-code file (or inline text containing a newline)")
    p.add_argument("-o", "--out", type=Path, default=None,
                   help="output .svg path (default: <gcode>.svg next to the input; "
                        "required when passing inline G-code text)")
    p.add_argument("--show-travel", action="store_true",
                   help="also draw G0 travel moves, dashed and grey (default: spray only)")
    p.add_argument("--no-flip-y", action="store_true",
                   help="don't flip Y back to SVG's Y-down convention (keep raw "
                        "G-code/machine coordinates)")
    p.add_argument("--margin-mm", type=float, default=DEFAULT_MARGIN_MM,
                   help=f"blank margin (mm) around the drawing (default: {DEFAULT_MARGIN_MM:g})")
    p.add_argument("--stroke-width-mm", type=float, default=DEFAULT_STROKE_WIDTH_MM,
                   help=f"stroke width (mm) of the spray path (default: {DEFAULT_STROKE_WIDTH_MM:g})")
    p.add_argument("--precision", type=int, default=3,
                   help="decimal places for emitted coordinates (default: 3)")
    args = p.parse_args()

    gcode_path = Path(args.gcode)
    is_file = "\n" not in args.gcode and gcode_path.exists()
    out_path = args.out
    if out_path is None:
        if not is_file:
            p.error("--out is required when passing inline G-code text")
        out_path = gcode_path.with_suffix(".svg")

    try:
        spray, travel = convert(
            args.gcode, out_path,
            show_travel=args.show_travel,
            flip_y=not args.no_flip_y,
            margin_mm=args.margin_mm,
            stroke_width_mm=args.stroke_width_mm,
            precision=args.precision,
        )
    except (ValueError, OSError) as exc:
        p.error(str(exc))

    n_points = sum(len(sp) for sp in spray) + sum(len(sp) for sp in travel)
    print(f"{len(spray)} spray subpath(s), {len(travel)} travel subpath(s), "
          f"{n_points} points -> {out_path}")


if __name__ == "__main__":
    main()
