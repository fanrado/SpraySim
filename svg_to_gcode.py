#!/usr/bin/env python3
"""Convert an SVG file's ``<path>`` elements into a G-code spray toolpath.

Parses a pragmatic subset of the SVG path mini-language (``M L H V C S Q T A
Z``, absolute and relative) and linearises everything — including cubic and
quadratic Beziers and elliptical arcs — into straight ``G1`` segments, since
:mod:`spraysim.gcode` (the consumer) rejects ``G2``/``G3`` arcs and expects a
pre-linearised path. Each SVG subpath (one ``M``/``m`` to the next) becomes a
``G0`` travel move to its start followed by ``G1`` moves along it — the same
travel/spray convention as a hand-written raster (see
``examples/raster.gcode``).

Only ``<path>`` elements are read; other drawable shapes (``<rect>``,
``<circle>``, ``<line>``, ...) are ignored with a warning. ``transform``
attributes on a ``<path>`` are **not** applied (also warned) — flatten
transforms in your SVG editor first (e.g. Inkscape's "Edit > Flatten
Transforms" or "Object to Path").

Units
-----
Path coordinates are taken directly as **millimetres** by default. If the
root ``<svg>`` has both a ``viewBox`` and a physical ``width``/``height``
(``mm``/``cm``/``in``/``pt``/``pc``), the mm-per-user-unit scale is inferred
from those — the one unambiguous case. Anything else (unitless or ``px``
width, no viewBox, ...) falls back to 1:1; pass ``--scale`` to override
either way.

Homing and Z are both optional
-------------------------------
Neither is required. Without ``--z-offset-mm``, no ``Z`` is written anywhere
in the output — the path runs entirely at whatever height the *simulator*
supplies via ``--standoff-mm`` / ``PathConfig.standoff`` (this is not a
special case: :mod:`spraysim.gcode` already carries a move's last-known Z
forward when a line doesn't set one, so an SVG-derived path with no Z at all
sprays at a single fixed, pre-set height with no leading move required).
Without ``--home``, the output skips straight to the first subpath's travel
+ spray moves; with it, a leading ``G0 X0 Y0`` (plus ``Z`` if given) is
emitted first.

Usage
-----
    python svg_to_gcode.py drawing.svg                       # -> drawing.gcode
    python svg_to_gcode.py drawing.svg -o out.gcode --feed 2000
    python svg_to_gcode.py drawing.svg --home --z-offset-mm 150
    python svg_to_gcode.py drawing.svg --scale 0.5 --tolerance-mm 0.05

(run from the repo root; see the module docstring of ``spraysim/gcode.py``
for what the resulting file means to the simulator.)
"""

from __future__ import annotations

import argparse
import math
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

Point = tuple[float, float]
Subpath = list[Point]

DEFAULT_FEED_MM_MIN = 3000.0   # mm/min, matches spraysim.gcode.DEFAULT_FEED
DEFAULT_TOLERANCE_MM = 0.2     # mm, max chord error when linearising curves/arcs

_UNIT_TO_MM = {
    "mm": 1.0, "cm": 10.0, "in": 25.4,
    "pt": 25.4 / 72.0, "pc": 25.4 / 6.0, "px": 25.4 / 96.0,
}
_LENGTH_RE = re.compile(r"^\s*([-+]?[0-9]*\.?[0-9]+)\s*([a-zA-Z%]*)\s*$")
_FLOAT_RE = re.compile(r"[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?")


# --------------------------------------------------------------------------- #
# SVG "d" attribute parsing
# --------------------------------------------------------------------------- #

class _PathScanner:
    """Cursor over an SVG path ``d`` string (commands, floats, arc flags)."""

    def __init__(self, text: str):
        self.s = text
        self.i = 0
        self.n = len(text)

    def _skip_sep(self) -> None:
        while self.i < self.n and self.s[self.i] in " \t\r\n,":
            self.i += 1

    def at_end(self) -> bool:
        self._skip_sep()
        return self.i >= self.n

    def peek_is_letter(self) -> bool:
        self._skip_sep()
        return self.i < self.n and self.s[self.i].isalpha()

    def next_letter(self) -> str:
        self._skip_sep()
        c = self.s[self.i]
        self.i += 1
        return c

    def next_float(self) -> float:
        self._skip_sep()
        m = _FLOAT_RE.match(self.s, self.i)
        if not m:
            raise ValueError(f"expected a number at offset {self.i} in path data")
        self.i = m.end()
        return float(m.group())

    def next_flag(self) -> bool:
        """A single unseparated SVG arc flag: exactly one '0' or '1' character."""
        self._skip_sep()
        if self.i >= self.n or self.s[self.i] not in "01":
            raise ValueError(f"expected a flag (0/1) at offset {self.i} in path data")
        v = self.s[self.i]
        self.i += 1
        return v == "1"


def _lerp(a: Point, b: Point, t: float) -> Point:
    return (a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t)


def _point_line_distance(p: Point, a: Point, b: Point) -> float:
    dx, dy = b[0] - a[0], b[1] - a[1]
    if dx == 0.0 and dy == 0.0:
        return math.hypot(p[0] - a[0], p[1] - a[1])
    t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy)
    proj = (a[0] + t * dx, a[1] + t * dy)
    return math.hypot(p[0] - proj[0], p[1] - proj[1])


def _cubic_flatten(p0: Point, p1: Point, p2: Point, p3: Point,
                    tol: float, depth: int = 0, max_depth: int = 24) -> list[Point]:
    """Recursive de Casteljau subdivision; returns points after p0, ending at p3."""
    flat = (_point_line_distance(p1, p0, p3) <= tol
            and _point_line_distance(p2, p0, p3) <= tol)
    if flat or depth >= max_depth:
        return [p3]
    p01, p12, p23 = _lerp(p0, p1, 0.5), _lerp(p1, p2, 0.5), _lerp(p2, p3, 0.5)
    p012, p123 = _lerp(p01, p12, 0.5), _lerp(p12, p23, 0.5)
    p0123 = _lerp(p012, p123, 0.5)
    left = _cubic_flatten(p0, p01, p012, p0123, tol, depth + 1, max_depth)
    right = _cubic_flatten(p0123, p123, p23, p3, tol, depth + 1, max_depth)
    return left + right


def _quadratic_to_cubic(p0: Point, q: Point, p2: Point) -> tuple[Point, Point]:
    c1 = (p0[0] + 2.0 / 3.0 * (q[0] - p0[0]), p0[1] + 2.0 / 3.0 * (q[1] - p0[1]))
    c2 = (p2[0] + 2.0 / 3.0 * (q[0] - p2[0]), p2[1] + 2.0 / 3.0 * (q[1] - p2[1]))
    return c1, c2


def _arc_to_points(p0: Point, rx: float, ry: float, rot_deg: float,
                    large_arc: bool, sweep: bool, p1: Point, tol: float) -> list[Point]:
    """SVG elliptical arc p0 -> p1, endpoint-to-center form (SVG 1.1 appendix F.6.5)."""
    if rx == 0.0 or ry == 0.0 or p0 == p1:
        return [p1]
    rx, ry = abs(rx), abs(ry)
    phi = math.radians(rot_deg % 360.0)
    cphi, sphi = math.cos(phi), math.sin(phi)

    dx2, dy2 = (p0[0] - p1[0]) / 2.0, (p0[1] - p1[1]) / 2.0
    x1p = cphi * dx2 + sphi * dy2
    y1p = -sphi * dx2 + cphi * dy2

    lam = (x1p ** 2) / (rx ** 2) + (y1p ** 2) / (ry ** 2)
    if lam > 1.0:
        s = math.sqrt(lam)
        rx, ry = rx * s, ry * s

    sign = -1.0 if large_arc == sweep else 1.0
    num = rx ** 2 * ry ** 2 - rx ** 2 * y1p ** 2 - ry ** 2 * x1p ** 2
    den = rx ** 2 * y1p ** 2 + ry ** 2 * x1p ** 2
    co = sign * math.sqrt(max(num / den, 0.0)) if den > 0.0 else 0.0
    cxp, cyp = co * rx * y1p / ry, -co * ry * x1p / rx

    cx = cphi * cxp - sphi * cyp + (p0[0] + p1[0]) / 2.0
    cy = sphi * cxp + cphi * cyp + (p0[1] + p1[1]) / 2.0

    def angle(ux: float, uy: float, vx: float, vy: float) -> float:
        dot = ux * vx + uy * vy
        length = math.hypot(ux, uy) * math.hypot(vx, vy)
        a = math.acos(max(-1.0, min(1.0, dot / length))) if length else 0.0
        return -a if ux * vy - uy * vx < 0 else a

    theta1 = angle(1.0, 0.0, (x1p - cxp) / rx, (y1p - cyp) / ry)
    dtheta = angle((x1p - cxp) / rx, (y1p - cyp) / ry, (-x1p - cxp) / rx, (-y1p - cyp) / ry)
    if not sweep and dtheta > 0.0:
        dtheta -= 2.0 * math.pi
    elif sweep and dtheta < 0.0:
        dtheta += 2.0 * math.pi

    r_max = max(rx, ry)
    ratio = 1.0 - tol / r_max if r_max > tol else 0.0
    max_step = 2.0 * math.acos(max(-1.0, min(1.0, ratio))) if ratio < 1.0 else math.pi / 2.0
    max_step = max(max_step, math.radians(1.0))
    steps = max(1, math.ceil(abs(dtheta) / max_step))

    points = []
    for i in range(1, steps + 1):
        t = theta1 + dtheta * i / steps
        ex = cx + rx * math.cos(t) * cphi - ry * math.sin(t) * sphi
        ey = cy + rx * math.cos(t) * sphi + ry * math.sin(t) * cphi
        points.append((ex, ey))
    points[-1] = p1  # snap out float drift at the exact endpoint
    return points


def parse_path_d(d: str, tolerance: float) -> list[Subpath]:
    """Parse an SVG path ``d`` string into subpaths of flattened (x, y) points.

    ``tolerance`` is in the same units as the coordinates in ``d`` (apply the
    unit scale to it *before* calling this, so curves flatten to a physical
    error bound rather than a raw-coordinate one).
    """
    sc = _PathScanner(d)
    subpaths: list[Subpath] = []
    cur: Subpath | None = None
    x = y = start_x = start_y = 0.0
    last_cmd: str | None = None
    last_ctrl: Point | None = None  # reflected control point, for S/T

    while not sc.at_end():
        if sc.peek_is_letter():
            cmd = sc.next_letter()
        elif last_cmd is not None:
            cmd = "L" if last_cmd == "M" else "l" if last_cmd == "m" else last_cmd
        else:
            raise ValueError("path data must start with a moveto (M/m) command")

        upper, rel = cmd.upper(), cmd.islower()

        if cur is None and upper != "M":
            raise ValueError(
                f"path data must start with a moveto (M/m) command, got {cmd!r}"
            )

        if upper == "M":
            nx = sc.next_float() + (x if rel else 0.0)
            ny = sc.next_float() + (y if rel else 0.0)
            x, y = nx, ny
            start_x, start_y = x, y
            cur = [(x, y)]
            subpaths.append(cur)
            last_ctrl = None
            last_cmd = cmd
            continue

        if upper == "Z":
            x, y = start_x, start_y
            cur.append((x, y))
            last_ctrl = None

        elif upper == "L":
            x = sc.next_float() + (x if rel else 0.0)
            y = sc.next_float() + (y if rel else 0.0)
            cur.append((x, y))
            last_ctrl = None

        elif upper == "H":
            x = sc.next_float() + (x if rel else 0.0)
            cur.append((x, y))
            last_ctrl = None

        elif upper == "V":
            y = sc.next_float() + (y if rel else 0.0)
            cur.append((x, y))
            last_ctrl = None

        elif upper == "C":
            x1 = sc.next_float() + (x if rel else 0.0)
            y1 = sc.next_float() + (y if rel else 0.0)
            x2 = sc.next_float() + (x if rel else 0.0)
            y2 = sc.next_float() + (y if rel else 0.0)
            nx = sc.next_float() + (x if rel else 0.0)
            ny = sc.next_float() + (y if rel else 0.0)
            cur.extend(_cubic_flatten((x, y), (x1, y1), (x2, y2), (nx, ny), tolerance))
            last_ctrl, x, y = (x2, y2), nx, ny

        elif upper == "S":
            if last_cmd and last_cmd.upper() in ("C", "S") and last_ctrl is not None:
                x1, y1 = 2 * x - last_ctrl[0], 2 * y - last_ctrl[1]
            else:
                x1, y1 = x, y
            x2 = sc.next_float() + (x if rel else 0.0)
            y2 = sc.next_float() + (y if rel else 0.0)
            nx = sc.next_float() + (x if rel else 0.0)
            ny = sc.next_float() + (y if rel else 0.0)
            cur.extend(_cubic_flatten((x, y), (x1, y1), (x2, y2), (nx, ny), tolerance))
            last_ctrl, x, y = (x2, y2), nx, ny

        elif upper == "Q":
            x1 = sc.next_float() + (x if rel else 0.0)
            y1 = sc.next_float() + (y if rel else 0.0)
            nx = sc.next_float() + (x if rel else 0.0)
            ny = sc.next_float() + (y if rel else 0.0)
            c1, c2 = _quadratic_to_cubic((x, y), (x1, y1), (nx, ny))
            cur.extend(_cubic_flatten((x, y), c1, c2, (nx, ny), tolerance))
            last_ctrl, x, y = (x1, y1), nx, ny

        elif upper == "T":
            if last_cmd and last_cmd.upper() in ("Q", "T") and last_ctrl is not None:
                x1, y1 = 2 * x - last_ctrl[0], 2 * y - last_ctrl[1]
            else:
                x1, y1 = x, y
            nx = sc.next_float() + (x if rel else 0.0)
            ny = sc.next_float() + (y if rel else 0.0)
            c1, c2 = _quadratic_to_cubic((x, y), (x1, y1), (nx, ny))
            cur.extend(_cubic_flatten((x, y), c1, c2, (nx, ny), tolerance))
            last_ctrl, x, y = (x1, y1), nx, ny

        elif upper == "A":
            rx, ry, rot = sc.next_float(), sc.next_float(), sc.next_float()
            large, sweep = sc.next_flag(), sc.next_flag()
            nx = sc.next_float() + (x if rel else 0.0)
            ny = sc.next_float() + (y if rel else 0.0)
            cur.extend(_arc_to_points((x, y), rx, ry, rot, large, sweep, (nx, ny), tolerance))
            last_ctrl, x, y = None, nx, ny

        else:
            raise ValueError(f"unsupported path command {cmd!r}")

        last_cmd = cmd

    return subpaths


# --------------------------------------------------------------------------- #
# SVG document handling: units, shape collection, geometric transforms
# --------------------------------------------------------------------------- #

def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _parse_length_mm(text: str | None) -> float | None:
    if not text:
        return None
    m = _LENGTH_RE.match(text)
    if not m or m.group(2) == "%":
        return None
    unit = m.group(2).lower()
    if unit not in _UNIT_TO_MM:
        return None
    return float(m.group(1)) * _UNIT_TO_MM[unit]


def infer_scale(root: ET.Element) -> float:
    """mm per SVG user unit, from ``<svg width height viewBox>`` (best-effort).

    Only the unambiguous case — an explicit ``viewBox`` plus a physical
    ``width``/``height`` (mm/cm/in/pt/pc) — is auto-detected; everything else
    (no viewBox, unitless/``px``/``%`` width, ...) falls back to 1.0, i.e.
    path coordinates are taken directly as millimetres. Pass ``--scale`` to
    override either way.
    """
    vb = root.get("viewBox")
    if not vb:
        return 1.0
    parts = re.split(r"[,\s]+", vb.strip())
    if len(parts) != 4:
        return 1.0
    vb_w = float(parts[2])
    width_mm = _parse_length_mm(root.get("width"))
    if width_mm and vb_w:
        return width_mm / vb_w
    return 1.0


def _bbox(subpaths: list[Subpath]) -> tuple[float, float, float, float]:
    xs = [x for sp in subpaths for x, _ in sp]
    ys = [y for sp in subpaths for _, y in sp]
    return min(xs), max(xs), min(ys), max(ys)


def _flip_y(subpaths: list[Subpath]) -> list[Subpath]:
    _, _, ymin, ymax = _bbox(subpaths)
    return [[(x, ymin + ymax - y) for x, y in sp] for sp in subpaths]


def _transform(subpaths: list[Subpath], sx: float, sy: float,
                dx: float, dy: float) -> list[Subpath]:
    return [[(x * sx + dx, y * sy + dy) for x, y in sp] for sp in subpaths]


def convert(
    svg_path: str | Path,
    *,
    scale: float | None = None,
    tolerance_mm: float = DEFAULT_TOLERANCE_MM,
    flip_y: bool = True,
    normalize: bool = True,
    origin_mm: tuple[float, float] = (0.0, 0.0),
) -> tuple[list[Subpath], list[str]]:
    """Read ``svg_path`` and return (subpaths in mm, warnings).

    Subpaths are lists of ``(x, y)`` points in millimetres, flip/scale/
    normalize already applied, ready for :func:`moves_to_gcode`.
    """
    root = ET.parse(svg_path).getroot()
    warnings: list[str] = []

    s = scale if scale is not None else infer_scale(root)
    if s <= 0.0:
        raise ValueError(f"scale must be positive, got {s!r}")
    tol_user = tolerance_mm / s

    subpaths: list[Subpath] = []
    other_shapes = 0
    for el in root.iter():
        tag = _local(el.tag)
        if tag == "path":
            d = el.get("d")
            if not d:
                continue
            if el.get("transform"):
                warnings.append(
                    f"<path id={el.get('id')!r}> has a transform attribute; "
                    "ignored (flatten transforms in your SVG editor first)"
                )
            subpaths.extend(parse_path_d(d, tol_user))
        elif tag in ("rect", "circle", "ellipse", "line", "polyline", "polygon"):
            other_shapes += 1

    if other_shapes:
        warnings.append(
            f"{other_shapes} non-<path> shape element(s) ignored "
            "(only <path> is supported)"
        )

    subpaths = [sp for sp in subpaths if len(sp) >= 2]
    if not subpaths:
        raise ValueError(f"no drawable <path> segments found in {svg_path}")

    if flip_y:
        subpaths = _flip_y(subpaths)
    subpaths = _transform(subpaths, s, s, 0.0, 0.0)
    if normalize:
        xmin, _, ymin, _ = _bbox(subpaths)
        subpaths = _transform(subpaths, 1.0, 1.0, -xmin, -ymin)
    ox, oy = origin_mm
    if ox or oy:
        subpaths = _transform(subpaths, 1.0, 1.0, ox, oy)

    return subpaths, warnings


# --------------------------------------------------------------------------- #
# G-code emission
# --------------------------------------------------------------------------- #

def moves_to_gcode(
    subpaths: list[Subpath],
    *,
    feed_mm_min: float = DEFAULT_FEED_MM_MIN,
    z_mm: float | None = None,
    home: bool = False,
    precision: int = 3,
) -> str:
    """Render subpaths (mm) as G-code: G0 travel to each subpath, G1 within it.

    ``z_mm`` is optional — omit it (``None``) to write no ``Z`` at all, so the
    simulator's ``--standoff-mm`` sets the height instead. ``home`` is also
    optional: it only adds a leading ``G0 X0 Y0`` before spraying starts. If
    both are given, ``Z`` is written exactly once, on whichever line comes
    first (the home move, or else the first subpath's travel move) — every
    later move carries it forward implicitly, the same way a hand-written
    raster does (see ``examples/raster.gcode``).
    """
    def fmt(v: float) -> str:
        return f"{v:.{precision}f}"

    lines = ["G21", "G90"]
    wrote_z = False

    def z_suffix() -> str:
        nonlocal wrote_z
        if z_mm is None or wrote_z:
            return ""
        wrote_z = True
        return f" Z{fmt(z_mm)}"

    if home:
        lines.append(f"G0 X{fmt(0.0)} Y{fmt(0.0)}{z_suffix()}")

    lines.append(f"G1 F{feed_mm_min:g}")

    for sp in subpaths:
        x0, y0 = sp[0]
        lines.append(f"G0 X{fmt(x0)} Y{fmt(y0)}{z_suffix()}")
        for x, y in sp[1:]:
            lines.append(f"G1 X{fmt(x)} Y{fmt(y)}")

    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def main() -> None:
    p = argparse.ArgumentParser(
        description="Convert an SVG file's <path> elements into a G-code spray toolpath."
    )
    p.add_argument("svg", help="input SVG file")
    p.add_argument("-o", "--out", type=Path, default=None,
                   help="output .gcode path (default: <svg>.gcode next to the input)")
    p.add_argument("--feed", type=float, default=DEFAULT_FEED_MM_MIN,
                   help=f"feed rate (mm/min) for the whole path "
                        f"(default: {DEFAULT_FEED_MM_MIN:g})")
    p.add_argument("--scale", type=float, default=None,
                   help="mm per SVG user unit; overrides auto-detection from "
                        "viewBox + width/height (default: auto, falling back to 1.0)")
    p.add_argument("--tolerance-mm", type=float, default=DEFAULT_TOLERANCE_MM,
                   help="max flatness error (mm) when linearising curves/arcs "
                        f"(default: {DEFAULT_TOLERANCE_MM:g})")
    p.add_argument("--precision", type=int, default=3,
                   help="decimal places for emitted coordinates (default: 3)")
    p.add_argument("--home", action="store_true",
                   help="emit a leading G0 X0 Y0 travel move before spraying starts "
                        "(optional; default: off, so the path sprays directly)")
    p.add_argument("--z-offset-mm", type=float, default=None,
                   help="nozzle Z height (mm) to write into the G-code (optional; "
                        "default: none written, so the simulator's --standoff-mm "
                        "controls the height instead)")
    p.add_argument("--no-flip-y", action="store_true",
                   help="don't flip Y (SVG's Y grows downward; flipping makes the "
                        "sprayed pattern match what you see on screen)")
    p.add_argument("--no-normalize", action="store_true",
                   help="don't shift the artwork so its bounding box starts at "
                        "(0, 0); keep raw (scaled) SVG coordinates")
    p.add_argument("--origin-x-mm", type=float, default=0.0,
                   help="X offset (mm) applied after normalisation (default: 0)")
    p.add_argument("--origin-y-mm", type=float, default=0.0,
                   help="Y offset (mm) applied after normalisation (default: 0)")
    args = p.parse_args()

    svg_path = Path(args.svg)
    out_path = args.out or svg_path.with_suffix(".gcode")

    try:
        subpaths, warnings = convert(
            svg_path,
            scale=args.scale,
            tolerance_mm=args.tolerance_mm,
            flip_y=not args.no_flip_y,
            normalize=not args.no_normalize,
            origin_mm=(args.origin_x_mm, args.origin_y_mm),
        )
    except (ET.ParseError, ValueError, OSError) as exc:
        p.error(str(exc))

    for w in warnings:
        print(f"warning: {w}", file=sys.stderr)

    text = moves_to_gcode(
        subpaths, feed_mm_min=args.feed, z_mm=args.z_offset_mm,
        home=args.home, precision=args.precision,
    )
    out_path.write_text(text)

    n_points = sum(len(sp) for sp in subpaths)
    print(f"{len(subpaths)} subpath(s), {n_points} points -> {out_path}")


if __name__ == "__main__":
    main()
