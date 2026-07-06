"""Minimal G-code parser for spray toolpaths.

Parses a pragmatic subset of G-code into a list of straight-line motion segments
(:class:`Move`). The nozzle sprays while it travels a **G1** (linear) move and is
off during a **G0** (rapid / travel) move.

Supported subset
----------------
- ``G0`` (travel, spray off) / ``G1`` (linear, spray on); modal (a bare
  coordinate line reuses the last motion mode).
- ``X Y Z`` coordinates and ``F`` feed rate.
- ``G90`` / ``G91`` absolute / relative positioning.
- ``G20`` / ``G21`` inch / millimetre units (**mm by default**).
- ``;`` and ``( ... )`` comments.

Everything is converted to **SI** internally: positions in metres, feed in m/s.
Arc moves ``G2``/``G3`` are **rejected** (linearise the path first); other
unrecognised codes (``M``/``T``/``S``/``E`` ...) are ignored.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

DEFAULT_FEED = 0.05        # m/s (= 3000 mm/min) when a program sets no feed
DEFAULT_STANDOFF = 0.15    # m, nozzle height used when the path has no Z

_WORD = re.compile(r"([A-Za-z])\s*([-+]?[0-9]*\.?[0-9]+)")


@dataclass(frozen=True)
class Move:
    """A straight-line nozzle motion segment (SI units)."""

    start: tuple[float, float, float]   # m
    end: tuple[float, float, float]     # m
    feed: float                          # m/s
    spray_on: bool

    @property
    def length(self) -> float:
        (x0, y0, z0), (x1, y1, z1) = self.start, self.end
        return math.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2)

    @property
    def duration(self) -> float:
        return self.length / self.feed if self.feed > 0.0 else 0.0


def _strip_comment(line: str) -> str:
    line = re.sub(r"\(.*?\)", "", line)   # ( inline ) comments
    semi = line.find(";")
    if semi >= 0:
        line = line[:semi]
    return line.strip()


def parse_gcode(
    text: str,
    *,
    default_feed: float = DEFAULT_FEED,
    standoff: float = DEFAULT_STANDOFF,
    feed_override: float | None = None,
) -> list[Move]:
    """Parse G-code ``text`` into a list of :class:`Move`.

    ``standoff`` (m) sets the initial nozzle height (used throughout unless the
    program moves in Z). ``feed_override`` (m/s), if given, forces the feed on
    every move and ignores ``F`` words.
    """
    absolute = True
    unit = 1.0e-3                      # mm -> m
    pos = [0.0, 0.0, float(standoff)]
    feed = float(feed_override if feed_override is not None else default_feed)
    motion: int | None = None          # modal motion mode: 0 or 1
    moves: list[Move] = []

    for raw in text.splitlines():
        line = _strip_comment(raw)
        if not line:
            continue

        target: dict[str, float] = {}
        cmd_motion: int | None = None
        for letter, value in _WORD.findall(line):
            L = letter.upper()
            num = float(value)
            if L == "G":
                g = int(round(num))
                if g in (0, 1):
                    cmd_motion = motion = g
                elif g == 90:
                    absolute = True
                elif g == 91:
                    absolute = False
                elif g == 20:
                    unit = 25.4e-3
                elif g == 21:
                    unit = 1.0e-3
                elif g in (2, 3):
                    raise ValueError(
                        "Arc moves (G2/G3) are not supported; linearise the path."
                    )
                # other G codes (G17, G54, ...) are ignored
            elif L == "F":
                if feed_override is None:
                    feed = num * unit / 60.0   # units/min -> m/s
            elif L in ("X", "Y", "Z"):
                target[L] = num * unit
            # M/T/S/E and other words are ignored

        if not target:
            continue
        mode = cmd_motion if cmd_motion is not None else motion
        if mode is None:
            continue  # coordinates before any G0/G1 — nothing to do

        new = list(pos)
        for i, ax in enumerate("XYZ"):
            if ax in target:
                new[i] = target[ax] if absolute else pos[i] + target[ax]
        moves.append(Move(tuple(pos), tuple(new), feed, spray_on=(mode == 1)))
        pos = new

    return moves


def load_moves(source: str, **kwargs) -> list[Move]:
    """Parse moves from inline G-code text (contains a newline) or a file path."""
    text = source if "\n" in source else Path(source).read_text()
    return parse_gcode(text, **kwargs)


def total_spray_time(moves: list[Move]) -> float:
    """Total time (s) spent spraying (G1 segments)."""
    return sum(m.duration for m in moves if m.spray_on)


def spray_length(moves: list[Move]) -> float:
    """Total path length (m) sprayed (G1 segments)."""
    return sum(m.length for m in moves if m.spray_on)


def bounds(moves: list[Move]) -> tuple[float, float, float, float]:
    """(xmin, xmax, ymin, ymax) over all move endpoints (m)."""
    xs, ys = [], []
    for m in moves:
        xs += [m.start[0], m.end[0]]
        ys += [m.start[1], m.end[1]]
    if not xs:
        return (0.0, 0.0, 0.0, 0.0)
    return (min(xs), max(xs), min(ys), max(ys))
