"""Tests for gcode_to_svg.py, including svg_to_gcode <-> gcode_to_svg round trips
(the reason this tool exists: validating svg_to_gcode.py end to end)."""

import math
import sys

import pytest

import gcode_to_svg as g2s
import svg_to_gcode as s2g
from spraysim import gcode


def _moves(text: str):
    return gcode.parse_gcode(text, standoff=0.15)


# --- moves_to_subpaths: geometry-based grouping ------------------------------ #

def test_chained_spray_moves_form_one_subpath():
    moves = _moves("G21 G90\nG1 F3000 X10 Y0\nG1 X10 Y10\n")
    subpaths = g2s.moves_to_subpaths(moves, lambda m: m.spray_on)
    assert subpaths == [[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0)]]


def test_travel_move_breaks_the_spray_chain():
    moves = _moves("G21 G90\nG1 F3000 X10 Y0\nG0 X10 Y10\nG1 X20 Y10\n")
    subpaths = g2s.moves_to_subpaths(moves, lambda m: m.spray_on)
    assert subpaths == [[(0.0, 0.0), (10.0, 0.0)], [(10.0, 10.0), (20.0, 10.0)]]


def test_predicate_selects_travel_only():
    moves = _moves("G21 G90\nG0 X10 Y0\nG1 F3000 X20 Y0\n")
    travel = g2s.moves_to_subpaths(moves, lambda m: not m.spray_on)
    spray = g2s.moves_to_subpaths(moves, lambda m: m.spray_on)
    assert travel == [[(0.0, 0.0), (10.0, 0.0)]]
    assert spray == [[(10.0, 0.0), (20.0, 0.0)]]


def test_disconnected_moves_never_merge_even_without_travel():
    """A coordinate jump with no explicit G0 in between still breaks the chain
    (grouping is geometry-based, not G-code-structure-based)."""
    moves = [
        gcode.Move((0.0, 0.0, 0.15), (0.01, 0.0, 0.15), 0.05, True),
        gcode.Move((0.05, 0.05, 0.15), (0.06, 0.05, 0.15), 0.05, True),  # doesn't connect
    ]
    subpaths = g2s.moves_to_subpaths(moves, lambda m: m.spray_on)
    assert len(subpaths) == 2


# --- prepare_subpaths: the flip is a well-defined, self-inverse reflection -- #

def test_flip_reflects_about_spray_bbox():
    # No leading "X0 Y0": the parser's default start position is already
    # (0, 0), so restating it would produce a zero-length phantom move.
    moves = _moves("G21 G90\nG1 F3000\nG1 X10 Y0\nG1 X10 Y10\n")
    spray, _ = g2s.prepare_subpaths(moves, flip_y=True)
    # y in [0, 10] before flip; after reflection about that range, 0 <-> 10.
    assert spray[0][0] == (0.0, 10.0)
    assert spray[0][2] == (10.0, 0.0)


def test_no_flip_keeps_raw_coordinates():
    moves = _moves("G21 G90\nG1 F3000\nG1 X10 Y10\n")
    spray, _ = g2s.prepare_subpaths(moves, flip_y=False)
    assert spray == [[(0.0, 0.0), (10.0, 10.0)]]


def test_double_flip_is_identity():
    """Flipping twice about the same bbox must recover the original points —
    this is exactly what makes gcode_to_svg the inverse of svg_to_gcode."""
    moves = _moves("G21 G90\nG1 F3000\nG1 X30 Y5\nG1 X10 Y20\n")
    once, _ = g2s.prepare_subpaths(moves, flip_y=True)
    # Flip again by treating `once`'s own bbox as the new basis.
    ys = [y for sp in once for _, y in sp]
    ymin, ymax = min(ys), max(ys)
    twice = [[(x, ymin + ymax - y) for x, y in sp] for sp in once]
    original, _ = g2s.prepare_subpaths(moves, flip_y=False)
    assert twice == original


# --- render_svg: output is well-formed and correctly sized ------------------ #

def test_render_svg_viewbox_covers_all_drawn_points(tmp_path):
    spray = [[(0.0, 0.0), (100.0, 0.0), (100.0, 50.0)]]
    out = tmp_path / "out.svg"
    g2s.render_svg(spray, [], out, margin_mm=5.0)
    text = out.read_text()
    assert 'viewBox="-5.000 -5.000 110.000 60.000"' in text
    assert "M0.000,0.000" in text and "L100.000,50.000" in text


def test_render_svg_travel_drawn_dashed_and_spray_solid(tmp_path):
    out = tmp_path / "out.svg"
    g2s.render_svg([[(0.0, 0.0), (10.0, 0.0)]], [[(10.0, 0.0), (10.0, 10.0)]], out)
    text = out.read_text()
    assert "stroke-dasharray" in text  # travel
    assert text.count("<path") == 2


def test_render_svg_raises_on_empty_input(tmp_path):
    with pytest.raises(ValueError, match="nothing to draw"):
        g2s.render_svg([], [], tmp_path / "empty.svg")


# --- convert(): CLI-level entry point --------------------------------------- #

def test_convert_writes_file_and_returns_written_points(tmp_path):
    gcode_path = tmp_path / "line.gcode"
    gcode_path.write_text("G21 G90\nG1 F3000 X10 Y0\n")
    out = tmp_path / "line.svg"

    spray, travel = g2s.convert(str(gcode_path), out, flip_y=False)
    assert spray == [[(0.0, 0.0), (10.0, 0.0)]]
    assert travel == []
    assert out.exists()
    assert "M0.000,0.000" in out.read_text()


def test_cli_end_to_end(tmp_path, capsys, monkeypatch):
    gcode_path = tmp_path / "line.gcode"
    gcode_path.write_text("G21 G90\nG1 F3000 X10 Y0\nG1 X10 Y10\n")

    monkeypatch.setattr(sys, "argv", ["gcode_to_svg.py", str(gcode_path)])
    g2s.main()
    out = capsys.readouterr().out

    default_out = gcode_path.with_suffix(".svg")
    assert default_out.exists()
    assert "1 spray subpath(s), 0 travel subpath(s)" in out


def test_cli_inline_text_requires_out(capsys, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["gcode_to_svg.py", "G21 G90\nG1 F3000 X10 Y0\n"])
    with pytest.raises(SystemExit):
        g2s.main()
    assert "--out is required" in capsys.readouterr().err


# --- round trip with svg_to_gcode.py: the actual validation use case -------- #

def test_round_trip_polyline_svg_is_exact(tmp_path):
    """A straight-edged SVG path survives svg_to_gcode -> gcode_to_svg exactly
    (no curve flattening involved, so there's no tolerance to account for)."""
    svg_path = tmp_path / "square.svg"
    d = "M10,10 L90,10 L90,90 L10,90 Z"
    svg_path.write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" width="100mm" height="100mm" '
        f'viewBox="0 0 100 100"><path d="{d}"/></svg>'
    )
    original = s2g.parse_path_d(d, tolerance=0.05)

    # normalize=False so the gcode keeps the same bbox as the original SVG —
    # otherwise the round trip is only exact up to a known translation.
    subpaths, _ = s2g.convert(svg_path, normalize=False)
    gcode_text = s2g.moves_to_gcode(subpaths, home=False, z_mm=None)
    gcode_path = tmp_path / "square.gcode"
    gcode_path.write_text(gcode_text)

    recovered, _ = g2s.convert(str(gcode_path), tmp_path / "square_rt.svg")
    assert recovered == original


def test_round_trip_curved_svg_stays_within_tolerance(tmp_path):
    """A curved path can't round-trip exactly (it's linearised into G-code),
    but every recovered point must lie within svg_to_gcode's flattening
    tolerance of the true analytic curve — proving the whole svg_to_gcode ->
    G-code -> gcode_to_svg pipeline preserves shape, not just endpoints."""
    svg_path = tmp_path / "curve.svg"
    d = "M10,50 C10,20 90,20 90,50"
    svg_path.write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" width="100mm" height="100mm" '
        f'viewBox="0 0 100 100"><path d="{d}"/></svg>'
    )
    tol = 0.05
    subpaths, _ = s2g.convert(svg_path, normalize=False, tolerance_mm=tol)
    gcode_text = s2g.moves_to_gcode(subpaths, home=False, z_mm=None)
    gcode_path = tmp_path / "curve.gcode"
    gcode_path.write_text(gcode_text)

    recovered, _ = g2s.convert(str(gcode_path), tmp_path / "curve_rt.svg")

    p0, p1, p2, p3 = (10.0, 50.0), (10.0, 20.0), (90.0, 20.0), (90.0, 50.0)

    def bezier(t):
        mt = 1 - t
        x = mt**3*p0[0] + 3*mt**2*t*p1[0] + 3*mt*t**2*p2[0] + t**3*p3[0]
        y = mt**3*p0[1] + 3*mt**2*t*p1[1] + 3*mt*t**2*p2[1] + t**3*p3[1]
        return x, y

    samples = [bezier(t / 400.0) for t in range(401)]
    assert recovered[0][0] == p0
    assert recovered[0][-1] == p3
    for px, py in recovered[0]:
        best = min(math.hypot(px - sx, py - sy) for sx, sy in samples)
        assert best < tol * 5  # generous vs. the subdivision tolerance


def test_round_trip_multi_subpath_svg_preserves_subpath_count(tmp_path):
    svg_path = tmp_path / "two_lines.svg"
    d = "M0,0 L100,0 M0,20 L100,20"
    svg_path.write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" width="100mm" height="100mm" '
        f'viewBox="0 0 100 100"><path d="{d}"/></svg>'
    )
    subpaths, _ = s2g.convert(svg_path, normalize=False)
    gcode_text = s2g.moves_to_gcode(subpaths, home=True, z_mm=150.0)
    gcode_path = tmp_path / "two_lines.gcode"
    gcode_path.write_text(gcode_text)

    recovered, _ = g2s.convert(str(gcode_path), tmp_path / "two_lines_rt.svg")
    assert len(recovered) == 2
    original = s2g.parse_path_d(d, tolerance=0.05)
    assert recovered == original
