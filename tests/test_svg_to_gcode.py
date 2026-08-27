"""Tests for svg_to_gcode.py: SVG path parsing, curve/arc flattening, and the
resulting G-code's compatibility with spraysim.gcode."""

import math
import sys

import pytest

import svg_to_gcode as s2g
from spraysim import gcode


def _svg(body: str, *, width="100mm", height="100mm", viewbox="0 0 100 100") -> str:
    vb = f' viewBox="{viewbox}"' if viewbox else ""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
            f'height="{height}"{vb}>{body}</svg>')


# --- path "d" parsing ------------------------------------------------------- #

def test_parse_absolute_lines():
    subpaths = s2g.parse_path_d("M0,0 L10,0 L10,10", tolerance=0.1)
    assert subpaths == [[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0)]]


def test_parse_relative_commands_match_absolute():
    rel = s2g.parse_path_d("m0,0 l10,0 l0,10", tolerance=0.1)
    absolute = s2g.parse_path_d("M0,0 L10,0 L10,10", tolerance=0.1)
    assert rel == absolute


def test_implicit_repeated_lineto():
    """Numbers after L without a new command letter repeat the lineto."""
    subpaths = s2g.parse_path_d("M0,0 L10,0 20,0 20,10", tolerance=0.1)
    assert subpaths == [[(0.0, 0.0), (10.0, 0.0), (20.0, 0.0), (20.0, 10.0)]]


def test_horizontal_vertical_shorthand():
    subpaths = s2g.parse_path_d("M0,0 H10 V10", tolerance=0.1)
    assert subpaths == [[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0)]]


def test_closepath_returns_to_subpath_start():
    subpaths = s2g.parse_path_d("M0,0 L10,0 L10,10 Z", tolerance=0.1)
    assert subpaths[0][-1] == (0.0, 0.0)


def test_multiple_subpaths_from_multiple_moveto():
    subpaths = s2g.parse_path_d("M0,0 L10,0 M5,5 L15,5", tolerance=0.1)
    assert len(subpaths) == 2
    assert subpaths[0] == [(0.0, 0.0), (10.0, 0.0)]
    assert subpaths[1] == [(5.0, 5.0), (15.0, 5.0)]


def test_cubic_bezier_endpoint_exact_and_interior_points_near_curve():
    p0, p1, p2, p3 = (0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0)
    tol = 0.05
    subpaths = s2g.parse_path_d(
        f"M{p0[0]},{p0[1]} C{p1[0]},{p1[1]} {p2[0]},{p2[1]} {p3[0]},{p3[1]}", tolerance=tol
    )
    pts = subpaths[0]
    assert pts[0] == p0
    assert pts[-1] == p3
    assert len(pts) > 4  # actually flattened, not just the 2 endpoints

    def bezier(t):
        mt = 1 - t
        x = mt**3*p0[0] + 3*mt**2*t*p1[0] + 3*mt*t**2*p2[0] + t**3*p3[0]
        y = mt**3*p0[1] + 3*mt**2*t*p1[1] + 3*mt*t**2*p2[1] + t**3*p3[1]
        return x, y

    # Every flattened point should land close to *some* point on the true
    # curve (loose sampling check, generous vs. the subdivision tolerance).
    samples = [bezier(t / 200.0) for t in range(201)]
    for px, py in pts:
        best = min(math.hypot(px - sx, py - sy) for sx, sy in samples)
        assert best < tol * 5


def test_quadratic_bezier_endpoint_exact():
    subpaths = s2g.parse_path_d("M0,0 Q5,10 10,0", tolerance=0.05)
    pts = subpaths[0]
    assert pts[0] == (0.0, 0.0)
    assert pts[-1] == (10.0, 0.0)
    # The curve should bulge toward the control point (peak y > 0).
    assert max(y for _, y in pts) > 3.0


def test_smooth_cubic_reflects_previous_control_point():
    """S after C should continue smoothly (S's implicit control point is the
    reflection of C's second control point)."""
    subpaths = s2g.parse_path_d("M0,0 C0,10 10,10 10,0 S20,-10 20,0", tolerance=0.05)
    assert len(subpaths) == 1
    assert subpaths[0][-1] == (20.0, 0.0)


def test_arc_quarter_circle_stays_on_radius():
    pts = s2g._arc_to_points((10.0, 0.0), 10.0, 10.0, 0.0, False, True, (0.0, 10.0), tol=0.01)
    assert pts[-1] == (0.0, 10.0)
    for x, y in pts:
        assert math.hypot(x, y) == pytest.approx(10.0, abs=1e-6)


def test_arc_degenerate_zero_radius_falls_back_to_line():
    pts = s2g._arc_to_points((0.0, 0.0), 0.0, 0.0, 0.0, False, True, (5.0, 5.0), tol=0.01)
    assert pts == [(5.0, 5.0)]


def test_unsupported_command_raises_clearly():
    with pytest.raises(ValueError, match="unsupported path command"):
        s2g.parse_path_d("M0,0 B10,10", tolerance=0.1)


def test_missing_moveto_raises_clearly():
    with pytest.raises(ValueError, match="moveto"):
        s2g.parse_path_d("L10,10", tolerance=0.1)


# --- unit / scale inference -------------------------------------------------- #

def test_infer_scale_from_viewbox_and_mm_width(tmp_path):
    import xml.etree.ElementTree as ET
    svg = _svg("", width="50mm", height="50mm", viewbox="0 0 500 500")
    root = ET.fromstring(svg)
    assert s2g.infer_scale(root) == pytest.approx(0.1)  # 50mm / 500 units


def test_infer_scale_defaults_to_one_without_viewbox():
    import xml.etree.ElementTree as ET
    root = ET.fromstring(_svg("", viewbox=None))
    assert s2g.infer_scale(root) == pytest.approx(1.0)


# --- convert(): full pipeline (flip, scale, normalize) ----------------------- #

def test_convert_flips_y_and_normalizes_to_origin(tmp_path):
    svg_path = tmp_path / "square.svg"
    svg_path.write_text(_svg('<path d="M10,10 L90,10 L90,90 L10,90 Z"/>'))

    subpaths, warnings = s2g.convert(svg_path, flip_y=True, normalize=True)
    assert not warnings
    xs = [x for sp in subpaths for x, _ in sp]
    ys = [y for sp in subpaths for _, y in sp]
    assert min(xs) == pytest.approx(0.0)
    assert min(ys) == pytest.approx(0.0)
    assert max(xs) == pytest.approx(80.0)  # 90-10 inset square, 1mm/unit scale
    assert max(ys) == pytest.approx(80.0)


def test_convert_no_flip_preserves_svg_y_orientation(tmp_path):
    svg_path = tmp_path / "tall.svg"
    svg_path.write_text(_svg('<path d="M0,0 L0,10 L10,10"/>'))

    flipped, _ = s2g.convert(svg_path, flip_y=True, normalize=False)
    unflipped, _ = s2g.convert(svg_path, flip_y=False, normalize=False)
    # Flipping should invert the relative y-order of the first two points.
    assert flipped[0][0][1] > flipped[0][1][1]
    assert unflipped[0][0][1] < unflipped[0][1][1]


def test_convert_warns_on_non_path_shapes_and_transform(tmp_path):
    svg_path = tmp_path / "mixed.svg"
    svg_path.write_text(_svg(
        '<rect x="0" y="0" width="10" height="10"/>'
        '<path id="p1" d="M0,0 L10,10" transform="translate(5,5)"/>'
    ))
    _, warnings = s2g.convert(svg_path)
    joined = " ".join(warnings)
    assert "non-<path>" in joined
    assert "transform" in joined


def test_convert_raises_on_no_paths(tmp_path):
    svg_path = tmp_path / "empty.svg"
    svg_path.write_text(_svg('<rect x="0" y="0" width="10" height="10"/>'))
    with pytest.raises(ValueError, match="no drawable"):
        s2g.convert(svg_path)


def test_convert_origin_offset_applied_after_normalize(tmp_path):
    svg_path = tmp_path / "square.svg"
    svg_path.write_text(_svg('<path d="M0,0 L10,0 L10,10 L0,10 Z"/>'))
    subpaths, _ = s2g.convert(svg_path, normalize=True, origin_mm=(5.0, 7.0))
    xs = [x for sp in subpaths for x, _ in sp]
    ys = [y for sp in subpaths for _, y in sp]
    assert min(xs) == pytest.approx(5.0)
    assert min(ys) == pytest.approx(7.0)


# --- moves_to_gcode(): homing / Z are independently optional ----------------- #

def test_no_home_no_z_writes_no_z_anywhere():
    text = s2g.moves_to_gcode([[(0.0, 0.0), (10.0, 0.0)]], home=False, z_mm=None)
    assert "Z" not in text


def test_home_without_z_has_no_z_either():
    text = s2g.moves_to_gcode([[(0.0, 0.0), (10.0, 0.0)]], home=True, z_mm=None)
    assert "Z" not in text
    assert any(line.startswith("G0 X0.000 Y0.000") for line in text.splitlines())


def test_z_without_home_written_once_on_first_travel_move():
    text = s2g.moves_to_gcode([[(0.0, 0.0), (10.0, 0.0)]], home=False, z_mm=150.0)
    z_lines = [l for l in text.splitlines() if "Z" in l]
    assert len(z_lines) == 1
    assert z_lines[0].startswith("G0") and "Z150.000" in z_lines[0]


def test_z_and_home_written_once_on_home_line_only():
    text = s2g.moves_to_gcode(
        [[(0.0, 0.0), (10.0, 0.0)], [(5.0, 5.0), (15.0, 5.0)]], home=True, z_mm=200.0
    )
    z_lines = [l for l in text.splitlines() if "Z" in l]
    assert len(z_lines) == 1
    assert z_lines[0].startswith("G0 X0.000 Y0.000")


def test_gcode_output_parses_with_spraysim_gcode_module(tmp_path):
    svg_path = tmp_path / "raster_like.svg"
    svg_path.write_text(_svg(
        '<path d="M0,0 L100,0 M0,20 L100,20"/>', width="100mm", height="100mm",
        viewbox="0 0 100 100",
    ))
    subpaths, _ = s2g.convert(svg_path)
    text = s2g.moves_to_gcode(subpaths, feed_mm_min=3000.0, z_mm=150.0, home=True)

    out = tmp_path / "raster_like.gcode"
    out.write_text(text)
    moves = gcode.load_moves(str(out))

    spray = [m for m in moves if m.spray_on]
    assert len(spray) == 2
    assert gcode.spray_length(moves) == pytest.approx(0.2)  # two 100mm passes, in metres
    xmin, xmax, ymin, ymax = gcode.bounds(moves)
    assert (xmin, xmax) == pytest.approx((0.0, 0.1))
    assert (ymin, ymax) == pytest.approx((0.0, 0.02))


def test_feed_override_shortens_spray_time(tmp_path):
    svg_path = tmp_path / "line.svg"
    svg_path.write_text(_svg('<path d="M0,0 L100,0"/>'))
    subpaths, _ = s2g.convert(svg_path)

    slow = s2g.moves_to_gcode(subpaths, feed_mm_min=1000.0)
    fast = s2g.moves_to_gcode(subpaths, feed_mm_min=5000.0)
    slow_moves = gcode.parse_gcode(slow)
    fast_moves = gcode.parse_gcode(fast)
    assert gcode.total_spray_time(fast_moves) < gcode.total_spray_time(slow_moves)


# --- moves_to_gcode(closed_loop=True): return pass ---------------------------- #

def test_closed_loop_false_is_no_regression():
    subpaths = [[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0)], [(5.0, 5.0), (15.0, 5.0)]]
    baseline = s2g.moves_to_gcode(subpaths, home=True, z_mm=150.0)
    explicit = s2g.moves_to_gcode(subpaths, home=True, z_mm=150.0, closed_loop=False)
    assert explicit == baseline


def test_closed_loop_return_waypoints_are_exact_reverse_of_forward_pass():
    subpaths = [[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0)]]
    text = s2g.moves_to_gcode(subpaths, home=False, closed_loop=True)
    lines = text.splitlines()

    return_start = next(i for i, l in enumerate(lines) if l.startswith("G0 F"))
    return_lines = lines[return_start + 1:]
    return_coords = [
        tuple(float(tok[1:]) for tok in l.split()[1:]) for l in return_lines
    ]
    # Forward waypoints, reversed, excluding the point already reached
    # (the end of the forward path).
    assert return_coords == [(10.0, 0.0), (0.0, 0.0)]


def test_closed_loop_return_moves_are_all_g0_never_g1():
    subpaths = [[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0)], [(5.0, 5.0), (15.0, 5.0)]]
    text = s2g.moves_to_gcode(subpaths, home=True, closed_loop=True)
    lines = text.splitlines()

    return_start = next(i for i, l in enumerate(lines) if l.startswith("G0 F"))
    return_lines = lines[return_start + 1:]
    assert return_lines  # non-empty for this multi-point input
    assert all(l.startswith("G0 ") for l in return_lines)
    assert not any(l.startswith("G1") for l in return_lines)


def test_closed_loop_return_section_starts_with_default_feed():
    text = s2g.moves_to_gcode([[(0.0, 0.0), (10.0, 0.0)]], home=False, closed_loop=True)
    assert f"G0 F{s2g.DEFAULT_RETURN_FEED_MM_MIN:g}" in text.splitlines()


def test_closed_loop_return_section_honors_custom_return_feed():
    text = s2g.moves_to_gcode(
        [[(0.0, 0.0), (10.0, 0.0)]], home=False, closed_loop=True, return_feed_mm_min=1500.0
    )
    assert "G0 F1500" in text.splitlines()
    assert f"G0 F{s2g.DEFAULT_RETURN_FEED_MM_MIN:g}" not in text.splitlines()


def test_closed_loop_ends_at_first_forward_waypoint_without_home():
    subpaths = [[(3.0, 4.0), (10.0, 0.0), (10.0, 10.0)]]
    text = s2g.moves_to_gcode(subpaths, home=False, closed_loop=True)
    last_line = text.splitlines()[-1]
    assert last_line == "G0 X3.000 Y4.000"


def test_closed_loop_ends_at_home_point_when_home_is_true():
    subpaths = [[(5.0, 5.0), (15.0, 5.0)]]
    text = s2g.moves_to_gcode(subpaths, home=True, closed_loop=True)
    last_line = text.splitlines()[-1]
    assert last_line == "G0 X0.000 Y0.000"


def test_closed_loop_single_point_forward_path_emits_no_return_moves():
    # Only one waypoint total -> waypoints[:-1] is empty -> nothing to retrace.
    text = s2g.moves_to_gcode([[(5.0, 5.0)]], home=False, closed_loop=True)
    lines = text.splitlines()
    return_start = next(i for i, l in enumerate(lines) if l.startswith("G0 F"))
    assert lines[return_start + 1:] == []


# --- convert(fit_box_mm=...) / _fit_to_box(): fit-and-center into a mm box -- #

def test_fit_box_wide_artwork_scales_to_width_and_centers_vertically():
    # 100x10, box is 50x50 -> width-limited scale 0.5 -> 50x5, centered on Y.
    subpaths = s2g._fit_to_box(
        [[(0.0, 0.0), (100.0, 0.0), (100.0, 10.0), (0.0, 10.0)]], (0.0, 0.0, 50.0, 50.0)
    )
    xs = [x for sp in subpaths for x, _ in sp]
    ys = [y for sp in subpaths for _, y in sp]
    assert (min(xs), max(xs)) == pytest.approx((0.0, 50.0))
    assert (min(ys), max(ys)) == pytest.approx((22.5, 27.5))


def test_fit_box_tall_artwork_scales_to_height_and_centers_horizontally():
    # 10x100, box is 50x50 -> height-limited scale 0.5 -> 5x50, centered on X.
    subpaths = s2g._fit_to_box(
        [[(0.0, 0.0), (10.0, 0.0), (10.0, 100.0), (0.0, 100.0)]], (0.0, 0.0, 50.0, 50.0)
    )
    xs = [x for sp in subpaths for x, _ in sp]
    ys = [y for sp in subpaths for _, y in sp]
    assert (min(xs), max(xs)) == pytest.approx((22.5, 27.5))
    assert (min(ys), max(ys)) == pytest.approx((0.0, 50.0))


def test_fit_box_offset_box_centers_relative_to_box_origin():
    subpaths = s2g._fit_to_box(
        [[(0.0, 0.0), (100.0, 0.0), (100.0, 10.0), (0.0, 10.0)]], (10.0, 20.0, 60.0, 70.0)
    )
    xs = [x for sp in subpaths for x, _ in sp]
    ys = [y for sp in subpaths for _, y in sp]
    assert (min(xs), max(xs)) == pytest.approx((10.0, 60.0))
    assert (min(ys), max(ys)) == pytest.approx((42.5, 47.5))


def test_fit_box_invalid_box_raises_value_error():
    with pytest.raises(ValueError, match="invalid fit box"):
        s2g._fit_to_box([[(0.0, 0.0), (10.0, 0.0)]], (10.0, 0.0, 10.0, 50.0))


def test_convert_fit_box_mm_handles_zero_width_source_without_zero_division(tmp_path):
    # A vertical line has a zero-width bbox; scale must fall back to the
    # height-only ratio instead of dividing by zero.
    svg_path = tmp_path / "vertical_line.svg"
    svg_path.write_text(_svg('<path d="M5,0 L5,10"/>'))
    subpaths, _ = s2g.convert(svg_path, fit_box_mm=(0.0, 0.0, 50.0, 50.0))
    xs = [x for sp in subpaths for x, _ in sp]
    ys = [y for sp in subpaths for _, y in sp]
    assert all(math.isfinite(v) for v in xs + ys)
    assert 0.0 <= min(xs) and max(xs) <= 50.0
    assert (min(ys), max(ys)) == pytest.approx((0.0, 50.0))


def test_convert_fit_box_mm_applied_after_normalize(tmp_path):
    svg_path = tmp_path / "wide.svg"
    svg_path.write_text(_svg('<path d="M0,0 L100,0 L100,10 L0,10 Z"/>'))
    subpaths, _ = s2g.convert(svg_path, fit_box_mm=(0.0, 0.0, 50.0, 50.0))
    xs = [x for sp in subpaths for x, _ in sp]
    ys = [y for sp in subpaths for _, y in sp]
    assert (min(xs), max(xs)) == pytest.approx((0.0, 50.0))
    assert (min(ys), max(ys)) == pytest.approx((22.5, 27.5))


# --- main(): --fit-box-mm CLI flag ------------------------------------------- #

def test_cli_fit_box_mm_fits_and_centers_artwork(tmp_path, capsys, monkeypatch):
    svg_path = tmp_path / "wide.svg"
    svg_path.write_text(_svg('<path d="M0,0 L100,0 L100,10 L0,10 Z"/>'))
    out_path = tmp_path / "wide.gcode"

    monkeypatch.setattr(sys, "argv", [
        "svg_to_gcode.py", str(svg_path), "-o", str(out_path),
        "--fit-box-mm", "0", "0", "50", "50",
    ])
    s2g.main()
    capsys.readouterr()

    # gcode.bounds() would also include the machine's implicit (0, 0) starting
    # position (the first G0 travel move's start point); restrict to the
    # sprayed (G1) segments, which is what --fit-box-mm actually controls.
    spray = [m for m in gcode.load_moves(str(out_path)) if m.spray_on]
    xs = [m.start[0] for m in spray] + [m.end[0] for m in spray]
    ys = [m.start[1] for m in spray] + [m.end[1] for m in spray]
    assert (min(xs), max(xs)) == pytest.approx((0.0, 0.05))
    assert (min(ys), max(ys)) == pytest.approx((0.0225, 0.0275))


@pytest.mark.parametrize("extra_args", [
    ["--scale", "2.0"],
    ["--no-normalize"],
    ["--origin-x-mm", "5"],
    ["--origin-y-mm", "5"],
])
def test_cli_fit_box_mm_rejects_combination_with_other_placement_flags(
    tmp_path, capsys, monkeypatch, extra_args
):
    svg_path = tmp_path / "square.svg"
    svg_path.write_text(_svg('<path d="M0,0 L10,0 L10,10 L0,10 Z"/>'))

    monkeypatch.setattr(sys, "argv", [
        "svg_to_gcode.py", str(svg_path),
        "--fit-box-mm", "0", "0", "50", "50",
        *extra_args,
    ])
    with pytest.raises(SystemExit):
        s2g.main()
    assert "cannot be combined" in capsys.readouterr().err


def test_cli_fit_box_mm_invalid_box_errors_clearly(tmp_path, capsys, monkeypatch):
    svg_path = tmp_path / "square.svg"
    svg_path.write_text(_svg('<path d="M0,0 L10,0 L10,10 L0,10 Z"/>'))

    monkeypatch.setattr(sys, "argv", [
        "svg_to_gcode.py", str(svg_path),
        "--fit-box-mm", "10", "0", "5", "50",
    ])
    with pytest.raises(SystemExit):
        s2g.main()
    assert "invalid fit box" in capsys.readouterr().err


# --- main(): --closed-loop / --return-feed CLI flags ------------------------- #

def test_cli_return_feed_without_closed_loop_errors_clearly(tmp_path, capsys, monkeypatch):
    svg_path = tmp_path / "square.svg"
    svg_path.write_text(_svg('<path d="M0,0 L10,0 L10,10 L0,10 Z"/>'))

    monkeypatch.setattr(sys, "argv", [
        "svg_to_gcode.py", str(svg_path), "--return-feed", "1500",
    ])
    with pytest.raises(SystemExit):
        s2g.main()
    assert "--return-feed requires --closed-loop" in capsys.readouterr().err


def test_cli_closed_loop_alone_uses_default_return_feed(tmp_path, capsys, monkeypatch):
    svg_path = tmp_path / "square.svg"
    svg_path.write_text(_svg('<path d="M0,0 L10,0 L10,10 L0,10 Z"/>'))
    out_path = tmp_path / "square.gcode"

    monkeypatch.setattr(sys, "argv", [
        "svg_to_gcode.py", str(svg_path), "-o", str(out_path), "--closed-loop",
    ])
    s2g.main()
    capsys.readouterr()

    text = out_path.read_text()
    assert f"G0 F{s2g.DEFAULT_RETURN_FEED_MM_MIN:g}" in text.splitlines()


def test_cli_closed_loop_with_return_feed_passes_custom_value_through(
    tmp_path, capsys, monkeypatch
):
    svg_path = tmp_path / "square.svg"
    svg_path.write_text(_svg('<path d="M0,0 L10,0 L10,10 L0,10 Z"/>'))
    out_path = tmp_path / "square.gcode"

    monkeypatch.setattr(sys, "argv", [
        "svg_to_gcode.py", str(svg_path), "-o", str(out_path),
        "--closed-loop", "--return-feed", "1200",
    ])
    s2g.main()
    capsys.readouterr()

    text = out_path.read_text()
    assert "G0 F1200" in text.splitlines()
    assert f"G0 F{s2g.DEFAULT_RETURN_FEED_MM_MIN:g}" not in text.splitlines()


def test_cli_without_closed_loop_omits_return_pass(tmp_path, capsys, monkeypatch):
    svg_path = tmp_path / "square.svg"
    svg_path.write_text(_svg('<path d="M0,0 L10,0 L10,10 L0,10 Z"/>'))
    out_path = tmp_path / "square.gcode"

    monkeypatch.setattr(sys, "argv", [
        "svg_to_gcode.py", str(svg_path), "-o", str(out_path),
    ])
    s2g.main()
    capsys.readouterr()

    text = out_path.read_text()
    assert not any(l.startswith("G0 F") for l in text.splitlines())
