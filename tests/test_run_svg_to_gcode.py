"""Tests for run_svg_to_gcode.sh's argument validation and how it invokes
svg_to_gcode.py's --fit-box-mm CLI.

These run the real script via subprocess with a fake `python` shim on PATH
(same pattern as test_main_launcher.py), so no actual conversion runs and
nothing is written under output/.
"""

import os
import re
import stat
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = "run_svg_to_gcode.sh"


def _fake_bin(tmp_path, name):
    """A directory containing one executable `name` that echoes its argv."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    fake = bin_dir / name
    fake.write_text("#!/bin/sh\necho FAKE_PYTHON_INVOKED \"$@\"\n")
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return bin_dir


def _run(args, env, cwd=REPO_ROOT):
    return subprocess.run(
        ["bash", str(REPO_ROOT / SCRIPT), *args],
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        timeout=15,
    )


def test_missing_argument_errors_with_usage(tmp_path):
    env = os.environ.copy()
    result = _run([], env)

    assert result.returncode == 1
    assert "missing required argument" in result.stderr
    assert f"usage: {REPO_ROOT / SCRIPT} <drawing.svg>" in result.stderr


def test_nonexistent_svg_errors_clearly(tmp_path):
    env = os.environ.copy()
    missing = tmp_path / "does_not_exist.svg"
    result = _run([str(missing)], env)

    assert result.returncode == 1
    assert f"no such file: {missing}" in result.stderr


def test_invokes_svg_to_gcode_with_fit_box_args(tmp_path):
    bin_dir = _fake_bin(tmp_path, "python")
    svg_path = tmp_path / "drawing.svg"
    svg_path.write_text("<svg/>")

    env = os.environ.copy()
    env.pop("PYTHON", None)
    env["PATH"] = f"{bin_dir}{os.pathsep}{env['PATH']}"

    result = _run([str(svg_path)], env)

    assert result.returncode == 0
    assert (
        f"FAKE_PYTHON_INVOKED svg_to_gcode.py {svg_path} "
        "-o output/paths_0_to_12cm.gcode "
        "--fit-box-mm 0 0 120 120"
    ) in result.stdout


def test_default_config_appends_closed_loop(tmp_path):
    # CLOSED_LOOP defaults to true in the committed script (continuous
    # multi-run spraying is the primary use case).
    bin_dir = _fake_bin(tmp_path, "python")
    svg_path = tmp_path / "drawing.svg"
    svg_path.write_text("<svg/>")

    env = os.environ.copy()
    env.pop("PYTHON", None)
    env["PATH"] = f"{bin_dir}{os.pathsep}{env['PATH']}"

    result = _run([str(svg_path)], env)

    assert result.returncode == 0
    assert result.stdout.rstrip("\n").endswith("--closed-loop")
    assert "--return-feed" not in result.stdout


def _copy_with_config(tmp_path, closed_loop, return_feed):
    """A copy of run_svg_to_gcode.sh with CLOSED_LOOP/RETURN_FEED overridden.

    The fake `python` shim never reads svg_to_gcode.py, so running from a
    copy (whose own SCRIPT_DIR the script cd's into) is safe.
    """
    text = (REPO_ROOT / SCRIPT).read_text()
    text = re.sub(r"^CLOSED_LOOP=.*$", f"CLOSED_LOOP={closed_loop}", text, count=1, flags=re.M)
    text = re.sub(r'^RETURN_FEED=".*"$', f'RETURN_FEED="{return_feed}"', text, count=1, flags=re.M)
    copy_path = tmp_path / SCRIPT
    copy_path.write_text(text)
    copy_path.chmod(copy_path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return copy_path


def test_closed_loop_true_appends_flag_without_return_feed(tmp_path):
    bin_dir = _fake_bin(tmp_path, "python")
    svg_path = tmp_path / "drawing.svg"
    svg_path.write_text("<svg/>")
    script = _copy_with_config(tmp_path, closed_loop="true", return_feed="")

    env = os.environ.copy()
    env.pop("PYTHON", None)
    env["PATH"] = f"{bin_dir}{os.pathsep}{env['PATH']}"

    result = subprocess.run(
        ["bash", str(script), str(svg_path)],
        env=env, capture_output=True, text=True, timeout=15,
    )

    assert result.returncode == 0
    assert result.stdout.rstrip("\n").endswith("--closed-loop")
    assert "--return-feed" not in result.stdout


def test_closed_loop_true_with_return_feed_appends_both(tmp_path):
    bin_dir = _fake_bin(tmp_path, "python")
    svg_path = tmp_path / "drawing.svg"
    svg_path.write_text("<svg/>")
    script = _copy_with_config(tmp_path, closed_loop="true", return_feed="1500")

    env = os.environ.copy()
    env.pop("PYTHON", None)
    env["PATH"] = f"{bin_dir}{os.pathsep}{env['PATH']}"

    result = subprocess.run(
        ["bash", str(script), str(svg_path)],
        env=env, capture_output=True, text=True, timeout=15,
    )

    assert result.returncode == 0
    assert result.stdout.rstrip("\n").endswith("--closed-loop --return-feed 1500")


def test_closed_loop_false_omits_return_feed_even_when_set(tmp_path):
    # RETURN_FEED alone (without CLOSED_LOOP=true) must not reach
    # svg_to_gcode.py, which rejects --return-feed without --closed-loop.
    bin_dir = _fake_bin(tmp_path, "python")
    svg_path = tmp_path / "drawing.svg"
    svg_path.write_text("<svg/>")
    script = _copy_with_config(tmp_path, closed_loop="false", return_feed="1500")

    env = os.environ.copy()
    env.pop("PYTHON", None)
    env["PATH"] = f"{bin_dir}{os.pathsep}{env['PATH']}"

    result = subprocess.run(
        ["bash", str(script), str(svg_path)],
        env=env, capture_output=True, text=True, timeout=15,
    )

    assert result.returncode == 0
    assert "--closed-loop" not in result.stdout
    assert "--return-feed" not in result.stdout


def test_python_env_override_used_instead_of_default(tmp_path):
    bin_dir = _fake_bin(tmp_path, "custom-python")
    svg_path = tmp_path / "drawing.svg"
    svg_path.write_text("<svg/>")

    env = os.environ.copy()
    env["PYTHON"] = str(bin_dir / "custom-python")

    result = _run([str(svg_path)], env)

    assert result.returncode == 0
    assert "FAKE_PYTHON_INVOKED svg_to_gcode.py" in result.stdout
