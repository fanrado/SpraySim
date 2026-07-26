"""Tests for main.sh's interpreter resolution (PYTHON default / override).

These run the real script via subprocess with a fake `python`/`custom-python`
shim on PATH, so no actual simulation runs and nothing is written under
output/ — main.sh's exec still fires, it just execs the shim instead of a
real interpreter.
"""

import os
import stat
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent


def _fake_bin(tmp_path, name):
    """A directory containing one executable `name` that echoes its argv."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    fake = bin_dir / name
    fake.write_text("#!/bin/sh\necho FAKE_PYTHON_INVOKED \"$@\"\n")
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return bin_dir


def _run_main_sh(env):
    return subprocess.run(
        ["bash", "main.sh"],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        timeout=15,
    )


def test_main_sh_defaults_to_python_not_python3(tmp_path):
    """With PYTHON unset, main.sh must invoke `python`, not `python3` (a
    Homebrew python3 on this machine lacks matplotlib; the project's deps
    live in the miniconda `python`)."""
    bin_dir = _fake_bin(tmp_path, "python")
    env = os.environ.copy()
    env.pop("PYTHON", None)
    env["PATH"] = f"{bin_dir}{os.pathsep}{env['PATH']}"

    result = _run_main_sh(env)

    assert "Running: python run.py" in result.stdout
    assert "python3" not in result.stdout.splitlines()[
        next(i for i, l in enumerate(result.stdout.splitlines()) if l.startswith("Running:"))
    ]
    assert "FAKE_PYTHON_INVOKED run.py" in result.stdout


def test_main_sh_python_env_override_still_works(tmp_path):
    """PYTHON=<explicit interpreter> must still take priority over the default."""
    bin_dir = _fake_bin(tmp_path, "custom-python")
    env = os.environ.copy()
    env["PYTHON"] = str(bin_dir / "custom-python")

    result = _run_main_sh(env)

    assert f"Running: {bin_dir / 'custom-python'} run.py" in result.stdout
    assert "FAKE_PYTHON_INVOKED run.py" in result.stdout


@pytest.mark.parametrize("missing_shim", [True])
def test_main_sh_errors_clearly_when_interpreter_missing(tmp_path, missing_shim):
    """Sanity check: if PYTHON points nowhere, main.sh fails at exec (not
    silently) rather than hanging or succeeding without running anything."""
    env = os.environ.copy()
    env["PYTHON"] = str(tmp_path / "does-not-exist")

    result = _run_main_sh(env)

    assert result.returncode != 0
    assert "Running:" in result.stdout  # got past config loading
    assert "FAKE_PYTHON_INVOKED" not in result.stdout
