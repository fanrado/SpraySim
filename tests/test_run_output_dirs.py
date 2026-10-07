"""Per-run output folders: output/<run-name>/ (default output/tmp/)."""

import os
import stat
import subprocess
import sys
from pathlib import Path

from spraysim import storage

REPO_ROOT = Path(__file__).resolve().parent.parent


def test_run_directory_defaults_to_tmp():
    assert storage.run_directory(None) == Path("output") / "tmp"
    assert storage.run_directory("") == Path("output") / "tmp"
    assert storage.run_directory("   ") == Path("output") / "tmp"


def test_run_directory_uses_name_and_root():
    assert storage.run_directory("trial_3") == Path("output") / "trial_3"
    assert storage.run_directory("x", "/elsewhere") == Path("/elsewhere") / "x"


def test_resolve_output_path_places_bare_names_in_run_dir():
    run_dir = Path("output/trial_3")
    assert storage.resolve_output_path("spray_data.npz", run_dir) == run_dir / "spray_data.npz"
    # An explicit directory wins.
    assert storage.resolve_output_path("custom/fig.png", run_dir) == Path("custom/fig.png")
    assert storage.resolve_output_path("/abs/fig.png", run_dir) == Path("/abs/fig.png")


def _run_py(tmp_path, *extra):
    cmd = [sys.executable, "run.py", "--droplets", "30", "--no-plot",
           "--output-root", str(tmp_path), *extra]
    return subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True, timeout=120)


def test_run_py_writes_into_tmp_folder_by_default(tmp_path):
    result = _run_py(tmp_path)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "tmp" / "spray_data.npz").exists()
    assert "Run folder:" in result.stdout


def test_run_py_writes_into_named_folder(tmp_path):
    result = _run_py(tmp_path, "--run-name", "trial_3", "--data", "mine.npz")
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "trial_3" / "mine.npz").exists()
    assert not (tmp_path / "tmp").exists()


def test_run_py_explicit_data_path_is_respected(tmp_path):
    target = tmp_path / "elsewhere" / "d.npz"
    result = _run_py(tmp_path, "--run-name", "ignored", "--data", str(target))
    assert result.returncode == 0, result.stderr
    assert target.exists()
    assert not (tmp_path / "ignored").exists()


def _fake_python(tmp_path):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    fake = bin_dir / "python"
    fake.write_text("#!/bin/sh\necho FAKE_PYTHON_INVOKED \"$@\"\n")
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    return bin_dir


def _main_sh(tmp_path, *args):
    env = os.environ.copy()
    env["PYTHON"] = str(_fake_python(tmp_path) / "python")
    return subprocess.run(["bash", "main.sh", *args], cwd=REPO_ROOT, env=env,
                          capture_output=True, text=True, timeout=15)


def test_main_sh_omits_run_name_when_config_leaves_it_empty(tmp_path):
    out = _main_sh(tmp_path).stdout
    assert "FAKE_PYTHON_INVOKED run.py" in out
    assert "--run-name" not in out
    assert "--out spray_summary.png" in out
    assert "--data spray_data.npz" in out


def test_main_sh_passes_run_name_from_config(tmp_path):
    conf = tmp_path / "named.conf"
    conf.write_text("RUN_NAME=trial_7\nNO_PLOT=true\n")
    out = _main_sh(tmp_path, str(conf)).stdout
    assert "--run-name trial_7" in out


def test_main_sh_passes_run_name_from_command_line(tmp_path):
    out = _main_sh(tmp_path, "default", "--run-name", "cli_run").stdout
    assert "--run-name cli_run" in out
