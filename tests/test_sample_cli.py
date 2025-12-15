from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


def _command_for(module: str, exe: str) -> list[str]:
    """
    Prefer the installed console script (true CLI integration), but fall back to `python -m ...`
    so the test still works in minimal environments.
    """
    resolved = shutil.which(exe)
    if resolved is not None:
        return [resolved]
    return [sys.executable, "-m", module]


def test_sample_cli_pipeline(tmp_path: Path) -> None:
    repo_root = Path(__file__).resolve().parents[1]
    matrix_path = repo_root / "data" / "sample" / "data.txt"
    geneidlist_path = repo_root / "data" / "sample" / "geneidlist.txt"

    assert matrix_path.exists()
    assert geneidlist_path.exists()

    # Write into a temp dir to avoid polluting data/output in the repo.
    out_prefix = tmp_path / "Sample"

    env = dict(os.environ)
    env.setdefault("MPLBACKEND", "Agg")

    eeisp_cmd = _command_for("eeisp.cli_eeisp", "eeisp")
    subprocess.run(
        [
            *eeisp_cmd,
            str(matrix_path),
            str(out_prefix),
            "--threCDI",
            "-1.0",
            "--threEEI",
            "-1.0",
            "-p",
            "1",
        ],
        cwd=str(repo_root),
        env=env,
        check=True,
    )

    cdi_txt = tmp_path / "Sample_CDI_score_data_thre-1.0.txt"
    eei_txt = tmp_path / "Sample_EEI_score_data_thre-1.0.txt"
    cdi_deg = tmp_path / "Sample_CDI_degree_distribution.tsv"
    eei_deg = tmp_path / "Sample_EEI_degree_distribution.tsv"
    assert cdi_txt.exists()
    assert eei_txt.exists()
    assert cdi_deg.exists()
    assert eei_deg.exists()

    add_names_cmd = _command_for("eeisp.cli_add_genename_from_geneid", "add-names")
    cdi_named = tmp_path / "Sample_CDI_score_data_thre-1.0.addgenename.txt"
    eei_named = tmp_path / "Sample_EEI_score_data_thre-1.0.addgenename.txt"
    subprocess.run(
        [*add_names_cmd, str(cdi_txt), str(cdi_named), str(geneidlist_path)],
        cwd=str(repo_root),
        env=env,
        check=True,
    )
    subprocess.run(
        [*add_names_cmd, str(eei_txt), str(eei_named), str(geneidlist_path)],
        cwd=str(repo_root),
        env=env,
        check=True,
    )

    for out in (cdi_named, eei_named):
        assert out.exists()
        first_line = out.read_text(encoding="utf-8").splitlines()[0]
        assert len(first_line.split("\t")) == 7


