import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import TypedDict

# Commits to benchmark (approx 12 selected from 42)


class Stats(TypedDict, total=False):
    total: int
    p_match: int
    spdi_match: int
    p_perc: float
    spdi_perc: float
    commit: str


COMMITS = [
    "08f9ef2",
    "8ed0107",
    "59292bd",
    "b9f618f",
    "1d04b5c",
    "8e92cf2",
    "4b164cc",
    "028489c",
    "b6888ea",
    "070cb95",
    "a551cff",
    "01372cc",
    "204362e",
    "e225ee6",
]

REPO_ROOT = Path(__file__).parent.parent.resolve()
WORKTREE_DIR = REPO_ROOT.parent / "hgvs-rs-bench-tmp"
INPUT_FILE = REPO_ROOT / "clinvar_baseline_validation.tsv"
RESULTS_DIR = REPO_ROOT / "benchmark_results"
RESULTS_DIR.mkdir(exist_ok=True)

# Data files to link into worktree
DATA_FILES = [
    "GRCh38_latest_genomic.gff.gz",
    "GCF_000001405.40_GRCh38.p14_genomic.fna",
    "clinvar_baseline_validation.tsv",
]


def run(
    cmd: list[str],
    cwd: Path = REPO_ROOT,
    check: bool = True,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    print(f"Running: {' '.join(cmd)}")
    res = subprocess.run(cmd, check=False, capture_output=True, text=True, cwd=cwd, env=env)  # noqa: S603
    if check and res.returncode != 0:
        print(f"ERROR: Command failed with code {res.returncode}")
        print(f"STDOUT: {res.stdout}")
        print(f"STDERR: {res.stderr}")
        raise subprocess.CalledProcessError(res.returncode, cmd, res.stdout, res.stderr)
    return res


def setup_worktree() -> None:
    if WORKTREE_DIR.exists():
        print(f"Cleaning up existing worktree at {WORKTREE_DIR}...")
        run(["git", "worktree", "remove", "-f", str(WORKTREE_DIR)])
        shutil.rmtree(WORKTREE_DIR, ignore_errors=True)

    print(f"Creating worktree at {WORKTREE_DIR}...")
    run(["git", "worktree", "add", "--detach", str(WORKTREE_DIR), "HEAD"])


def link_data(target_dir: Path) -> None:
    for f in DATA_FILES:
        src = REPO_ROOT / f
        dst = target_dir / f
        if src.exists() and not dst.exists():
            print(f"Linking {f}...")
            os.link(src, dst)

    # Symlink .venv to make maturin happy
    venv_src = REPO_ROOT / ".venv"
    venv_dst = target_dir / ".venv"
    if venv_src.exists() and not venv_dst.exists():
        print("Symlinking .venv...")
        os.symlink(venv_src, venv_dst)


def build_weaver(target_dir: Path) -> None:
    print(f"Building weaver in {target_dir}...")
    # Use the venv in the root
    venv_bin = REPO_ROOT / ".venv" / "bin"
    maturin_exe = str(venv_bin / "maturin")
    if not os.path.exists(maturin_exe):
        maturin_exe = "maturin"  # Fallback

    env = os.environ.copy()
    env["VIRTUAL_ENV"] = str(REPO_ROOT / ".venv")
    env["PATH"] = f"{venv_bin}:{env.get('PATH', '')}"

    run([maturin_exe, "develop"], cwd=target_dir, env=env)


def validate(target_dir: Path, commit_hash: str, max_variants: int = 100000) -> Path | None:
    output_file = RESULTS_DIR / f"results_{commit_hash}.tsv"
    print(f"Validating {commit_hash} in worktree (max={max_variants})...")

    venv_bin = REPO_ROOT / ".venv" / "bin"
    python_exe = str(venv_bin / "python")
    if not os.path.exists(python_exe):
        python_exe = sys.executable

    env = os.environ.copy()
    env["VIRTUAL_ENV"] = str(REPO_ROOT / ".venv")
    env["PATH"] = f"{venv_bin}:{env.get('PATH', '')}"

    cmd = [
        str(python_exe),
        "-m",
        "weaver.cli.validate",
        "clinvar_baseline_validation.tsv",
        "--max-variants",
        str(max_variants),
        "--output-file",
        str(output_file),
        "--workers",
        "8",
    ]
    try:
        run(cmd, cwd=target_dir, env=env)
        return output_file
    except subprocess.CalledProcessError:
        print(f"Validation failed for {commit_hash}")
        return None


def standardize_p(s: str) -> str:
    if not s:
        return ""
    # Remove accession prefix
    if ":" in s:
        s = s.split(":")[-1]
    # Remove parentheses
    s = s.replace("(", "").replace(")", "")
    # Standardize Ter/*
    return s.replace("Ter", "*")


def analyze(results_file: Path | None) -> Stats | None:
    if not results_file or not results_file.exists():
        return None

    stats: Stats = {"total": 0, "p_match": 0, "spdi_match": 0}
    try:
        with open(results_file, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                stats["total"] += 1
                # Protein Match normalization
                gp = standardize_p(row.get("variant_prot", ""))
                rp = standardize_p(row.get("rs_p", ""))
                if gp and gp == rp:
                    stats["p_match"] += 1

                # SPDI Match
                if row.get("spdi") and row.get("spdi") == row.get("rs_spdi"):
                    stats["spdi_match"] += 1
    except (OSError, csv.Error) as e:
        print(f"Analysis failed: {e}")
        return None

    if stats["total"] > 0:
        stats["p_perc"] = (stats["p_match"] / stats["total"]) * 100
        stats["spdi_perc"] = (stats["spdi_match"] / stats["total"]) * 100
    else:
        stats["p_perc"] = stats["spdi_perc"] = 0.0
    return stats


def main() -> None:
    history: list[Stats] = []
    try:
        setup_worktree()
        link_data(WORKTREE_DIR)

        for i, commit in enumerate(COMMITS):
            print(f"\n=== Benchmark {i + 1}/{len(COMMITS)}: {commit} ===")
            run(["git", "checkout", commit], cwd=WORKTREE_DIR)

            try:
                build_weaver(WORKTREE_DIR)
            except subprocess.CalledProcessError:
                print(f"Skipping {commit} due to build failure.")
                continue

            # Always use 100k set
            count = 100000
            res_file = validate(WORKTREE_DIR, commit, max_variants=count)
            stats = analyze(res_file)
            if stats:
                stats["commit"] = commit
                history.append(stats)
                print(f"Result: P={stats['p_perc']:.2f}%, SPDI={stats['spdi_perc']:.2f}%")

            with open(RESULTS_DIR / "history.json", "w") as f:
                json.dump(history, f, indent=2)

    finally:
        if WORKTREE_DIR.exists():
            print("Cleaning up worktree...")
            run(["git", "worktree", "remove", "-f", str(WORKTREE_DIR)])
            shutil.rmtree(WORKTREE_DIR, ignore_errors=True)


if __name__ == "__main__":
    main()
