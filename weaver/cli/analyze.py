"""Analysis of HGVS validation results with integrated contingency reporting."""

import argparse
import csv
import io
import json
import re
import subprocess
from pathlib import Path


def clean_hgvs(s: str) -> str:
    if not s:
        return ""
    # Remove accession prefix
    if ":" in s:
        s = s.split(":")[-1]
    # Remove parentheses
    s = s.replace("(", "").replace(")", "")
    # Standardize Ter/*
    return s.replace("Ter", "*")


def is_p_match(pred: str, truth: str) -> bool:
    """Checks if predicted protein matches truth."""
    if not pred or pred.startswith("ERR:"):
        return False
    p = clean_hgvs(pred)
    t = clean_hgvs(truth)
    return p == t


def get_repo_root() -> Path:
    return Path(__file__).parent.parent.parent.resolve()


def get_tags(repo_root: Path) -> dict[str, str]:
    """Returns a mapping of commit hash to tag name."""
    tags = {}
    try:
        lines = (
            subprocess.check_output(  # noqa: S603
                ["git", "show-ref", "--tags"],
                text=True,
                cwd=repo_root,
                shell=False,
            )
            .strip()
            .split("\n")
        )
        for line in lines:
            if not line:
                continue
            h, ref = line.split()
            tag = ref.split("/")[-1]
            tags[h[:7]] = tag
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass
    return tags


def get_current_version(repo_root: Path) -> str:
    pyproject_file = repo_root / "pyproject.toml"
    if not pyproject_file.exists():
        return "unknown"
    content = pyproject_file.read_text()
    match = re.search(r'version\s*=\s*"([^"]+)"', content)
    return match.group(1) if match else "unknown"


def generate_svg(data_points: list[dict], mode: str = "light") -> str:
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns

    # Prepare data for plotting
    results_df = pd.DataFrame(data_points)

    # Set style based on mode
    if mode == "dark":
        # GitHub dark mode colors: bg=#0d1117, grid=#30363d
        sns.set_theme(
            style="darkgrid",
            context="talk",
            rc={
                "axes.facecolor": "#0d1117",
                "figure.facecolor": "#0d1117",
                "text.color": "#e6edf3",
                "axes.labelcolor": "#e6edf3",
                "xtick.color": "#e6edf3",
                "ytick.color": "#e6edf3",
                "grid.color": "#30363d",
                "patch.edgecolor": "#30363d",
            },
        )
    else:
        sns.set_theme(style="whitegrid", context="talk")

    _fig, ax = plt.subplots(figsize=(12, 6))

    # Standardize data: Identity and Analogous as percentages
    results_df["Identity %"] = (results_df["identity"] / results_df["total"]) * 100
    results_df["Analogous %"] = (results_df["analogous"] / results_df["total"]) * 100

    versions = results_df["version"].unique()
    tools = ["Weaver", "Ref-HGVS"]

    x = range(len(versions))
    width = 0.35  # width of bars

    # Colors
    # Blue for Weaver, Green for Ref
    if mode == "dark":
        colors = {
            ("Weaver", "Identity"): "#2980b9",
            ("Weaver", "Analogous"): "#5dade2",
            ("Ref-HGVS", "Identity"): "#27ae60",
            ("Ref-HGVS", "Analogous"): "#52be80",
        }
    else:
        colors = {
            ("Weaver", "Identity"): "#3498db",
            ("Weaver", "Analogous"): "#85c1e9",
            ("Ref-HGVS", "Identity"): "#27ae60",
            ("Ref-HGVS", "Analogous"): "#7dcea0",
        }

    for i, tool in enumerate(tools):
        tool_df = results_df[results_df["tool"] == tool]
        offset = (i - 0.5) * width

        ax.bar(
            [pos + offset for pos in x],
            tool_df["Identity %"],
            width,
            label=f"{tool} Identity",
            color=colors[(tool, "Identity")],
            edgecolor="#ffffff" if mode == "light" else "#0d1117",
        )

        bottom = tool_df["Identity %"].values
        ax.bar(
            [pos + offset for pos in x],
            tool_df["Analogous %"],
            width,
            bottom=bottom,
            label=f"{tool} Analogous",
            color=colors[(tool, "Analogous")],
            edgecolor="#ffffff" if mode == "light" else "#0d1117",
        )

    plt.title("Protein Projection Performance (100k ClinVar Variants)", fontsize=18, pad=20)
    plt.xlabel("Release", fontsize=14)
    plt.ylabel("Match %", fontsize=14)
    plt.ylim(85, 100)
    plt.xticks(x, versions)

    legend = plt.legend(title=None, bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.0)
    if mode == "dark":
        plt.setp(legend.get_texts(), color="#e6edf3")

    plt.tight_layout()
    img_data = io.StringIO()
    plt.savefig(img_data, format="svg", bbox_inches="tight", transparent=True)
    plt.close()

    svg_val = img_data.getvalue()
    return svg_val[svg_val.find("<svg") :]


def update_performance_graphs(repo_root: Path) -> None:
    history_file = repo_root / "benchmark_results" / "history.json"
    readme_file = repo_root / "README.md"

    if not history_file.exists():
        print(f"Warning: {history_file} not found. Skipping graph update.")
        return

    try:
        import pandas as pd  # noqa: F401
    except ImportError:
        print("Warning: pandas/seaborn/matplotlib not found. Skipping graph update.")
        return

    with open(history_file, encoding="utf-8") as f:
        history = json.load(f)

    tags = get_tags(repo_root)
    current_version = get_current_version(repo_root)

    seen_versions = set()
    data_points = []

    for entry in history:
        commit = entry["commit"][:7]
        version = tags.get(commit)

        if not version and entry == history[0] and current_version not in tags.values():
            version = f"{current_version} (dev)"

        if version and version not in seen_versions:
            seen_versions.add(version)
            data_points.append(
                {
                    "version": version,
                    "tool": "Weaver",
                    "identity": entry.get("w_identity", entry.get("p_match", 0)),
                    "analogous": entry.get("w_analogous", 0),
                    "total": entry["total"],
                }
            )
            data_points.append(
                {
                    "version": version,
                    "tool": "Ref-HGVS",
                    "identity": entry.get("ref_identity", 0),
                    "analogous": entry.get("ref_analogous", 0),
                    "total": entry["total"],
                }
            )

    data_points.reverse()
    if not data_points:
        return

    svg_light = generate_svg(data_points, mode="light")
    svg_dark = generate_svg(data_points, mode="dark")

    (repo_root / "benchmark_results" / "performance_light.svg").write_text(svg_light)
    (repo_root / "benchmark_results" / "performance_dark.svg").write_text(svg_dark)

    if readme_file.exists():
        content = readme_file.read_text()
        start_marker = "<!-- PERFORMANCE_GRAPH_START -->"
        end_marker = "<!-- PERFORMANCE_GRAPH_END -->"

        svg_tag = f"""{start_marker}
<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="benchmark_results/performance_dark.svg">
    <source media="(prefers-color-scheme: light)" srcset="benchmark_results/performance_light.svg">
    <img alt="Performance Graph" src="benchmark_results/performance_light.svg" width="800">
  </picture>
</p>
{end_marker}"""

        pattern = re.compile(f"{start_marker}.*?{end_marker}", re.DOTALL)
        if start_marker in content and end_marker in content:
            new_content = pattern.sub(svg_tag, content)
            readme_file.write_text(new_content)
            print("README.md updated with performance graphs.")


def main() -> None:
    """Main analysis entry point."""
    parser = argparse.ArgumentParser(description="Analyze full HGVS validation results.")
    parser.add_argument("input_file", help="Input validation TSV file.")
    parser.add_argument(
        "--update-readme",
        action="store_true",
        help="Update the project README.md with the latest results and performance graphs.",
    )
    args = parser.parse_args()

    total = 0
    rs_p_match = 0
    ref_p_match = 0
    rs_spdi_match = 0
    ref_spdi_match = 0

    rs_parse_err = 0
    ref_parse_err = 0
    rs_ref_mismatch = 0

    p_stats = {"both": 0, "rs_only": 0, "ref_only": 0, "neither": 0}
    spdi_stats = {"both": 0, "rs_only": 0, "ref_only": 0, "neither": 0}

    with open(args.input_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            total += 1
            cv_p = row["variant_prot"]
            cv_spdi = row["spdi"]
            rs_p_raw = row["rs_p"]
            ref_p_raw = row["ref_p"]

            if rs_p_raw.startswith("ERR:RefMismatch"):
                rs_ref_mismatch += 1
            if rs_p_raw.startswith("ERR:Parse"):
                rs_parse_err += 1
            if ref_p_raw.startswith("ERR:Parse"):
                ref_parse_err += 1

            rs_p_ok = is_p_match(rs_p_raw, cv_p)
            ref_p_ok = is_p_match(ref_p_raw, cv_p)

            if rs_p_ok:
                rs_p_match += 1
            if ref_p_ok:
                ref_p_match += 1

            if rs_p_ok and ref_p_ok:
                p_stats["both"] += 1
            elif rs_p_ok:
                p_stats["rs_only"] += 1
            elif ref_p_ok:
                p_stats["ref_only"] += 1
            else:
                p_stats["neither"] += 1

            rs_spdi_ok = row["rs_spdi"] == cv_spdi
            ref_spdi_ok = row["ref_spdi"] == cv_spdi

            if rs_spdi_ok:
                rs_spdi_match += 1
            if ref_spdi_ok:
                ref_spdi_match += 1

            if rs_spdi_ok and ref_spdi_ok:
                spdi_stats["both"] += 1
            elif rs_spdi_ok:
                spdi_stats["rs_only"] += 1
            elif ref_spdi_ok:
                spdi_stats["ref_only"] += 1
            else:
                spdi_stats["neither"] += 1

    if total == 0:
        print("No variants processed.")
        return

    rs_p_pct = rs_p_match / total * 100
    ref_p_pct = ref_p_match / total * 100
    rs_spdi_pct = rs_spdi_match / total * 100
    ref_spdi_pct = ref_spdi_match / total * 100

    rs_analogous_pct = rs_spdi_pct - rs_p_pct
    ref_analogous_pct = ref_spdi_pct - ref_p_pct

    def fmt_pct(pct: float, bold: bool = False) -> str:
        s = f"{pct:.3f}%"
        return f"**{s}**" if bold else s

    rs_p_str = fmt_pct(rs_p_pct, rs_p_pct > ref_p_pct)
    ref_p_str = fmt_pct(ref_p_pct, ref_p_pct > rs_p_pct)

    rs_ana_str = fmt_pct(rs_analogous_pct, rs_analogous_pct > ref_analogous_pct)
    ref_ana_str = fmt_pct(ref_analogous_pct, ref_analogous_pct > rs_analogous_pct)

    rs_spdi_str = fmt_pct(rs_spdi_pct, rs_spdi_pct > ref_spdi_pct)
    ref_spdi_str = fmt_pct(ref_spdi_pct, ref_spdi_pct > rs_spdi_pct)

    rs_err_str = f"{rs_parse_err:,}"
    ref_err_str = f"{ref_parse_err:,}"
    if rs_parse_err < ref_parse_err:
        rs_err_str = f"**{rs_err_str}**"
    elif ref_parse_err < rs_parse_err:
        ref_err_str = f"**{ref_err_str}**"

    report = [
        f"### Validation Results ({total:,} variants)",
        "",
        "Summary of results comparing `weaver` and `ref-hgvs` against ClinVar ground truth:",
        "",
        "| Implementation | Protein Identity | Protein Analogous | SPDI (Genomic) | Parse Errors |",
        "| :------------- | :--------------: | :---------------: | :------------: | :----------: |",
        f"| weaver         |  {rs_p_str}  | {rs_ana_str} | {rs_spdi_str} | {rs_err_str} |",
        f"| ref-hgvs       |  {ref_p_str}  | {ref_ana_str} | {ref_spdi_str} | {ref_err_str} |",
        "",
        f"RefSeq Data Mismatches: {rs_ref_mismatch:,} ({rs_ref_mismatch / total * 100:.1f}%)",
        "",
        "#### Protein Translation Agreement",
        "",
        "|                     | ref-hgvs Match | ref-hgvs Mismatch |",
        "| :------------------ | :------------: | :---------------: |",
        f"| **weaver Match**    |     {p_stats['both']:,}     |     {p_stats['rs_only']:,}     |",
        f"| **weaver Mismatch** |     {p_stats['ref_only']:,}     |     {p_stats['neither']:,}     |",
        "",
        "#### SPDI Mapping Agreement",
        "",
        "|                     | ref-hgvs Match | ref-hgvs Mismatch |",
        "| :------------------ | :------------: | :---------------: |",
        f"| **weaver Match**    |     {spdi_stats['both']:,}     |     {spdi_stats['rs_only']:,}     |",
        f"| **weaver Mismatch** |     {spdi_stats['ref_only']:,}     |     {spdi_stats['neither']:,}     |",
    ]

    out_text = "\n".join(report)
    print(out_text)

    if args.update_readme:
        repo_root = get_repo_root()
        readme_path = repo_root / "README.md"
        if readme_path.exists():
            content = readme_path.read_text()
            pattern = re.compile(
                r"### Validation Results \(.*?\).*?(?=\n- \*\*Variant Equivalence\*\*)",
                re.DOTALL,
            )
            if pattern.search(content):
                new_content = pattern.sub(out_text.replace("\\", "\\\\"), content)
                readme_path.write_text(new_content)
                print(f"\n[Updated {readme_path}]")
                update_performance_graphs(repo_root)
            else:
                print("\n[Error: Could not find Validation Results section in README.md]")
        else:
            print("\n[Error: README.md not found]")


if __name__ == "__main__":
    main()
