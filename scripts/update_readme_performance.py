#!/usr/bin/env python3
# /// script
# dependencies = [
#   "seaborn",
#   "matplotlib",
#   "pandas",
# ]
# ///

import io
import json
import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

REPO_ROOT = Path(__file__).parent.parent.resolve()
HISTORY_FILE = REPO_ROOT / "benchmark_results" / "history.json"
README_FILE = REPO_ROOT / "README.md"
PYPROJECT_FILE = REPO_ROOT / "pyproject.toml"


def run(cmd: list[str]) -> str:
    return subprocess.check_output(cmd, text=True, cwd=REPO_ROOT).strip()


def get_tags() -> dict[str, str]:
    """Returns a mapping of commit hash to tag name."""
    tags = {}
    try:
        lines = run(["git", "show-ref", "--tags"]).split("\n")
        for line in lines:
            if not line:
                continue
            h, ref = line.split()
            tag = ref.split("/")[-1]
            tags[h[:7]] = tag
    except subprocess.CalledProcessError:
        pass
    return tags


def get_current_version() -> str:
    content = PYPROJECT_FILE.read_text()
    match = re.search(r'version\s*=\s*"([^"]+)"', content)
    return match.group(1) if match else "unknown"


def generate_svg(data_points: list[dict], mode: str = "light") -> str:
    # Prepare data for plotting
    df = pd.DataFrame(data_points)

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
    df["Identity %"] = (df["identity"] / df["total"]) * 100
    df["Analogous %"] = (df["analogous"] / df["total"]) * 100

    versions = df["version"].unique()
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
        tool_df = df[df["tool"] == tool]

        # Calculate offsets for grouped bars
        offset = (i - 0.5) * width

        # Identity bar
        _rects1 = ax.bar(
            [pos + offset for pos in x],
            tool_df["Identity %"],
            width,
            label=f"{tool} Identity",
            color=colors[(tool, "Identity")],
            edgecolor="#ffffff" if mode == "light" else "#0d1117",
        )

        # Analogous bar (stacked)
        bottom = tool_df["Identity %"].values
        _rects2 = ax.bar(
            [pos + offset for pos in x],
            tool_df["Analogous %"],
            width,
            bottom=bottom,
            label=f"{tool} Analogous",
            color=colors[(tool, "Analogous")],
            edgecolor="#ffffff" if mode == "light" else "#0d1117",
        )

    # Customize labels and title
    plt.title("Protein Projection Performance (100k ClinVar Variants)", fontsize=18, pad=20)
    plt.xlabel("Release", fontsize=14)
    plt.ylabel("Match %", fontsize=14)
    plt.ylim(85, 100)  # Zoom in on the high performance range
    plt.xticks(x, versions)

    # Legend cleanup
    legend = plt.legend(title=None, bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.0)
    if mode == "dark":
        plt.setp(legend.get_texts(), color="#e6edf3")

    plt.tight_layout()

    # Save to SVG string
    img_data = io.StringIO()
    plt.savefig(img_data, format="svg", bbox_inches="tight", transparent=True)
    plt.close()

    svg_val = img_data.getvalue()
    return svg_val[svg_val.find("<svg") :]


def update_readme() -> None:
    if not HISTORY_FILE.exists():
        print(f"Error: {HISTORY_FILE} not found.")
        return

    with open(HISTORY_FILE, encoding="utf-8") as f:
        history = json.load(f)

    tags = get_tags()
    current_version = get_current_version()

    seen_versions = set()
    data_points = []

    for entry in history:
        commit = entry["commit"][:7]
        version = tags.get(commit)

        if not version and entry == history[0] and current_version not in tags.values():
            version = f"{current_version} (dev)"

        if version and version not in seen_versions:
            seen_versions.add(version)

            # Add Weaver data point
            data_points.append(
                {
                    "version": version,
                    "tool": "Weaver",
                    "identity": entry.get("w_identity", entry.get("p_match", 0)),
                    "analogous": entry.get("w_analogous", 0),
                    "total": entry["total"],
                },
            )

            # Add Ref data point
            data_points.append(
                {
                    "version": version,
                    "tool": "Ref-HGVS",
                    "identity": entry.get("ref_identity", 0),
                    "analogous": entry.get("ref_analogous", 0),
                    "total": entry["total"],
                },
            )

    data_points.reverse()

    if not data_points:
        print("No valid data points found for graph.")
        return

    # Generate light and dark versions
    svg_light = generate_svg(data_points, mode="light")
    svg_dark = generate_svg(data_points, mode="dark")

    # Save to files
    (REPO_ROOT / "benchmark_results" / "performance_light.svg").write_text(svg_light)
    (REPO_ROOT / "benchmark_results" / "performance_dark.svg").write_text(svg_dark)

    # Inject into README using <picture> for theme awareness
    content = README_FILE.read_text()
    start_marker = "<!-- PERFORMANCE_GRAPH_START -->"
    end_marker = "<!-- PERFORMANCE_GRAPH_END -->"

    svg_tag = """
<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="benchmark_results/performance_dark.svg">
    <source media="(prefers-color-scheme: light)" srcset="benchmark_results/performance_light.svg">
    <img alt="Performance Graph" src="benchmark_results/performance_light.svg" width="800">
  </picture>
</p>
"""

    pattern = re.compile(f"{start_marker}.*?{end_marker}", re.DOTALL)
    replacement = f"{start_marker}\n{svg_tag}\n{end_marker}"

    if start_marker in content and end_marker in content:
        new_content = pattern.sub(replacement, content)
        README_FILE.write_text(new_content)
        print("README.md updated with theme-aware performance graph.")
    else:
        print("Error: Performance markers not found in README.md")


if __name__ == "__main__":
    update_readme()
