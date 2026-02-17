#!/usr/bin/env python3
# /// script
# dependencies = [
#   "seaborn",
#   "matplotlib",
#   "pandas",
# ]
# ///

import json
import re
import subprocess
import io
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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


def generate_svg(data_points: list[dict]):
    # Prepare data for plotting
    df = pd.DataFrame(data_points)
    df_melted = df.melt(id_vars=["version"], value_vars=["p", "spdi"], var_name="Metric", value_name="Match %")

    # Rename metrics for better legend
    df_melted["Metric"] = df_melted["Metric"].map({"p": "Protein Match %", "spdi": "SPDI Match %"})

    # Set style
    sns.set_theme(style="whitegrid", context="talk")
    plt.figure(figsize=(10, 5))

    # Create plot
    palette = {"Protein Match %": "#3498db", "SPDI Match %": "#2ecc71"}
    ax = sns.lineplot(data=df_melted, x="version", y="Match %", hue="Metric", marker="o", palette=palette, linewidth=3)

    # Customize labels and title
    plt.title("Performance Trend (100k Variants)", fontsize=16, pad=20)
    plt.xlabel("Version", fontsize=14)
    plt.ylabel("Match %", fontsize=14)
    plt.ylim(80, 100)

    # Adjust legend
    plt.legend(title=None, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    # Save to SVG string
    img_data = io.StringIO()
    plt.savefig(img_data, format="svg", bbox_inches="tight")
    plt.close()

    svg_str = img_data.getvalue()
    # Strip everything before <svg
    svg_str = svg_str[svg_str.find("<svg") :]
    return svg_str


def update_readme():
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

        if not version and entry == history[0]:
            if current_version not in tags.values():
                version = f"{current_version} (dev)"

        if version and version not in seen_versions:
            seen_versions.add(version)
            data_points.append(
                {"version": version, "p": round(entry["p_perc"], 2), "spdi": round(entry["spdi_perc"], 2)}
            )

    data_points.reverse()

    if not data_points:
        print("No valid data points found for graph.")
        return

    svg_content = generate_svg(data_points)

    # Inject into README
    content = README_FILE.read_text()
    start_marker = "<!-- PERFORMANCE_GRAPH_START -->"
    end_marker = "<!-- PERFORMANCE_GRAPH_END -->"

    pattern = re.compile(f"{start_marker}.*?{end_marker}", re.DOTALL)
    replacement = f'{start_marker}\n\n<p align="center">\n{svg_content}\n</p>\n\n{end_marker}'

    if start_marker in content and end_marker in content:
        new_content = pattern.sub(replacement, content)
        README_FILE.write_text(new_content)
        print("README.md updated with SVG performance graph.")
    else:
        print("Error: Performance markers not found in README.md")


if __name__ == "__main__":
    update_readme()
