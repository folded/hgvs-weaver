import re
import os


def main():
    if not os.path.exists("tests/data"):
        os.makedirs("tests/data")

    with open("full_regression_data_v2.txt", "r") as f:
        content = f.read()

    # Find all let seq_NM_... = "...";
    matches = re.findall(r'let seq_(\w+) = "([^"]+)";', content)

    for name, seq in matches:
        # name looks like NM_000038_6, convert back to NM_000038.6
        # Actually, let's just use the name as is if that's what we need.
        # Regressions.rs uses include_str!("data/NM_000038.6.seq")
        # So name should be reconstructed or matched.

        # Let's find the // NM_... lines instead
        pass

    # Better approach: split by // NM
    sections = re.split(r"// (NM_\d+\.\d+)", content)
    # sections[0] is header
    # sections[1] is ac, sections[2] is data
    for i in range(1, len(sections), 2):
        ac = sections[i]
        data = sections[i + 1]
        seq_match = re.search(r'let seq_\w+ = "([^"]+)";', data)
        if seq_match:
            seq = seq_match.group(1)
            filename = f"tests/data/{ac}.seq"
            with open(filename, "w") as out:
                out.write(seq)
            print(f"Wrote {filename}")


if __name__ == "__main__":
    main()
