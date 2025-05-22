import glob
import re
import sys
import argparse

KEYS = [
    "Protein Backbone RMSD (Direct)",
    "Protein Backbone RMSD (MDTraj)",
    "Ligand RMSD (Direct)",
    "Ligand RMSD (MDTraj)",
    "Protein Backbone Pocket RMSD (Direct)",
    "Protein Backbone Pocket RMSD (MDTraj)"
]

def extract_min_and_array(filepath, num_per_array):
    results = {}
    for key in KEYS:
        results[key + " Min"] = "None"
        results[key + " Array"] = ["None"] * num_per_array

    try:
        with open(filepath, 'r') as file:
            for line in file:
                for key in KEYS:
                    if key in line:
                        min_match = re.search(r"Min\s*=\s*([\d.]+)", line)
                        array_match = re.search(r"\[([^\]]+)\]", line)

                        if min_match:
                            results[key + " Min"] = min_match.group(1)

                        if array_match:
                            array_values = array_match.group(1).split()
                            padded_array = array_values + ["None"] * (num_per_array - len(array_values))
                            results[key + " Array"] = padded_array[:num_per_array]
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)

    output = [filepath]
    for key in KEYS:
        output.append(results[key + " Min"])
        output.extend(results[key + " Array"])
    return output

def main():
    parser = argparse.ArgumentParser(description="Extract min and array RMSD values from file(s).")
    parser.add_argument("pattern", type=str, help="File pattern to match (e.g., '*/af3_model.rmsd')")
    parser.add_argument("--num_per_array", type=int, default=6, help="Number of values per array (default: 6)")

    args = parser.parse_args()
    files = glob.glob(args.pattern)

    if not files:
        print("No files matched the given pattern.")
        sys.exit(0)

    # Print header
    header = ["#File"]
    for key in KEYS:
        key_base = key.replace(" ", "_").replace("(", "").replace(")", "")
        header.append(f"{key_base}_Min")
        for i in range(1, args.num_per_array + 1):
            header.append(f"{key_base}_Array_{i}")
    print("\t".join(header))

    # Print data
    for f in files:
        result = extract_min_and_array(f, args.num_per_array)
        print("\t".join(result))

if __name__ == "__main__":
    main()

