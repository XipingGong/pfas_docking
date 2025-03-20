import argparse
import numpy as np

def extract_column_values(filename, column_indices):
    values = [[] for _ in column_indices]
    with open(filename, 'r') as file:
        for line in file:
            columns = line.strip().split()
            for i, col_idx in enumerate(column_indices):
                if len(columns) > col_idx:  # Ensure column exists
                    try:
                        values[i].append(float(columns[col_idx]))
                    except ValueError:
                        pass  # Skip invalid entries
    return values

def compute_histogram(values):
    all_histograms = []
    bin_edges = np.arange(0, 0.01 + np.max(values), 0.01)
    bin_mids = (bin_edges[:-1] + bin_edges[1:]) / 2  # Compute midpoints of bins
    for value_set in values:
        hist, _ = np.histogram(value_set, bins=bin_edges, density=True)  # Use density=True for normalized histogram
        all_histograms.append(hist)
    return bin_mids, all_histograms

def main():
    parser = argparse.ArgumentParser(description="Extract multiple column values and compute histogram bins")
    parser.add_argument("filename", type=str, help="Path to the input text file")
    parser.add_argument("--columns", type=str, required=True, help="Indices of the columns to analyze (e.g., '3 4 5')")
    args = parser.parse_args()

    column_indices = list(map(int, args.columns.split()))
    values = extract_column_values(args.filename, column_indices)
    
    if all(values):
        bin_mids, hist_values = compute_histogram(values)
        print("# x_bins", " ".join([f"x_values_col{col}" for col in column_indices]))  # Comment line
        for i in range(len(hist_values[0])):
            hist_values_str = " ".join(str(hist_values[j][i]) for j in range(len(column_indices)))
            print(f"{bin_mids[i]:.6f} {hist_values_str}")
    else:
        print("No valid numerical values found in the specified columns.")

if __name__ == "__main__":
    main()

