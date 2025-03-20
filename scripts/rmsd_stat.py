import numpy as np
import pandas as pd
import argparse
import re

# Parse command-line arguments
def get_arguments():
    parser = argparse.ArgumentParser(description="Process data with a range for specified columns.")
    parser.add_argument("filename", type=str, help="Path to the data file")
    parser.add_argument("--range", type=str, default="[0.0, 0.2]", help="Range [a, b] or [a, b) to filter values (default: [0.0, 0.2])")
    parser.add_argument("--columns", type=str, required=True, help="Specify columns in the format '4 and 5', '7', or '7 or 9'")
    return parser.parse_args()

args = get_arguments()
filename = args.filename

# Parse range input with inclusive/exclusive handling
range_match = re.match(r"(\[|\()\s*(-?\d*\.?\d+)\s*,\s*(-?\d*\.?\d+)\s*(\]|\))", args.range)
if not range_match:
    print("Error: Invalid range format. Use '[a, b]' or '[a, b)'.")
    exit(1)

range_min, range_max = float(range_match.group(2)), float(range_match.group(3))
inclusive_min = range_match.group(1) == "["
inclusive_max = range_match.group(4) == "]"

# Parse columns input and determine logic
if " and " in args.columns:
    logical_op = "and"
    column_parts = args.columns.replace(" and ", " ").split()
elif " or " in args.columns:
    logical_op = "or"
    column_parts = args.columns.replace(" or ", " ").split()
else:
    logical_op = "single"
    column_parts = [args.columns]

column_indices = [int(col) - 1 for col in column_parts]  # Convert to 0-based index

# Load the data from the text file
data = []
binary_data = []

# Read the file and process specified columns
try:
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            selected_values = []
            for col in column_indices:
                if col >= len(parts):
                    print(f"Error: Column {col + 1} out of range in line: {line.strip()}")
                    exit(1)
                try:
                    value = float(parts[col])
                    selected_values.append(value)
                except ValueError:
                    print(f"Error: Non-numeric value '{parts[col]}' found in column {col + 1} in line: {line.strip()}")
                    exit(1)
            
            data.append(selected_values)
            
            # Apply logical AND/OR conditions
            if logical_op == "and":
                binary_values = all(
                    (inclusive_min and range_min <= val or not inclusive_min and range_min < val) and \
                    (inclusive_max and val <= range_max or not inclusive_max and val < range_max)
                    for val in selected_values
                )
            elif logical_op == "or":
                binary_values = any(
                    (inclusive_min and range_min <= val or not inclusive_min and range_min < val) and \
                    (inclusive_max and val <= range_max or not inclusive_max and val < range_max)
                    for val in selected_values
                )
            else:
                binary_values = (inclusive_min and range_min <= selected_values[0] or not inclusive_min and range_min < selected_values[0]) and \
                                (inclusive_max and selected_values[0] <= range_max or not inclusive_max and selected_values[0] < range_max)

            binary_data.append(1 if binary_values else 0)
except FileNotFoundError:
    print(f"Error: File '{filename}' not found.")
    exit(1)

# Convert data to NumPy arrays
if data:
    data_array = np.array(data, dtype=float)
    binary_array = np.array(binary_data, dtype=int)
else:
    data_array = np.array([], dtype=float)
    binary_array = np.array([], dtype=int)

# Compute statistics based on binary values
stats = {
    f"Proportion in range {args.range}": np.mean(binary_array) if binary_array.size > 0 else 0.0,
    f"Count (total = {binary_array.shape[0]})": np.sum(binary_array) if binary_array.size > 0 else 0
}

# Perform bootstrap resampling to compute 95% confidence intervals
if binary_array.size > 0:
    num_resamples = 10000
    bootstrap_samples = np.random.choice(binary_array, size=(num_resamples, len(binary_array)), replace=True)
    bootstrap_means = np.mean(bootstrap_samples, axis=1)
    lower_bound = np.percentile(bootstrap_means, 2.5)
    upper_bound = np.percentile(bootstrap_means, 97.5)
else:
    lower_bound, upper_bound = 0.0, 0.0  # If no data, set confidence intervals to zero

# Convert results into a DataFrame for display
stats_df = pd.DataFrame(stats, index=[f"Column {args.columns}"])
stats_df["95% CI Lower"], stats_df["95% CI Upper"] = lower_bound, upper_bound

# Print the DataFrame
print(stats_df.to_string(index=True))

