import os

# Define input and output files
input_file = "ana.dat"
filtered_file = "x1.log"
no_none_file = "x2.log"
sequence_filtered_file = "x3.log"
before_set_file = "before_set.dat"
after_set_file = "after_set.dat"

# Define the columns to extract (Convert from 1-based to 0-based indexing)
columns_to_extract = [x - 1 for x in [1, 2, 3, 5, 6, 7, 11, 15, 19, 23, 27, 17, 21, 36, 38, 35, 37, 48, 50, 47, 49]]

# Define sequence length and date thresholds
sequence_length_threshold = 50
date_threshold = "2021-09-30"

def safe_min(a, b):
    """Returns the minimum of two values, handling 'None' or missing values."""
    try:
        return str(min(float(a), float(b)))
    except (ValueError, TypeError):
        return "None"

def filter_comments_and_extract_columns(input_path, output_path, columns):
    """ Step 1: Remove comment lines, extract specific columns, and add computed values. """
    with open(input_path, "r") as infile, open(output_path, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                continue  # Skip comment lines

            columns_split = line.strip().split()
            
            # Extract only selected columns
            selected_columns = [columns_split[i] for i in columns if i < len(columns_split)]

            # Compute additional columns based on minimum values
            min_15_48 = safe_min(columns_split[14], columns_split[47]) if len(columns_split) > 47 else "None"
            min_19_50 = safe_min(columns_split[18], columns_split[49]) if len(columns_split) > 49 else "None"
            min_17_47 = safe_min(columns_split[16], columns_split[46]) if len(columns_split) > 46 else "None"
            min_21_49 = safe_min(columns_split[20], columns_split[48]) if len(columns_split) > 48 else "None"

            # Append computed columns
            selected_columns.extend([min_15_48, min_19_50, min_17_47, min_21_49])

            outfile.write(" ".join(selected_columns) + "\n")

def remove_failed_jobs(input_path, output_path):
    """ Step 2: Remove jobs that contain 'None' """
    with open(input_path, "r") as infile, open(output_path, "w") as outfile:
        outfile.writelines(line for line in infile if "None" not in line)

def filter_by_sequence_length(input_path, output_path, threshold):
    """ Step 3: Remove jobs where sequence length (3rd column) is less than the threshold """
    with open(input_path, "r") as infile, open(output_path, "w") as outfile:
        for line in infile:
            columns = line.strip().split()
            try:
                if int(columns[2]) > threshold:  # Ensure sequence length > threshold
                    outfile.write(line)
            except (IndexError, ValueError):  # Handle cases where conversion fails
                continue

def split_by_date(input_path, before_path, after_path, threshold):
    """ Step 4: Split dataset into 'Before Set' and 'After Set' based on date """
    with open(input_path, "r") as infile, open(before_path, "w") as before_file, open(after_path, "w") as after_file:
        for line in infile:
            columns = line.strip().split()
            if len(columns) > 1:  # Ensure there's a second column
                job_date = columns[1]
                if job_date <= threshold:
                    before_file.write(line)
                else:
                    after_file.write(line)

def main():
    """ Main function to process the dataset """
    print("Step 1: Filtering comments, extracting columns, and computing additional values...")
    filter_comments_and_extract_columns(input_file, filtered_file, columns_to_extract)

    print("Step 2: Removing failed jobs (None values)...")
    remove_failed_jobs(filtered_file, no_none_file)

    print("Step 3: Filtering by sequence length...")
    filter_by_sequence_length(no_none_file, sequence_filtered_file, sequence_length_threshold)

    print("Step 4: Splitting into Before and After sets...")
    split_by_date(sequence_filtered_file, before_set_file, after_set_file, date_threshold)

    print("✅ Processing completed successfully: before_set.dat and after_set.dat saved!")

if __name__ == "__main__":
    main()

