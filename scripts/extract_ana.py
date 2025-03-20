import re
import argparse

# Define RMSD fields globally so they are accessible in all functions
rmsd_fields = [
    "Protein Backbone RMSD (Direct)",
    "Protein Backbone RMSD (MDTraj)",
    "Ligand RMSD (Direct)",
    "Ligand RMSD (MDTraj)",
    "Protein Pocket RMSD (Direct)",
    "Protein Pocket RMSD (MDTraj)"
]

def extract_data(file_path):
    """
    Extracts relevant data for multiple jobs, including:
    - PDB Released Date
    - Sequence Residues Count
    - Elapsed Time
    - Processing Chain Time
    - Total Residue Number
    - RMSD Values for `<< AF3 >>`, `<< Vina >>`, `<< AF3Pocket-Vina >>`
    """
    job_data = {}
    current_job = None
    current_section = None  # Track the current section (AF3, Vina, AF3Pocket-Vina)

    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            # Detect Job Name
            job_match = re.search(r"Processing:\s*(\S+)", line)
            if job_match:
                current_job = job_match.group(1)
                job_data[current_job] = {
                    "PDB Released Date": None,
                    "Sequence Residues Count": None,
                    "Total Residues": None,
                    "Elapsed Time": None,
                    "Processing Chain Time": None,
                    "AF3": {field: [] for field in rmsd_fields},  # Ensure list for 2 trials
                    "Vina": {field: (None, None) for field in rmsd_fields},
                    "AF3Pocket-Vina": {field: (None, None) for field in rmsd_fields}
                }

            if not current_job:
                continue  # Skip lines before detecting a job

            # Detect section
            section_match = re.search(r"<<\s*(.*?)\s*>>", line)
            if section_match:
                current_section = section_match.group(1).strip()

            # Extract PDB Released Date
            pdb_release_match = re.search(r"Initially released date:\s*([\d-]+)", line)
            if pdb_release_match:
                job_data[current_job]["PDB Released Date"] = pdb_release_match.group(1)

            # Extract Sequence Residues Count
            seq_match = re.search(r'"sequence":\s*"([A-Z]+)"', line)
            if seq_match:
                job_data[current_job]["Sequence Residues Count"] = len(seq_match.group(1))

            # Extract Elapsed Time
            elapsed_match = re.search(r"# Elapsed time:\s*([\d.]+)\s*seconds", line)
            if elapsed_match:
                job_data[current_job]["Elapsed Time"] = elapsed_match.group(1)

            # Extract Processing Chain Time
            processing_match = re.search(r"Processing chain.*?([\d.]+)\s*seconds", line)
            if processing_match:
                job_data[current_job]["Processing Chain Time"] = processing_match.group(1)

            # Extract Total Residues from mdtraj line
            residues_match = re.search(r"<mdtraj\.Trajectory.*?(\d+)\s*residues", line)
            if residues_match:
                job_data[current_job]["Total Residues"] = residues_match.group(1)

            # Extract RMSD values based on section
            if current_section in ["AF3", "Vina", "AF3Pocket-Vina"]:
                for field in rmsd_fields:
                    pattern = rf"{re.escape(field)}:\s*Min\s*=\s*([\d.]+|None)\s*nm\s*;\s*\[([\d.]+|None)"
                    match = re.search(pattern, line)
                    if match:
                        min_value = match.group(1)
                        first_bracket_value = match.group(2)

                        if current_section == "AF3":
                            # Store multiple trials in a list (ensuring 24 values)
                            job_data[current_job]["AF3"][field].append((min_value, first_bracket_value))
                        else:
                            # Store a single tuple for Vina & AF3Pocket-Vina
                            job_data[current_job][current_section][field] = (min_value, first_bracket_value)

    # Ensure AF3 has exactly **two trials** (even if missing, add None)
    for job in job_data:
        for field in rmsd_fields:
            while len(job_data[job]["AF3"][field]) < 2:
                job_data[job]["AF3"][field].append((None, None))  # Fill missing trial

    return job_data


def main():
    parser = argparse.ArgumentParser(description="Extract RMSD data from a job processing log file.")
    parser.add_argument("file_path", help="Path to the text file containing the job data.")
    args = parser.parse_args()

    extracted_data = extract_data(args.file_path)

    # Print header description
    print("# 1: Job Name")
    print("# 2: PDB Released Date")
    print("# 3: Sequence Residues Count")
    print("# 4: Total Residues")
    print("# 5: Elapsed Time")
    print("# 6: Processing Chain Time")
    
    column_number = 7
    for section in ["AF3", "Vina", "AF3Pocket-Vina"]:
        for field in rmsd_fields:
            if section == "AF3":
                print(f"# {column_number}: {section} {field} Trial 1 Min")
                print(f"# {column_number + 1}: {section} {field} Trial 1 First Value")
                print(f"# {column_number + 2}: {section} {field} Trial 2 Min")
                print(f"# {column_number + 3}: {section} {field} Trial 2 First Value")
                column_number += 4
            else:
                print(f"# {column_number}: {section} {field} Min")
                print(f"# {column_number + 1}: {section} {field} First Value")
                column_number += 2

    # Print job data
    for job, data in extracted_data.items():
        row = [
            job,
            data["PDB Released Date"],
            data["Sequence Residues Count"],
            data["Total Residues"],
            data["Elapsed Time"],
            data["Processing Chain Time"]
        ]

        for section in ["AF3", "Vina", "AF3Pocket-Vina"]:
            for field in rmsd_fields:
                if section == "AF3":
                    trial_1 = data[section][field][0] if len(data[section][field]) > 0 else (None, None)
                    trial_2 = data[section][field][1] if len(data[section][field]) > 1 else (None, None)
                    row.extend([trial_1[0], trial_1[1], trial_2[0], trial_2[1]])
                else:
                    row.extend([data[section][field][0], data[section][field][1]])

        print(" ".join(map(str, row)))  # Print values in space-separated format

if __name__ == "__main__":
    main()

