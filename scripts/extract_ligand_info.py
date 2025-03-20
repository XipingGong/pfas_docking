import argparse

def extract_ligand_info(ccdcif_file, ligand_file):
    # Read ligand IDs from the second file
    with open(ligand_file, 'r') as file:
        ligand_ids = [line.strip() for line in file if line.strip()]  # Store as a list to loop in order
    
    # Read the first file into memory
    with open(ccdcif_file, 'r') as file:
        ccdcif_lines = file.readlines()
    
    # Loop through each ligand ID and extract matching lines
    for ligand_id in ligand_ids:
        for line in ccdcif_lines:
            if f'Ligand ID: {ligand_id} ' in line:
                print(f"{ligand_id} ; {line.strip()}")
                break  # Stop searching once a match is found for this ligand

def main():
    parser = argparse.ArgumentParser(description="Extract ligand info based on ligand IDs from another file")
    parser.add_argument("ccdcif_file", type=str, help="Path to the ccdcif_info.log file")
    parser.add_argument("ligand_file", type=str, help="Path to the file containing ligand IDs")
    args = parser.parse_args()
    
    extract_ligand_info(args.ccdcif_file, args.ligand_file)

if __name__ == "__main__":
    main()

