#!/usr/bin/env python3

import requests
import sys
import argparse

def search_pdb_for_ligand(ligand_id):
    """
    Searches the RCSB PDB database for PDB entries that contain a given ligand as a free ligand.

    Parameters:
        ligand_id (str): The ligand ID to search for in the PDB database.

    Returns:
        list: A list of PDB IDs containing the specified ligand as a free ligand.
    """
    url = "https://search.rcsb.org/rcsbsearch/v2/query"

    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_nonpolymer_instance_annotation.comp_id",
                        "operator": "exact_match",
                        "value": ligand_id.upper()  # Ensure ligand ID is uppercase
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_nonpolymer_instance_annotation.type",
                        "operator": "exact_match",
                        "value": "HAS_NO_COVALENT_LINKAGE"
                    }
                }
            ]
        },
        "return_type": "entry",
        "request_options": {
            "results_verbosity": "compact"
        }
    }

    try:
        response = requests.post(url, json=query)
        response.raise_for_status()  # Raises an HTTPError for bad responses
        results = response.json()

        if "result_set" in results and isinstance(results["result_set"], list):
            return results["result_set"]  # Return the list of PDB IDs
        else:
            return []

    except requests.exceptions.HTTPError as http_err:
        print(f"❌ HTTP error occurred: {http_err}")
    except requests.exceptions.RequestException as req_err:
        print(f"❌ Request error occurred: {req_err}")
    except Exception as err:
        print(f"❌ An error occurred: {err}")

    return []

def main():
    """
    Main execution function.

    Usage:
        python search_pdb_by_ligand_id.py <ligand_id>

    Arguments:
        ligand_id : The ligand ID to search for in the PDB database.

    Output:
        Prints the list of PDB entries that contain the specified ligand as a free ligand.
    """
    parser = argparse.ArgumentParser(description="Search for PDB entries containing a given ligand as a free ligand.")
    parser.add_argument("ligand_id", type=str, help="Ligand ID to search in the PDB database.")
    args = parser.parse_args()

    pdb_entries = search_pdb_for_ligand(args.ligand_id.upper())

    if pdb_entries:
        print(f"Ligand ID: {args.ligand_id.upper()} ; PDBID: {' '.join(pdb_entries)} ;")
    else:
        print(f"Ligand ID: {args.ligand_id.upper()} ; PDBID: not_available ;")

if __name__ == "__main__":
    main()
