#!/usr/bin/env python3

import requests
import argparse
import os

def search_pdb_for_ligand(ligand_id):
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
                        "value": ligand_id.upper()
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

        if response.status_code != 200:
            print(f"❌ HTTP error for ligand {ligand_id}: Status code {response.status_code}")
            return []

        # Check if response body is empty
        if not response.text.strip():
            print(f"❌ Empty response for ligand {ligand_id}")
            return []

        results = response.json()

        if "result_set" in results:
            return [entry.get("identifier", "") if isinstance(entry, dict) else entry for entry in results["result_set"]]
        else:
            return []

    except requests.exceptions.RequestException as err:
        print(f"❌ Network error for ligand {ligand_id}: {err}")
        return []

    except ValueError as err:
        print(f"❌ JSON parsing error for ligand {ligand_id}: {err}")
        return []


def main():
    parser = argparse.ArgumentParser(description="Search PDB for free ligands by ID or list file.")
    parser.add_argument("input", help="Ligand ID or a text file of ligand IDs (one per line)")
    args = parser.parse_args()

    if os.path.isfile(args.input):
        with open(args.input, 'r') as f:
            ligand_ids = [line.strip() for line in f if line.strip()]
    else:
        ligand_ids = [args.input.strip()]

    for ligand_id in ligand_ids:
        pdb_entries = search_pdb_for_ligand(ligand_id)
        if pdb_entries:
            print(f"Ligand ID: {ligand_id.upper()} ; PDBID: {' '.join(pdb_entries)} ;")
        else:
            print(f"Ligand ID: {ligand_id.upper()} ; PDBID: not_available ;")

if __name__ == "__main__":
    main()

