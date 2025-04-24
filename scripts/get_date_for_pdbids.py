#!/usr/bin/env python3
import argparse
import requests
from datetime import datetime
import sys
import os

def parse_date(date_str):
    try:
        dt = datetime.strptime(date_str, "%Y-%m-%dT%H:%M:%S%z")
        return dt.strftime("%Y-%m-%d")
    except ValueError:
        return date_str

def get_pdb_dates(pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.lower()}"
    response = requests.get(url)
    try:
        response.raise_for_status()
    except requests.exceptions.HTTPError as err:
        raise Exception(f"HTTP error occurred: {err}")

    data = response.json()

    try:
        deposition_date_raw = data["rcsb_accession_info"]["deposit_date"]
        release_date_raw = data["rcsb_accession_info"]["initial_release_date"]
    except KeyError as e:
        raise Exception("Expected date information not found in the API response.") from e

    deposition_date = parse_date(deposition_date_raw)
    release_date = parse_date(release_date_raw)

    return deposition_date, release_date

def main():
    parser = argparse.ArgumentParser(
        description="Get deposition and release dates for a PDB ID or a file of IDs."
    )
    parser.add_argument("input", help="A PDB ID or a file containing multiple PDB IDs (one per line)")
    args = parser.parse_args()

    # Read one or more PDB IDs
    if os.path.isfile(args.input):
        with open(args.input, "r") as f:
            pdb_ids = [line.strip() for line in f if line.strip()]
    else:
        pdb_ids = [args.input.strip()]

    for pdb_id in pdb_ids:
        try:
            deposition_date, release_date = get_pdb_dates(pdb_id)
            print(f"PDBID: {pdb_id.upper()} ; Deposited date: {deposition_date} ; Initially released date: {release_date} ;")
        except Exception as e:
            print(f"PDBID: {pdb_id.upper()} ; Error: {e}")

if __name__ == "__main__":
    main()

