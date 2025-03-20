#!/usr/bin/env python3
import argparse
import requests
from datetime import datetime
import sys

def parse_date(date_str):
    """
    Parse an ISO datetime string and return it formatted as YYYY-MM-DD.

    Args:
        date_str (str): The date string from the API, e.g. "2011-08-10T00:00:00+0000".

    Returns:
        str: The formatted date "YYYY-MM-DD" or the original string if parsing fails.
    """
    try:
        dt = datetime.strptime(date_str, "%Y-%m-%dT%H:%M:%S%z")
        return dt.strftime("%Y-%m-%d")
    except ValueError:
        return date_str

def get_pdb_dates(pdb_id):
    """
    Query the RCSB PDB API to retrieve deposition and initial release dates for a given PDB ID.

    Args:
        pdb_id (str): The PDB identifier.

    Returns:
        tuple: (deposition_date, release_date) formatted as YYYY-MM-DD.

    Raises:
        Exception: If the API call fails or the expected data is not found.
    """
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.lower()}"
    response = requests.get(url)
    try:
        response.raise_for_status()  # Raise an error for bad HTTP responses
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
        description="Query the RCSB PDB API to get deposition and release dates for a given PDB ID."
    )
    parser.add_argument("pdbid", help="The PDB ID for which to retrieve the dates")
    args = parser.parse_args()

    try:
        deposition_date, release_date = get_pdb_dates(args.pdbid)
        print(f"PDBID: {args.pdbid.upper()} ; Deposited date: {deposition_date} ;")
        print(f"PDBID: {args.pdbid.upper()} ; Initially released date: {release_date} ;")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

