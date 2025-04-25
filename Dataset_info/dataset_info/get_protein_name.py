#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script reads a protein information file (e.g. "protein_info.txt") where the first
column is the pdb id. For each pdb id, the script retrieves the corresponding protein title
from the RCSB PDB REST API and writes the results (pdb id and protein title) into a new file.
"""

import requests


INFO_FILE = "benchmark_info.txt"  # Input info file containing protein data
OUTPUT_FILE = "pdb_names.txt"  # Output file for pdb id and protein title pairs


# ==============================================================

def get_protein_title(pdb_id):
    """
    Retrieve the protein title for the given pdb id from the RCSB PDB database.

    Args:
        pdb_id (str): The PDB ID (case-insensitive).

    Returns:
        str: The protein title if found, otherwise None.
    """
    base_url = "https://data.rcsb.org/rest/v1/core/entry"
    pdb_id_lower = pdb_id.lower()  # API expects lowercase pdb id
    url = f"{base_url}/{pdb_id_lower}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        # Protein title is located under "struct" -> "title"
        title = data.get("struct", {}).get("title", None)
        return title
    except requests.exceptions.RequestException as e:
        print(f"Error retrieving data for pdb id {pdb_id}: {e}")
        return None


def read_protein_ids(info_file):
    """
    Read protein ids from the info file.
    Assumes the file is a tab or space-separated file where the first column is pdb id.

    Args:
        info_file (str): Path to the input info file.

    Returns:
        list: Sorted list of unique pdb ids.
    """
    protein_ids = set()
    with open(info_file, "r", encoding="utf-8") as f:
        header = f.readline()  # Skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            pdb_id = parts[0]
            protein_ids.add(pdb_id)
    return sorted(protein_ids)


def write_protein_titles(protein_titles, output_file):
    """
    Write the protein titles to an output file.

    Args:
        protein_titles (dict): Dictionary mapping pdb id to protein title.
        output_file (str): Path to the output file.
    """
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("pdb_id\tprotein_title\n")
        for pdb_id, title in protein_titles.items():
            if title is None:
                title = "Not Found"
            f.write(f"{pdb_id}\t{title}\n")


def main():
    # Read protein ids from the info file.
    protein_ids = read_protein_ids(INFO_FILE)
    print(f"Found {len(protein_ids)} unique protein ids.")

    # Retrieve protein title for each pdb id.
    protein_titles = {}
    for pdb_id in protein_ids:
        title = get_protein_title(pdb_id)
        protein_titles[pdb_id] = title
        print(f"{pdb_id}: {title}")

    # Write the results to the output file.
    write_protein_titles(protein_titles, OUTPUT_FILE)
    print(f"Protein titles saved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
