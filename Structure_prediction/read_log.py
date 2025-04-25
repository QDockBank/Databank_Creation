# --*-- conding:utf-8 --*--
# @time:4/9/25 9:43 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : benchmark.py


"""
This script reads a folder containing both .log and .out files. Each file may contain
multiple fragments for different proteins, with lines like:

--- Fragment 1/25: PDB=6g3a ---
=== Processing protein 6g3a ===
Residue sequence: KIGEGVFGEVFQ
Sequence length: 12
Number of qubits: 82
Iters. done: 1 [Current cost: 14201.490132392033]
...
Top 1 best energy = 164.81302034505325, xyz saved: ...
Finished processing: 6g3a

It extracts these fields for each pdb_id:
1) pdb_id,
2) Residue sequence,
3) Sequence length,
4) Number of qubits,
5) Top 1 best energy,
6) The maximum 'Current cost' observed during the optimization,
7) The range of 'Current cost' (max_cost - min_cost)

Finally, it writes a summary file, one line per pdb_id:

pdb_id    sequence    seq_length    qubits    top1_energy    cost_max    cost_range
"""

import os
import re

LOG_FOLDER = "log"  # Folder containing .log and .out files
OUTPUT_FILE = "log/log_summary.txt"  # Final summary output

# ======================================================

def parse_log_file(filepath):
    """
    Parse a single .log or .out file, which may contain multiple protein fragments.
    Return a list of dictionaries, each containing:
      {
        'pdb_id': str,
        'sequence': str,
        'seq_length': int,
        'qubits': int,
        'top1_energy': float or None,
        'cost_values': list of float
      }
    """
    results = []
    current_data = None  # used to accumulate data for current fragment

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line_s = line.strip()
            if not line_s:
                continue

            # Detect start of a new protein fragment
            # e.g.: === Processing protein 6g3a ===
            m_start = re.match(r'^=== Processing protein (\S+) ===', line_s)
            if m_start:
                # If we already had a current_data for previous protein, finalize it
                if current_data:
                    results.append(current_data)
                pdb_id = m_start.group(1)
                current_data = {
                    'pdb_id': pdb_id,
                    'sequence': "",
                    'seq_length': None,
                    'qubits': None,
                    'top1_energy': None,
                    'cost_values': []
                }
                continue

            # Residue sequence line: "Residue sequence: KIGEGVFGEVFQ"
            m_seq = re.match(r'^Residue sequence:\s+(.+)$', line_s)
            if m_seq and current_data:
                current_data['sequence'] = m_seq.group(1).strip()
                continue

            # Sequence length line: "Sequence length: 12"
            m_len = re.match(r'^Sequence length:\s+(\d+)$', line_s)
            if m_len and current_data:
                current_data['seq_length'] = int(m_len.group(1))
                continue

            # Number of qubits line: "Number of qubits: 82"
            m_qubits = re.match(r'^Number of qubits:\s+(\d+)$', line_s)
            if m_qubits and current_data:
                current_data['qubits'] = int(m_qubits.group(1))
                continue

            # Current cost lines: "Iters. done: N [Current cost: 14201.49]"
            # We'll store these cost values
            m_cost = re.search(r'Current cost:\s*([\-\d\.]+)', line_s)
            if m_cost and current_data:
                try:
                    val = float(m_cost.group(1))
                    current_data['cost_values'].append(val)
                except ValueError:
                    pass
                continue

            # Top 1 best energy line: "Top 1 best energy = 164.81302034505325, xyz saved: 6ugp_top_1.xyz"
            m_top1 = re.match(r'^Top 1 best energy\s*=\s*([\-\d\.]+)', line_s)
            if m_top1 and current_data:
                try:
                    top1_val = float(m_top1.group(1))
                    current_data['top1_energy'] = top1_val
                except ValueError:
                    pass
                continue

            # Finished processing line: "Finished processing: 6g3a"
            # This indicates the end of a protein block
            m_end = re.match(r'^Finished processing:\s*(\S+)', line_s)
            if m_end and current_data:
                # If there's a mismatch, assume it's the correct end
                if m_end.group(1) == current_data['pdb_id']:
                    # finalize
                    results.append(current_data)
                    current_data = None
                else:
                    # If the file is not well-structured, we still finalize anyway
                    results.append(current_data)
                    current_data = None
                continue

    # If file ended but we still have a current_data
    if current_data:
        results.append(current_data)

    return results


def main():
    # We'll collect a dictionary keyed by pdb_id to store the final combined info
    summary_dict = {}

    # Traverse .log or .out files in LOG_FOLDER
    for fname in os.listdir(LOG_FOLDER):
        if not (fname.endswith(".log") or fname.endswith(".out")):
            continue
        filepath = os.path.join(LOG_FOLDER, fname)
        fragments = parse_log_file(filepath)
        # Merge them into summary_dict
        for frag in fragments:
            pdb_id = frag['pdb_id']
            if pdb_id not in summary_dict:
                summary_dict[pdb_id] = frag
            else:
                # If the same pdb_id appears multiple times, decide how to handle conflict
                # For simplicity, we can overwrite or unify data. We'll do an overwrite if needed.
                # You may also want to check if e.g. cost_values are bigger, etc.
                summary_dict[pdb_id] = frag

    # Now we build the final results
    # For each pdb_id, we want: pdb_id, residue_seq, seq_length, qubits, top1_best_energy, cost_max, cost_range
    final_records = []
    for pdb_id in summary_dict:
        data = summary_dict[pdb_id]
        seq = data['sequence']
        seq_len = data['seq_length']
        qubits = data['qubits']
        top1 = data['top1_energy']
        cost_vals = data['cost_values']
        if cost_vals:
            cost_max = max(cost_vals)
            cost_min = min(cost_vals)
            cost_range = cost_max - cost_min
        else:
            cost_max = None
            cost_range = None
        final_records.append((pdb_id, seq, seq_len, qubits, top1, cost_max, cost_range))

    # Sort final records by pdb_id
    final_records.sort(key=lambda x: x[0])

    # Write output
    with open(OUTPUT_FILE, "w", encoding="utf-8") as out_f:
        # header
        out_f.write(
            "pdb_id\tResidue_sequence\tSequence_length\tNumber_of_qubits\tTop1_best_energy\tCost_max\tCost_range\n")
        for (pdb_id, seq, seq_len, qubits, top1, cost_max, cost_range) in final_records:
            # Convert None to empty or something
            seq_len_str = str(seq_len) if seq_len is not None else ""
            qubits_str = str(qubits) if qubits is not None else ""
            top1_str = f"{top1:.6f}" if top1 is not None else ""
            cost_max_str = f"{cost_max:.6f}" if cost_max is not None else ""
            cost_range_str = f"{cost_range:.6f}" if cost_range is not None else ""
            out_line = f"{pdb_id}\t{seq}\t{seq_len_str}\t{qubits_str}\t{top1_str}\t{cost_max_str}\t{cost_range_str}\n"
            out_f.write(out_line)

    print(f"Done! Processed {len(final_records)} proteins. Summary written to {OUTPUT_FILE}.")


if __name__ == "__main__":
    main()
