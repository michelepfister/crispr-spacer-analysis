import csv
import sys
from typing import Dict, List, Set, Any
from parser import parse_spacer_fasta, collect_spacers_dict

def write_csv_dict(parsed_data: List[Dict[str, Any]], spacers_dict: Dict[str, List[str]], output_file: str) -> None:
    """
    Write a CSV with spacer IDs and their corresponding sequences.
    Includes both grouped spacers and ungrouped sequences from parsed_data.

    Args:
        parsed_data (List[Dict[str, Any]]): The parsed FASTA data
        spacers_dict (Dict[str, List[str]]): Dictionary mapping counter IDs to grouped spacer sequences
        output_file (str): Output file path, or None to print to stdout
    """
    # Create header row
    header = ["Id", "Sequence"]

    # Prepare rows - sort by ID numerically
    rows = []
    
    # First, add all sequences from spacers_dict
    # Sort spacer IDs numerically and format as S001, S002, etc.
    spacer_ids = sorted(spacers_dict.keys(), key=lambda x: int(x[1:]))
    spacer_ids = [f"S{int(sid[1:]):03d}" for sid in spacer_ids]
    
    for spacer_id in spacer_ids:
        sequences = spacers_dict[spacer_id]
        rows.append([spacer_id, "; ".join(sequences)])

    # Collect all sequences that are in spacers_dict
    grouped_sequences: Set[str] = set()
    for sequence_list in spacers_dict.values():
        grouped_sequences.update(sequence_list)
    
    # Find sequences in parsed_data that are NOT in spacers_dict
    ungrouped_sequences: Set[str] = set()
    for entry in parsed_data:
        for array_data in entry["arrays"]:
            for spacer in array_data["spacers"]:
                if spacer not in grouped_sequences:
                    ungrouped_sequences.add(spacer)
    
    # Add ungrouped sequences with U001, U002, etc.
    ungrouped_counter = 1
    for sequence in sorted(ungrouped_sequences):  # Sort for consistent output
        u_id = f"U{ungrouped_counter:03d}"  # Format as U001, U002, etc.
        rows.append([u_id, sequence])
        ungrouped_counter += 1

    # Write CSV
    if output_file:
        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(header)
            writer.writerows(rows)
        print(f"Spacers CSV written to: {output_file}")
    else:
        # Print to stdout
        writer = csv.writer(sys.stdout)
        writer.writerow(header)
        writer.writerows(rows)


def process(args):
    # Parse the FASTA file
    parsed_data = parse_spacer_fasta(args.file_path)
    # Get the dict of id + spacers
    spacers_dict = collect_spacers_dict(parsed_data)

    write_csv_result(parsed_data, spacers_dict, args.output + ".csv")
    write_csv_dict(parsed_data, spacers_dict, args.output + "_definition.csv")

def write_csv_result(parsed_data: List[Dict[str, Any]], spacers_dict: Dict[str, List[str]], output_file: str = None) -> None:
    """
    Write a CSV with presence/absence matrix of spacers.
    
    Args:
        parsed_data (List[Dict[str, Any]]): The parsed FASTA data
        spacers_dict (Dict[str, List[str]]): Dictionary mapping counter IDs to spacer sequences
        output_file (str): Output file path, or None to print to stdout
    """
    spacer_ids = sorted(spacers_dict.keys(), key=lambda x: int(x[1:]))  # Sort numerically

    # Create header row
    header = ["Name/Id"] + spacer_ids
    
    # Prepare rows
    rows = []
    
    for entry in parsed_data:
        row = [entry["id"]]  # Start with the ID
        
        # Collect all spacers from all arrays for this entry
        entry_spacers = set()
        for array_data in entry["arrays"]:
            entry_spacers.update(array_data["spacers"])
        
        # Check presence/absence for each spacer in spacers_dict
        for spacer_id in spacer_ids:
            spacer_group = spacers_dict[spacer_id]
            
            # Check if any spacer from the entry matches any sequence in this group
            is_present = False
            
            for entry_spacer in entry_spacers:
                if entry_spacer in spacer_group:
                    is_present = True
                    break
            
            row.append("x" if is_present else "-")
        
        rows.append(row)
    
    # Write CSV
    if output_file:
        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(header)
            writer.writerows(rows)
        print(f"CSV written to: {output_file}")
    else:
        # Print to stdout
        writer = csv.writer(sys.stdout)
        writer.writerow(header)
        writer.writerows(rows)