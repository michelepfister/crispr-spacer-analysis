from typing import List, Dict, Set, Union, Tuple, Any

def parse_spacer_fasta(file_path: str) -> List[Dict[str, Any]]:
    """
    Parse a FASTA file containing spacer sequences and organize them by ID and arrays.
    
    Args:
        file_path (str): Path to the FASTA file
        
    Returns:
        List[Dict[str, Any]]: List of dictionaries with structure:
        [
            {
                "id": "ID",
                "arrays": [
                    {
                        "spacers": ["sequence1", "sequence2", ...]
                    },
                    ...
                ]
            },
            ...
        ]
    """
    with open(file_path, 'r') as file:
        content = file.read().strip()
    
    if not content:
        return []
    
    # Split content into lines and filter out empty lines
    lines = [line.strip() for line in content.split('\n') if line.strip()]
    
    if not lines:
        return []
    
    result: List[Dict[str, Any]] = []
    current_id: Union[str, None] = None
    current_array: List[Dict[str, List[str]]] = []
    current_spacers: List[str] = []
    
    i = 0
    while i < len(lines):
        line = lines[i]
        
        if line.startswith('>'):
            # This is a header line
            header: str = line[1:]  # Remove the '>' character
            
            # Extract ID and spacer info
            if '_spacer' in header:
                # This is a header with ID (e.g., "2203436_spacer1")
                parts: List[str] = header.split('_spacer')
                if len(parts) != 2:
                    raise ValueError(f"Malformed header: {header}")
                
                new_id: str = parts[0]
                spacer_num: str = parts[1]
                
                # If we encounter spacer1 and we already have data, it means a new array
                if spacer_num == '1' and current_spacers:
                    # Save current array
                    if current_spacers:
                        current_array.append({"spacers": current_spacers})
                        current_spacers = []
                
                # If this is a completely new ID
                if new_id != current_id:
                    # Save previous ID's data if it exists
                    if current_id is not None:
                        if current_spacers:
                            current_array.append({"spacers": current_spacers})
                        if current_array:
                            result.append({
                                "id": current_id,
                                "arrays": current_array
                            })
                    
                    # Start new ID
                    current_id = new_id
                    current_array = []
                    current_spacers = []
                
            else:
                # This is a spacer without ID prefix (e.g., "spacer2")
                if not header.startswith('spacer'):
                    raise ValueError(f"Malformed header: {header}")
                
                if current_id is None:
                    raise ValueError(f"Spacer without ID: {header}")
            
            # Get the sequence on the next line
            if i + 1 >= len(lines):
                raise ValueError(f"Header without sequence: {header}")
            
            sequence: str = lines[i + 1].upper()
            if sequence.startswith('>'):
                raise ValueError(f"Expected sequence after header: {header}")
            
            current_spacers.append(sequence)
            i += 2  # Skip the sequence line
            
        else:
            # This shouldn't happen if file is properly formatted
            raise ValueError(f"Unexpected line: {line}")
    
    # Save the last ID's data
    if current_id is not None:
        if current_spacers:
            current_array.append({"spacers": current_spacers})
        if current_array:
            result.append({
                "id": current_id,
                "arrays": current_array
            })
    
    return result

def hamming_distance(seq1: str, seq2: str) -> Union[int, float]:
    """
    Calculate the Hamming distance (number of differing positions) between two sequences.
    
    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence
        
    Returns:
        Union[int, float]: Number of positions where sequences differ, or float('inf') if lengths differ
    """
    if len(seq1) != len(seq2):
        return float('inf')  # Can't compare sequences of different lengths
    
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def group_similar_sequences(sequences: List[str], max_distance: int = 2) -> List[List[str]]:
    """
    Group sequences that have less than max_distance+1 differences between them.
    
    Args:
        sequences (List[str]): List of sequences to group
        max_distance (int): Maximum allowed differences (default: 2, meaning <3 differences)
        
    Returns:
        List[List[str]]: List of groups, where each group is a list of similar sequences
    """
    if not sequences:
        return []
    
    groups: List[List[str]] = []
    used: Set[str] = set()
    
    for i, seq1 in enumerate(sequences):
        if seq1 in used:
            continue
            
        # Start a new group with this sequence
        current_group: List[str] = [seq1]
        used.add(seq1)
        
        # Find all sequences similar to this one
        for j, seq2 in enumerate(sequences):
            if i != j and seq2 not in used:
                distance = hamming_distance(seq1, seq2)
                if isinstance(distance, int) and distance <= max_distance:
                    current_group.append(seq2)
                    used.add(seq2)
        
        groups.append(current_group)
    
    return groups

def collect_spacers_dict(parsed_data: List[Dict[str, Any]]) -> Dict[str, List[str]]:
    """
    Collect spacers into groups of similar sequences, keeping groups that appear at least twice total.
    Group similar sequences first, then filter based on total group frequency.
    
    Args:
        parsed_data (List[Dict[str, Any]]): The parsed FASTA data from parse_spacer_fasta
        
    Returns:
        Dict[str, List[str]]: A mapping of incrementing counter to arrays of similar spacer sequences
        {
            "S001": ["ACGTACC", "ACGTGCC"],  # Each appears once, but group total = 2
            "S002": ["CGGCTTA"],             # Appears twice
            "S003": ["TTAG", "TTGG", "CTAG"] # Different frequencies, but group total >= 2
        }
    """
    # First pass: count frequency of each spacer across all isolates
    spacer_counts: Dict[str, int] = {}
    
    for entry in parsed_data:
        # Collect all spacers from all arrays for this ID
        for array_data in entry["arrays"]:
            for spacer in array_data["spacers"]:
                spacer_counts[spacer] = spacer_counts.get(spacer, 0) + 1
    
    # Get all unique spacers (don't filter yet)
    all_spacers: List[str] = list(spacer_counts.keys())
    
    if not all_spacers:
        return {}
    
    # Group similar sequences (less than 3 differences)
    sequence_groups: List[List[str]] = group_similar_sequences(all_spacers, max_distance=2)
    
    # Now filter groups based on total frequency of the group
    filtered_groups: List[List[str]] = []
    
    for group in sequence_groups:
        # Calculate total frequency for this group
        group_total_frequency = sum(spacer_counts[spacer] for spacer in group)
        
        # Keep group if total frequency is at least 2
        if group_total_frequency >= 2:
            filtered_groups.append(group)
    
    # Create the result dictionary
    spacers_dict: Dict[str, List[str]] = {}
    counter: int = 1
    
    for group in filtered_groups:
        spacers_dict[f"S{counter:03d}"] = group
        counter += 1
    
    return spacers_dict