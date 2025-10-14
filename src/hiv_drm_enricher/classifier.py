import pysam

# --- Configuration ---
MAPQ_THRESHOLD = 20

# --- Placeholder DRM Database ---
# This remains the same as before.
DRM_DATABASE = {
    244: { 
        'A': { 
            'gene': 'protease',
            'aa_change': 'D30N',
            'drugs': [{'name': 'Nelfinavir', 'level': 5}]
        }
    },
    254: {
        'A': {
            'gene': 'protease',
            'aa_change': 'I85V',
            'drugs': [{'name': 'Atazanavir', 'level': 2}]
        }
    }
}

# This function remains the same as it's a helper utility.
def parse_mismatches(read: pysam.AlignedSegment) -> dict[int, str]:
    """Identifies all single-base mismatches from an aligned read."""
    mismatches = {}
    aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=True)
    for read_pos, ref_pos, ref_base in aligned_pairs:
        if read_pos is not None and ref_pos is not None and ref_base is not None:
            read_base = read.query_sequence[read_pos]
            if read_base.upper() != ref_base.upper():
                mismatches[ref_pos] = read_base.upper()
    return mismatches


def classify_read_verbose(read: pysam.AlignedSegment) -> str:
    """
    Executes and explains the full classification logic for a single read.

    This version prints its decision-making process at each step for clarity.

    Args:
        read (pysam.AlignedSegment): A single read from the BAM file.

    Returns:
        str: The final classification label for the read.
    """
    print(f"\nProcessing Read ID: {read.query_name}")

    # Step 1: Alignment Validation
    if read.is_unmapped:
        print("  └─ Step 1: FAIL - Read is unmapped.")
        return "Non-HIV"
    print("  ├─ Step 1: PASS - Read is mapped to the reference.")

    # Step 2: Mapping Quality Check
    mapq = read.mapping_quality
    if mapq < MAPQ_THRESHOLD:
        print(f"  └─ Step 2: FAIL - MAPQ ({mapq}) is below threshold ({MAPQ_THRESHOLD}).")
        return "Low-Quality"
    print(f"  ├─ Step 2: PASS - MAPQ ({mapq}) is acceptable.")

    # Step 3: Variant Identification
    mismatches = parse_mismatches(read)
    if not mismatches:
        print("  └─ Step 3: INFO - No mismatches found.")
        return "DRM-Negative"
    print(f"  ├─ Step 3: INFO - Found {len(mismatches)} mismatch(es): {mismatches}")

    # Step 4: DRM Lookup
    print("  └─ Step 4: Checking mismatches against DRM database...")
    for position, base in mismatches.items():
        # Check if the position and base combination matches a known DRM
        if position in DRM_DATABASE and base in DRM_DATABASE[position]:
            print(f"    └─> MATCH! Position {position+1} with base '{base}' is a known DRM.")
            return "DRM-Positive"
    
    # If the loop finishes without finding a match
    print("      └─> No matches found in DRM database.")
    return "DRM-Negative"