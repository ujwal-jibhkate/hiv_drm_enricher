import pysam
import json

# --- Configuration ---
MAPQ_THRESHOLD = 20

# --- Standard Genetic Code (Codon -> Amino Acid) ---
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def load_drm_database(json_path: str) -> dict:
    """Loads the DRM database from the specified JSON file."""
    print(f"INFO: Loading DRM database from '{json_path}'...")
    try:
        with open(json_path, 'r') as f:
            # The JSON keys are strings, so we must convert them to integers
            db_str_keys = json.load(f)
            return {int(k): v for k, v in db_str_keys.items()}
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"FATAL: Could not load or parse DRM database. Error: {e}", file=sys.stderr)
        raise SystemExit(1)

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

def get_amino_acid_change(mismatch_pos: int, new_base: str, ref_seq: str) -> tuple[str, str, int] | None:
    """
    Translates a nucleotide mismatch into an amino acid change.
    
    Returns:
        A tuple of (original_aa, new_aa, aa_position), or None if translation fails.
    """
    # 1. Determine the position within the codon (0, 1, or 2)
    codon_pos = mismatch_pos % 3
    
    # 2. Find the start of the codon
    codon_start = mismatch_pos - codon_pos
    
    # 3. Extract the original codon from the reference
    ref_codon = ref_seq[codon_start : codon_start + 3]
    if len(ref_codon) < 3: return None # Can't translate incomplete codon

    # 4. Construct the new, mutated codon
    mut_codon_list = list(ref_codon)
    mut_codon_list[codon_pos] = new_base
    mut_codon = "".join(mut_codon_list)

    # 5. Translate both codons to amino acids
    ref_aa = CODON_TABLE.get(ref_codon.upper(), '?')
    mut_aa = CODON_TABLE.get(mut_codon.upper(), '?')
    
    # 6. Calculate the amino acid position (1-based)
    aa_position = (codon_start // 3) + 1

    # Only return a result if the amino acid has actually changed
    if ref_aa != mut_aa:
        return ref_aa, mut_aa, aa_position
    return None

def classify_read_verbose(read: pysam.AlignedSegment, drm_db: dict, ref_seq: str) -> str:
    """
    Executes the full, biologically-aware classification logic for a single read.
    """
    print(f"\nProcessing Read ID: {read.query_name}")
    
    if read.is_unmapped:
        print("  └─ Step 1: FAIL - Read is unmapped.")
        return "Non-HIV"
    print("  ├─ Step 1: PASS - Read is mapped.")

    mapq = read.mapping_quality
    if mapq < MAPQ_THRESHOLD:
        print(f"  └─ Step 2: FAIL - MAPQ ({mapq}) is below threshold ({MAPQ_THRESHOLD}).")
        return "Low-Quality"
    print(f"  ├─ Step 2: PASS - MAPQ ({mapq}) is acceptable.")

    mismatches = parse_mismatches(read)
    if not mismatches:
        print("  └─ Step 3: INFO - No nucleotide mismatches found.")
        return "DRM-Negative"
    print(f"  ├─ Step 3: INFO - Found {len(mismatches)} mismatch(es): {mismatches}")

    print("  └─ Step 4: Translating mismatches and checking against DRM database...")
    for pos, base in mismatches.items():
        # Translate the nucleotide change to an amino acid change
        aa_change_info = get_amino_acid_change(pos, base, ref_seq)
        
        if aa_change_info:
            ref_aa, mut_aa, aa_pos = aa_change_info
            print(f"    ├─> Mismatch at nt {pos+1} ('{base}') causes AA change: {ref_aa}{aa_pos}{mut_aa}")

            # Now, look up this change in our loaded database
            # The database keys are 0-based nucleotide positions
            if pos in drm_db and mut_aa in drm_db[pos]:
                print(f"    └─> MATCH! AA change {mut_aa} at nt position {pos+1} is a known DRM.")
                return "DRM-Positive"

    print("      └─> No matches found in DRM database.")
    return "DRM-Negative"