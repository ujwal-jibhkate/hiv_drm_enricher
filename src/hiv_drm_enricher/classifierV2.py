import pysam
import json

GENE_COORDINATES = {
    "PR": {"start": 2253, "end": 2550},
    "RT": {"start": 2550, "end": 4230},
    "IN": {"start": 4230, "end": 5096}
}
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def load_drm_database(json_path: str) -> dict:
    """ This is the corrected function that loads the gene-centric database. """
    print(f"INFO: Loading gene-centric DRM database from '{json_path}'...")
    try:
        with open(json_path, 'r') as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"FATAL: Could not load or parse DRM database. Error: {e}", file=sys.stderr)
        raise SystemExit(1)

def parse_mismatches(read: pysam.AlignedSegment) -> dict[int, str]:
    mismatches = {}
    aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=True)
    for read_pos, ref_pos, ref_base in aligned_pairs:
        if read_pos is not None and ref_pos is not None and ref_base is not None:
            read_base = read.query_sequence[read_pos]
            if read_base.upper() != ref_base.upper():
                mismatches[ref_pos] = read_base.upper()
    return mismatches

def get_amino_acid_change(mismatch_pos: int, new_base: str, ref_seq: str) -> tuple[str, str, int, str] | None:
    for gene, coords in GENE_COORDINATES.items():
        if coords['start'] - 1 <= mismatch_pos < coords['end']:
            gene_start = coords['start'] - 1
            pos_in_gene = mismatch_pos - gene_start
            codon_pos_in_gene = pos_in_gene % 3
            codon_start_in_gene = pos_in_gene - codon_pos_in_gene
            codon_start_abs = gene_start + codon_start_in_gene
            ref_codon = ref_seq[codon_start_abs : codon_start_abs + 3]
            if len(ref_codon) < 3: return None
            mut_codon_list = list(ref_codon)
            mut_codon_list[codon_pos_in_gene] = new_base
            mut_codon = "".join(mut_codon_list)
            ref_aa = CODON_TABLE.get(ref_codon.upper(), '?')
            mut_aa = CODON_TABLE.get(mut_codon.upper(), '?')
            aa_position = (codon_start_in_gene // 3) + 1
            if ref_aa != mut_aa:
                return ref_aa, mut_aa, aa_position, gene
            return None
    return None

def classify_read_verbose(read: pysam.AlignedSegment, drm_db: dict, ref_sequences: dict) -> str:
    """ This is the final classifier that works with the gene-centric DB. """
    if read.is_unmapped: return "Non-HIV"
    ref_name = read.reference_name
    ref_seq = ref_sequences.get(ref_name)
    if not ref_seq:
        if len(ref_sequences) == 1:
            _, ref_seq = next(iter(ref_sequences.items()))
        else: return "Reference-Not-Found"
    if read.mapping_quality < 20: return "Low-Quality"
    mismatches = parse_mismatches(read)
    if not mismatches: return "DRM-Negative"
    
    print(f"\nProcessing Read ID: {read.query_name}")
    print(f"  └─ INFO: Found {len(mismatches)} mismatches. Translating and checking...")

    for pos, base in mismatches.items():
        aa_change_info = get_amino_acid_change(pos, base, ref_seq)
        if aa_change_info:
            ref_aa, mut_aa, aa_pos, gene = aa_change_info
            
            # Construct the string to look up in the database (e.g., "M184V")
            aa_change_str = f"{ref_aa}{aa_pos}{mut_aa}"
            
            # --- The Corrected, Direct Lookup Logic ---
            if gene in drm_db and aa_change_str in drm_db[gene]:
                print(f"    └─> MATCH! Found DRM '{aa_change_str}' in Gene '{gene}'.")
                return "DRM-Positive"

    print("      └─> No known DRMs found in this read.")
    return "DRM-Negative"