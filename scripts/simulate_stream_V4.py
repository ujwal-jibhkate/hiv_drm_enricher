import argparse
import sys
import json
import pysam

# --- Configuration & Helper Functions ---
GENE_COORDINATES = {"PR": {"start": 2253, "end": 2550}, "RT": {"start": 2550, "end": 4230}}
CODON_TABLE = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

def load_drm_database_v4(path: str) -> dict:
    print(f"INFO: Loading V4 DRM database from '{path}'...")
    with open(path, 'r') as f:
        db_str_keys = json.load(f)
        return {gene: {int(pos): data for pos, data in muts.items()} for gene, muts in db_str_keys.items()}

def load_reference_fasta(path: str) -> dict:
    sequences = {}
    current_seq_name = ""
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                current_seq_name = line[1:].split()[0]
                sequences[current_seq_name] = ""
            elif current_seq_name:
                sequences[current_seq_name] += line.strip()
    return sequences

def stream_bam_as_reads(path: str):
    with pysam.AlignmentFile(path, "rb") as bamfile:
        for read in bamfile:
            yield read

def get_amino_acid_change(pos: int, base: str, ref: str) -> tuple | None:
    for gene, coords in GENE_COORDINATES.items():
        if coords['start'] - 1 <= pos < coords['end']:
            gene_start = coords['start'] - 1
            pos_in_gene = pos - gene_start
            codon_start = pos_in_gene - (pos_in_gene % 3)
            ref_codon_str = ref[gene_start + codon_start : gene_start + codon_start + 3]
            if len(ref_codon_str) < 3: return None
            mut_codon = list(ref_codon_str)
            mut_codon[pos_in_gene % 3] = base
            ref_aa = CODON_TABLE.get(ref_codon_str, '?')
            mut_aa = CODON_TABLE.get("".join(mut_codon), '?')
            aa_pos = (codon_start // 3) + 1
            if ref_aa != mut_aa:
                return ref_aa, mut_aa, aa_pos, gene
    return None

def classify_read_final(read, drm_db, ref_sequences) -> str:
    if read.is_unmapped: return "Non-HIV"
    ref_seq = ref_sequences.get(read.reference_name)
    if not ref_seq and len(ref_sequences) == 1:
        _, ref_seq = next(iter(ref_sequences.items()))
    elif not ref_seq: return "Reference-Not-Found"
    if read.mapping_quality < 20: return "Low-Quality"
    
    aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=True)
    # --- THIS IS THE CORRECTED LINE ---
    for read_pos, ref_pos, ref_base in aligned_pairs:
        if read_pos is not None and ref_pos is not None and ref_base is not None:
            read_base = read.query_sequence[read_pos]
            if read_base.upper() != ref_base.upper():
                aa_change = get_amino_acid_change(ref_pos, read_base.upper(), ref_seq)
                if aa_change:
                    ref_aa, mut_aa, aa_pos, gene = aa_change
                    if gene in drm_db and aa_pos in drm_db[gene] and mut_aa in drm_db[gene][aa_pos]:
                        print(f"\nProcessing Read ID: {read.query_name}")
                        print(f"  └─> MATCH! Found DRM {ref_aa}{aa_pos}{mut_aa} in Gene '{gene}'.")
                        return "DRM-Positive"
    return "DRM-Negative"

def main():
    parser = argparse.ArgumentParser(description="Final HIV DRM Simulation (V4).")
    parser.add_argument("-i", "--input_bam", required=True)
    parser.add_argument("--ref_fasta", required=True)
    parser.add_argument("--db_json", required=True)
    parser.add_argument("-n", "--num_reads", type=int, default=1000)
    args = parser.parse_args()

    print("--- Running Final HIV DRM Simulation (V4) ---")
    drm_db = load_drm_database_v4(args.db_json)
    ref_sequences = load_reference_fasta(args.ref_fasta)
    print("-" * 50)

    classifications = {"DRM-Positive": 0, "DRM-Negative": 0, "Low-Quality": 0, "Non-HIV": 0, "Reference-Not-Found": 0}
    read_iterator = stream_bam_as_reads(args.input_bam)
    for i, read in enumerate(read_iterator):
        if i >= args.num_reads: break
        classification = classify_read_final(read, drm_db, ref_sequences)
        classifications[classification] += 1
    
    print("\n" + "="*25 + " FINAL REPORT " + "="*25)
    print(f"Total Reads Processed: {sum(classifications.values())}")
    for category, count in classifications.items():
        print(f"  - {category}: {count}")
    print("="*64)

if __name__ == "__main__":
    main()