import pandas as pd
import json
import sys
import re
import argparse
import pysam

# --- FINAL CONFIGURATION & DATABASE BUILDER ---
GENE_COORDINATES = {
    "PR": {"start": 1385, "end": 1682},
    "RT": {"start": 1683, "end": 3363}
}
CODON_TABLE = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

def build_definitive_database(csv_path: str) -> dict:
    print(f"INFO: Building definitive database from '{csv_path}'...")
    df = pd.read_csv(csv_path, dtype=str)
    df_long = pd.melt(df, id_vars=['mutations'], var_name='Drug', value_name='Score')
    df_long = df_long[df_long['Score'] != 'S'].dropna(subset=['Score'])
    drm_db = {"PR": {}, "RT": {}}
    mutation_pattern = re.compile(r'([A-Z]?)(\d+)([A-Z])', re.IGNORECASE)

    for _, row in df_long.iterrows():
        cleaned_str = re.sub(r"['\"\[\]\s]", "", str(row['mutations']))
        individual_mutations = cleaned_str.split(',')
        for mut in individual_mutations:
            if not mut: continue
            match = mutation_pattern.match(mut)
            if match:
                _, position_str, alt_aa = match.groups()
                position = int(position_str)
                alt_aa = alt_aa.upper()
                gene = "PR" if 1 <= position <= 99 else "RT" if 100 <= position <= 560 else None
                if not gene: continue
                drug_info = {'name': row['Drug'], 'score': row['Score']}
                if position not in drm_db[gene]:
                    drm_db[gene][position] = {}
                if alt_aa not in drm_db[gene][position]:
                    drm_db[gene][position][alt_aa] = {'drugs': []}
                drm_db[gene][position][alt_aa]['drugs'].append(drug_info)
    print("INFO: Database build complete.")
    return drm_db

# --- FINAL SIMULATION LOGIC ---

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
    for read_pos, ref_pos, ref_base in aligned_pairs:
        if read_pos is not None and ref_pos is not None and ref_base is not None:
            read_base = read.query_sequence[read_pos]
            if read_base.upper() != ref_base.upper():
                aa_change = get_amino_acid_change(ref_pos, read_base.upper(), ref_seq)
                if aa_change:
                    ref_aa, mut_aa, aa_pos, gene = aa_change
                    if gene in drm_db and aa_pos in drm_db[gene] and mut_aa in drm_db[gene][aa_pos]:
                        print(f"\n  └─> MATCH! Found DRM {ref_aa}{aa_pos}{mut_aa} in Gene '{gene}' on Read '{read.query_name}'.")
                        return "DRM-Positive"
    return "DRM-Negative"

def main():
    parser = argparse.ArgumentParser(description="Final, All-in-One HIV DRM Simulation Pipeline.")
    parser.add_argument("-i", "--input_bam", required=True)
    parser.add_argument("--ref_fasta", required=True)
    parser.add_argument("--profile_csv", required=True, help="Path to the synthetic_profiles.csv file.")
    parser.add_argument("-n", "--num_reads", type=int, default=5000)
    args = parser.parse_args()

    # --- Run Everything in One Go ---
    drm_db = build_definitive_database(args.profile_csv)
    ref_sequences = load_reference_fasta(args.ref_fasta)
    print("-" * 50)
    print("--- Running Final Simulation ---")

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
    if classifications["DRM-Positive"] == 0:
        print("\nNOTE: No DRM-Positive reads were found. This is expected if the input BAM file does not contain reads covering the PR or RT gene regions.")

if __name__ == "__main__":
    main()