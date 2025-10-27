import argparse
import sys
import pysam

# --- THE CORRECTED GENE COORDINATES for AF033819.3 ---
# Based on public database records for this specific HIV-1 strain.
GENE_COORDINATES = {
    "PR": {"start": 1385, "end": 1682},  # Protease (99 aa * 3 bp/aa = 297 bp)
    "RT": {"start": 1683, "end": 3363}   # Reverse Transcriptase (560 aa * 3 bp/aa = 1680 bp)
}
# ----------------------------------------------------

CODON_TABLE = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

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

def run_diagnostics(read, ref_sequences):
    if read.is_unmapped or read.mapping_quality < 20:
        return
    ref_seq = ref_sequences.get(read.reference_name)
    if not ref_seq and len(ref_sequences) == 1:
        _, ref_seq = next(iter(ref_sequences.items()))
    elif not ref_seq:
        return

    aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=True)
    for read_pos, ref_pos, ref_base in aligned_pairs:
        if read_pos is not None and ref_pos is not None and ref_base is not None:
            read_base = read.query_sequence[read_pos]
            if read_base.upper() != ref_base.upper():
                aa_change = get_amino_acid_change(ref_pos, read_base.upper(), ref_seq)
                if aa_change and aa_change[3] == "RT":
                    ref_aa, mut_aa, aa_pos, gene = aa_change
                    print("\n" + "="*20 + " RT GENE MUTATION FOUND " + "="*20)
                    print(f"  - Read ID:          {read.query_name}")
                    print(f"  - Nucleotide Pos:   {ref_pos + 1}")
                    print(f"  - Gene:             {gene}")
                    print(f"  - Calculated Change:  {ref_aa}{aa_pos}{mut_aa}")
                    print("="*62)

def main():
    parser = argparse.ArgumentParser(description="Final Diagnostic Tool for the HIV DRM pipeline.")
    parser.add_argument("-i", "--input_bam", required=True)
    parser.add_argument("--ref_fasta", required=True)
    parser.add_argument("-n", "--num_reads", type=int, default=1000)
    args = parser.parse_args()

    print("--- Running Final Diagnostic with Corrected Coordinates ---")
    ref_sequences = load_reference_fasta(args.ref_fasta)
    print(f"--- Analyzing the first {args.num_reads} reads for any RT mutations... ---")

    read_iterator = stream_bam_as_reads(args.input_bam)
    for i, read in enumerate(read_iterator):
        if i >= args.num_reads:
            break
        run_diagnostics(read, ref_sequences)
    
    print("\n--- Final Diagnostic Run Complete ---")

if __name__ == "__main__":
    main()