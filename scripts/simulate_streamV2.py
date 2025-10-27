import argparse
import sys
sys.path.append('src')

from hiv_drm_enricher.streaming import stream_bam_as_reads
# These imports are now correct for the final version
from hiv_drm_enricher.classifierV2 import classify_read_verbose, load_drm_database

def load_reference_fasta(fasta_path: str) -> dict[str, str]:
    sequences = {}
    current_seq_name = ""
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_seq_name = line[1:].split()[0]
                sequences[current_seq_name] = ""
            elif current_seq_name:
                sequences[current_seq_name] += line
    return sequences

def main():
    parser = argparse.ArgumentParser(description="Final HIV DRM Simulation Pipeline.")
    parser.add_argument("-i", "--input_bam", required=True, help="Path to input BAM file.")
    parser.add_argument("--ref_fasta", required=True, help="Path to reference FASTA file.")
    parser.add_argument("--db_json", required=True, help="Path to the final DRM database JSON file.")
    parser.add_argument("-n", "--num_reads", type=int, default=500, help="Number of reads to process.")
    args = parser.parse_args()

    print("--- Running Final HIV DRM Simulation ---")
    
    drm_db = load_drm_database(args.db_json)
    ref_sequences = load_reference_fasta(args.ref_fasta)
    print("-" * 50)

    # A dictionary to count the final classifications
    classifications = {"DRM-Positive": 0, "DRM-Negative": 0, "Low-Quality": 0, "Non-HIV": 0, "Reference-Not-Found": 0}
    
    read_iterator = stream_bam_as_reads(args.input_bam, delay_ms=0)
    for i, read in enumerate(read_iterator):
        if i >= args.num_reads:
            break
        classification = classify_read_verbose(read, drm_db, ref_sequences)
        classifications[classification] += 1
    
    # --- FINAL REPORT ---
    print("\n" + "="*25 + " FINAL REPORT " + "="*25)
    print(f"Total Reads Processed: {sum(classifications.values())}")
    for category, count in classifications.items():
        print(f"  - {category}: {count}")
    print("="*64)

if __name__ == "__main__":
    main()