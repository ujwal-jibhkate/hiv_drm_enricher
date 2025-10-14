import argparse
import sys
sys.path.append('src')

from hiv_drm_enricher.streaming import stream_bam_as_reads
# 1. UPDATE THE IMPORT
from hiv_drm_enricher.classifierV2 import classify_read_verbose, load_drm_database

def load_reference_fasta(fasta_path: str) -> str:
    """Loads the first sequence from a FASTA file."""
    print(f"INFO: Loading reference genome from '{fasta_path}'...")
    try:
        with open(fasta_path, 'r') as f:
            lines = f.readlines()
            # Filter out header and join sequence lines
            sequence = "".join([line.strip() for line in lines if not line.startswith('>')])
        print("INFO: Reference genome loaded successfully.")
        return sequence
    except FileNotFoundError:
        print(f"FATAL: Reference FASTA file not found at '{fasta_path}'.", file=sys.stderr)
        raise SystemExit(1)

def main():
    """Main function to run the full, integrated simulation pipeline."""
    parser = argparse.ArgumentParser(description="Simulate a real-time sequencing stream from a BAM file.")
    parser.add_argument("-i", "--input_bam", required=True, help="Path to the input BAM file.")
    # 2. ADD ARGUMENTS FOR DATABASE AND FASTA
    parser.add_argument("--db_json", default="data/drm_database.json", help="Path to the DRM database JSON file.")
    parser.add_argument("--ref_fasta", default="data/HIV_reference_dummy.fasta", help="Path to the reference FASTA file.")
    parser.add_argument("-d", "--delay", type=int, default=10, help="Delay in ms between reads. Default: 10ms.")
    parser.add_argument("-n", "--num_reads", type=int, default=50, help="Number of reads to process. Default: 50.")
    args = parser.parse_args()

    print("--- HIV DRM Enrichment Pipeline: Full Simulation ---")
    
    # 3. LOAD THE DATABASE AND REFERENCE
    drm_db = load_drm_database(args.db_json)
    ref_sequence = load_reference_fasta(args.ref_fasta)

    print("-" * 50)

    read_counter = 0
    try:
        for read in stream_bam_as_reads(args.input_bam, args.delay):
            if read_counter >= args.num_reads:
                print(f"\nINFO: Reached max reads ({args.num_reads}). Stopping.")
                break

            # 4. CALL THE CLASSIFIER WITH THE NEW ARGUMENTS
            classification = classify_read_verbose(read, drm_db, ref_sequence) 
            
            read_counter += 1
            
    except Exception as e:
        print(f"\nA critical error occurred: {e}")
        sys.exit(1)

    print("\n--- Simulation Complete ---")
    print(f"Total reads processed: {read_counter}")

if __name__ == "__main__":
    main()