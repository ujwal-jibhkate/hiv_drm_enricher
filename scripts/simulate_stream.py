import argparse
import sys
# We need to add the src directory to the path to allow our script to find the package
sys.path.append('src')

from hiv_drm_enricher.streaming import stream_bam_as_reads
from hiv_drm_enricher.classifier import classify_read_verbose

def main():
    """
    Main function to run the BAM streaming simulation.
    Parses command-line arguments and prints information for each streamed read.
    """
    parser = argparse.ArgumentParser(
        description="Simulate a real-time sequencing stream from a BAM file."
    )
    parser.add_argument(
        "-i", "--input_bam",
        required=True,
        help="Path to the input BAM file."
    )
    parser.add_argument(
        "-d", "--delay",
        type=int,
        default=10,
        help="Delay in milliseconds between reads to simulate streaming. Default: 10ms."
    )
    parser.add_argument(
        "-n", "--num_reads",
        type=int,
        default=50,
        help="Number of reads to process before stopping. Default: 50."
    )
    args = parser.parse_args()

    print("--- HIV DRM Enrichment Pipeline: Streaming Simulator ---")
    print(f"Configuration:")
    print(f"  - Input BAM: {args.input_bam}")
    print(f"  - Delay: {args.delay} ms")
    print(f"  - Max Reads: {args.num_reads}")
    print("-" * 50)

    read_counter = 0
    try:
        # Use the generator to get reads one by one
        for read in stream_bam_as_reads(args.input_bam, args.delay):
            if read_counter >= args.num_reads:
                print(f"\nINFO: Reached the maximum number of reads ({args.num_reads}). Stopping.")
                break
            
            # For now, just print basic info to prove it's working
            classification = classify_read_verbose(read) 
            
            read_counter += 1
            
    except Exception as e:
        print(f"\nA critical error occurred: {e}")
        sys.exit(1)

    print("\n--- Simulation Complete ---")
    print(f"Total reads processed: {read_counter}")

if __name__ == "__main__":
    main()