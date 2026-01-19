# ============================================================================
# FASTQ INTEGRATION MODULE FOR HIV DRM PIPELINE
# ============================================================================
#
# This module adds FASTQ input support to the existing BAM-based pipeline.
# It provides alignment functionality for Oxford Nanopore, PacBio, and Illumina data.
#
# INTEGRATION INSTRUCTIONS:
# 1. Add this entire file content to the TOP of run_entire_pipeline_improved.py
#    (right after the imports and before PART 1)
# 2. Modify the main() function as shown at the end of this file
# 3. The pipeline will then accept both --input_fastq and --input_bam
#
# ============================================================================

import subprocess
import gzip
import time
from pathlib import Path

# ============================================================================
# DEPENDENCY CHECKING
# ============================================================================

def check_tool_available(tool: str) -> bool:
    """Check if a command-line tool is available."""
    try:
        result = subprocess.run([tool, '--version'], 
                              capture_output=True, 
                              timeout=5,
                              text=True)
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
        return False


def check_dependencies(required_tools: List[str]) -> Dict[str, bool]:
    """Check if all required tools are installed."""
    print("INFO: Checking dependencies...")
    deps = {}
    
    for tool in required_tools:
        available = check_tool_available(tool)
        deps[tool] = available
        status = "✓" if available else "✗"
        print(f"  {status} {tool}: {'Available' if available else 'NOT FOUND'}")
    
    return deps


# ============================================================================
# FASTQ VALIDATION
# ============================================================================

def validate_fastq(fastq_path: str) -> bool:
    """Validate FASTQ file format and integrity."""
    try:
        # Determine opener based on file extension
        if fastq_path.endswith('.gz'):
            opener = gzip.open
            mode = 'rt'
        else:
            opener = open
            mode = 'r'
        
        with opener(fastq_path, mode) as f:
            # Read first 4 lines (one FASTQ record)
            lines = []
            for i in range(4):
                line = f.readline()
                if not line:
                    raise ValueError("File too short (less than 4 lines)")
                lines.append(line.strip())
            
            # Validate FASTQ format
            if not lines[0].startswith('@'):
                raise ValueError(f"Invalid header line: {lines[0][:50]}")
            
            if not lines[2].startswith('+'):
                raise ValueError(f"Invalid separator line: {lines[2][:50]}")
            
            if len(lines[1]) != len(lines[3]):
                raise ValueError(f"Sequence length != Quality length")
        
        print(f"INFO: FASTQ validation passed: {fastq_path}")
        return True
        
    except Exception as e:
        print(f"ERROR: FASTQ validation failed: {e}")
        return False


# ============================================================================
# PLATFORM DETECTION
# ============================================================================

def detect_sequencing_platform(fastq_path: str, sample_size: int = 100) -> str:
    """
    Detect sequencing platform from read lengths.
    
    Returns: 'long' for Nanopore/PacBio, 'short' for Illumina
    """
    try:
        if fastq_path.endswith('.gz'):
            opener = gzip.open
            mode = 'rt'
        else:
            opener = open
            mode = 'r'
        
        read_lengths = []
        
        with opener(fastq_path, mode) as f:
            read_count = 0
            for line in f:
                if line.startswith('@'):
                    # Next line is sequence
                    seq_line = f.readline().strip()
                    if seq_line:
                        read_lengths.append(len(seq_line))
                        read_count += 1
                    
                    if read_count >= sample_size:
                        break
        
        if read_lengths:
            median_length = np.median(read_lengths)
            platform = 'long' if median_length > 500 else 'short'
            
            print(f"INFO: Detected sequencing platform: {platform}-read")
            print(f"      Median read length: {median_length:.0f} bp")
            
            return platform
        
        # Default to long reads (Nanopore) if can't determine
        print("WARNING: Could not determine platform, defaulting to long-read (Nanopore)")
        return 'long'
        
    except Exception as e:
        print(f"WARNING: Platform detection failed: {e}. Defaulting to long-read")
        return 'long'


# ============================================================================
# CORE ALIGNMENT FUNCTION (PIPED FOR EFFICIENCY)
# ============================================================================

def align_fastq_to_bam_piped(
    fastq_path: str,
    ref_fasta: str,
    output_bam: str,
    aligner: str = 'minimap2',
    threads: int = 4,
    preset: str = 'map-ont'  # map-ont for Nanopore, map-pb for PacBio
) -> bool:
    """
    Align FASTQ to reference and directly create sorted BAM (piped).
    
    This is optimized for Oxford Nanopore data (map-ont preset).
    Pipeline: minimap2 | samtools view | samtools sort | samtools index
    
    Args:
        fastq_path: Path to input FASTQ file (can be .gz)
        ref_fasta: Path to reference FASTA
        output_bam: Path for output sorted BAM
        aligner: 'minimap2' (recommended for Nanopore) or 'bwa'
        threads: Number of CPU threads (4-8 recommended for Jetson AGX Orin)
        preset: minimap2 preset (map-ont=Nanopore, map-pb=PacBio, sr=Illumina)
    
    Returns:
        True if successful, False otherwise
    """
    print("="*80)
    print(" ALIGNMENT STEP")
    print("="*80)
    print(f"INFO: Aligner: {aligner}")
    print(f"INFO: Preset: {preset}")
    print(f"INFO: Input FASTQ: {fastq_path}")
    print(f"INFO: Reference: {ref_fasta}")
    print(f"INFO: Output BAM: {output_bam}")
    print(f"INFO: Threads: {threads}")
    print("="*80)
    
    start_time = time.time()
    
    try:
        # Build alignment command
        if aligner == 'minimap2':
            # Optimized for Nanopore (map-ont) or PacBio (map-pb)
            align_cmd = [
                'minimap2',
                '-ax', preset,  # Sequencing platform preset
                '-t', str(threads),
                '--secondary=no',  # No secondary alignments (cleaner results)
                '-L',  # Write CIGAR with >65535 ops at CG tag
                ref_fasta,
                fastq_path
            ]
        elif aligner == 'bwa':
            print("ERROR: BWA not yet implemented in this version.")
            print("       Use minimap2 for Nanopore data.")
            return False
        else:
            print(f"ERROR: Unknown aligner: {aligner}")
            return False
        
        # Build SAM to BAM conversion command
        view_cmd = [
            'samtools', 'view',
            '-bS',  # BAM output, SAM input
            '-@', str(threads),
            '-'  # Read from stdin
        ]
        
        # Build sort command
        sort_cmd = [
            'samtools', 'sort',
            '-@', str(threads),
            '-o', output_bam,
            '-'  # Read from stdin
        ]
        
        # Create pipeline: aligner | view | sort
        print("INFO: Running alignment pipeline (this may take 10-20 minutes)...")
        print("      minimap2 | samtools view | samtools sort")
        
        p1 = subprocess.Popen(align_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.Popen(view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p1.stdout.close()  # Allow p1 to receive SIGPIPE if p2 exits
        p3 = subprocess.Popen(sort_cmd, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2.stdout.close()  # Allow p2 to receive SIGPIPE if p3 exits
        
        # Wait for completion
        stdout3, stderr3 = p3.communicate()
        
        # Check return codes
        if p1.returncode != 0:
            stderr1 = p1.stderr.read().decode('utf-8') if p1.stderr else ''
            print("="*80)
            print(" ALIGNMENT ERROR")
            print("="*80)
            print("Minimap2 alignment failed. Common causes:")
            print("  1. Wrong reference genome (make sure it's HXB2 for HIV)")
            print("  2. Corrupted FASTQ file")
            print("  3. Out of memory (try reducing --threads)")
            print("="*80)
            print(f"Error details: {stderr1[:500]}")
            return False
        
        if p2.returncode != 0:
            stderr2 = p2.stderr.read().decode('utf-8') if p2.stderr else ''
            print(f"ERROR: SAM to BAM conversion failed: {stderr2[:500]}")
            return False
        
        if p3.returncode != 0:
            print(f"ERROR: BAM sorting failed: {stderr3.decode('utf-8')[:500]}")
            return False
        
        elapsed = time.time() - start_time
        print(f"INFO: Alignment completed in {elapsed/60:.1f} minutes ({elapsed:.1f} seconds)")
        
        # Index the BAM file
        print("INFO: Indexing BAM file...")
        index_cmd = ['samtools', 'index', output_bam]
        result = subprocess.run(index_cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"ERROR: BAM indexing failed: {result.stderr}")
            return False
        
        print("INFO: BAM file indexed successfully.")
        
        # Verify BAM file was created
        if not os.path.exists(output_bam):
            print("ERROR: Output BAM file was not created!")
            return False
        
        bam_size = os.path.getsize(output_bam) / (1024**2)  # Size in MB
        print(f"INFO: Created BAM file: {output_bam} ({bam_size:.1f} MB)")
        print("="*80)
        
        return True
        
    except Exception as e:
        print("="*80)
        print(" CRITICAL ERROR")
        print("="*80)
        print(f"Alignment pipeline failed: {e}")
        import traceback
        traceback.print_exc()
        print("="*80)
        return False


# ============================================================================
# ALIGNMENT QUALITY CHECKS
# ============================================================================

def get_alignment_statistics(bam_path: str, sample_size: int = 10000) -> Dict:
    """
    Calculate alignment statistics from BAM file.
    
    Args:
        bam_path: Path to BAM file
        sample_size: Number of reads to sample for stats
    
    Returns:
        Dictionary with statistics
    """
    stats = {
        'total_reads': 0,
        'mapped_reads': 0,
        'unmapped_reads': 0,
        'mapping_rate': 0.0,
        'avg_mapping_quality': 0.0
    }
    
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            read_count = 0
            quality_sum = 0
            
            for read in bam:
                stats['total_reads'] += 1
                read_count += 1
                
                if read.is_unmapped:
                    stats['unmapped_reads'] += 1
                else:
                    stats['mapped_reads'] += 1
                    quality_sum += read.mapping_quality
                
                # Sample for large files
                if read_count >= sample_size:
                    # Estimate total
                    try:
                        total_in_bam = bam.count(until_eof=True)
                        scaling_factor = total_in_bam / sample_size
                        stats['total_reads'] = int(stats['total_reads'] * scaling_factor)
                        stats['mapped_reads'] = int(stats['mapped_reads'] * scaling_factor)
                        stats['unmapped_reads'] = int(stats['unmapped_reads'] * scaling_factor)
                    except:
                        pass
                    break
            
            # Calculate rates
            if stats['total_reads'] > 0:
                stats['mapping_rate'] = (stats['mapped_reads'] / stats['total_reads']) * 100
            
            if stats['mapped_reads'] > 0:
                stats['avg_mapping_quality'] = quality_sum / min(stats['mapped_reads'], read_count)
        
        return stats
        
    except Exception as e:
        print(f"WARNING: Could not calculate alignment statistics: {e}")
        return stats


# ============================================================================
# MAIN FUNCTION MODIFICATIONS
# ============================================================================

def create_hybrid_argument_parser():
    """
    Create argument parser that accepts both FASTQ and BAM inputs.
    
    USAGE:
      # Option 1: From FASTQ (automatic alignment)
      python script.py --input_fastq reads.fastq.gz --ref_fasta ref.fa ...
      
      # Option 2: From BAM (skip alignment)
      python script.py --input_bam aligned.bam --ref_fasta ref.fa ...
    """
    parser = argparse.ArgumentParser(
        description=(
            "HIV Drug Resistance Mutation Analysis Pipeline v3.0 - HYBRID VERSION\n"
            "  ✓ Accepts FASTQ (raw reads) or BAM (aligned reads)\n"
            "  ✓ Optimized for Oxford Nanopore sequencing\n"
            "  ✓ Indel Detection (Insertions + Deletions)\n"
            "  ✓ Complete Gene Support (PR, RT, IN)\n"
            "  ✓ Drug-Specific Clinical Thresholds\n"
            "  ✓ Publication-Quality Statistics"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Input options (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input_fastq",
        help="Path to input FASTQ file (.fastq or .fastq.gz). Pipeline will perform alignment."
    )
    input_group.add_argument(
        "--input_bam",
        help="Path to input BAM file (pre-aligned). Pipeline will skip alignment step."
    )
    
    # Required for both modes
    parser.add_argument(
        "--ref_fasta",
        required=True,
        help="Path to reference FASTA (HXB2 for HIV)"
    )
    parser.add_argument(
        "--profile_csv",
        required=True,
        help="Path to drug resistance profile CSV"
    )
    
    # Optional parameters
    parser.add_argument(
        "-n", "--num_reads",
        type=int,
        default=0,
        help="Number of reads to process (0=all, default: 0)"
    )
    parser.add_argument(
        "-o", "--output_dir",
        default="data",
        help="Output directory (default: data)"
    )
    
    # Alignment-specific parameters (only used if --input_fastq provided)
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of CPU threads for alignment (default: 4, recommend 4-8 for Jetson AGX Orin)"
    )
    parser.add_argument(
        "--aligner",
        choices=['minimap2', 'auto'],
        default='minimap2',
        help="Alignment tool (default: minimap2 for Nanopore)"
    )
    parser.add_argument(
        "--preset",
        choices=['map-ont', 'map-pb', 'sr'],
        default='map-ont',
        help="Minimap2 preset: map-ont (Nanopore), map-pb (PacBio), sr (Illumina)"
    )
    parser.add_argument(
        "--keep-bam",
        action='store_true',
        help="Keep intermediate BAM file after analysis (useful for debugging)"
    )
    
    return parser


# ============================================================================
# EXAMPLE: HOW TO MODIFY YOUR EXISTING main() FUNCTION
# ============================================================================

"""
In your existing run_entire_pipeline_improved.py, replace the main() function with:

def main():
    # Use the new argument parser
    parser = create_hybrid_argument_parser()
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    start_time = datetime.datetime.now()
    timestamp_file = start_time.strftime("%Y%m%d_%H%M%S")
    timestamp_report = start_time.strftime("%Y-%m-%d %H:%M:%S")
    
    print("="*80)
    print(" HIV DRUG RESISTANCE MUTATION ANALYSIS PIPELINE v3.0")
    print(" HYBRID VERSION: FASTQ + BAM Support")
    print("="*80)
    print(f"Started: {timestamp_report}")
    
    # Determine if we need to align
    if args.input_fastq:
        print("INFO: Input mode: FASTQ (will perform alignment)")
        print(f"INFO: Input FASTQ: {args.input_fastq}")
        
        # Validate FASTQ
        if not os.path.exists(args.input_fastq):
            print(f"ERROR: FASTQ file not found: {args.input_fastq}")
            sys.exit(1)
        
        if not validate_fastq(args.input_fastq):
            print("ERROR: FASTQ validation failed!")
            sys.exit(1)
        
        # Check dependencies
        required_tools = ['minimap2', 'samtools']
        deps = check_dependencies(required_tools)
        
        missing = [tool for tool, available in deps.items() if not available]
        if missing:
            print("="*80)
            print(" MISSING DEPENDENCIES")
            print("="*80)
            print(f"The following tools are required but not installed: {', '.join(missing)}")
            print("\nInstall with:")
            print("  Ubuntu/Debian: sudo apt-get install minimap2 samtools")
            print("  macOS: brew install minimap2 samtools")
            print("  Conda: conda install -c bioconda minimap2 samtools")
            print("="*80)
            sys.exit(1)
        
        # Detect platform (optional)
        platform = detect_sequencing_platform(args.input_fastq)
        
        # NEW STEP: Alignment
        print("\nSTEP 1/6: Aligning FASTQ to reference genome...")
        
        base_name = os.path.basename(args.input_fastq).replace('.fastq.gz', '').replace('.fastq', '')
        temp_bam = os.path.join(args.output_dir, f"{base_name}_aligned.sorted.bam")
        
        success = align_fastq_to_bam_piped(
            args.input_fastq,
            args.ref_fasta,
            temp_bam,
            args.aligner,
            args.threads,
            args.preset
        )
        
        if not success:
            print("FATAL: Alignment failed. Cannot continue.")
            sys.exit(1)
        
        # Check alignment quality
        print("\nSTEP 2/6: Checking alignment quality...")
        stats = get_alignment_statistics(temp_bam)
        print(f"  Total reads: {stats['total_reads']:,}")
        print(f"  Mapped reads: {stats['mapped_reads']:,} ({stats['mapping_rate']:.1f}%)")
        print(f"  Avg mapping quality: {stats['avg_mapping_quality']:.1f}")
        
        if stats['mapping_rate'] < 50:
            print("="*80)
            print(" WARNING: LOW MAPPING RATE")
            print("="*80)
            print("Less than 50% of reads mapped to the reference.")
            print("Possible causes:")
            print("  1. Wrong reference genome (ensure it's HXB2 for HIV)")
            print("  2. Low quality sequencing data")
            print("  3. Sample contamination (non-HIV reads)")
            print("="*80)
            response = input("Continue anyway? (y/n): ")
            if response.lower() != 'y':
                sys.exit(1)
        
        # Use the aligned BAM for downstream analysis
        input_bam_path = temp_bam
        step_offset = 2  # We've done 2 extra steps
        
    else:
        print("INFO: Input mode: BAM (skipping alignment)")
        print(f"INFO: Input BAM: {args.input_bam}")
        
        # Validate BAM exists
        if not os.path.exists(args.input_bam):
            print(f"ERROR: BAM file not found: {args.input_bam}")
            sys.exit(1)
        
        input_bam_path = args.input_bam
        step_offset = 0
    
    print(f"INFO: Reference: {args.ref_fasta}")
    print(f"INFO: Profile DB: {args.profile_csv}")
    print("="*80)
    print()
    
    # From here, everything is the same as the original pipeline
    # Just update step numbers (add step_offset)
    
    print(f"STEP {3+step_offset}/6: Building drug resistance database...")
    drm_db = build_definitive_database(args.profile_csv)
    print()
    
    print(f"STEP {4+step_offset}/6: Loading reference genome...")
    ref_sequences = load_reference_fasta(args.ref_fasta)
    print(f"INFO: Loaded {len(ref_sequences)} reference sequence(s).")
    print()
    
    print(f"STEP {5+step_offset}/6: Processing reads for mutations...")
    
    # ... rest of your existing processing code ...
    # Use input_bam_path instead of args.input_bam
    
    # ... existing code for processing, analysis, plots, PDF ...
    
    print(f"STEP {6+step_offset}/6: Generating report...")
    
    # ... existing reporting code ...
    
    print("\n" + "="*80)
    print(" PIPELINE COMPLETED SUCCESSFULLY")
    print("="*80)
    
    # Cleanup
    if args.input_fastq and not args.keep_bam:
        print("INFO: Cleaning up intermediate BAM (use --keep-bam to preserve)")
        # Optionally delete temp_bam here


if __name__ == "__main__":
    main()
"""