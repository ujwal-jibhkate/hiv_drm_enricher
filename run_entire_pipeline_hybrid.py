import pandas as pd
import json
import sys
import re
import argparse
import pysam
import os
import numpy as np
from scipy.stats import norm, binomtest
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_effectsize
from statsmodels.stats.power import zt_ind_solve_power
from fpdf import FPDF, XPos, YPos
import matplotlib.pyplot as plt
import matplotlib
import datetime
import itertools
from typing import Dict, Set, Tuple, List, Optional
from collections import defaultdict, Counter

# --- Configuration & Constants ---
matplotlib.use('Agg')

# HXB2 Reference Coordinates (K03455.1) - 0-based, half-open intervals
GENE_COORDINATES = {
    "PR": {"start": 2252, "end": 2549, "frame": 0},
    "RT": {"start": 2549, "end": 4229, "frame": 0},
    "IN": {"start": 4229, "end": 5096, "frame": 0}
}

# Standard Genetic Code
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}

# Analysis parameters
NOISE_THRESHOLD = 30
MAX_COMBO_SIZE = 4

# ============================================================================
# IMPROVEMENT #3: DRUG-SPECIFIC CLINICAL THRESHOLDS
# ============================================================================
# Different drug classes have different barriers to resistance
# NNRTIs: Low barrier (single mutation can cause resistance)
# NRTIs: Variable barrier (some high, some low)
# PIs: High barrier (multiple mutations usually needed)
# INSTIs: Variable barrier

DRUG_SPECIFIC_THRESHOLDS = {
    # NNRTIs (Non-Nucleoside Reverse Transcriptase Inhibitors) - Low barrier
    'EFV': 0.01,   # Efavirenz
    'NVP': 0.01,   # Nevirapine
    'ETR': 0.01,   # Etravirine
    'RPV': 0.01,   # Rilpivirine
    'DOR': 0.01,   # Doravirine
    
    # NRTIs (Nucleoside Reverse Transcriptase Inhibitors) - Variable
    '3TC': 0.20,   # Lamivudine - high barrier (M184V needed)
    'FTC': 0.20,   # Emtricitabine - high barrier
    'ABC': 0.10,   # Abacavir - medium barrier
    'AZT': 0.10,   # Zidovudine - medium barrier
    'D4T': 0.10,   # Stavudine - medium barrier
    'DDI': 0.10,   # Didanosine - medium barrier
    'TDF': 0.05,   # Tenofovir - low-medium barrier
    'TAF': 0.05,   # Tenofovir alafenamide
    
    # PIs (Protease Inhibitors) - High barrier
    'ATV/r': 0.10, # Atazanavir/ritonavir
    'DRV/r': 0.15, # Darunavir/ritonavir - very high barrier
    'LPV/r': 0.10, # Lopinavir/ritonavir
    'FPV': 0.10,   # Fosamprenavir
    'IDV': 0.10,   # Indinavir
    'NFV': 0.10,   # Nelfinavir
    'SQV': 0.10,   # Saquinavir
    'TPV': 0.10,   # Tipranavir
    'ATV': 0.10,   # Atazanavir (unboosted)
    'DRV': 0.15,   # Darunavir (unboosted)
    'LPV': 0.10,   # Lopinavir (unboosted)
    
    # INSTIs (Integrase Strand Transfer Inhibitors) - Variable
    'RAL': 0.05,   # Raltegravir - medium barrier
    'EVG': 0.05,   # Elvitegravir - medium barrier
    'DTG': 0.10,   # Dolutegravir - high barrier
    'BIC': 0.10,   # Bictegravir - high barrier
    'CAB': 0.10,   # Cabotegravir - high barrier
}

# Default threshold for drugs not in the specific list
DEFAULT_CLINICAL_THRESHOLD = 0.01

def get_drug_threshold(drug_name: str) -> float:
    """
    Returns the clinical threshold for a given drug.
    Uses drug-specific threshold if available, otherwise default.
    """
    return DRUG_SPECIFIC_THRESHOLDS.get(drug_name, DEFAULT_CLINICAL_THRESHOLD)


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

# ============================================================================
# PART 1: GENE ASSIGNMENT & DATABASE BUILDER
# ============================================================================

def determine_gene_from_position(aa_position: int) -> Optional[str]:
    """
    Determines gene from amino acid position.
    
    IMPROVEMENT #2: Now properly handles IN gene (integrase).
    
    Gene boundaries (amino acid positions):
    - PR (Protease): 1-99
    - RT (Reverse Transcriptase): 100-560  
    - IN (Integrase): 561-849
    """
    if 1 <= aa_position <= 99:
        return "PR"
    elif 100 <= aa_position <= 560:
        return "RT"
    elif 561 <= aa_position <= 849:
        return "IN"
    else:
        return None


def parse_mutation_list_from_csv(mutation_cell: str) -> List[str]:
    """Parses mutation column from CSV."""
    cleaned = mutation_cell.strip().strip('"').strip("'")
    
    if cleaned.startswith('[') and cleaned.endswith(']'):
        inner = cleaned[1:-1]
        mutations = [m.strip().strip('"').strip("'") for m in inner.split(',')]
        return [m for m in mutations if m]
    
    return [m.strip() for m in cleaned.split(',') if m.strip()]


def build_definitive_database(csv_path: str) -> Dict[str, Dict[str, Dict]]:
    """Builds comprehensive drug resistance database from CSV."""
    print(f"INFO: Building drug resistance database from '{csv_path}'...")
    
    try:
        df = pd.read_csv(csv_path, dtype=str)
    except FileNotFoundError:
        print(f"FATAL: CSV file not found at '{csv_path}'.", file=sys.stderr)
        sys.exit(1)
    
    df_long = pd.melt(df, id_vars=['mutations'], var_name='Drug', value_name='Score')
    df_long = df_long[df_long['Score'] != 'S'].dropna(subset=['Score'])
    
    drm_db = {"PR": {}, "RT": {}, "IN": {}, "CROSS_GENE": {}}
    mutation_pattern = re.compile(r'([A-Z]?)(\d+)([A-Z*])', re.IGNORECASE)
    
    for idx, row in df_long.iterrows():
        raw_mutations = parse_mutation_list_from_csv(row['mutations'])
        if not raw_mutations:
            continue
        
        parsed_mutations = []
        genes_involved = set()
        
        for mut in raw_mutations:
            match = mutation_pattern.match(mut.strip())
            if not match:
                continue
            
            ref_aa, position_str, alt_aa = match.groups()
            position = int(position_str)
            gene = determine_gene_from_position(position)
            
            if not gene:
                continue
            
            genes_involved.add(gene)
            mutation_canonical = f"?{position}{alt_aa.upper()}"
            parsed_mutations.append((gene, mutation_canonical))
        
        if not parsed_mutations:
            continue
        
        is_cross_gene = len(genes_involved) > 1
        drug_info = {'name': row['Drug'], 'score': row['Score']}
        
        if is_cross_gene:
            mutation_key = ",".join(sorted([f"{gene}:{mut}" for gene, mut in parsed_mutations]))
            
            if mutation_key not in drm_db["CROSS_GENE"]:
                drm_db["CROSS_GENE"][mutation_key] = {'drugs': [], 'genes': list(genes_involved)}
            
            if drug_info not in drm_db["CROSS_GENE"][mutation_key]['drugs']:
                drm_db["CROSS_GENE"][mutation_key]['drugs'].append(drug_info)
        else:
            gene = list(genes_involved)[0]
            mutation_list = sorted([mut for g, mut in parsed_mutations])
            mutation_key = ",".join(mutation_list)
            
            if mutation_key not in drm_db[gene]:
                drm_db[gene][mutation_key] = {'drugs': []}
            
            if drug_info not in drm_db[gene][mutation_key]['drugs']:
                drm_db[gene][mutation_key]['drugs'].append(drug_info)
    
    total_entries = sum(len(drm_db[gene]) for gene in drm_db)
    print(f"INFO: Database built with {total_entries} unique mutation patterns:")
    print(f"      - PR: {len(drm_db['PR'])} patterns")
    print(f"      - RT: {len(drm_db['RT'])} patterns")
    print(f"      - IN: {len(drm_db['IN'])} patterns")
    print(f"      - Cross-gene: {len(drm_db['CROSS_GENE'])} patterns")
    
    return drm_db


def get_mutation_subsets(mutation_set: Set[str], max_size: int = MAX_COMBO_SIZE) -> Set[str]:
    """
    Generates all non-empty subsets up to max_size.
    
    Safety features:
    - Caps combination size at max_size (default: 4)
    - If mutation_set is too large, only generates up to max_size subsets
    - Prevents combinatorial explosion
    
    Complexity: O(n choose k) where k = min(max_size, len(mutation_set))
    Worst case: C(30,4) = 27,405 combinations (acceptable)
    """
    if not mutation_set:
        return set()
    
    subsets_as_keys = set()
    mut_list = sorted(list(mutation_set))
    
    # Safety check: warn if mutation count is suspiciously high
    if len(mut_list) > 20:
        print(f"WARNING: Large mutation set ({len(mut_list)} mutations). "
              f"Generating combinations up to size {max_size}.")
    
    # Cap at max_size
    max_k = min(len(mut_list), max_size)
    
    # Generate all combinations from size 1 to max_k
    for k in range(1, max_k + 1):
        for subset_tuple in itertools.combinations(mut_list, k):
            canonical_key = ",".join(subset_tuple)
            subsets_as_keys.add(canonical_key)
    
    return subsets_as_keys


# ============================================================================
# PART 2: REFERENCE GENOME & READ PROCESSING (WITH INDEL SUPPORT)
# ============================================================================

def load_reference_fasta(path: str) -> Dict[str, str]:
    """Loads reference sequences from FASTA file."""
    sequences = {}
    current_seq_name = ""
    
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_seq_name = line[1:].split()[0]
                sequences[current_seq_name] = ""
            elif current_seq_name:
                sequences[current_seq_name] += line
    
    return sequences


def get_amino_acid_change(ref_pos: int, alt_base: str, ref_seq: str) -> Optional[Tuple[str, str, int, str]]:
    """Determines amino acid change from SNP."""
    for gene, coords in GENE_COORDINATES.items():
        if coords['start'] <= ref_pos < coords['end']:
            gene_start = coords['start']
            pos_in_gene = ref_pos - gene_start
            codon_start_in_gene = pos_in_gene - (pos_in_gene % 3)
            codon_start_in_ref = gene_start + codon_start_in_gene
            
            ref_codon = ref_seq[codon_start_in_ref:codon_start_in_ref + 3].upper()
            if len(ref_codon) < 3:
                return None
            
            mut_codon = list(ref_codon)
            position_in_codon = pos_in_gene % 3
            mut_codon[position_in_codon] = alt_base.upper()
            mut_codon_str = "".join(mut_codon)
            
            ref_aa = CODON_TABLE.get(ref_codon, '?')
            alt_aa = CODON_TABLE.get(mut_codon_str, '?')
            
            if ref_aa != alt_aa:
                aa_position = (codon_start_in_gene // 3) + 1
                return (ref_aa, alt_aa, aa_position, gene)
    
    return None


# ============================================================================
# IMPROVEMENT #1: INDEL DETECTION AND PROCESSING
# ============================================================================

def get_amino_acid_deletion(ref_pos: int, deletion_length: int, ref_seq: str) -> Optional[Tuple[str, int, str, str]]:
    """
    Determines amino acid impact of a deletion.
    
    Returns: (ref_aa, aa_position, deletion_type, gene)
    deletion_type can be: 'single_aa', 'frameshift', 'multi_aa', 'partial_codon'
    """
    for gene, coords in GENE_COORDINATES.items():
        if coords['start'] <= ref_pos < coords['end']:
            gene_start = coords['start']
            pos_in_gene = ref_pos - gene_start
            
            # Determine if deletion affects whole codons or causes frameshift
            if deletion_length % 3 == 0:
                # In-frame deletion (removes whole amino acids)
                num_aa_deleted = deletion_length // 3
                codon_start_in_gene = pos_in_gene - (pos_in_gene % 3)
                aa_position = (codon_start_in_gene // 3) + 1
                
                # Get reference amino acid(s) that will be deleted
                codon_start_in_ref = gene_start + codon_start_in_gene
                deleted_codons = ref_seq[codon_start_in_ref:codon_start_in_ref + deletion_length].upper()
                ref_aa = ''.join([CODON_TABLE.get(deleted_codons[i:i+3], '?') for i in range(0, len(deleted_codons), 3)])
                
                deletion_type = 'single_aa' if num_aa_deleted == 1 else 'multi_aa'
                return (ref_aa, aa_position, deletion_type, gene)
            else:
                # Frameshift deletion (not multiple of 3)
                codon_start_in_gene = pos_in_gene - (pos_in_gene % 3)
                aa_position = (codon_start_in_gene // 3) + 1
                return ('?', aa_position, 'frameshift', gene)
    
    return None


def get_amino_acid_insertion(ref_pos: int, insertion_seq: str, ref_seq: str) -> Optional[Tuple[str, int, str, str]]:
    """
    Determines amino acid impact of an insertion.
    
    Returns: (inserted_aa, aa_position, insertion_type, gene)
    insertion_type can be: 'single_aa', 'frameshift', 'multi_aa'
    """
    for gene, coords in GENE_COORDINATES.items():
        if coords['start'] <= ref_pos < coords['end']:
            gene_start = coords['start']
            pos_in_gene = ref_pos - gene_start
            
            insertion_length = len(insertion_seq)
            
            # Determine if insertion is in-frame or frameshift
            if insertion_length % 3 == 0:
                # In-frame insertion (adds whole amino acids)
                num_aa_inserted = insertion_length // 3
                codon_start_in_gene = pos_in_gene - (pos_in_gene % 3)
                aa_position = (codon_start_in_gene // 3) + 1
                
                # Translate inserted sequence
                inserted_aa = ''.join([CODON_TABLE.get(insertion_seq[i:i+3].upper(), '?') 
                                      for i in range(0, len(insertion_seq), 3)])
                
                insertion_type = 'single_aa' if num_aa_inserted == 1 else 'multi_aa'
                return (inserted_aa, aa_position, insertion_type, gene)
            else:
                # Frameshift insertion
                codon_start_in_gene = pos_in_gene - (pos_in_gene % 3)
                aa_position = (codon_start_in_gene // 3) + 1
                return ('?', aa_position, 'frameshift', gene)
    
    return None


def process_read_for_mutations(read, ref_sequences: Dict[str, str], drm_db: Dict[str, Dict]) -> Dict:
    """
    Processes a single read to identify mutations (SNPs + INDELs) and classify as DRM-positive.
    
    IMPROVEMENT #1: Now detects insertions and deletions in addition to SNPs.
    
    CRITICAL: Returns consistent structure for deterministic results.
    """
    # Initialize result with all required fields
    result = {
        "read_id": read.query_name,
        "status": "Unknown",
        "mutations_by_gene": {},
        "indels_by_gene": {},  # NEW: Track indels separately
        "quality_score": None,
        "alignment_start": None,
        "alignment_end": None,
        "covered_positions": []
    }
    
    # Quick rejections
    if read.is_unmapped:
        result["status"] = "Unmapped"
        return result
    
    if read.mapping_quality < 20:
        result["status"] = "Low-Quality"
        result["quality_score"] = read.mapping_quality
        return result
    
    # Get reference sequence
    ref_seq = ref_sequences.get(read.reference_name)
    if not ref_seq and len(ref_sequences) == 1:
        _, ref_seq = next(iter(ref_sequences.items()))
    elif not ref_seq:
        result["status"] = "Reference-Not-Found"
        return result
    
    # Store quality and alignment info
    result["quality_score"] = read.mapping_quality
    result["alignment_start"] = read.reference_start
    result["alignment_end"] = read.reference_end
    
    # Track mutations and coverage
    found_mutations = {"PR": set(), "RT": set(), "IN": set()}
    found_indels = {"PR": [], "RT": [], "IN": []}  # NEW: Store indel details
    mutation_count = 0
    covered_positions = []
    
    # Process aligned pairs for SNPs and coverage
    aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=True)
    
    for read_pos, ref_pos, ref_base in aligned_pairs:
        # Track coverage
        if ref_pos is not None:
            covered_positions.append(ref_pos)
        
        # Only process matches/mismatches (SNPs)
        if read_pos is not None and ref_pos is not None and ref_base is not None:
            read_base = read.query_sequence[read_pos]
            
            # Check for mismatch
            if read_base.upper() != ref_base.upper():
                aa_change = get_amino_acid_change(ref_pos, read_base, ref_seq)
                
                if aa_change:
                    ref_aa, alt_aa, aa_pos, gene = aa_change
                    mutation_str = f"?{aa_pos}{alt_aa}"
                    found_mutations[gene].add(mutation_str)
                    mutation_count += 1
                    
                    # EARLY EXIT for noisy reads
                    if mutation_count > NOISE_THRESHOLD:
                        result["status"] = "Noisy-Read-Filtered"
                        result["mutations_by_gene"] = {}
                        result["indels_by_gene"] = {}
                        result["covered_positions"] = []
                        return result
    
    # NEW: Process CIGAR string for indels
    cigar_tuples = read.cigartuples
    if cigar_tuples:
        current_ref_pos = read.reference_start
        current_read_pos = 0
        
        for cigar_op, cigar_len in cigar_tuples:
            # CIGAR operations:
            # 0 = M (match/mismatch)
            # 1 = I (insertion to reference)
            # 2 = D (deletion from reference)
            # 4 = S (soft clipping)
            # 5 = H (hard clipping)
            
            if cigar_op == 0:  # Match/mismatch (M)
                current_ref_pos += cigar_len
                current_read_pos += cigar_len
                
            elif cigar_op == 1:  # Insertion (I)
                # Bases present in read but not in reference
                inserted_seq = read.query_sequence[current_read_pos:current_read_pos + cigar_len]
                
                insertion_info = get_amino_acid_insertion(current_ref_pos, inserted_seq, ref_seq)
                if insertion_info:
                    inserted_aa, aa_pos, ins_type, gene = insertion_info
                    
                    # Format: position_ins_AA (e.g., "69_ins_SS" for T69SS insertion)
                    indel_notation = f"{aa_pos}_ins_{inserted_aa}"
                    found_indels[gene].append({
                        'type': 'insertion',
                        'position': aa_pos,
                        'notation': indel_notation,
                        'inserted_aa': inserted_aa,
                        'insertion_type': ins_type
                    })
                    
                    # Also add to mutations set for DRM matching
                    found_mutations[gene].add(indel_notation)
                    mutation_count += 1
                    
                    if mutation_count > NOISE_THRESHOLD:
                        result["status"] = "Noisy-Read-Filtered"
                        result["mutations_by_gene"] = {}
                        result["indels_by_gene"] = {}
                        result["covered_positions"] = []
                        return result
                
                current_read_pos += cigar_len
                
            elif cigar_op == 2:  # Deletion (D)
                # Bases present in reference but not in read
                deletion_info = get_amino_acid_deletion(current_ref_pos, cigar_len, ref_seq)
                if deletion_info:
                    ref_aa, aa_pos, del_type, gene = deletion_info
                    
                    # Format: position_del (e.g., "69_del" for deletion at position 69)
                    indel_notation = f"{aa_pos}_del"
                    found_indels[gene].append({
                        'type': 'deletion',
                        'position': aa_pos,
                        'notation': indel_notation,
                        'deleted_aa': ref_aa,
                        'deletion_type': del_type
                    })
                    
                    # Also add to mutations set for DRM matching
                    found_mutations[gene].add(indel_notation)
                    mutation_count += 1
                    
                    if mutation_count > NOISE_THRESHOLD:
                        result["status"] = "Noisy-Read-Filtered"
                        result["mutations_by_gene"] = {}
                        result["indels_by_gene"] = {}
                        result["covered_positions"] = []
                        return result
                
                current_ref_pos += cigar_len
                
            elif cigar_op == 4:  # Soft clip (S)
                current_read_pos += cigar_len
                
            elif cigar_op == 5:  # Hard clip (H)
                pass  # No change to positions
    
    # Store covered positions
    result["covered_positions"] = covered_positions
    
    # Check if any mutations are known DRMs (check ALL combinations up to size 4)
    has_known_drm = False
    if mutation_count > 0:
        # Safety: If too many mutations, only check singles (performance protection)
        # This should never happen due to NOISE_THRESHOLD, but added as defense-in-depth
        effective_max_size = 1 if mutation_count > 20 else MAX_COMBO_SIZE
        
        for gene, mut_set in found_mutations.items():
            if not mut_set or gene not in drm_db:
                continue
            
            # Generate all subsets up to effective_max_size
            subsets = get_mutation_subsets(mut_set, max_size=effective_max_size)
            
            for subset_key in subsets:
                if subset_key in drm_db[gene]:
                    has_known_drm = True
                    break
            
            if has_known_drm:
                break
        
        # Also check cross-gene combinations if multiple genes have mutations
        if not has_known_drm and len(found_mutations) > 1 and drm_db.get("CROSS_GENE"):
            all_cross_mutations = []
            for gene, muts in found_mutations.items():
                if muts:  # Only include genes with mutations
                    for mut in muts:
                        all_cross_mutations.append(f"{gene}:{mut}")
            
            if all_cross_mutations and len(all_cross_mutations) <= 20:  # Safety check
                # Generate cross-gene combinations (also capped at size 4)
                cross_subsets = get_mutation_subsets(set(all_cross_mutations), max_size=effective_max_size)
                for cross_key in cross_subsets:
                    if cross_key in drm_db["CROSS_GENE"]:
                        has_known_drm = True
                        break
    
    # Determine final status
    if mutation_count == 0:
        result["status"] = "Wild-Type"
    elif has_known_drm:
        result["status"] = "DRM-Positive"
    else:
        result["status"] = "Polymorphic"
    
    # Convert sets to sorted lists for JSON serialization
    result["mutations_by_gene"] = {
        gene: sorted(list(mut_set))
        for gene, mut_set in found_mutations.items()
        if mut_set
    }
    
    # Store indel information
    result["indels_by_gene"] = {
        gene: indel_list
        for gene, indel_list in found_indels.items()
        if indel_list
    }
    
    return result


# ============================================================================
# PART 3: ENHANCED STATISTICAL FUNCTIONS
# ============================================================================

def calculate_wilson_ci(count: int, total: int, confidence: float = 0.95) -> Tuple[float, float]:
    """
    Wilson score confidence interval for proportions.
    More accurate than normal approximation for small counts.
    """
    if total == 0:
        return (0.0, 0.0)
    
    p_hat = count / total
    z = norm.ppf((1 + confidence) / 2)
    
    denominator = 1 + z**2 / total
    center = (p_hat + z**2 / (2 * total)) / denominator
    margin = z * np.sqrt(p_hat * (1 - p_hat) / total + z**2 / (4 * total**2)) / denominator
    
    lower = max(0, center - margin)
    upper = min(1, center + margin)
    
    return (lower * 100, upper * 100)


def binomial_exact_test(count: int, total: int, threshold: float) -> float:
    """Exact binomial test for rare variants."""
    if total == 0:
        return 1.0
    result = binomtest(count, total, threshold, alternative='greater')
    return result.pvalue


def calculate_statistical_power(n_reads: int, threshold: float, 
                                alpha: float = 0.05, detect_at: float = None) -> float:
    """Calculate statistical power to detect resistance."""
    if detect_at is None:
        detect_at = threshold * 1.5
    
    if n_reads < 10:
        return 0.0
    
    try:
        effect_size = proportion_effectsize(detect_at, threshold)
        power = zt_ind_solve_power(effect_size=effect_size, nobs1=n_reads, 
                                   alpha=alpha, alternative='larger')
        return max(0.0, min(1.0, power))
    except:
        return 0.0


def calculate_coverage_by_position(reads_data: List[Dict]) -> Dict[int, int]:
    """Calculate coverage depth at each genomic position."""
    coverage = defaultdict(int)
    
    for read in reads_data:
        for pos in read.get('covered_positions', []):
            coverage[pos] += 1
    
    return dict(coverage)


def calculate_mutation_counts_by_gene(reads_data: List[Dict]) -> Dict[str, int]:
    """Count total mutations detected in each gene region."""
    gene_mutation_counts = {"PR": 0, "RT": 0, "IN": 0}
    
    for read in reads_data:
        for gene, mutations in read.get('mutations_by_gene', {}).items():
            if gene in gene_mutation_counts:
                gene_mutation_counts[gene] += len(mutations)
    
    return gene_mutation_counts


def calculate_indel_counts_by_gene(reads_data: List[Dict]) -> Dict[str, Dict[str, int]]:
    """Count insertions and deletions in each gene region."""
    indel_counts = {
        "PR": {"insertions": 0, "deletions": 0},
        "RT": {"insertions": 0, "deletions": 0},
        "IN": {"insertions": 0, "deletions": 0}
    }
    
    for read in reads_data:
        for gene, indels in read.get('indels_by_gene', {}).items():
            if gene in indel_counts:
                for indel in indels:
                    if indel['type'] == 'insertion':
                        indel_counts[gene]['insertions'] += 1
                    elif indel['type'] == 'deletion':
                        indel_counts[gene]['deletions'] += 1
    
    return indel_counts


# ============================================================================
# PART 4: COMPREHENSIVE STATISTICAL ANALYSIS (WITH DRUG-SPECIFIC THRESHOLDS)
# ============================================================================

def perform_statistical_analysis(
    json_path: str,
    drm_db: Dict[str, Dict]
) -> Optional[Dict]:
    """
    Performs comprehensive statistical analysis.
    
    IMPROVEMENT #3: Uses drug-specific clinical thresholds.
    
    CRITICAL: Deterministic - no random operations, no read dropping.
    """
    print(f"INFO: Loading results from '{json_path}'...")
    
    try:
        with open(json_path, 'r') as f:
            all_reads = json.load(f)
    except Exception as e:
        print(f"ERROR: Could not load JSON: {e}")
        return None
    
    if not all_reads:
        print("WARNING: Empty results file.")
        return None
    
    # Convert to DataFrame for analysis
    df = pd.DataFrame(all_reads)
    
    total_reads = len(df)
    
    # Filter for usable reads - DETERMINISTIC
    usable_reads_df = df[df['status'].isin(['DRM-Positive', 'Wild-Type', 'Polymorphic'])].copy()
    usable_reads_count = len(usable_reads_df)
    
    print(f"INFO: Total reads: {total_reads}")
    print(f"INFO: Usable reads: {usable_reads_count}")
    
    if usable_reads_count == 0:
        return {
            "total_reads": total_reads,
            "usable_reads_count": 0,
            "classification_counts": df['status'].value_counts().to_dict(),
            "drug_stats": {},
            "statistical_power_avg": 0.0,
            "quality_scores": [],
            "coverage_data": {},
            "gene_mutation_counts": {},
            "indel_counts_by_gene": {}
        }
    
    # Calculate average statistical power (will vary by drug threshold)
    # Use median threshold for overall power estimate
    all_thresholds = list(DRUG_SPECIFIC_THRESHOLDS.values())
    median_threshold = np.median(all_thresholds) if all_thresholds else DEFAULT_CLINICAL_THRESHOLD
    print(f"INFO: Calculating statistical power (median threshold: {median_threshold*100:.1f}%)...")
    statistical_power_avg = calculate_statistical_power(usable_reads_count, median_threshold)
    print(f"INFO: Average statistical power = {statistical_power_avg:.2f}")
    
    # Extract quality scores for plotting
    quality_scores = df[df['quality_score'].notna()]['quality_score'].tolist()
    
    # Calculate coverage
    print("INFO: Calculating genome coverage...")
    coverage_data = calculate_coverage_by_position(all_reads)
    
    # Calculate mutations by gene
    gene_mutation_counts = calculate_mutation_counts_by_gene(all_reads)
    
    # NEW: Calculate indel counts
    indel_counts_by_gene = calculate_indel_counts_by_gene(all_reads)
    
    # Drug resistance analysis
    print(f"INFO: Analyzing drug resistance in {usable_reads_count} reads...")
    drug_counts = {}
    
    for _, row in usable_reads_df.iterrows():
        mutations_by_gene = row['mutations_by_gene']
        if not mutations_by_gene:
            continue
        
        drugs_for_this_read = set()
        
        # Check single-gene mutations
        for gene, mutation_list in mutations_by_gene.items():
            if gene not in drm_db or not mutation_list:
                continue
            
            mutation_set = set(mutation_list)
            subsets_to_check = get_mutation_subsets(mutation_set)
            
            for subset_key in subsets_to_check:
                if subset_key in drm_db[gene]:
                    for drug_entry in drm_db[gene][subset_key]['drugs']:
                        drugs_for_this_read.add(drug_entry['name'])
        
        # Check cross-gene combinations
        if len(mutations_by_gene) > 1 and drm_db.get("CROSS_GENE"):
            all_cross_mutations = []
            for gene, muts in mutations_by_gene.items():
                for mut in muts:
                    all_cross_mutations.append(f"{gene}:{mut}")
            
            cross_subsets = get_mutation_subsets(set(all_cross_mutations))
            for cross_key in cross_subsets:
                if cross_key in drm_db["CROSS_GENE"]:
                    for drug_entry in drm_db["CROSS_GENE"][cross_key]['drugs']:
                        drugs_for_this_read.add(drug_entry['name'])
        
        # Increment counts
        for drug_name in drugs_for_this_read:
            drug_counts[drug_name] = drug_counts.get(drug_name, 0) + 1
    
    # Calculate prevalence
    drug_prevalence = {drug: count / usable_reads_count for drug, count in drug_counts.items()}
    
    # IMPROVEMENT #3: Perform statistical tests with drug-specific thresholds
    drug_list = sorted(drug_prevalence.keys())
    p_values_raw = []
    drug_stats_temp = {}
    
    for drug in drug_list:
        p_obs = drug_prevalence[drug]
        count = drug_counts[drug]
        
        # Get drug-specific threshold
        drug_threshold = get_drug_threshold(drug)
        
        # Calculate drug-specific statistical power
        drug_power = calculate_statistical_power(usable_reads_count, drug_threshold)
        
        # Use exact test for low counts, Z-test otherwise
        if count < 30:
            p_value = binomial_exact_test(count, usable_reads_count, drug_threshold)
            test_used = "Exact"
        else:
            p_hyp = drug_threshold
            variance = (p_hyp * (1 - p_hyp)) / usable_reads_count
            z_score = (p_obs - p_hyp) / np.sqrt(variance)
            p_value = norm.sf(z_score)
            test_used = "Z-test"
        
        p_values_raw.append(p_value)
        
        # Calculate confidence intervals
        ci_lower, ci_upper = calculate_wilson_ci(count, usable_reads_count)
        
        drug_stats_temp[drug] = {
            "count": count,
            "prevalence": p_obs,
            "ci_lower": ci_lower,
            "ci_upper": ci_upper,
            "p_value_raw": p_value,
            "test_used": test_used,
            "threshold": drug_threshold,  # NEW: Store drug-specific threshold
            "power": drug_power  # NEW: Store drug-specific power
        }
    
    # FDR correction for multiple testing
    print("INFO: Applying FDR correction...")
    if p_values_raw:
        reject, p_adjusted, _, _ = multipletests(p_values_raw, alpha=0.05, method='fdr_bh')
        
        for i, drug in enumerate(drug_list):
            drug_stats_temp[drug]['p_value_adjusted'] = p_adjusted[i]
            # IMPROVEMENT #3: Use drug-specific threshold for significance
            drug_stats_temp[drug]['is_significant'] = (
                reject[i] and 
                (drug_stats_temp[drug]['prevalence'] > drug_stats_temp[drug]['threshold'])
            )
    
    analysis_results = {
        "total_reads": total_reads,
        "usable_reads_count": usable_reads_count,
        "statistical_power_avg": statistical_power_avg,
        "classification_counts": df['status'].value_counts().to_dict(),
        "drug_stats": drug_stats_temp,
        "quality_scores": quality_scores,
        "coverage_data": coverage_data,
        "gene_mutation_counts": gene_mutation_counts,
        "indel_counts_by_gene": indel_counts_by_gene
    }
    
    print(f"INFO: Analysis complete. Found resistance to {len(drug_stats_temp)} drugs.")
    
    return analysis_results


# ============================================================================
# PART 5: ENHANCED PLOTTING FUNCTIONS
# ============================================================================

def generate_plots(results: Dict, base_name: str) -> Dict[str, str]:
    """Generates all visualization plots including new diagnostic plots."""
    print("INFO: Generating visualization plots...")
    plot_paths = {}
    
    # Plot 1: Read Classification Pie Chart
    plt.figure(figsize=(10, 7))
    counts = results['classification_counts']
    filtered_counts = {k: v for k, v in counts.items() if v > 0}
    
    if filtered_counts:
        labels = list(filtered_counts.keys())
        sizes = list(filtered_counts.values())
        colors = plt.cm.Set3(range(len(labels)))
        
        plt.pie(sizes, labels=labels, autopct='%1.1f%%', colors=colors,
                startangle=90, pctdistance=0.85)
        plt.title("Read Classification Breakdown", fontsize=16, fontweight='bold')
        plt.axis('equal')
    
    pie_path = f"{base_name}_classification.png"
    plt.savefig(pie_path, dpi=150, bbox_inches='tight')
    plt.close()
    plot_paths['pie_chart'] = pie_path
    
    # Plot 2: Drug Resistance with Error Bars (color-coded by threshold)
    plt.figure(figsize=(12, 10))
    drug_stats = results['drug_stats']
    
    if drug_stats:
        sorted_drugs = sorted(drug_stats.items(), key=lambda x: x[1]['prevalence'], reverse=True)
        drug_names = [d[0] for d in sorted_drugs]
        prevalences = [d[1]['prevalence'] * 100 for d in sorted_drugs]
        ci_lowers = [d[1]['ci_lower'] for d in sorted_drugs]
        ci_uppers = [d[1]['ci_upper'] for d in sorted_drugs]
        is_significant = [d[1]['is_significant'] for d in sorted_drugs]
        thresholds = [d[1]['threshold'] * 100 for d in sorted_drugs]  # NEW
        
        yerr_lower = [prevalences[i] - ci_lowers[i] for i in range(len(prevalences))]
        yerr_upper = [ci_uppers[i] - prevalences[i] for i in range(len(prevalences))]
        
        colors = ['#d62728' if sig else '#1f77b4' for sig in is_significant]
        
        bars = plt.barh(drug_names, prevalences, color=colors, xerr=[yerr_lower, yerr_upper],
                       capsize=3, error_kw={'elinewidth': 1, 'alpha': 0.7})
        
        # Draw drug-specific threshold lines
        for i, (drug_name, threshold) in enumerate(zip(drug_names, thresholds)):
            plt.plot([threshold, threshold], [i-0.4, i+0.4], 
                    color='black', linestyle='--', linewidth=1.5, alpha=0.6)
        
        plt.xlabel("Prevalence (%) with 95% CI", fontsize=12, fontweight='bold')
        plt.ylabel("Antiretroviral Drug", fontsize=12, fontweight='bold')
        plt.title("Drug Resistance Prevalence (with Drug-Specific Thresholds)", 
                 fontsize=16, fontweight='bold')
        plt.grid(axis='x', alpha=0.3)
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#d62728', label='Significant'),
            Patch(facecolor='#1f77b4', label='Not Significant'),
            plt.Line2D([0], [0], color='black', linestyle='--', label='Drug Threshold')
        ]
        plt.legend(handles=legend_elements, loc='lower right')
        
        plt.tight_layout()
    
    drug_path = f"{base_name}_drug_prevalence.png"
    plt.savefig(drug_path, dpi=150, bbox_inches='tight')
    plt.close()
    plot_paths['drug_chart'] = drug_path
    
    # Plot 3: Q-Score Distribution
    plt.figure(figsize=(10, 6))
    quality_scores = results.get('quality_scores', [])
    
    if quality_scores:
        plt.hist(quality_scores, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
        plt.axvline(x=20, color='red', linestyle='--', linewidth=2, label='Min Quality Threshold (Q20)')
        plt.xlabel("Mapping Quality Score (MAPQ)", fontsize=12, fontweight='bold')
        plt.ylabel("Read Count", fontsize=12, fontweight='bold')
        plt.title("Q-Score Distribution", fontsize=16, fontweight='bold')
        plt.legend()
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
    else:
        plt.text(0.5, 0.5, "No Quality Score Data Available", 
                ha='center', va='center', fontsize=14)
    
    qscore_path = f"{base_name}_qscore_dist.png"
    plt.savefig(qscore_path, dpi=150, bbox_inches='tight')
    plt.close()
    plot_paths['qscore_chart'] = qscore_path
    
    # Plot 4: Coverage Depth
    plt.figure(figsize=(14, 6))
    coverage_data = results.get('coverage_data', {})
    
    if coverage_data:
        positions = sorted(coverage_data.keys())
        depths = [coverage_data[pos] for pos in positions]
        
        plt.fill_between(positions, depths, alpha=0.4, color='green')
        plt.plot(positions, depths, color='darkgreen', linewidth=0.5)
        
        # Mark gene boundaries
        for gene, coords in GENE_COORDINATES.items():
            plt.axvline(x=coords['start'], color='red', linestyle=':', alpha=0.5)
            plt.axvline(x=coords['end'], color='red', linestyle=':', alpha=0.5)
            mid_point = (coords['start'] + coords['end']) / 2
            plt.text(mid_point, max(depths) * 0.9, gene, ha='center', fontsize=10, fontweight='bold')
        
        plt.xlabel("Genome Position (bp)", fontsize=12, fontweight='bold')
        plt.ylabel("Coverage Depth", fontsize=12, fontweight='bold')
        plt.title("Sequencing Coverage Across HIV Genome", fontsize=16, fontweight='bold')
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
    else:
        plt.text(0.5, 0.5, "No Coverage Data Available", 
                ha='center', va='center', fontsize=14)
    
    coverage_path = f"{base_name}_coverage.png"
    plt.savefig(coverage_path, dpi=150, bbox_inches='tight')
    plt.close()
    plot_paths['coverage_chart'] = coverage_path
    
    # Plot 5: Mutation Count by Gene Region (with indels)
    plt.figure(figsize=(12, 6))
    gene_mutation_counts = results.get('gene_mutation_counts', {})
    indel_counts = results.get('indel_counts_by_gene', {})
    
    if gene_mutation_counts and any(gene_mutation_counts.values()):
        genes = list(gene_mutation_counts.keys())
        snp_counts = [gene_mutation_counts[g] for g in genes]
        insertion_counts = [indel_counts.get(g, {}).get('insertions', 0) for g in genes]
        deletion_counts = [indel_counts.get(g, {}).get('deletions', 0) for g in genes]
        
        x = np.arange(len(genes))
        width = 0.25
        
        plt.bar(x - width, snp_counts, width, label='SNPs', color='#3498db', edgecolor='black')
        plt.bar(x, insertion_counts, width, label='Insertions', color='#2ecc71', edgecolor='black')
        plt.bar(x + width, deletion_counts, width, label='Deletions', color='#e74c3c', edgecolor='black')
        
        plt.xlabel("Gene Region", fontsize=12, fontweight='bold')
        plt.ylabel("Mutation Count", fontsize=12, fontweight='bold')
        plt.title("Detected Mutations by Gene Region (SNPs + Indels)", fontsize=16, fontweight='bold')
        plt.xticks(x, genes)
        plt.legend()
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
    else:
        plt.text(0.5, 0.5, "No Mutations Detected", 
                ha='center', va='center', fontsize=14)
    
    gene_mut_path = f"{base_name}_gene_mutations.png"
    plt.savefig(gene_mut_path, dpi=150, bbox_inches='tight')
    plt.close()
    plot_paths['gene_mutations_chart'] = gene_mut_path
    
    print("INFO: All plots generated successfully.")
    return plot_paths


# ============================================================================
# PART 6: PDF REPORT WITH DETAILED EXPLANATIONS
# ============================================================================

class PDF(FPDF):
    """Custom PDF with enhanced reporting."""
    
    def __init__(self, *args, report_time="N/A", **kwargs):
        super().__init__(*args, **kwargs)
        self.report_time = report_time
    
    def header(self):
        self.set_font("Helvetica", "B", 12)
        self.cell(0, 10, "HIV Drug Resistance Mutation Analysis Report", align="C", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 8)
        self.cell(0, 6, f"Generated: {self.report_time}", align="C", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.ln(10)
    
    def footer(self):
        self.set_y(-15)
        self.set_font("Helvetica", "I", 8)
        self.cell(0, 10, f"Page {self.page_no()}", align="C")
    
    def chapter_title(self, title: str):
        self.set_font("Helvetica", "B", 14)
        self.cell(0, 10, title, new_x=XPos.LMARGIN, new_y=YPos.NEXT, border="B")
        self.ln(5)
    
    def add_methodology_page(self, results: Dict):
        """Enhanced methodology with detailed explanations."""
        self.add_page()
        
        self.set_font("Helvetica", "B", 18)
        self.cell(0, 12, "Statistical Methodology Guide", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
        self.ln(8)
        
        # Introduction
        self.set_font("Helvetica", "", 11)
        self.multi_cell(0, 6,
            "This report uses four advanced statistical methods combined with drug-specific clinical "
            "thresholds to ensure accurate, clinically reliable results. This page explains each method "
            "in plain English, why we use it, and how to interpret the results."
        )
        self.ln(5)
        
        # Method 1: Confidence Intervals
        self.set_font("Helvetica", "B", 13)
        self.cell(0, 8, "1. Confidence Intervals (95% CI)", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "B", 11)
        self.cell(0, 7, "What it is:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 10)
        self.multi_cell(0, 5,
            "A confidence interval gives us a range where the TRUE resistance prevalence likely falls. "
            "Think of it as a 'margin of error' around our measurement."
        )
        self.ln(2)
        
        self.set_font("Helvetica", "B", 11)
        self.cell(0, 7, "Why we use it:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 10)
        self.multi_cell(0, 5,
            "A single percentage (e.g., '2.5% resistance') doesn't tell us how CONFIDENT we are. "
            "The CI tells us: 'We're 95% sure the true value is between X% and Y%'."
        )
        self.ln(2)
        
        self.set_font("Helvetica", "B", 11)
        self.cell(0, 7, "How to interpret it:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 10)
        self.multi_cell(0, 5,
            "Example: NVP shows 2.5% resistance with CI [1.8% - 3.4%]\n"
            "  - NARROW interval (1.8-3.4): High confidence, good data\n"
            "  - If CI crosses threshold: Result is UNCERTAIN\n"
            "  - If entire CI is above threshold: Strong evidence of resistance"
        )
        self.ln(4)
        
        # Method 2: Multiple Testing Correction
        self.set_font("Helvetica", "B", 13)
        self.cell(0, 8, "2. Multiple Testing Correction (FDR)", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "B", 11)
        self.cell(0, 7, "What it is:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 10)
        self.multi_cell(0, 5,
            "When we test 27 drugs, we're doing 27 separate tests. This increases the chance "
            "of false positives. FDR correction adjusts p-values to control the overall error rate."
        )
        self.ln(2)
        
        self.set_font("Helvetica", "B", 11)
        self.cell(0, 7, "Why we use it:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 10)
        self.multi_cell(0, 5,
            "Without correction: Testing 27 drugs at 5% error rate each gives ~75% chance of at least "
            "ONE false positive! With FDR, we control the OVERALL error rate at 5%."
        )
        self.ln(2)
        
        self.set_font("Helvetica", "B", 11)
        self.cell(0, 7, "How to interpret it:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 10)
        self.multi_cell(0, 5,
            "Example: Drug shows p=0.03 (raw) but p=0.08 (FDR-adjusted)\n"
            "  - Raw test says 'significant' but FDR says 'not significant after correction'\n"
            "  - ALWAYS use the FDR-adjusted p-value for clinical decisions\n"
            f"  - This report tested {len(results.get('drug_stats', {}))} drugs simultaneously"
        )
        self.ln(4)
        
        # Method 3: Statistical Power
        self.set_font("Helvetica", "B", 13)
        self.cell(0, 8, "3. Statistical Power Analysis", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "B", 11)
        self.cell(0, 7, "What it is:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 10)
        self.multi_cell(0, 5,
            "Power tells us: IF resistance exists at the threshold level, what's the probability we'll "
            "detect it? It measures whether we have ENOUGH data to find resistance if it's there."
        )
        self.ln(2)
        
        self.set_font("Helvetica", "B", 11)
        self.cell(0, 7, "Why we use it:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 10)
        self.multi_cell(0, 5,
            "If we say 'no resistance found,' we need to prove we had enough data. Low power "
            "means: 'We might have missed it due to insufficient sequencing depth.'"
        )
        self.ln(2)
        
        self.set_font("Helvetica", "B", 11)
        self.cell(0, 7, "How to interpret it:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 10)
        power = results.get('statistical_power_avg', 0.0)
        power_interpretation = "EXCELLENT - Highly reliable" if power >= 0.8 else \
                              "MODERATE - Acceptable but borderline" if power >= 0.6 else \
                              "LOW - Need more sequencing"
        
        self.multi_cell(0, 5,
            f"This analysis: Average Power = {power:.2f} ({power*100:.0f}%) - {power_interpretation}\n"
            "  - Power >= 0.80: Gold standard, results trustworthy\n"
            "  - Power 0.60-0.79: Acceptable, but borderline\n"
            "  - Power < 0.60: CAUTION - May miss real resistance\n"
            "Note: Power varies by drug based on its specific threshold"
        )
        self.ln(4)
        
        # Method 4: Binomial Exact Test
        self.set_font("Helvetica", "B", 13)
        self.cell(0, 8, "4. Binomial Exact Test (for Rare Variants)", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "B", 11)
        self.cell(0, 7, "What it is:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 10)
        self.multi_cell(0, 5,
            "When a drug has fewer than 30 resistant reads, the standard Z-test becomes unreliable. "
            "The exact test uses precise mathematical formulas instead of approximations."
        )
        self.ln(2)
        
        self.set_font("Helvetica", "B", 11)
        self.cell(0, 7, "Why we use it:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 10)
        self.multi_cell(0, 5,
            "For rare variants (e.g., 5 reads out of 1,000), the Z-test can give WRONG p-values. "
            "The exact test is mathematically correct for all sample sizes."
        )
        self.ln(2)
        
        self.set_font("Helvetica", "B", 11)
        self.cell(0, 7, "How to interpret it:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 10)
        self.multi_cell(0, 5,
            "In the results table, you'll see which test was used:\n"
            "  - 'Exact': For counts < 30 (more accurate for rare variants)\n"
            "  - 'Z-test': For counts >= 30 (equally accurate, faster)\n"
            "Both test: 'Is prevalence significantly above the drug-specific threshold?'"
        )
        self.ln(6)
        
        # NEW: Drug-Specific Thresholds
        self.set_font("Helvetica", "B", 13)
        self.cell(0, 8, "5. Drug-Specific Clinical Thresholds (NEW)", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "B", 11)
        self.cell(0, 7, "What it is:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 10)
        self.multi_cell(0, 5,
            "Different drugs have different 'barriers to resistance.' This pipeline uses drug-specific "
            "thresholds based on clinical evidence rather than a one-size-fits-all cutoff."
        )
        self.ln(2)
        
        self.set_font("Helvetica", "B", 11)
        self.cell(0, 7, "Why we use it:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 10)
        self.multi_cell(0, 5,
            "Example thresholds:\n"
            "  - NNRTIs (e.g., Efavirenz): 1% (single mutation can cause resistance)\n"
            "  - NRTIs (e.g., 3TC): 20% (M184V highly prevalent but tolerable)\n"
            "  - PIs (e.g., Darunavir): 15% (requires multiple mutations)\n"
            "Using the correct threshold for each drug improves clinical relevance."
        )
        self.ln(2)
        
        self.set_font("Helvetica", "B", 11)
        self.cell(0, 7, "How to interpret it:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 10)
        self.multi_cell(0, 5,
            "In the results:\n"
            "  - Each drug's threshold is shown in the detailed table\n"
            "  - Drug-specific threshold lines appear on the prevalence plot\n"
            "  - Significance is determined relative to EACH drug's threshold"
        )
        self.ln(6)
        
        # Summary
        self.set_font("Helvetica", "B", 13)
        self.cell(0, 8, "Summary: Ensuring Clinical Accuracy", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 10)
        self.multi_cell(0, 5,
            "These five methods work together:\n"
            "  1. CIs quantify CONFIDENCE in each measurement\n"
            "  2. FDR correction prevents FALSE POSITIVES from multiple testing\n"
            "  3. Power analysis validates NEGATIVE results\n"
            "  4. Exact tests ensure ACCURACY for rare variants\n"
            "  5. Drug-specific thresholds provide CLINICAL RELEVANCE\n\n"
            "This approach meets publication standards and ensures results are clinically "
            "actionable and scientifically defensible."
        )


def create_pdf_report(
    results: Dict,
    output_path: str,
    plot_paths: Dict[str, str],
    report_time: str,
    duration_str: str
):
    """Creates comprehensive PDF report with all visualizations."""
    print(f"INFO: Generating PDF report...")
    
    pdf = PDF(report_time=report_time)
    pdf.set_auto_page_break(auto=True, margin=15)
    
    # Page 1: Methodology
    pdf.add_methodology_page(results)
    
    # Page 2: Executive Summary
    pdf.add_page()
    pdf.chapter_title("Executive Summary")
    
    pdf.set_font("Helvetica", "", 11)
    pdf.cell(0, 8, f"Total Reads Processed: {results['total_reads']:,}", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.cell(0, 8, f"Usable High-Quality Reads: {results['usable_reads_count']:,}", 
            new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.cell(0, 8, f"Analysis Duration: {duration_str}", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    
    usable_pct = (results['usable_reads_count'] / results['total_reads'] * 100) if results['total_reads'] > 0 else 0
    pdf.cell(0, 8, f"Usable Read Percentage: {usable_pct:.1f}%", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    
    # Statistical Power
    power = results.get('statistical_power_avg', 0.0)
    pdf.set_font("Helvetica", "B", 11)
    pdf.cell(0, 8, f"Average Statistical Power: {power:.2f} ({power*100:.0f}%)", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font("Helvetica", "I", 9)
    
    if power >= 0.8:
        pdf.cell(0, 6, "   Excellent - High confidence in detecting resistance at threshold levels", 
                new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    elif power >= 0.6:
        pdf.cell(0, 6, "   Moderate - Reasonable confidence, but borderline", 
                new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    else:
        pdf.cell(0, 6, "   Low - May miss low-frequency resistance. Consider more sequencing.", 
                new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    
    pdf.set_font("Helvetica", "", 11)
    pdf.ln(5)
    
    # Classification breakdown
    pdf.set_font("Helvetica", "B", 12)
    pdf.cell(0, 8, "Read Classification:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font("Helvetica", "", 11)
    
    for status, count in sorted(results['classification_counts'].items(), key=lambda x: x[1], reverse=True):
        pct = (count / results['total_reads'] * 100) if results['total_reads'] > 0 else 0
        pdf.cell(0, 6, f"  - {status}: {count:,} ({pct:.1f}%)", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    
    pdf.ln(5)
    
    if 'pie_chart' in plot_paths:
        pdf.image(plot_paths['pie_chart'], w=170)
    
    # Page 3: Sequencing Quality Summary
    pdf.add_page()
    pdf.chapter_title("Sequencing Quality Summary")
    
    pdf.set_font("Helvetica", "", 11)
    pdf.multi_cell(0, 6,
        "This section shows the quality metrics of the sequencing run, including mapping quality "
        "distribution and genome coverage depth. High-quality data ensures reliable mutation calling."
    )
    pdf.ln(5)
    
    # Q-Score Distribution
    if 'qscore_chart' in plot_paths:
        pdf.set_font("Helvetica", "B", 12)
        pdf.cell(0, 8, "Q-Score Distribution", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        pdf.set_font("Helvetica", "", 10)
        pdf.multi_cell(0, 5,
            "This plot shows the distribution of mapping quality scores (MAPQ). Higher scores indicate "
            "more confident alignments. Reads with MAPQ < 20 are filtered out."
        )
        pdf.ln(3)
        pdf.image(plot_paths['qscore_chart'], w=170)
    
    # Coverage Depth
    pdf.add_page()
    if 'coverage_chart' in plot_paths:
        pdf.set_font("Helvetica", "B", 12)
        pdf.cell(0, 8, "Coverage Depth Across Genome", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        pdf.set_font("Helvetica", "", 10)
        pdf.multi_cell(0, 5,
            "This plot shows sequencing depth at each position in the HIV genome. Vertical lines "
            "mark gene boundaries (PR, RT, IN). Uniform, high coverage ensures reliable mutation detection."
        )
        pdf.ln(3)
        pdf.image(plot_paths['coverage_chart'], w=180)
    
    # Page 4: Detected Mutations Summary
    pdf.add_page()
    pdf.chapter_title("Detected Mutations")
    
    if 'gene_mutations_chart' in plot_paths:
        pdf.set_font("Helvetica", "B", 12)
        pdf.cell(0, 8, "Mutation Count by Gene Region (SNPs + Indels)", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        pdf.set_font("Helvetica", "", 10)
        pdf.multi_cell(0, 5,
            "This chart shows the total number of mutations detected in each gene region, including "
            "single nucleotide polymorphisms (SNPs), insertions, and deletions. These include both "
            "drug-resistance mutations (DRMs) and natural polymorphisms."
        )
        pdf.ln(3)
        pdf.image(plot_paths['gene_mutations_chart'], w=170)
    
    # NEW: Indel Summary
    indel_counts = results.get('indel_counts_by_gene', {})
    if indel_counts and any(indel_counts.get(g, {}).get('insertions', 0) + 
                            indel_counts.get(g, {}).get('deletions', 0) > 0 
                            for g in ['PR', 'RT', 'IN']):
        pdf.ln(5)
        pdf.set_font("Helvetica", "B", 12)
        pdf.cell(0, 8, "Indel Detection Summary", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        pdf.set_font("Helvetica", "", 10)
        pdf.multi_cell(0, 5,
            "Insertions and deletions (indels) are important for HIV drug resistance. Notable examples "
            "include the T69 insertion complex (NRTI resistance) and various deletions in protease."
        )
        pdf.ln(2)
        
        pdf.set_font("Helvetica", "", 10)
        for gene in ['PR', 'RT', 'IN']:
            ins = indel_counts.get(gene, {}).get('insertions', 0)
            dels = indel_counts.get(gene, {}).get('deletions', 0)
            if ins > 0 or dels > 0:
                pdf.cell(0, 6, f"  {gene}: {ins} insertions, {dels} deletions", 
                        new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    
    # Page 5: Drug Resistance Analysis
    if results['drug_stats']:
        pdf.add_page()
        pdf.chapter_title("Drug Resistance Analysis")
        
        pdf.set_font("Helvetica", "", 11)
        pdf.multi_cell(0, 6,
            "The following analysis shows resistance prevalence with 95% confidence intervals. "
            "P-values are FDR-corrected for multiple testing. Each drug uses its specific clinical "
            "threshold (shown in the table). Drugs marked as 'Significant' have FDR-adjusted p < 0.05 "
            "AND prevalence above their drug-specific threshold."
        )
        pdf.ln(5)
        
        if 'drug_chart' in plot_paths:
            pdf.image(plot_paths['drug_chart'], w=170)
        
        # Detailed table
        pdf.add_page()
        pdf.chapter_title("Detailed Statistical Results")
        
        pdf.set_font("Helvetica", "B", 7)
        pdf.cell(30, 8, "Drug", border=1, align="C")
        pdf.cell(15, 8, "Count", border=1, align="C")
        pdf.cell(20, 8, "Prevalence", border=1, align="C")
        pdf.cell(25, 8, "95% CI", border=1, align="C")
        pdf.cell(20, 8, "Threshold", border=1, align="C")
        pdf.cell(20, 8, "Power", border=1, align="C")
        pdf.cell(20, 8, "P (FDR)", border=1, align="C")
        pdf.cell(15, 8, "Test", border=1, align="C")
        pdf.cell(20, 8, "Significant", border=1, align="C")
        pdf.ln()
        
        pdf.set_font("Helvetica", "", 7)
        
        for drug, stats in sorted(results['drug_stats'].items(), 
                                  key=lambda x: x[1]['prevalence'], reverse=True):
            pdf.cell(30, 6, drug, border=1)
            pdf.cell(15, 6, str(stats['count']), border=1, align="C")
            pdf.cell(20, 6, f"{stats['prevalence']*100:.2f}%", border=1, align="C")
            pdf.cell(25, 6, f"[{stats['ci_lower']:.1f}-{stats['ci_upper']:.1f}]", border=1, align="C")
            pdf.cell(20, 6, f"{stats['threshold']*100:.1f}%", border=1, align="C")
            pdf.cell(20, 6, f"{stats['power']:.2f}", border=1, align="C")
            pdf.cell(20, 6, f"{stats['p_value_adjusted']:.2e}", border=1, align="C")
            pdf.cell(15, 6, stats['test_used'], border=1, align="C")
            pdf.cell(20, 6, "YES" if stats['is_significant'] else "NO", border=1, align="C")
            pdf.ln()
        
        # Clinical Interpretation
        pdf.add_page()
        pdf.chapter_title("Clinical Interpretation")
        
        significant_drugs = [drug for drug, stats in results['drug_stats'].items() 
                           if stats['is_significant']]
        
        if significant_drugs:
            pdf.set_font("Helvetica", "B", 12)
            pdf.cell(0, 8, f"Significant Resistance Detected: {len(significant_drugs)} Drug(s)", 
                    new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            pdf.ln(3)
            
            pdf.set_font("Helvetica", "", 11)
            pdf.multi_cell(0, 6,
                "The following drugs showed statistically significant resistance after FDR correction, "
                "with prevalence exceeding their drug-specific clinical thresholds. These drugs are "
                "likely to have reduced clinical efficacy:"
            )
            pdf.ln(3)
            
            for drug in sorted(significant_drugs):
                stats = results['drug_stats'][drug]
                pdf.set_font("Helvetica", "B", 11)
                pdf.cell(0, 7, f"  * {drug}", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
                pdf.set_font("Helvetica", "", 10)
                pdf.cell(0, 6, 
                        f"      Prevalence: {stats['prevalence']*100:.2f}% "
                        f"[CI: {stats['ci_lower']:.1f}-{stats['ci_upper']:.1f}%]", 
                        new_x=XPos.LMARGIN, new_y=YPos.NEXT)
                pdf.cell(0, 6, 
                        f"      Threshold: {stats['threshold']*100:.1f}% | "
                        f"Power: {stats['power']:.2f} | "
                        f"FDR p: {stats['p_value_adjusted']:.2e}", 
                        new_x=XPos.LMARGIN, new_y=YPos.NEXT)
                pdf.ln(2)
        else:
            pdf.set_font("Helvetica", "B", 12)
            pdf.cell(0, 8, "No Significant Resistance Detected", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            pdf.ln(3)
            pdf.set_font("Helvetica", "", 11)
            pdf.multi_cell(0, 6,
                "No drug resistance mutations were found at statistically significant levels after "
                "FDR correction. All detected mutations are below their drug-specific thresholds or "
                "not statistically significant.\n\n"
                f"Note: With average statistical power of {results.get('statistical_power_avg', 0)*100:.0f}%, "
                "this analysis " + ("IS" if results.get('statistical_power_avg', 0) >= 0.8 else "may NOT be") + 
                " adequately powered to detect low-frequency resistance."
            )
    else:
        pdf.add_page()
        pdf.chapter_title("Drug Resistance Analysis")
        pdf.set_font("Helvetica", "", 11)
        pdf.multi_cell(0, 6,
            "No drug resistance mutations were detected in the analyzed reads."
        )
    
    pdf.output(output_path)
    print(f"INFO: PDF report saved to '{output_path}'.")


# ============================================================================
# PART 7: MAIN PIPELINE EXECUTION
# ============================================================================

def main():
    """Main pipeline orchestrator with validation checks."""
    parser = argparse.ArgumentParser(
        description="HIV Drug Resistance Mutation Analysis Pipeline - IMPROVED VERSION\n"
                    "  * Now detects insertions and deletions (indels)\n"
                    "  * Fixed IN gene support (positions 561-849)\n"
                    "  * Drug-specific clinical thresholds\n"
                    "  * Publication-quality statistical analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("-i", "--input_bam", required=True, help="Path to input BAM file")
    parser.add_argument("--ref_fasta", required=True, help="Path to reference FASTA (HXB2)")
    parser.add_argument("--profile_csv", required=True, help="Path to resistance profile CSV")
    parser.add_argument("-n", "--num_reads", type=int, default=0, help="Number of reads (0=all)")
    parser.add_argument("-o", "--output_dir", default="data", help="Output directory")
    
    args = parser.parse_args()
    
    # Validate inputs
    for path, name in [(args.input_bam, "BAM"), (args.ref_fasta, "FASTA"), (args.profile_csv, "CSV")]:
        if not os.path.exists(path):
            print(f"ERROR: {name} file not found: {path}", file=sys.stderr)
            sys.exit(1)
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    start_time = datetime.datetime.now()
    timestamp_file = start_time.strftime("%Y%m%d_%H%M%S")
    timestamp_report = start_time.strftime("%Y-%m-%d %H:%M:%S")
    
    print("="*80)
    print(" HIV DRUG RESISTANCE MUTATION ANALYSIS PIPELINE - IMPROVED VERSION")
    print(" ✓ Indel Detection (Insertions + Deletions)")
    print(" ✓ Complete Gene Support (PR, RT, IN)")
    print(" ✓ Drug-Specific Clinical Thresholds")
    print(" ✓ Publication-Quality Statistics")
    print("="*80)
    print(f"Started: {timestamp_report}")
    print(f"Input BAM: {args.input_bam}")
    print(f"Reference: {args.ref_fasta}")
    print(f"Profile DB: {args.profile_csv}")
    print("="*80)
    print()
    
    # STEP 1: Build database
    print("STEP 1/4: Building drug resistance database...")
    drm_db = build_definitive_database(args.profile_csv)
    print()
    
    # STEP 2: Load reference
    print("STEP 2/4: Loading reference genome...")
    ref_sequences = load_reference_fasta(args.ref_fasta)
    print(f"INFO: Loaded {len(ref_sequences)} reference sequence(s).")
    print()
    
    # STEP 3: Process reads
    print("STEP 3/4: Processing reads (with indel detection)...")
    
    try:
        with pysam.AlignmentFile(args.input_bam, "rb") as bam_count:
            total_reads_in_bam = bam_count.count(until_eof=True)
    except Exception as e:
        print(f"ERROR: Could not read BAM: {e}", file=sys.stderr)
        sys.exit(1)
    
    num_to_process = total_reads_in_bam if args.num_reads == 0 else min(args.num_reads, total_reads_in_bam)
    print(f"INFO: Processing {num_to_process:,} reads...")
    
    all_results = []
    progress_interval = max(1, num_to_process // 20)
    
    try:
        with pysam.AlignmentFile(args.input_bam, "rb") as bamfile:
            for i, read in enumerate(bamfile):
                if i >= num_to_process:
                    break
                
                result = process_read_for_mutations(read, ref_sequences, drm_db)
                all_results.append(result)
                
                if (i + 1) % progress_interval == 0:
                    progress_pct = ((i + 1) / num_to_process) * 100
                    print(f"  Progress: {i+1:,}/{num_to_process:,} reads ({progress_pct:.1f}%)")
    except Exception as e:
        print(f"ERROR: Failed during processing: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"INFO: Completed processing {len(all_results):,} reads.")
    print()
    
    # Save results 
    base_name = os.path.basename(args.input_bam).replace(".sorted.bam", "").replace(".bam", "")
    file_base = f"{base_name}_{timestamp_file}"
    json_path = os.path.join(args.output_dir, f"{file_base}_results.json")
    pdf_path = os.path.join(args.output_dir, f"{file_base}_report.pdf")
    plot_base = os.path.join(args.output_dir, file_base)
    
    print(f"INFO: Saving results to '{json_path}'...")
    with open(json_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print()
    
    # STEP 4: Statistical analysis
    print("STEP 4/4: Performing enhanced statistical analysis...")
    analysis_results = perform_statistical_analysis(json_path, drm_db)
    
    if not analysis_results or analysis_results['usable_reads_count'] == 0:
        print("WARNING: No usable data for report.")
        print("\nPipeline completed with no actionable results.")
        return
    
    print()
    
    # Validation check
    print("INFO: Validation check - comparing intermediate results...")
    usable_count = analysis_results['usable_reads_count']
    drm_positive = analysis_results['classification_counts'].get('DRM-Positive', 0)
    polymorphic = analysis_results['classification_counts'].get('Polymorphic', 0)
    print(f"  - Usable reads: {usable_count}")
    print(f"  - DRM-Positive: {drm_positive}")
    print(f"  - Polymorphic: {polymorphic}")
    print("  - These counts should be IDENTICAL across runs with same input")
    
    # NEW: Show indel counts
    indel_counts = analysis_results.get('indel_counts_by_gene', {})
    total_insertions = sum(indel_counts.get(g, {}).get('insertions', 0) for g in ['PR', 'RT', 'IN'])
    total_deletions = sum(indel_counts.get(g, {}).get('deletions', 0) for g in ['PR', 'RT', 'IN'])
    print(f"  - Total insertions detected: {total_insertions}")
    print(f"  - Total deletions detected: {total_deletions}")
    print()
    
    # Calculate duration
    end_time = datetime.datetime.now()
    duration = end_time - start_time
    minutes = int(duration.total_seconds() // 60)
    seconds = int(duration.total_seconds() % 60)
    duration_str = f"{minutes}m {seconds}s"
    
    # Generate plots
    print("INFO: Generating all diagnostic and analysis plots...")
    plot_paths = generate_plots(analysis_results, plot_base)
    
    # Generate PDF
    create_pdf_report(analysis_results, pdf_path, plot_paths, timestamp_report, duration_str)
    
    # Cleanup temporary plots
    for plot_path in plot_paths.values():
        if os.path.exists(plot_path):
            os.remove(plot_path)
    
    # Final summary
    print()
    print("="*80)
    print(" PIPELINE COMPLETED SUCCESSFULLY")
    print("="*80)
    print(f"Total Runtime: {duration_str}")
    print(f"Results JSON: {json_path}")
    print(f"Final Report: {pdf_path}")
    
    if analysis_results['drug_stats']:
        significant_count = sum(1 for stats in analysis_results['drug_stats'].values() 
                               if stats['is_significant'])
        print(f"\nDrugs Analyzed: {len(analysis_results['drug_stats'])}")
        print(f"Significant Resistance Found: {significant_count}")
        
        if significant_count > 0:
            print("\nSignificant Drugs (with drug-specific thresholds):")
            for drug, stats in sorted(analysis_results['drug_stats'].items(), 
                                     key=lambda x: x[1]['prevalence'], reverse=True):
                if stats['is_significant']:
                    print(f"  - {drug}: {stats['prevalence']*100:.2f}% "
                          f"[{stats['ci_lower']:.1f}-{stats['ci_upper']:.1f}%], "
                          f"threshold={stats['threshold']*100:.1f}%, "
                          f"p={stats['p_value_adjusted']:.2e}")
    
    print(f"\nAverage Statistical Power: {analysis_results['statistical_power_avg']:.2f}")
    print("="*80)


if __name__ == "__main__":
    main()