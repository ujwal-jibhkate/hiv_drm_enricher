import pandas as pd
import json
import sys
import re
import argparse
import pysam
import os
import numpy as np
from scipy.stats import norm
from fpdf import FPDF, XPos, YPos
import matplotlib.pyplot as plt
import matplotlib
import datetime
import itertools

# --- Configuration & Helper Functions ---
matplotlib.use('Agg')

# These are the CORRECT HXB2 coordinates (K03455.1)
GENE_COORDINATES = {
    "PR": {"start": 2253, "end": 2550},
    "RT": {"start": 2550, "end": 4230},
    "IN": {"start": 4230, "end": 5096}
}

CODON_TABLE = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

# --- PART 1: DATABASE BUILDER ---

def get_mutation_subsets(mutation_set: set) -> set:
    """
    Generates all non-empty, canonical-keyed subsets of mutations
    UP TO A CLINICALLY-RELEVANT LIMIT.
    """
    subsets_as_keys = set()
    mut_list = sorted(list(mutation_set)) 
    
    # --- THIS IS THE FIX ---
    # We will only check for combinations of size 1, 2, 3, and 4.
    # This prevents combinatorial explosion from noisy reads.
    MAX_COMBO_SIZE = 4 
    
    # We take the smaller of the list size or our cap
    max_k_to_check = min(len(mut_list), MAX_COMBO_SIZE)
    
    for k in range(1, max_k_to_check + 1):
    # --- END OF FIX ---
        for subset_tuple in itertools.combinations(mut_list, k):
            canonical_key = ",".join(subset_tuple)
            subsets_as_keys.add(canonical_key)
            
    return subsets_as_keys

def build_definitive_database(csv_path: str) -> dict:
    """
    Builds a definitive database from BOTH single mutations and combinations.
    It validates all mutations in a combo and stores them under a
    single, sorted, canonical key.
    """
    print(f"INFO: Building definitive database from '{csv_path}'...")
    try:
        df = pd.read_csv(csv_path, dtype=str)
    except FileNotFoundError:
        print(f"FATAL: Source CSV not found at '{csv_path}'.", file=sys.stderr)
        sys.exit(1)
        
    df_long = pd.melt(df, id_vars=['mutations'], var_name='Drug', value_name='Score')
    df_long = df_long[df_long['Score'] != 'S'].dropna(subset=['Score'])
    
    drm_db = {"PR": {}, "RT": {}, "IN": {}}
    mutation_pattern = re.compile(r'([A-Z]?)(\d+)([A-Z])', re.IGNORECASE)

    for _, row in df_long.iterrows():
        cleaned_str = re.sub(r"['\"\[\]\s]", "", str(row['mutations']))
        individual_mutations_raw = cleaned_str.split(',')
        
        if not individual_mutations_raw: continue
        
        parsed_muts = []
        gene_for_combo = None
        is_valid_combo = True
        
        for mut in individual_mutations_raw:
            if not mut: continue
            match = mutation_pattern.match(mut)
            if not match:
                is_valid_combo = False; break
            
            ref_aa, position_str, alt_aa = match.groups()
            position = int(position_str)
            
            # Determine gene
            current_gene = "PR" if 1 <= position <= 99 else "RT" if 100 <= position <= 560 else None
            if not current_gene:
                is_valid_combo = False; break
            
            # Check for cross-gene combos (e.g., PR + RT), which are invalid
            if gene_for_combo is None:
                gene_for_combo = current_gene
            elif gene_for_combo != current_gene:
                is_valid_combo = False; break # This combo mixes genes
                
            # Reconstruct the canonical mutation string (e.g., ?103N)
            ref_aa_final = ref_aa.upper() if ref_aa else '?'
            aa_change_str = f"{ref_aa_final}{position}{alt_aa.upper()}"
            parsed_muts.append(aa_change_str)

        if is_valid_combo and gene_for_combo and parsed_muts:
            # Create the final canonical key, e.g., "K103N,M184V"
            parsed_muts.sort()
            mutation_key = ",".join(parsed_muts)
            
            drug_info = {'name': row['Drug'], 'score': row['Score']}
            if mutation_key not in drm_db[gene_for_combo]:
                drm_db[gene_for_combo][mutation_key] = {'drugs': []}
            
            # De-duplicate drugs at the source
            if drug_info not in drm_db[gene_for_combo][mutation_key]['drugs']:
                drm_db[gene_for_combo][mutation_key]['drugs'].append(drug_info)

    print("INFO: Database build complete.")
    return drm_db
# --- PART 2: SIMULATION PIPELINE ---

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
            aa_position = (codon_start // 3) + 1
            if ref_aa != mut_aa:
                return ref_aa, mut_aa, aa_position, gene
    return None

def process_read_for_mutations(read, ref_sequences) -> dict:
    """
    Finds ALL amino acid mutations on a single read.
    It no longer checks the DRM database.
    Returns a dict with the read's status and a list of found mutations.
    """
    if read.is_unmapped: 
        return {"read_id": read.query_name, "status": "Non-HIV", "mutations_by_gene": {}}
    
    ref_seq = ref_sequences.get(read.reference_name)
    if not ref_seq and len(ref_sequences) == 1:
        _, ref_seq = next(iter(ref_sequences.items()))
    elif not ref_seq: 
        return {"read_id": read.query_name, "status": "Reference-Not-Found", "mutations_by_gene": {}}

    if read.mapping_quality < 20: 
        return {"read_id": read.query_name, "status": "Low-Quality", "mutations_by_gene": {}}
    
    found_aa_changes = {"PR": set(), "RT": set(), "IN": set()}
    
    aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=True)
    for read_pos, ref_pos, ref_base in aligned_pairs:
        if read_pos is not None and ref_pos is not None and ref_base is not None:
            read_base = read.query_sequence[read_pos]
            if read_base.upper() != ref_base.upper():
                aa_change = get_amino_acid_change(ref_pos, read_base.upper(), ref_seq)
                
                if aa_change:
                    ref_aa, mut_aa, aa_pos, gene = aa_change
                    aa_change_full = f"{ref_aa}{aa_pos}{mut_aa}"
                    aa_change_partial = f"?{aa_pos}{mut_aa}"
                    
                    found_aa_changes[gene].add(aa_change_full)
                    found_aa_changes[gene].add(aa_change_partial)
                    
    # --- THIS IS THE NEW NOISE-FILTERING LOGIC ---
    # A generous cap. A real read should not have this many mutations.
    NOISE_THRESHOLD = 30 
    total_mutations_found = sum(len(s) for s in found_aa_changes.values())

    if total_mutations_found > NOISE_THRESHOLD:
        return {"read_id": read.query_name, "status": "Noisy-Read-Filtered", "mutations_by_gene": {}}
    # --- END OF NEW LOGIC ---

    # Now, determine the read's final status
    status = "DRM-Negative"
    if total_mutations_found > 0:
        status = "DRM-Positive" # This just means "has a mutation"

    # Convert sets to lists for JSON serialization
    final_mutations = {
        gene: list(mut_set)
        for gene, mut_set in found_aa_changes.items()
        if mut_set # Only include genes that have mutations
    }

    return {"read_id": read.query_name, "status": status, "mutations_by_gene": final_mutations}

# --- PART 3: ANALYSIS & PDF REPORTING ---

def perform_statistical_analysis(json_path: str, drm_db: dict, threshold: float = 0.01):
    """
    Loads the JSON results, performs statistical analysis, and returns a results dictionary.
    
    This version correctly handles the new data structure by:
    1. Getting all mutations from a read (e.g., {"K103N", "M184V"})
    2. Generating all subsets (e.g., "K103N", "M184V", "K103N,M184V")
    3. Checking all subsets against the combination-aware database.
    """
    print(f"INFO: Loading results from {json_path} for analysis...")
    df = pd.read_json(json_path)
    
    if df.empty:
        print("WARNING: The results JSON is empty. No analysis can be performed.")
        return None

    # 1. Overall Metrics
    total_reads = len(df)
    usable_reads_df = df[df['status'].isin(['DRM-Positive', 'DRM-Negative'])]
    usable_reads_count = len(usable_reads_df)
    
    if usable_reads_count == 0:
        print("WARNING: No usable reads found. Cannot perform statistical analysis.")
        return {
            "total_reads": total_reads, 
            "usable_reads_count": 0, 
            "classification_counts": df['status'].value_counts().to_dict(),
            "drug_stats": {}
        }

    # --- THIS IS THE NEW, COMBINATION-AWARE ANALYSIS LOGIC ---
    
    # 3. Analyze Drug Resistance Prevalence
    drug_counts = {}
    
    # We iterate over *all* usable reads
    for _, row in usable_reads_df.iterrows():
        mutations_by_gene = row['mutations_by_gene']
        
        # This set ensures we count a drug *only once* for this read,
        # even if it has 5 mutations that all cause NVP resistance.
        read_drug_set = set()
        
        # 'mutations_by_gene' looks like: {"RT": ["K103N", "?103N", "M184V", "?184V"]}
        for gene, mutation_list in mutations_by_gene.items():
            if gene not in drm_db:
                continue
                
            # 1. Get the set of unique mutations found on this read for this gene
            #    e.g., {"K103N", "?103N", "M184V", "?184V"}
            mutation_set = set(mutation_list)
            
            # 2. Generate all possible subsets (powerset) to check against the DB
            #    e.g., {"K103N", "?103N", "M184V", ..., "K103N,M184V", ...}
            subsets_to_check = get_mutation_subsets(mutation_set)
            
            # 3. Check each subset key against the database
            for subset_key in subsets_to_check:
                if subset_key in drm_db[gene]:
                    # MATCH! This read has a single or combo resistance
                    for drug in drm_db[gene][subset_key]['drugs']:
                        read_drug_set.add(drug['name'])
        
        # 4. Now that we've processed all mutations on this read,
        #    increment the master counts for the unique drugs found.
        for drug_name in read_drug_set:
            drug_counts[drug_name] = drug_counts.get(drug_name, 0) + 1
    # --- END OF NEW LOGIC ---

    drug_prevalence = {drug: (count / usable_reads_count) for drug, count in drug_counts.items()}
    
    # 4. Perform Hypothesis Testing (This logic is statistically sound and unchanged)
    drug_stats = {}
    for drug, p_obs in drug_prevalence.items():
        count = drug_counts[drug]
        p_hyp = threshold
        
        # Add a check for 0 or negative variance
        variance = (p_hyp * (1 - p_hyp)) / usable_reads_count
        if variance <= 0:
            z_score = 0
            p_value = 1.0
        else:
            z_score = (p_obs - p_hyp) / np.sqrt(variance)
            p_value = norm.sf(z_score) # sf() is correct for one-tailed (p > p_hyp)
        
        drug_stats[drug] = {
            "count": count,
            "prevalence": p_obs,
            "p_value": p_value,
            "is_significant": p_value < 0.05 and p_obs > p_hyp
        }

    # 5. Collate all results
    analysis_results = {
        "total_reads": total_reads,
        "usable_reads_count": usable_reads_count,
        "classification_counts": df['status'].value_counts().to_dict(),
        "drug_stats": drug_stats
    }
    
    print("INFO: Statistical analysis complete.")
    return analysis_results


class PDF(FPDF):
    def __init__(self, *args, report_time="N/A", **kwargs):
        """
        Custom PDF class to hold the report generation time.
        """
        super().__init__(*args, **kwargs)
        self.report_time = report_time

    def header(self):
        self.set_font("Helvetica", "B", 12)
        self.cell(0, 10, "HIV DRM Analysis Report", align="C", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 8)
        # Use the passed-in, human-readable report time
        self.cell(0, 6, f"Report Date: {self.report_time}", align="C", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.ln(10)

    def footer(self):
        self.set_y(-15)
        self.set_font("Helvetica", "I", 8)
        self.cell(0, 10, f"Page {self.page_no()}", align="C")

    def chapter_title(self, title):
        self.set_font("Helvetica", "B", 16)
        self.cell(0, 10, title, new_x=XPos.LMARGIN, new_y=YPos.NEXT, border="B")
        self.ln(5)

    def chapter_body(self, text):
        self.set_font("Helvetica", "", 11)
        self.multi_cell(0, 6, text)
        self.ln()


    def add_methodology_page(self):
        """
        Adds the introductory and methodology pages to the PDF report.
        """
        self.add_page()
        
        # Title
        self.set_font("Helvetica", "B", 20)
        self.cell(0, 10, "HIV Drug Resistance Analysis Report", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
        self.ln(10)

        # --- Section 1: What is this Report? ---
        self.set_font("Helvetica", "B", 16)
        self.cell(0, 10, "1. What is this Report?", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 11)
        self.multi_cell(0, 6, 
            "This report analyzes sequencing data from an HIV sample to identify genetic mutations "
            "known to cause resistance to antiretroviral drugs. It provides a summary of all reads processed, "
            "a breakdown of the specific mutations found, and a statistical analysis of which drugs are likely "
            "to be ineffective for this sample."
        )
        self.ln(5)

        # --- Section 2: Why do we perform this analysis? ---
        self.set_font("Helvetica", "B", 16)
        self.cell(0, 10, "2. Why We Perform This Analysis", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 11)
        self.multi_cell(0, 6, 
            "When sequencing a viral population, we must distinguish between true, low-frequency mutations "
            "and simple, random sequencing errors. A mutation must be present at a high enough prevalence "
            "to be considered 'real' and clinically relevant.\n\n"
            "For this analysis, we use a standard clinical research threshold of 1.0%. Any mutation or drug "
            "resistance profile found in more than 1% of the viral population is considered a significant "
            "finding. This report tests our observed results against this 1.0% threshold."
        )
        self.ln(5)

        # --- Section 3: How do we test for significance? (The Statistics) ---
        self.set_font("Helvetica", "B", 16)
        self.cell(0, 10, "3. How We Test for Significance: The One-Proportion Z-Test", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 11)
        self.multi_cell(0, 6, 
            "To determine if an observed resistance prevalence (e.g., 1.5%) is 'statistically significant' or "
            "just a random chance event, we use a standard statistical tool called the one-proportion z-test. "
            "This test helps us decide between two competing claims:"
        )
        self.ln(4)

        # The Hypotheses (Corrected: H0, Ha)
        self.set_font("Helvetica", "B", 12)
        self.cell(0, 8, "The Null Hypothesis (H0): The 'Default' Assumption", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 11)
        self.multi_cell(0, 6, 
            "The null hypothesis states that our observation is not special. The true prevalence of the drug "
            "resistance in the entire viral population (p) is less than or equal to the 1.0% threshold."
        )
        self.set_font("Helvetica", "I", 11)
        
        # --- FIX 1: H0 Formula ---
        self.cell(0, 8, "H0: p <= 0.01", align="C", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.ln(4)

        self.set_font("Helvetica", "B", 12)
        self.cell(0, 8, "The Alternative Hypothesis (Ha): What We Want to Prove", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 11)
        self.multi_cell(0, 6, 
            "The alternative hypothesis states that our observation is significant. The true prevalence of "
            "the resistance (p) is statistically greater than the 1.0% threshold."
        )
        self.set_font("Helvetica", "I", 11)
        
        # --- FIX 2: Ha Formula ---
        self.cell(0, 8, "Ha: p > 0.01", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
        self.ln(4)

        # The Formula
        self.set_font("Helvetica", "B", 12)
        self.cell(0, 8, "The Z-Score Formula", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 11)
        self.multi_cell(0, 6, 
            "We calculate a 'Z-score' that measures how many standard deviations our observed result is from "
            "the 1.0% threshold. A large Z-score means our result is 'surprising'."
        )
        self.set_font("Helvetica", "I", 11)
        
        # --- FIX 3: Z-Score Formula ---
        self.cell(0, 8, "Z = (p_obs - p_hyp) / sqrt( (p_hyp * (1 - p_hyp)) / n )")
        
        self.set_font("Helvetica", "", 10)
        self.cell(0, 6, "     - p_obs = Observed prevalence (e.g., 20 resistant reads / 1000 total reads = 0.02)", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.cell(0, 6, "     - p_hyp = The hypothesized prevalence (0.01)", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.cell(0, 6, "     - n = The total number of usable reads", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.ln(4)

        # The Conclusion
        self.set_font("Helvetica", "B", 12)
        self.cell(0, 8, "The p-value and Our Decision", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 11)
        self.multi_cell(0, 6, 
            "The Z-score is converted into a 'p-value'. The p-value is the probability of seeing our result (or "
            "one more extreme) *if the null hypothesis were true*.\n"
            "We use a standard significance level (alpha) of 0.05:\n"
            "   - If p-value < 0.05 AND p_obs > p_hyp: The result is too 'weird' to be random. We reject the null hypothesis and "
            "conclude the drug resistance is STATISTICALLY SIGNIFICANT.\n"
            "   - Otherwise: The result is not surprising. We fail to reject the null hypothesis and "
            "conclude the finding is NOT statistically significant."
        )


def create_pdf_report(results: dict, output_path: str, plot_paths: dict, report_time: str, duration_str: str):
    """
    Generates a multi-page PDF report from the analysis results.
    """
    print(f"INFO: Generating PDF report at {output_path}...")
    
    pdf = PDF(report_time=report_time)
    pdf.set_auto_page_break(auto=True, margin=15)
    
    # --- ADD THE NEW METHODOLOGY PAGE FIRST ---
    pdf.add_methodology_page()

    # --- ADD THE RESULTS PAGE ---
    pdf.add_page()
    
    # --- Section 1: Run Summary ---
    pdf.chapter_title("Section 1: Run Summary & Quality Metrics")
    pdf.set_font("Helvetica", "", 12)
    pdf.cell(0, 8, f"Total Reads Processed: {results['total_reads']}", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.cell(0, 8, f"Total Usable Reads (High-Quality, Mapped): {results['usable_reads_count']}", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    
    # --- THIS IS THE NEW LINE ---
    pdf.cell(0, 8, f"Total Analysis Time: {duration_str}", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    # --- END OF UPDATE ---
    
    pdf.ln(5)
    
    # --- Section 2: Classification Breakdown ---
    pdf.chapter_title("Section 2: Overall Read Classification")
    pdf.image(plot_paths['pie_chart'], w=180)
    pdf.ln(5)

    # --- Section 3: Drug Resistance Prevalence ---
    if results['drug_stats']:
        pdf.add_page()
        pdf.chapter_title("Section 3: Drug Resistance Prevalence")
        pdf.set_font("Helvetica", "", 11)
        pdf.multi_cell(0, 6, "This plot shows the percentage of usable reads containing mutations that confer resistance to a specific drug. Any bar that crosses the red 1.0% threshold is considered a statistically significant finding.", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        pdf.ln(5)
        pdf.image(plot_paths['drug_chart'], w=180)
        pdf.ln(5)

        pdf.add_page()

        # --- Section 4: Detailed Statistical Table ---
        pdf.chapter_title("Section 4: Detailed Resistance Statistics")
        pdf.set_font("Helvetica", "B", 10)
        
        # Table Header
        pdf.cell(60, 8, "Drug", border=1)
        pdf.cell(30, 8, "Resist. Count", border=1, align="C")
        pdf.cell(35, 8, "Prevalence (%)", border=1, align="C")
        pdf.cell(35, 8, "p-value", border=1, align="C")
        pdf.cell(30, 8, "Significant?", border=1, align="C")
        pdf.set_font("Helvetica", "", 10)
        
        # Table Body
        for drug, stats in sorted(results['drug_stats'].items(), key=lambda item: item[1]['prevalence'], reverse=True):
            pdf.ln()
            pdf.cell(60, 8, drug, border=1)
            pdf.cell(30, 8, str(stats['count']), border=1, align="C")
            pdf.cell(35, 8, f"{stats['prevalence'] * 100:.2f}%", border=1, align="C")
            pdf.cell(35, 8, f"{stats['p_value']:.2e}", border=1, align="C")
            pdf.cell(30, 8, "Yes" if stats['is_significant'] else "No", border=1, align="C")
    
    pdf.output(output_path)
    print(f"INFO: PDF report generation complete at '{output_path}'.")

def generate_plots(results: dict, base_name: str) -> dict:
    """
    Generates all plots and saves them to disk.
    Returns a dictionary of plot file paths.
    """
    print("INFO: Generating plots...")
    plot_paths = {}

    # Plot 1: Pie Chart
    plt.figure(figsize=(10, 7))
    counts = results['classification_counts']
    filtered_counts = {k: v for k, v in counts.items() if v > 0}
    labels = filtered_counts.keys()
    sizes = filtered_counts.values()
    if not sizes:
        print("WARNING: No data to plot for pie chart.")
    else:
        plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, pctdistance=0.85)
        plt.title("Overall Read Classification", fontsize=16)
        plt.axis('equal')
    pie_path = f"{base_name}_plot_pie.png"
    plt.savefig(pie_path)
    plot_paths['pie_chart'] = pie_path

    # Plot 2: Drug Prevalence Bar Chart
    plt.figure(figsize=(10, 10))
    drug_stats = results['drug_stats']
    if drug_stats:
        sorted_drugs = sorted(drug_stats.items(), key=lambda item: item[1]['prevalence'])
        drug_names = [d[0] for d in sorted_drugs]
        drug_prevalence_pct = [d[1]['prevalence'] * 100 for d in sorted_drugs]
        
        plt.barh(drug_names, drug_prevalence_pct)
        plt.axvline(x=1.0, color='r', linestyle='--', label='1% Significance Threshold')
        plt.xlabel("Prevalence of Resistant Reads (%)")
        plt.title("Drug Resistance Prevalence", fontsize=16)
        plt.legend()
        plt.tight_layout()
    else:
        print("INFO: No drug resistance data to plot.")
        plt.text(0.5, 0.5, "No Drug Resistance Detected", horizontalalignment='center', verticalalignment='center', fontsize=14)
    
    drug_chart_path = f"{base_name}_plot_drugs.png"
    plt.savefig(drug_chart_path)
    
    plot_paths['drug_chart'] = drug_chart_path

    print("INFO: Plots generated and saved.")
    return plot_paths

# --- PART 4: MAIN EXECUTION ---

def main():
    parser = argparse.ArgumentParser(description="Final HIV DRM Analysis Pipeline (Build -> Run -> Report).")
    parser.add_argument("-i", "--input_bam", required=True)
    parser.add_argument("--ref_fasta", required=True)
    parser.add_argument("--profile_csv", required=True, help="Path to synthetic_profiles.csv")
    parser.add_argument("-n", "--num_reads", type=int, default=0, help="Number of reads to process. Set to 0 to process all reads.")
    args = parser.parse_args()

    # --- TIMESTAMP & START TIME ---
    start_time = datetime.datetime.now()
    timestamp_file = start_time.strftime("%Y%m%d_%H%M%S")
    timestamp_report = start_time.strftime("%Y-%m-%d %H:%M")
    # --- END OF UPDATE ---

    # --- 1. BUILD THE DATABASE ---
    drm_db = build_definitive_database(args.profile_csv)
    
    # --- 2. RUN THE SIMULATION ---
    print("\n--- Running Final Simulation ---")
    ref_sequences = load_reference_fasta(args.ref_fasta)
    print("-" * 50)
    
    all_read_results = []
    
    try:
        with pysam.AlignmentFile(args.input_bam, "rb") as bam_file_for_count:
            total_reads_in_bam = bam_file_for_count.count(until_eof=True)
    except Exception as e:
        print(f"FATAL: Could not read BAM file: {e}", file=sys.stderr)
        sys.exit(1)

    if args.num_reads == 0:
        num_to_process = total_reads_in_bam
    else:
        num_to_process = min(args.num_reads, total_reads_in_bam)

    print(f"INFO: Processing {num_to_process} reads from {args.input_bam}...")
    
    read_iterator = stream_bam_as_reads(args.input_bam)
    for i, read in enumerate(read_iterator):
        if i >= num_to_process:
            print(f"\nINFO: Reached read limit ({num_to_process}).")
            break
        classification_result = process_read_for_mutations(read, ref_sequences)
        all_read_results.append(classification_result)
        
        if (i + 1) % 5000 == 0:
            print(f"  ...processed {i+1} reads...")
    
    print("INFO: Simulation run complete.")
    
    # --- 3. SAVE THE JSON RESULTS (WITH TIMESTAMP) ---
    base_name = os.path.basename(args.input_bam).replace(".sorted.bam", "")
    output_dir = "data"
    os.makedirs(output_dir, exist_ok=True)
    
    file_base_name = f"{base_name}_{timestamp_file}"
    json_output_path = f"{output_dir}/{file_base_name}_results.json"
    plot_base_path = f"{output_dir}/{file_base_name}"
    pdf_output_path = f"{output_dir}/{file_base_name}_report.pdf"
    
    print(f"INFO: Saving raw simulation results to {json_output_path}...")
    with open(json_output_path, 'w') as f:
        json.dump(all_read_results, f, indent=2)
    
    # --- 4. GENERATE THE FINAL REPORT (WITH TIMESTAMP & DURATION) ---
    analysis_results = perform_statistical_analysis(json_output_path, drm_db)
    
    if analysis_results and analysis_results["usable_reads_count"] > 0:
        
        # --- CALCULATE DURATION ---
        end_time = datetime.datetime.now()
        duration = end_time - start_time
        total_seconds = duration.total_seconds()
        minutes = int(total_seconds // 60)
        seconds = int(total_seconds % 60)
        duration_str = f"{minutes} minute(s), {seconds} second(s)"
        # --- END OF UPDATE ---

        plot_paths = generate_plots(analysis_results, plot_base_path)
        
        # Pass the human-readable timestamp and duration string
        create_pdf_report(analysis_results, pdf_output_path, plot_paths, timestamp_report, duration_str)
        
        for path in plot_paths.values():
            if os.path.exists(path):
                os.remove(path)
                
        print("\n--- Pipeline Complete. Final report generated. ---")
        print(f"  - Results JSON: {json_output_path}")
        print(f"  - Final Report: {pdf_output_path}")
        print(f"  - Total Run Time: {duration_str}") # Added console output
    else:
        print("\n--- Pipeline Complete. No usable data to report on. ---")

if __name__ == "__main__":
    main()