<p align="center">
# HIV DRM Enricher 🧬
</p>

<p align="center">
<strong>A high-precision bioinformatics pipeline for detecting HIV drug resistance mutations from Nanopore sequencing data, featuring a hybrid alignment engine and rigorous clinical statistical validation.</strong>
</p>

<p align="center">
<img src="[https://img.shields.io/badge/Status-Research_Beta-yellow?style=flat-square](https://www.google.com/search?q=https://img.shields.io/badge/Status-Research_Beta-yellow%3Fstyle%3Dflat-square)" alt="Status">
<img src="[https://img.shields.io/badge/Python-3.10+-blue?style=flat-square&logo=python](https://www.google.com/search?q=https://img.shields.io/badge/Python-3.10%2B-blue%3Fstyle%3Dflat-square%26logo%3Dpython)" alt="Python">
<img src="[https://img.shields.io/badge/Platform-Linux_%7C_WSL-orange?style=flat-square&logo=linux](https://www.google.com/search?q=https://img.shields.io/badge/Platform-Linux_%257C_WSL-orange%3Fstyle%3Dflat-square%26logo%3Dlinux)" alt="Platform">
<img src="[https://img.shields.io/badge/Bioinformatics-Minimap2_%7C_Samtools-green?style=flat-square](https://www.google.com/search?q=https://img.shields.io/badge/Bioinformatics-Minimap2_%257C_Samtools-green%3Fstyle%3Dflat-square)" alt="Bioinformatics">
<a href="/LICENSE">
<img src="[https://img.shields.io/badge/License-MIT-lightgrey?style=flat-square](https://www.google.com/search?q=https://img.shields.io/badge/License-MIT-lightgrey%3Fstyle%3Dflat-square)" alt="License">
</a>
</p>

---

## 📖 Overview

**HIV DRM Enricher** is a specialized analysis system designed to bridge the gap between Third-Generation Sequencing (TGS) and clinical diagnostics. Targeted for deployment on edge devices like the **NVIDIA Jetson AGX Orin**, this pipeline processes Oxford Nanopore data to identify Drug Resistance Mutations (DRMs) in the HIV-1 *pol* gene (Protease, Reverse Transcriptase, and Integrase).

Unlike general-purpose variant callers, this system implements a **hybrid analysis engine** capable of handling both raw `FASTQ` reads and pre-aligned `BAM` files. It applies advanced statistical safeguards—including **FDR correction** and **statistical power analysis**—to distinguish genuine resistance from sequencing noise, generating publication-ready clinical reports.

## 🚀 Key Features

* **Hybrid Input Processing:** auto-detects sequencing platforms (Nanopore/PacBio/Illumina) and manages alignment via `minimap2` and `samtools` or accepts pre-aligned BAMs.
* **Deep Mutation Analysis:** Detects both **SNPs** and complex **Indels** (Insertions/Deletions) across PR, RT, and IN genes, essential for identifying complex resistance patterns like the T69 insertion complex.
* **Clinical Statistical Rigor:**
* **Drug-Specific Thresholds:** Applies distinct resistance barriers (e.g., 1% for NNRTIs vs. 15% for PIs) based on clinical evidence.
* **Power Analysis:** Calculates the statistical power (`1 - β`) to validate negative results.
* **FDR Correction:** Uses Benjamini-Hochberg correction to control false discovery rates across multi-drug testing.


* **Automated Reporting:** Generates a PDF clinical report with antibiograms, coverage plots, and statistical methodology.

## 🛠️ Installation

The project includes a **self-bootstrapping** setup script. It automatically detects, creates, and manages the required Conda environment.

### Prerequisites

* **Linux** or **Windows Subsystem for Linux (WSL2)**
* **Conda** (Miniconda or Mambaforge)

### Quick Start

1. Clone the repository:
```bash
git clone https://github.com/ujwal-jibhkate/hiv_drm_enricher.git
cd hiv_drm_enricher

```


2. Make the script executable and run:
```bash
chmod +x run_pipeline.sh
./run_pipeline.sh

```


*The script will automatically build the `hiv_drm_enricher` environment defined in `environment.yml` on the first run.*

## 💻 Usage

### 1. Interactive Mode (Drag & Drop)

Simply run the script without arguments to enter the interactive wizard:

```bash
./run_pipeline.sh

```

The system will prompt you to drag and drop your input file (`.fastq.gz`, `.fastq`, or `.bam`) directly into the terminal window.

### 2. CLI Mode (Batch Processing)

For server integration or batch jobs, use command-line arguments:

```bash
./run_pipeline.sh -i <input_file> [-n <reads>] [-p <preset>]

```

| Flag | Description | Default |
| --- | --- | --- |
| `-i` | Path to input file (`.fastq`, `.gz`, `.bam`) | **Required** |
| `-n` | Number of reads to process (0 for all) | `1600` |
| `-p` | Minimap2 preset (`map-ont`, `sr`, `map-pb`) | `sr` |

**Example:**

```bash
./run_pipeline.sh -i ./data/sample_01.fastq.gz -n 5000 -p map-ont

```

## 📊 Output & Artifacts

All results are generated in the `results/` directory:

1. **Clinical Report (`.pdf`)**: A comprehensive document containing executive summaries, drug resistance antibiograms with error bars, and quality control metrics.
2. **Raw Results (`.json`)**: Machine-readable classification data for every processed read.
3. **Visualization Plots (`.png`)**:
* **Drug Prevalence:** Resistance levels vs. clinical thresholds.
* **Genome Coverage:** Sequencing depth across HXB2 reference.
* **Mutation Breakdown:** SNP and Indel counts per gene.



## 🔬 Scientific Methodology

### Gene-Centric Database

The pipeline utilizes a compiled JSON database derived from the Stanford HIV Drug Resistance Database. It maps mutations using an  lookup strategy for Protease (AA 1-99), Reverse Transcriptase (AA 100-560), and Integrase (AA 561-849).

### Statistical Framework

To ensure clinical validity, the pipeline employs a multi-tiered statistical approach:

* **Binomial Exact Test:** Used for rare variants (count < 30) where normal approximation is insufficient.
* **Wilson Score Intervals:** Calculates 95% Confidence Intervals for all prevalence estimates.
* **Benjamini-Hochberg:** Adjusts P-values to correct for multiple hypothesis testing errors when screening 20+ antiretroviral drugs.

## 👤 Author

**Ujwal Jibhkate**

* *Master of Science in Data Science, Indiana University*
* *Research Assistant, Point-of-Care Bioinformatics & Genomic Analysis, Guan lab, IU*

## 📄 License

This project is licensed under the **MIT License**.