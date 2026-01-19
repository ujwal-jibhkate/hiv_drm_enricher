#!/bin/bash
set -e

# --- 1. Dynamic Setup ---
# Finds the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
REPO_ROOT="$SCRIPT_DIR"

# Paths
ENV_YML="$REPO_ROOT/environment.yml"
REF_FASTA="$REPO_ROOT/data/public/HXB2_reference.fasta"
PROFILE_CSV="$REPO_ROOT/data/synthetic_profiles.csv"
PYTHON_SCRIPT="$REPO_ROOT/scripts/run_entire_pipeline.py"
ENV_NAME="hiv_drm_enricher"

# --- 2. Environment Auto-Loader ---
echo "🔍 Checking environment..."

# Detect Conda or Mamba
if command -v mamba &> /dev/null; then
    CONDA_EXE="mamba"
elif command -v conda &> /dev/null; then
    CONDA_EXE="conda"
else
    echo "❌ Error: Neither 'conda' nor 'mamba' was found."
    echo "   To run this pipeline, please install Miniconda or Mambaforge first."
    echo "   Download: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Initialize Conda for this script session
# This allows 'conda activate' to work inside the script
eval "$($CONDA_EXE shell.bash hook)"

# Check if environment exists; if not, create it
if $CONDA_EXE env list | grep -q "$ENV_NAME"; then
    echo "✅ Environment '$ENV_NAME' found."
else
    echo "⚙️  Environment '$ENV_NAME' not found. Creating it now..."
    echo "   (This may take a few minutes, but only needs to happen once)"
    $CONDA_EXE env create -f "$ENV_YML"
fi

# Activate the environment
echo "🔌 Activating environment: $ENV_NAME"
conda activate "$ENV_NAME"

# --- 3. Interactive Input (The User Experience) ---
# Defaults
PRESET="sr"
N_VALUE=1600

usage() {
    echo "Usage: $0 -i <input_file> [-n <n_value>] [-p <preset>]"
    echo "  -i  Input file (Accepts .fastq, .fastq.gz, .fq, or .bam)"
    exit 1
}

# Check if arguments are passed
if [ $# -eq 0 ]; then
    # Interactive Mode
    echo " "
    echo "👋 Welcome to the HIV DRM Enricher Pipeline!"
    echo "-------------------------------------------"
    read -p "📂 Drag and drop your input file here: " INPUT_FILE
    # Remove quotes from drag-and-drop
    INPUT_FILE=$(echo "$INPUT_FILE" | tr -d "'\"")
    
    read -p "⚙️  Enter N Value (Default 1600): " USER_N
    if [ ! -z "$USER_N" ]; then N_VALUE=$USER_N; fi

    read -p "⚙️  Enter Minimap Preset (Default 'sr'): " USER_P
    if [ ! -z "$USER_P" ]; then PRESET=$USER_P; fi
else
    # Argument Mode
    while getopts "i:n:p:h" opt; do
        case $opt in
            i) INPUT_FILE="$OPTARG" ;;
            n) N_VALUE="$OPTARG" ;;
            p) PRESET="$OPTARG" ;;
            h) usage ;;
            *) usage ;;
        esac
    done
fi

# --- 4. Validation & Execution ---
if [ ! -f "$INPUT_FILE" ]; then
    echo "❌ Error: File '$INPUT_FILE' not found."
    exit 1
fi

FILENAME=$(basename -- "$INPUT_FILE")
BASENAME="${FILENAME%%.*}"
OUTPUT_DIR="$REPO_ROOT/results/${BASENAME}_output"
mkdir -p "$OUTPUT_DIR"

echo " "
echo "🚀 Running analysis for: $BASENAME"
echo "📂 Output Directory: $OUTPUT_DIR"
echo "-----------------------------------"

EXT="${INPUT_FILE##*.}"

if [[ "$EXT" == "bam" ]]; then
    echo "✅ Input is BAM. Skipping alignment."
    SORTED_BAM="$INPUT_FILE"
else 
    SORTED_BAM="$OUTPUT_DIR/${BASENAME}.sorted.bam"
    if [ -f "$SORTED_BAM" ]; then
        echo "⚠️  Sorted BAM exists. Skipping alignment."
    else
        echo "⚡ FASTQ detected. Starting alignment..."
        TEMP_BAM="$OUTPUT_DIR/${BASENAME}.temp.bam"
        
        # Note: We don't need absolute paths to minimap2/samtools 
        # because the conda environment is active!
        minimap2 -ax $PRESET "$REF_FASTA" "$INPUT_FILE" | \
        samtools view -bS - | \
        samtools sort -o "$TEMP_BAM" -

        samtools calmd -b "$TEMP_BAM" "$REF_FASTA" > "$SORTED_BAM"
        samtools index "$SORTED_BAM"
        rm "$TEMP_BAM"
    fi
fi

echo "🐍 Running Python analysis..."
python "$PYTHON_SCRIPT" \
  --input_bam "$SORTED_BAM" \
  --ref_fasta "$REF_FASTA" \
  --profile_csv "$PROFILE_CSV" \
  -n $N_VALUE \
  -o "$OUTPUT_DIR"

echo " "
echo "✅ Done! Check results in: $OUTPUT_DIR"