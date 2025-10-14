import pandas as pd
import json
import sys
import re

# --- Configuration ---
INPUT_CSV_PATH = "data/synthetic_df.csv"
OUTPUT_PATH = "data/drm_database.json"

# Set a row limit for testing. Set to None to process the entire file.
TEST_MODE_ROW_LIMIT = None 

def parse_and_transform_csv(csv_path: str, row_limit: int = None) -> dict:
    """
    Reads the CSV, transforms it, and uses a robust parsing
    strategy to build the final nested dictionary.
    """
    print(f"INFO: Reading DRM data from '{csv_path}'...")
    if row_limit:
        print(f"--- RUNNING IN TEST MODE: Processing only the first {row_limit} rows. ---")
    try:
        # The 'nrows' parameter limits how many lines of the CSV are read.
        df = pd.read_csv(csv_path, nrows=row_limit)
    except FileNotFoundError:
        print(f"FATAL: File not found at '{csv_path}'.", file=sys.stderr)
        raise SystemExit(1)

    print("INFO: Transforming data from wide to long format...")
    id_vars = ['mutations']
    value_vars = df.columns.drop(id_vars)
    df_long = pd.melt(df, id_vars=id_vars, value_vars=value_vars, var_name='Drug', value_name='Score')
    df_long = df_long[df_long['Score'] != 'S'].dropna(subset=['Score'])

    print("INFO: Parsing complex mutation strings with robust logic...")
    drm_records = []
    for _, row in df_long.iterrows():
        mutation_groups = re.findall(r"'([^']*)'", row['mutations'])
        for group in mutation_groups:
            individual_mutations = re.split(r'\s+AND\s+', group)
            for mut in individual_mutations:
                mut = mut.strip()
                if not mut: continue
                try:
                    match = re.match(r'([A-Z]?)(\d+)([A-Z])$', mut, re.IGNORECASE)
                    if not match: continue
                    position = int(match.group(2))
                    alt_aa = match.group(3).upper()
                    drm_records.append({
                        'Position': position,
                        'AltAA': alt_aa,
                        'Drug': row['Drug'],
                        'Score': row['Score']
                    })
                except (ValueError, IndexError):
                    continue
    
    if not drm_records:
        print("FATAL: No valid DRM records could be parsed.", file=sys.stderr)
        raise SystemExit(1)
        
    parsed_df = pd.DataFrame(drm_records)
    
    print("INFO: Building the final nested dictionary...")
    drm_db = {}
    for _, row in parsed_df.iterrows():
        position_0_based = int(row['Position']) - 1
        alt_aa = row['AltAA']
        drug_info = {'name': row['Drug'], 'score': row['Score']}
        if position_0_based not in drm_db:
            drm_db[position_0_based] = {}
        if alt_aa not in drm_db[position_0_based]:
            drm_db[position_0_based][alt_aa] = {'aa_change': f"{row['Position']}{row['AltAA']}",'drugs': []}
        drm_db[position_0_based][alt_aa]['drugs'].append(drug_info)
        
    print(f"INFO: Parsing complete. Found DRM information for {len(drm_db)} unique positions.")
    return drm_db

def validate_database(db: dict):
    """
    Performs basic checks on the generated database to ensure it's valid.
    """
    print("\n--- Running Validation ---")
    # 1. Check if the database is not empty
    if not db:
        print("❌ VALIDATION FAILED: The generated database is empty.")
        return False
    print("✅ Validation 1: Database is not empty.")

    # 2. Check if there are a reasonable number of positions
    if len(db) < 10:
        print(f"⚠️  VALIDATION WARNING: Only {len(db)} unique positions found. This seems low.")
    else:
        print(f"✅ Validation 2: Found {len(db)} unique positions, which is a reasonable number.")

    # 3. Print a sample record for visual inspection
    try:
        # Get the first position and the first mutation within it
        first_pos = next(iter(db))
        first_alt = next(iter(db[first_pos]))
        sample_record = db[first_pos][first_alt]
        
        print("\n--- Sample Record for Visual Inspection ---")
        print(f"Position (0-based): {first_pos}")
        print(f"Alternate AA: {first_alt}")
        print("Data:")
        # Pretty print the sample record
        print(json.dumps(sample_record, indent=2))
        print("-----------------------------------------")
        print("✅ Validation 3: Sample record structure looks correct.")
        
    except (StopIteration, KeyError):
        print("❌ VALIDATION FAILED: Could not retrieve a sample record.")
        return False
        
    print("\n--- VALIDATION SUCCESSFUL ---")
    return True

def save_database(db: dict, path: str):
    """Saves the processed database dictionary to a JSON file."""
    print(f"\nINFO: Saving processed database to {path}...")
    with open(path, 'w') as f:
        json.dump(db, f, indent=2)
    print("INFO: Database saved successfully.")

def main():
    """Main function to build the DRM database from the local CSV file."""
    print("--- Building HIV DRM Database from Local CSV (Test & Validate) ---")
    drm_database = parse_and_transform_csv(INPUT_CSV_PATH, row_limit=TEST_MODE_ROW_LIMIT)
    
    if validate_database(drm_database):
        save_database(drm_database, OUTPUT_PATH)
        print("\n--- Database Build Complete ---")
    else:
        print("\n--- Database Build FAILED due to validation errors. ---", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()