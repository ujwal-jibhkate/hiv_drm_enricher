import pandas as pd
import json
import sys
import re

# --- Configuration ---
INPUT_CSV_PATH = "data/synthetic_profiles.csv"
OUTPUT_PATH = "data/drm_database_definitive.json"

def build_definitive_database(csv_path: str) -> dict:
    """
    Reads the synthetic_profiles.csv and uses a robust, multi-step cleaning
    function to correctly parse the 'mutations' column.
    """
    print(f"INFO: Reading source data from '{csv_path}'...")
    try:
        # Read the CSV, treating everything as a string to avoid parsing errors
        df = pd.read_csv(csv_path, dtype=str)
    except FileNotFoundError:
        print(f"FATAL: The file was not found at '{csv_path}'.", file=sys.stderr)
        raise SystemExit(1)

    print("INFO: Transforming data from wide to long format...")
    df_long = pd.melt(df, id_vars=['mutations'], var_name='Drug', value_name='Score')
    df_long = df_long[df_long['Score'] != 'S'].dropna(subset=['Score'])

    print("INFO: Parsing complex mutation strings with definitive multi-step cleaner...")
    drm_db = {"PR": {}, "RT": {}, "IN": {}}
    mutation_pattern = re.compile(r'([A-Z])(\d+)([A-Z])', re.IGNORECASE)

    for _, row in df_long.iterrows():
        # --- The Definitive Parsing Logic ---
        # 1. Start with the raw string: "['75T', '54A', '70T', '44A']"
        mutations_str = str(row['mutations'])
        
        # 2. Remove all non-mutation characters: quotes, brackets, spaces
        cleaned_str = re.sub(r"['\"\[\]\s]", "", mutations_str)
        # Result: "75T,54A,70T,44A"
        
        # 3. Split by the comma to get individual mutations
        individual_mutations = cleaned_str.split(',')

        for mut in individual_mutations:
            if not mut: continue

            # 4. Now, match the clean string against our simple pattern
            match = mutation_pattern.match(mut)
            if match:
                ref_aa, position_str, alt_aa = match.groups()
                position = int(position_str)
                aa_change_str = f"{ref_aa.upper()}{position}{alt_aa.upper()}"

                gene = None
                if 1 <= position <= 99: gene = "PR"
                elif 100 <= position <= 560: gene = "RT"
                
                if not gene: continue

                drug_info = {'name': row['Drug'], 'score': row['Score']}
                if aa_change_str not in drm_db[gene]:
                    drm_db[gene][aa_change_str] = {'drugs': []}
                drm_db[gene][aa_change_str]['drugs'].append(drug_info)

    print(f"INFO: Parsing complete. Found DRM info for {len(drm_db['PR']) + len(drm_db['RT']) + len(drm_db['IN'])} unique mutations.")
    return drm_db

def save_database(db: dict, path: str):
    """Saves the processed database to a JSON file."""
    print(f"INFO: Saving definitive database to {path}...")
    with open(path, 'w') as f:
        json.dump(db, f, indent=2)
    print("INFO: Database saved successfully.")

def main():
    """Main function to build the definitive DRM database."""
    print("--- Building Definitive HIV DRM Database ---")
    drm_database = build_definitive_database(INPUT_CSV_PATH)
    save_database(drm_database, OUTPUT_PATH)
    print("\n--- Definitive Database Build Complete ---")

if __name__ == "__main__":
    main()