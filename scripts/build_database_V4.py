import pandas as pd
import json
import sys
import re

INPUT_CSV_PATH = "data/synthetic_profiles.csv"
OUTPUT_PATH = "data/drm_database_V4.json"

def build_final_database(csv_path: str) -> dict:
    print(f"INFO: Reading source data from '{csv_path}'...")
    df = pd.read_csv(csv_path, dtype=str)
    df_long = pd.melt(df, id_vars=['mutations'], var_name='Drug', value_name='Score')
    df_long = df_long[df_long['Score'] != 'S'].dropna(subset=['Score'])

    print("INFO: Parsing mutations with final, robust logic...")
    # New Structure: { GENE: { POSITION_INT: { ALT_AMINO_ACID: { drugs: [...] } } } }
    drm_db = {"PR": {}, "RT": {}, "IN": {}}
    mutation_pattern = re.compile(r'([A-Z]?)(\d+)([A-Z])', re.IGNORECASE)

    for _, row in df_long.iterrows():
        cleaned_str = re.sub(r"['\"\[\]\s]", "", str(row['mutations']))
        individual_mutations = cleaned_str.split(',')
        for mut in individual_mutations:
            if not mut: continue
            match = mutation_pattern.match(mut)
            if match:
                _, position_str, alt_aa = match.groups()
                position = int(position_str)
                alt_aa = alt_aa.upper()
                
                gene = None
                if 1 <= position <= 99: gene = "PR"
                elif 100 <= position <= 560: gene = "RT"
                else: continue

                drug_info = {'name': row['Drug'], 'score': row['Score']}
                
                # Create nested dictionaries if they don't exist
                if position not in drm_db[gene]:
                    drm_db[gene][position] = {}
                if alt_aa not in drm_db[gene][position]:
                    drm_db[gene][position][alt_aa] = {'drugs': []}
                
                drm_db[gene][position][alt_aa]['drugs'].append(drug_info)

    print(f"INFO: Parsing complete.")
    return drm_db

def save_database(db: dict, path: str):
    print(f"INFO: Saving final V4 database to {path}...")
    with open(path, 'w') as f:
        json.dump(db, f, indent=2)
    print("INFO: Database saved successfully.")

def main():
    print("--- Building Final, Correct HIV DRM Database (V4) ---")
    drm_database = build_final_database(INPUT_CSV_PATH)
    save_database(drm_database, OUTPUT_PATH)

if __name__ == "__main__":
    main()