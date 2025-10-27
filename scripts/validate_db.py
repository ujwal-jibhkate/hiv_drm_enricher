import json
import sys

def validate_database_structure(db_path: str):
    """
    Performs a series of checks on the specified JSON database file to
    ensure its structure is correct for the simulation pipeline.
    """
    print(f"\n--- Validating DRM Database: '{db_path}' ---")
    
    try:
        with open(db_path, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"❌ VALIDATION FAILED: The file was not found at '{db_path}'.")
        return
    except json.JSONDecodeError as e:
        print(f"❌ VALIDATION FAILED: The file is not a valid JSON. Error: {e}")
        return

    # --- Structure Checks ---
    checks_passed = True

    # Check 1: Top level is a dictionary
    if isinstance(data, dict):
        print("✅ Check 1: Top-level structure is a dictionary.")
    else:
        print("❌ Check 1: Top-level structure is NOT a dictionary.")
        checks_passed = False

    # Check 2: Expected gene keys are present
    expected_genes = ["PR", "RT", "IN"]
    if all(gene in data for gene in expected_genes):
        print("✅ Check 2: Contains the expected gene keys (PR, RT, IN).")
    else:
        print(f"❌ Check 2: Missing one or more gene keys. Expected {expected_genes}, but found {list(data.keys())}.")
        checks_passed = False

    # Check 3: Check a sample gene value
    if checks_passed and isinstance(data["RT"], dict):
        print("✅ Check 3: The value for 'RT' is a dictionary (correct).")
    else:
        print("❌ Check 3: The value for 'RT' is NOT a dictionary.")
        checks_passed = False

    # Check 4 & 5: Check a sample mutation entry
    try:
        if checks_passed:
            # Get the first mutation in the RT gene to use as a sample
            sample_mutation_key = next(iter(data["RT"]))
            sample_mutation_value = data["RT"][sample_mutation_key]
            
            if isinstance(sample_mutation_value, dict):
                print(f"✅ Check 4: Sample mutation entry ('{sample_mutation_key}') is a dictionary.")
            else:
                print(f"❌ Check 4: Sample mutation entry ('{sample_mutation_key}') is NOT a dictionary.")
                checks_passed = False

            if "drugs" in sample_mutation_value:
                print("✅ Check 5: Sample mutation contains the 'drugs' key.")
            else:
                print("❌ Check 5: Sample mutation is missing the 'drugs' key.")
                checks_passed = False
            
            # Check 6, 7, 8: Check the 'drugs' list
            drugs_list = sample_mutation_value.get("drugs")
            if isinstance(drugs_list, list):
                print("✅ Check 6: The 'drugs' value is a list.")
            else:
                print("❌ Check 6: The 'drugs' value is NOT a list.")
                checks_passed = False

            if checks_passed and drugs_list: # Check if the list is not empty
                 print("✅ Check 7: The 'drugs' list is not empty.")
                 first_drug = drugs_list[0]
                 if isinstance(first_drug, dict) and "name" in first_drug and "score" in first_drug:
                     print("✅ Check 8: Drug entry has the correct structure (contains 'name' and 'score').")
                 else:
                     print("❌ Check 8: Drug entry has an incorrect structure.")
                     checks_passed = False
            elif checks_passed:
                 print("⚠️  Check 7: The 'drugs' list is empty. This might be valid but is worth noting.")

    except (StopIteration, KeyError, TypeError):
        # This can happen if the 'RT' section is empty
        print("⚠️  Could not perform deep checks on a mutation entry because the 'RT' section appears to be empty.")
        # We don't fail the check, as an empty section could be valid.

    # --- Final Summary ---
    print("\n--- Validation Summary ---")
    if checks_passed:
        print("🎉 SUCCESS: The database file appears to be correctly structured for the pipeline.")
    else:
        print("🔥 FAILURE: The database file has one or more structure errors. Please review the checks above.")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        file_path = sys.argv[1]
        validate_database_structure(file_path)
    else:
        print("Usage: python validate_db.py <path_to_your_json_file>")
        print("Example: python scripts/validate_db.py data/resistance_data.json")