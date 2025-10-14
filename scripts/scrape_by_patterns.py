import pandas as pd
import json
import sys
import time

from selenium import webdriver
from selenium.webdriver.chrome.service import Service as ChromeService
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager

# --- Configuration ---
STANFORD_URL = "https://hivdb.stanford.edu/hivdb/by-patterns/"
OUTPUT_PATH = "data/drm_database.json"

# We will provide the major known mutations for each gene.
# The scraper will paste these into the form to get the resistance data.
MUTATIONS_TO_QUERY = {
    "PR": [
        "L10F", "V11I", "G16E", "K20R", "L24I", "D30N", "V32I", "L33F",
        "E35D", "M36I", "K43T", "M46I", "I47V", "G48V", "I50V", "F53L",
        "I54V", "Q58E", "D60E", "I62V", "L63P", "I64M", "A71V", "G73S",
        "V77I", "V82A", "I84V", "I85V", "N88D", "L90M"
    ],
    "RT": [
        "M41L", "A62V", "K65R", "D67N", "T69D", "K70R", "L74V", "V75I",
        "F77L", "Y115F", "F116Y", "Q151M", "M184V", "L210W", "T215Y", "K219Q",
        "K103N", "V106M", "E138K", "V179D", "Y181C", "Y188L", "G190A", "P225H", "M230L"
    ],
    "IN": [
        "T66I", "E92Q", "G118R", "Y143R", "S147G", "Q148H", "N155H", "R263K"
    ]
}

def setup_driver():
    """Initializes a Selenium Chrome WebDriver."""
    print("INFO: Setting up Selenium WebDriver...")
    # These options are crucial for running Chrome in a headless (no GUI)
    # environment like WSL.
    chrome_options = webdriver.ChromeOptions()
    chrome_options.add_argument("--headless")
    chrome_options.add_argument("--no-sandbox")
    chrome_options.add_argument("--disable-dev-shm-usage")
    
    # webdriver-manager will automatically download and manage the correct
    # driver for your installed version of Chrome.
    service = ChromeService(ChromeDriverManager().install())
    driver = webdriver.Chrome(service=service, options=chrome_options)
    print("INFO: WebDriver is ready.")
    return driver

def scrape_gene_data(driver, gene: str, mutations: list) -> pd.DataFrame:
    """
    Navigates the Stanford page, enters mutations, and scrapes the results.
    """
    print(f"\n--- Scraping data for gene: {gene} ---")
    try:
        driver.get(STANFORD_URL)

        # 1. Find the textarea and enter the mutations
        textarea = driver.find_element(By.ID, "patterns")
        textarea.clear()
        mutations_text = "\\n".join(mutations) # Use \n for newlines
        textarea.send_keys(mutations_text)
        print(f"INFO: Entered {len(mutations)} mutations into the form.")

        # 2. Find and click the "Analyze" button
        analyze_button = driver.find_element(By.XPATH, "//input[@value='Analyze']")
        analyze_button.click()
        print("INFO: Clicked 'Analyze'. Waiting for results...")

        # 3. Wait intelligently for the results table to appear
        # This is the most important part. We wait up to 60 seconds
        # for an element with the ID 'resultTable' to be visible.
        wait = WebDriverWait(driver, 60)
        wait.until(EC.visibility_of_element_located((By.ID, "resultTable")))
        print("INFO: Results table is now visible.")

        # 4. Scrape the table using pandas
        # We pass the page source of the now-updated page to pandas.
        page_source = driver.page_source
        tables = pd.read_html(StringIO(page_source), attrs={'id': 'resultTable'})
        
        if not tables:
            raise ValueError("Could not find the resultTable even after waiting.")
            
        df = tables[0]
        df['Gene'] = gene # Add gene info to each row
        print(f"INFO: Successfully scraped {len(df)} rows of data.")
        return df

    except Exception as e:
        print(f"FATAL: An error occurred during scraping for gene '{gene}'. Error: {e}", file=sys.stderr)
        # Take a screenshot for debugging if something goes wrong
        driver.save_screenshot("error_screenshot.png")
        print("INFO: A screenshot has been saved as error_screenshot.png")
        raise

def parse_scraped_data(df: pd.DataFrame) -> dict:
    """Transforms the scraped DataFrame into our final nested dictionary."""
    print("\nINFO: Transforming all scraped data into final JSON database...")
    drm_db = {}
    
    # Melt the DataFrame to convert drug columns into rows
    id_vars = ['Mutation', 'Gene']
    value_vars = df.columns.drop(id_vars)
    df_long = pd.melt(df, id_vars=id_vars, value_vars=value_vars, var_name='Drug', value_name='Score')
    df_long = df_long.dropna(subset=['Score']) # Drop rows with no score

    for _, row in df_long.iterrows():
        mutation = row['Mutation']
        try:
            position = int(''.join(filter(str.isdigit, mutation)))
            ref_aa = mutation[0]
            alt_aa = mutation[-1]
            position_0_based = position - 1

            drug_info = {'name': row['Drug'], 'score': str(row['Score'])}

            if position_0_based not in drm_db:
                drm_db[position_0_based] = {}
            if alt_aa not in drm_db[position_0_based]:
                drm_db[position_0_based][alt_aa] = {
                    'gene': row['Gene'],
                    'aa_change': mutation,
                    'drugs': []
                }
            drm_db[position_0_based][alt_aa]['drugs'].append(drug_info)
        except (ValueError, IndexError):
            continue

    print(f"INFO: Parsing complete. Found DRM info for {len(drm_db)} unique positions.")
    return drm_db

def save_database(db: dict, path: str):
    """Saves the processed database to a JSON file."""
    print(f"INFO: Saving processed database to {path}...")
    with open(path, 'w') as f:
        json.dump(db, f, indent=2)
    print("INFO: Database saved successfully.")

def main():
    driver = setup_driver()
    all_results = []
    try:
        for gene, mutations in MUTATIONS_TO_QUERY.items():
            gene_df = scrape_gene_data(driver, gene, mutations)
            all_results.append(gene_df)
    finally:
        # It's crucial to close the browser session, even if errors occur
        print("INFO: Closing WebDriver.")
        driver.quit()

    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        final_db = parse_scraped_data(combined_df)
        save_database(final_db, OUTPUT_PATH)
        print("\n--- Database Build Complete ---")
    else:
        print("ERROR: No data was scraped. The database was not built.", file=sys.stderr)

if __name__ == "__main__":
    from io import StringIO
    main()