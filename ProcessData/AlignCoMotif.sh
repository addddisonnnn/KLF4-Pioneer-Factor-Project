import pandas as pd
import subprocess
import os
import numpy as np # For handling NaN values

def run_command(command, error_message):
    """Helper function to run shell commands and check for errors."""
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error: {error_message}")
        print(f"Command: {e.cmd}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        raise

def align_fimo_results_to_csv(
    input_csv_file,
    input_delimiter,
    fimo_result_files,
    output_csv_file
):
    """
    Aligns multiple FIMO motif results to a base CSV file, adding boolean columns
    indicating motif overlap.

    Args:
        input_csv_file (str): Path to the base CSV file (e.g., Human_histone.csv).
        input_delimiter (str): Delimiter for the input CSV (e.g., ',', '\t').
        fimo_result_files (list): List of paths to FIMO TSV result files.
        output_csv_file (str): Path for the final output CSV file.
    """
    print("===================================================")
    print("    Aligning Multiple FIMO Results to CSV File     ")
    print("===================================================")

    # --- Configuration ---
    temp_dir = f"temp_fimo_align_{os.getpid()}" # Unique temp directory
    os.makedirs(temp_dir, exist_ok=True)

    print(f"Input CSV: '{input_csv_file}'")
    print(f"Input Delimiter: '{input_delimiter}'")
    print(f"Output CSV: '{output_csv_file}'")
    print(f"Temporary Directory: '{temp_dir}'")
    print("---------------------------------------------------")

    try:
        # --- Input Validation ---
        if not os.path.exists(input_csv_file):
            raise FileNotFoundError(f"Input CSV file '{input_csv_file}' not found.")

        # Check if bedtools is available
        run_command("command -v bedtools", "bedtools is not found. Please install it and ensure it's in your PATH.")

        for fimo_file in fimo_result_files:
            if not os.path.exists(fimo_file):
                raise FileNotFoundError(f"FIMO result file '{fimo_file}' not found.")
            # Check if FIMO file has data lines beyond header
            with open(fimo_file, 'r') as f:
                if len(f.readlines()) <= 1:
                    print(f"Warning: FIMO result file '{fimo_file}' contains only a header or is empty. No motif occurrences will be processed for this file.")

        # --- Step 1: Load base CSV and prepare as BED with unique IDs ---
        print("\nStep 1: Loading base CSV and preparing as BED with unique IDs...")
        try:
            # Read the input CSV, assuming first three columns are chrom, start, end
            # We don't need to specify column names here, pandas will infer or use header
            base_df = pd.read_csv(input_csv_file, sep=input_delimiter, dtype={'start': int, 'end': int})
        except Exception as e:
            raise ValueError(f"Error reading input CSV '{input_csv_file}' with delimiter '{input_delimiter}': {e}")

        # Add a unique ID column to the DataFrame
        # This ID will link back to the original rows after bedtools operations
        base_df['unique_id'] = base_df['chrom'].astype(str) + ':' + \
                               base_df['start'].astype(str) + '-' + \
                               base_df['end'].astype(str) + ':' + \
                               base_df.index.astype(str) # Use pandas index for uniqueness

        # Create a temporary BED file for bedtools intersection
        # Only chrom, start, end, and unique_id are needed
        temp_base_bed_path = os.path.join(temp_dir, "input_base_with_ids.bed")
        base_df[['chrom', 'start', 'end', 'unique_id']].to_csv(
            temp_base_bed_path, sep='\t', header=False, index=False
        )
        
        if os.path.getsize(temp_base_bed_path) == 0:
            raise ValueError(f"Failed to create temporary base BED file or it's empty. Check '{input_csv_file}' and its delimiter.")

        # --- Step 2: Process each FIMO result file and generate overlap IDs ---
        print("\nStep 2: Processing each FIMO result file for overlaps...")
        all_motif_overlap_ids = {} # Dictionary to store sets of overlapping IDs for each motif

        for fimo_file in fimo_result_files:
            column_name = os.path.basename(fimo_file).replace("fimo.tsv", "")
            print(f"  - Processing '{fimo_file}' (new column: '{column_name}')...")

            # Parse FIMO TSV robustly
            try:
                fimo_df = pd.read_csv(fimo_file, sep='\t')
            except Exception as e:
                print(f"Warning: Error reading FIMO TSV '{fimo_file}': {e}. Skipping this motif.")
                all_motif_overlap_ids[column_name] = set() # Add empty set
                continue

            # Validate required FIMO columns
            required_fimo_cols = ['sequence_name', 'start', 'stop']
            if not all(col in fimo_df.columns for col in required_fimo_cols):
                print(f"Warning: FIMO TSV '{fimo_file}' is missing one or more required columns: {required_fimo_cols}. Skipping this motif.")
                all_motif_overlap_ids[column_name] = set() # Add empty set
                continue

            # Filter out rows with NaN in critical columns (start, stop)
            original_fimo_rows = len(fimo_df)
            fimo_df.dropna(subset=['start', 'stop'], inplace=True)
            if len(fimo_df) < original_fimo_rows:
                print(f"    Warning: Removed {original_fimo_rows - len(fimo_df)} rows from FIMO file '{fimo_file}' due to missing start/stop coordinates.")

            # Convert FIMO 1-based coordinates to 0-based BED coordinates
            # FIMO 'start' is 1-based inclusive, 'stop' is 1-based inclusive.
            # BED 'start' is 0-based inclusive, 'end' is 0-based exclusive.
            fimo_df['bed_start'] = fimo_df['start'].astype(int) - 1
            fimo_df['bed_end'] = fimo_df['stop'].astype(int) # FIMO stop is inclusive, so it's the 0-based exclusive end

            # Ensure coordinates are valid (start < end) and non-negative
            fimo_df = fimo_df[fimo_df['bed_start'] >= 0]
            fimo_df = fimo_df[fimo_df['bed_start'] < fimo_df['bed_end']]
            if len(fimo_df) == 0:
                print(f"    Warning: No valid motif hits found in '{fimo_file}' after parsing and validation. No overlaps will be found for this motif.")
                all_motif_overlap_ids[column_name] = set() # Add empty set
                continue

            # Create a temporary BED file for the current FIMO results
            temp_fimo_bed_path = os.path.join(temp_dir, f"{column_name}_fimo.bed")
            fimo_df[['sequence_name', 'bed_start', 'bed_end']].to_csv(
                temp_fimo_bed_path, sep='\t', header=False, index=False
            )

            # Perform bedtools intersect
            # -a: base_df (Human_histone peaks with unique_id)
            # -b: current FIMO motif hits
            # -u: Report only unique entries in -a that overlap -b
            # cut -f4: Extract only the unique_id from the overlapping base peaks
            print(f"    Running bedtools intersect for {column_name}...")
            bedtools_cmd = f"bedtools intersect -a {temp_base_bed_path} -b {temp_fimo_bed_path} -u | cut -f4"
            overlapping_ids_str = run_command(bedtools_cmd, f"Error running bedtools for {column_name}")
            
            # Store overlapping IDs in a set for efficient lookup
            all_motif_overlap_ids[column_name] = set(overlapping_ids_str.strip().split('\n'))
            if '' in all_motif_overlap_ids[column_name]: # Remove potential empty string from split
                all_motif_overlap_ids[column_name].remove('')
            
            print(f"    Found {len(all_motif_overlap_ids[column_name])} overlaps for {column_name}.")

        # --- Step 3: Add new boolean columns to the original DataFrame ---
        print("\nStep 3: Adding new boolean columns to the original CSV...")
        for motif_name, overlapping_ids_set in all_motif_overlap_ids.items():
            # Create the boolean column based on whether the unique_id is in the set of overlaps
            base_df[motif_name] = base_df['unique_id'].apply(lambda x: 1 if x in overlapping_ids_set else 0)
            print(f"  - Added column '{motif_name}' with {base_df[motif_name].sum()} '1's.")

        # Drop the temporary unique_id column before saving
        final_df = base_df.drop(columns=['unique_id'])

        # --- Step 4: Save the final DataFrame to CSV ---
        print("\nStep 4: Saving the final DataFrame to CSV...")
        final_df.to_csv(output_csv_file, sep=input_delimiter, index=False)
        print(f"Analysis complete. Output saved to '{output_csv_file}'.")

    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")
        print("Script failed.")
    finally:
        # --- Cleanup Temporary Files ---
        print(f"\nCleaning up temporary directory '{temp_dir}'...")
        if os.path.exists(temp_dir):
            import shutil
            shutil.rmtree(temp_dir)
        print("Cleanup complete.")
        print("===================================================")


if __name__ == "__main__":
    # --- Define your input and output file paths here ---
    input_csv = "Human_histone.csv" # Your base CSV file
    input_sep = "," # Delimiter of your Human_histone.csv (e.g., ',' or '\t')

    fimo_files = [
        "FOSL2fimo.tsv",
        "FOXC2fimo.tsv",
        "OTX2fimo.tsv",
        "POU6F2fimo.tsv",
        "TP53fimo.tsv",
        "ZBED4fimo.tsv",
        "ZBTB6fimo.tsv",
        "ZNF213fimo.tsv",
        "ZNF384fimo.tsv",
        "ZNF460fimo.tsv"
    ]

    output_csv = "Human_histone_with_motifs.csv"

    align_fimo_results_to_csv(input_csv, input_sep, fimo_files, output_csv)
