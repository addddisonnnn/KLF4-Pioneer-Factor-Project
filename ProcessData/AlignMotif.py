# Script to align the BED contents with the KLF4 motif information in the form of a csv file
# Author: Addison Yam

import pandas as pd
from collections import defaultdict
import numpy as np # Import numpy for NaN checks

def parse_fasta(fasta_file_path):
    """
    Parses a FASTA file and returns a dictionary mapping sequence IDs to their sequences.
    Assumes standard FASTA format where header lines start with '>'
    and sequence lines follow.
    """
    sequences = {}
    current_seq_id = None
    current_sequence = []

    with open(fasta_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq_id:
                    sequences[current_seq_id] = "".join(current_sequence)
                current_seq_id = line[1:].split()[0] # Take only the first word as ID
                current_sequence = []
            else:
                current_sequence.append(line)
        if current_seq_id: # Add the last sequence
            sequences[current_seq_id] = "".join(current_sequence)
    return sequences

def parse_fimo_tsv(fimo_tsv_file_path):
    """
    Parses a FIMO TSV file and organizes motif hits by chromosome.
    Returns a dictionary where keys are chromosome names and values are lists of
    dictionaries, each representing a motif hit with 'start', 'stop', 'strand', 'score'.
    """
    fimo_motifs = defaultdict(list)

    # Read the FIMO TSV directly with pandas to handle headers and types
    try:
        fimo_df = pd.read_csv(fimo_tsv_file_path, sep='\t')
    except FileNotFoundError:
        print(f"Error: FIMO TSV file not found at {fimo_tsv_file_path}")
        return fimo_motifs
    except Exception as e:
        print(f"Error reading FIMO TSV: {e}")
        return fimo_motifs

    # Ensure required columns exist
    required_fimo_cols = ['sequence_name', 'start', 'stop', 'strand', 'score']
    if not all(col in fimo_df.columns for col in required_fimo_cols):
        print(f"Error: FIMO TSV is missing one or more required columns: {required_fimo_cols}")
        print(f"Available columns: {fimo_df.columns.tolist()}")
        return fimo_motifs

    # Iterate through rows and handle potential NaN values in 'start' and 'stop'
    for idx, row in fimo_df.iterrows():
        # FIX: Check for NaN values before attempting conversion to int
        if pd.isna(row['start']) or pd.isna(row['stop']):
            print(f"Warning: Skipping row {idx + 1} in FIMO TSV due to NaN in 'start' or 'stop' column: {row.to_dict()}")
            continue

        # Ensure 'score' is not NaN either, as it's used for max()
        if pd.isna(row['score']):
            print(f"Warning: Skipping row {idx + 1} in FIMO TSV due to NaN in 'score' column: {row.to_dict()}")
            continue

        chrom = row['sequence_name']
        # FIMO coordinates are 1-based, BED are 0-based. Convert FIMO to 0-based for consistency.
        # FIMO 'start' is 1-based inclusive, 'stop' is 1-based inclusive.
        # So, 0-based start is row['start'] - 1.
        # 0-based end is row['stop'] (since it's exclusive for 0-based, and FIMO stop is inclusive).
        fimo_start = int(row['start']) - 1
        fimo_stop = int(row['stop']) # FIMO stop is inclusive, so it's the 0-based exclusive end
        strand = str(row['strand'])
        score = float(row['score'])

        fimo_motifs[chrom].append({
            'start': fimo_start,
            'stop': fimo_stop,
            'strand': strand,
            'score': score
        })
    return fimo_motifs

def analyze_klf4_binding(bed_file_path, fasta_file_path, fimo_tsv_file_path, output_csv_path):
    """
    Analyzes KLF4 binding by combining information from BED, FASTA, and FIMO TSV files.
    """
    # --- 1. Load BED file ---
    # The BED file is assumed to have 10 columns as specified.
    # We'll assign the column names directly during loading.
    # UPDATED: Column names for the BED file
    bed_column_names = [
        "chrom", "start", "end", "BJ", "H9_ESC", "iPSC",
        "LIS49_hESC", "MCF7", "U87", "HN_SCC"
    ]
    try:
        # Read the BED file without a header, specifying column names
        bed_df = pd.read_csv(bed_file_path, sep='\t', header=None, names=bed_column_names)
        # Ensure start and end are integers
        bed_df['start'] = bed_df['start'].astype(int)
        bed_df['end'] = bed_df['end'].astype(int)
    except FileNotFoundError:
        print(f"Error: BED file not found at {bed_file_path}")
        return
    except Exception as e:
        print(f"Error reading BED file: {e}")
        return

    print(f"Loaded {len(bed_df)} peaks from BED file.")

    # --- 2. Parse FIMO TSV ---
    print("Parsing FIMO TSV file...")
    fimo_motifs_by_chrom = parse_fimo_tsv(fimo_tsv_file_path)
    if not fimo_motifs_by_chrom:
        print("No FIMO motifs loaded or an error occurred during FIMO parsing. Proceeding without motif data.")
        # Initialize new columns with default values if FIMO parsing failed
        bed_df['Klf4 motif'] = 0
        bed_df['Klf4 motif score'] = np.nan # Use np.nan for consistency with pandas
        bed_df['Number of Klf4 motifs'] = 0
        bed_df['Strand'] = np.nan # Use np.nan for consistency with pandas
        bed_df['Total motif score'] = np.nan
        bed_df.to_csv(output_csv_path, index=False)
        print(f"Output saved to {output_csv_path} (without motif details due to FIMO error/absence).")
        return

    print(f"Loaded FIMO motifs for {len(fimo_motifs_by_chrom)} chromosomes.")

    # --- 3. Parse FASTA file ---
    # Although the FASTA is not strictly used for motif finding in this script (FIMO is),
    # we'll parse it as requested, in case it's needed for future extensions or validation.
    print("Parsing FASTA file...")
    fasta_sequences = parse_fasta(fasta_file_path)
    if not fasta_sequences:
        print("Warning: No sequences loaded from FASTA file. FASTA-based motif check not possible if implemented.")
    else:
        print(f"Loaded {len(fasta_sequences)} sequences from FASTA file.")


    # --- 4. Add new columns and populate data ---
    klf4_motif_found = []
    klf4_motif_score = []
    num_klf4_motifs = []
    klf4_motif_strand = []
    total_motif_score = []

    print("Analyzing each BED peak for KLF4 motifs...")
    for index, row in bed_df.iterrows():
        chrom = row['chrom']
        bed_start = row['start']
        bed_end = row['end']

        overlapping_motifs = []
        if chrom in fimo_motifs_by_chrom:
            for motif in fimo_motifs_by_chrom[chrom]:
                fimo_start = motif['start']
                fimo_stop = motif['stop']

                # Check for overlap between BED peak and FIMO motif
                # Overlap exists if max(start1, start2) < min(end1, end2)
                if max(bed_start, fimo_start) < min(bed_end, fimo_stop):
                    overlapping_motifs.append(motif)

        if overlapping_motifs:
            klf4_motif_found.append(1)
            num_klf4_motifs.append(len(overlapping_motifs))

            # Find the motif with the highest score among overlapping ones
            best_motif = max(overlapping_motifs, key=lambda x: x['score'])
            klf4_motif_score.append(best_motif['score'])
            klf4_motif_strand.append(best_motif['strand'])

            total_score = sum(m['score'] for m in overlapping_motifs)
            total_motif_score.append(total_score)
        else:
            klf4_motif_found.append(0)
            klf4_motif_score.append(np.nan) # Use np.nan for missing score
            num_klf4_motifs.append(0)
            klf4_motif_strand.append(np.nan) # Use np.nan for consistency with pandas
            total_motif_score.append(np.nan)
    # Add the new columns to the DataFrame
    bed_df['Klf4 motif'] = klf4_motif_found
    bed_df['Klf4 motif score'] = klf4_motif_score
    bed_df['Number of Klf4 motifs'] = num_klf4_motifs
    bed_df['Strand'] = klf4_motif_strand
    bed_df['Total motif score'] = total_motif_score

    # --- 5. Save to CSV ---
    try:
        bed_df.to_csv(output_csv_path, index=False)
        print(f"Analysis complete. Output saved to {output_csv_path}")
    except Exception as e:
        print(f"Error saving output CSV: {e}")

# --- Main execution block ---
if __name__ == "__main__":
    # Define your input and output file paths here
    bed_file_path = "humanmaster.bed" # Path to your BED file (output from Step 3 Bash script)
    fasta_file_path = "humansequence.fa" # Path to your FASTA genome file (e.g., mm10.fasta)
    fimo_tsv_file_path = "humanfimo.tsv" # Path to your FIMO TSV output for KLF4 motif
    output_csv_path = "Human_motif.csv" # Desired output CSV file name

    analyze_klf4_binding(bed_file_path, fasta_file_path, fimo_tsv_file_path, output_csv_path)
