#!/bin/bash

# Script to add an 'H3K4ME3' boolean column (0 or 1) to a KLF4 binding CSV file by checking for overlaps between KLF4 binding sites and H3K4ME3 regions using bedtools.
# Author: Addison Yam

# Inputs
KLF4_CSV_FILE="Human_motif.csv"
H3K4ME3_BED_FILE="H3K4me3.bed"
INPUT_DELIMITER=','

# Output File
OUTPUT_CSV_FILE="Human_histone.csv"

# Temporary Files
TEMP_KLF4_BED_WITH_IDS="temp_klf4_with_ids.bed"
TEMP_OVERLAPPING_IDS="temp_overlapping_h3k4me3_ids.txt"

# Creates IDs for each line
awk -F"$(printf "%q" "$INPUT_DELIMITER")" '
    NR==1 {next} # Skip header
    {
        # Generate a unique ID for each line (chrom:start-end:line_num)
        # Using NR (record number) as part of the ID makes it unique even for identical coordinates
        id = $1 ":" $2 "-" $3 ":" NR;
        print $1 "\t" $2 "\t" $3 "\t" id;
    }
' "$KLF4_CSV_FILE" > "$TEMP_KLF4_BED_WITH_IDS"

# Identifies KLF4 peaks and histone methylation regions
bedtools intersect -a "$TEMP_KLF4_BED_WITH_IDS" -b "$H3K4ME3_BED_FILE" -u | cut -f4 | sort -u > "$TEMP_OVERLAPPING_IDS"
echo "Done identifying KLF4 and histone methylation overlaps
date"

# Adds H3K4ME3 column to the original KLF4 CSV
awk -F"$(printf "%q" "$INPUT_DELIMITER")" -v OFS="$INPUT_DELIMITER" \
    -v overlapping_ids_file="$TEMP_OVERLAPPING_IDS" '
    BEGIN {
        # Load overlapping IDs into an associative array for quick lookup
        while ((getline id < overlapping_ids_file) > 0) {
            overlapping[id] = 1;
        }
        close(overlapping_ids_file);
    }
    NR==1 {
        print $0, "H3K4ME3"; # Prints header and adds new column name
        next;
    }
    {
        current_id = $1 ":" $2 "-" $3 ":" NR; # Reconstruct the unique ID for the current line
        
        if (current_id in overlapping) { # Check if the current ID is in our list of overlapping IDs
            print $0, 1; # Overlaps, so add 1
        } else {
            print $0, 0; # No overlap, so add 0
        }
    }
' "$KLF4_CSV_FILE" > "$OUTPUT_CSV_FILE"

rm -f "$TEMP_KLF4_BED_WITH_IDS" "$TEMP_OVERLAPPING_IDS"

echo "Done. Output saved to: '$OUTPUT_CSV_FILE'"
date
