#!/bin/bash

# Script to add boolean columns (0 or 1) for multiple histone methylations
# to a KLF4 binding CSV file by checking for overlaps between KLF4 binding sites
# and histone methylation regions using bedtools.
# This version handles gzipped (.bed.gz) histone methylation files.

# --- Configuration ---

# Input KLF4 Binding CSV File (expected format: chrom,start,end,...)
KLF4_CSV_FILE="Human_histone_with_motifs.csv"

# Directory containing the gzipped histone methylation BED files
HISTONE_BED_DIR="Histone_BED_Files" # <--- IMPORTANT: Adjust this path!

# List of histone methylation marks to check (must match BED file names without '.bed.gz')
HISTONE_MARKS=(
    "H2AFZ" "H3F3A" "H3K27ac" "H3K27me3" "H3K36me3" "H3K4me1"
    "H3K4me2" "H3K79me2" "H3K9ac" "H3K9me2" "H3K9me3" "H4K20me1"
)

# Delimiter for the input and output CSV files
INPUT_DELIMITER=','

# Output File
OUTPUT_CSV_FILE="Human_Histone_Overlaps.csv" # A more descriptive name

# --- Temporary Files ---
TEMP_KLF4_BED_WITH_IDS="temp_klf4_with_ids.bed"
TEMP_CURRENT_OUTPUT="temp_current_output.csv"
TEMP_OVERLAPPING_IDS_PREFIX="temp_overlapping_ids" # Prefix for multiple temp files
TEMP_UNZIPPED_BEDS_DIR="temp_unzipped_beds" # Directory for unzipped BED files

# --- Main Script ---

echo "Starting histone methylation overlap analysis: $(date)"

# Create temporary directory for unzipped files
mkdir -p "$TEMP_UNZIPPED_BEDS_DIR"
echo "Created temporary directory for unzipped BED files: '$TEMP_UNZIPPED_BEDS_DIR'"

# 1. Create a temporary BED file for KLF4 sites with unique IDs
echo "1. Generating KLF4 BED file with unique IDs..."
awk -F"$(printf "%q" "$INPUT_DELIMITER")" '
    NR==1 {next} # Skip header
    {
        # Generate a unique ID for each line (chrom:start-end:line_num)
        # Using NR (record number) as part of the ID makes it unique even for identical coordinates
        id = $1 ":" $2 "-" $3 ":" NR;
        print $1 "\t" $2 "\t" $3 "\t" id;
    }
' "$KLF4_CSV_FILE" > "$TEMP_KLF4_BED_WITH_IDS"
echo "   KLF4 BED file created: '$TEMP_KLF4_BED_WITH_IDS'"

# Initialize the output file with the original KLF4 data.
# This file will be progressively updated with new columns.
cp "$KLF4_CSV_FILE" "$TEMP_CURRENT_OUTPUT"

# 2. Iterate through each histone methylation mark and add a column
for MARK in "${HISTONE_MARKS[@]}"; do
    GZIPPED_HISTONE_BED_FILE="${HISTONE_BED_DIR}/${MARK}.bed.gz"
    UNZIPPED_HISTONE_BED_FILE="${TEMP_UNZIPPED_BEDS_DIR}/${MARK}.bed"
    TEMP_OVERLAPPING_IDS="${TEMP_OVERLAPPING_IDS_PREFIX}_${MARK}.txt"

    echo "2. Processing histone mark: ${MARK}"

    if [ ! -f "$GZIPPED_HISTONE_BED_FILE" ]; then
        echo "   WARNING: Gzipped histone BED file not found: '$GZIPPED_HISTONE_BED_FILE'. Skipping this mark and adding a column of zeros."
        # If the gzipped file is missing, we still need to add a column of zeros for consistency
        awk -F"$(printf "%q" "$INPUT_DELIMITER")" -v OFS="$INPUT_DELIMITER" \
            -v mark_name="$MARK" '
            NR==1 { print $0, mark_name; next; }
            { print $0, 0; }
        ' "$TEMP_CURRENT_OUTPUT" > "${TEMP_CURRENT_OUTPUT}.new"
        mv "${TEMP_CURRENT_OUTPUT}.new" "$TEMP_CURRENT_OUTPUT"
        continue
    fi

    echo "   Unzipping '$GZIPPED_HISTONE_BED_FILE' to '$UNZIPPED_HISTONE_BED_FILE'..."
    gunzip -c "$GZIPPED_HISTONE_BED_FILE" > "$UNZIPPED_HISTONE_BED_FILE"
    if [ $? -ne 0 ]; then
        echo "   ERROR: Failed to unzip '$GZIPPED_HISTONE_BED_FILE'. Skipping this mark."
        # Clean up if gunzip failed
        rm -f "$UNZIPPED_HISTONE_BED_FILE"
        awk -F"$(printf "%q" "$INPUT_DELIMITER")" -v OFS="$INPUT_DELIMITER" \
            -v mark_name="$MARK" '
            NR==1 { print $0, mark_name; next; }
            { print $0, 0; }
        ' "$TEMP_CURRENT_OUTPUT" > "${TEMP_CURRENT_OUTPUT}.new"
        mv "${TEMP_CURRENT_OUTPUT}.new" "$TEMP_CURRENT_OUTPUT"
        continue
    fi

    # Identify KLF4 peaks that overlap with the current histone methylation regions
    echo "   Running bedtools intersect for ${MARK}..."
    bedtools intersect -a "$TEMP_KLF4_BED_WITH_IDS" -b "$UNZIPPED_HISTONE_BED_FILE" -u | cut -f4 | sort -u > "$TEMP_OVERLAPPING_IDS"
    echo "   Identified overlaps for ${MARK}."

    # Add the current histone methylation column to the growing CSV file
    awk -F"$(printf "%q" "$INPUT_DELIMITER")" -v OFS="$INPUT_DELIMITER" \
        -v overlapping_ids_file="$TEMP_OVERLAPPING_IDS" \
        -v mark_name="$MARK" '
        BEGIN {
            # Load overlapping IDs into an associative array for quick lookup
            while ((getline id < overlapping_ids_file) > 0) {
                overlapping[id] = 1;
            }
            close(overlapping_ids_file);
        }
        NR==1 {
            print $0, mark_name; # Print header and add new column name
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
    ' "$TEMP_CURRENT_OUTPUT" > "${TEMP_CURRENT_OUTPUT}.new"

    # Replace the old temporary output with the new one containing the added column
    mv "${TEMP_CURRENT_OUTPUT}.new" "$TEMP_CURRENT_OUTPUT"

    # Clean up temporary files for the current mark
    rm -f "$TEMP_OVERLAPPING_IDS" "$UNZIPPED_HISTONE_BED_FILE" # Clean up unzipped BED file here
    echo "   Column '${MARK}' added and temporary files cleaned."
done

# 3. Finalize output and cleanup
echo "3. Finalizing output and cleaning up..."
mv "$TEMP_CURRENT_OUTPUT" "$OUTPUT_CSV_FILE"
rm -f "$TEMP_KLF4_BED_WITH_IDS" "${TEMP_OVERLAPPING_IDS_PREFIX}"*.txt # Clean up any lingering temp ID files
rmdir "$TEMP_UNZIPPED_BEDS_DIR" # Remove the temporary unzipped directory

echo "Done. Final output saved to: '$OUTPUT_CSV_FILE'"
echo "Script finished: $(date)"
