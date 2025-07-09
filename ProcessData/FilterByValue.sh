#!/bin/bash

# Script to filter a delimited file (CSV) based on a specified column containing a '1' (boolean True).
# The output will now be a standard 3-column BED file (chrom, start, end) and be TAB-DELIMITED.
# This version has hardcoded input parameters for automation.

echo "Starting to filter data
date

# Specify the input file, the delimiter, and column name to filter by
input_file="oldHuman_motif.csv"
INPUT_DELIMITER="," # Set to comma (',') as your file is comma-separated
column_name="MCF7" # As per your example

# Create an output file name based on the input file prefix and column name
base_name=$(basename "$input_file")
file_prefix="${base_name%.*}" # Remove extension
# Output file will now be a standard 3-column BED file
output_file="${file_prefix}_filtered_by_${column_name// /_}.bed" # MODIFIED: Removed _all_cols

echo "Input File: '$input_file'"
echo "Input Delimiter: '$(printf "%q" "$INPUT_DELIMITER")'" # Print delimiter for confirmation
echo "Filtering Column: '$column_name'"
echo "Output File: '$output_file' (3-column Tab-delimited BED)" # MODIFIED: Updated description

# Filtering Logic

# 1. Get the header line from the input file (used only for finding column index)
header_line=$(head -n 1 "$input_file")

# Check if the header line is empty (e.g., if file is empty)
if [ -z "$header_line" ]; then
    echo "Error: Input file '$input_file' appears to be empty or has no header."
    exit 1
fi

# 2. Find the column index (1/True-based) of the specified column name
#    We use 'awk' to split the header by the specified INPUT_DELIMITER and find the position.
column_index=$(echo -e "$header_line" | awk -F"$(printf "%q" "$INPUT_DELIMITER")" -v col_name="$column_name" '{
    for (i=1; i<=NF; i++) {
        # Using a regex match (~=) and trimming spaces for robustness
        if (gensub(/^[[:space:]]+|[[:space:]]+$/, "", "g", $i) ~ "^" col_name "$") {
            print i;
            exit;
        }
    }
}')

# Check if the column name was found
if [ -z "$column_index" ]; then
    echo "Error: Column '$column_name' not found in the header of '$input_file'."
    echo "Available columns are: $(echo -e "$header_line" | tr "$INPUT_DELIMITER" '\n' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | paste -s -d',')"
    exit 1
fi

echo "Column '$column_name' found at index: $column_index"
echo "Filtering rows where '$column_name' is '1' and outputting only chrom, start, end..."

# 3. Filter the file and output only the first three columns (tab-delimited)
#    - Print the new 3-column header.
#    - Then, process the rest of the file using awk:
#      - Set input field separator (-F) to the INPUT_DELIMITER.
#      - Set output field separator (OFS) to a tab.
#      - Skip the header (NR==1{next}).
#      - If the value in the target column ($col_idx) is '1':
#        - Print the first three fields ($1, $2, $3) separated by OFS (tab).
{
    echo -e "chrom\tstart\tend" # MODIFIED: Explicitly print the 3-column header
    awk -F"$(printf "%q" "$INPUT_DELIMITER")" -v col_idx="$column_index" '
        NR==1 {next} # Skip the original header
        BEGIN{OFS="\t"}
        $col_idx == 1 {
            print $1, $2, $3 # MODIFIED: Print only the first three fields
        }
    ' "$input_file"
} > "$output_file"

# --- Completion Message ---
if [ $? -eq 0 ]; then
    echo "Filtering complete!"
    echo "Filtered 3-column BED file saved to: '$output_file'" # MODIFIED: Updated description
    echo "Number of filtered rows (including header): $(wc -l < "$output_file")"
else
    echo "An error occurred during filtering."
fi

echo ""
