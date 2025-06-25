import os
import sys

def expand_bed_file(input_filename):
    """
    Reads a 3-column BED file and outputs a 7-column BED file by adding
    a peak ID, two '0' fields, and a '+' field.

    Args:
        input_filename (str): The path to the input 3-column BED file.
    """

    # Generate an output filename by appending "_expanded"
    base_name = os.path.splitext(input_filename)[0]
    output_filename = base_name + "_expanded.bed"

    lines_processed = 0
    lines_written = 0

    print(f"Processing input file: {input_filename}")
    print(f"Output will be saved to: {output_filename}")

    try:
        with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
            for line_num, line in enumerate(infile, 1):
                lines_processed += 1
                parts = line.strip().split('\t')

                # Ensure the input line has at least 3 columns (chrom, chromStart, chromEnd)
                if len(parts) < 3:
                    print(f"Warning: Line {line_num} in '{input_filename}' is malformed (expected at least 3 columns, found {len(parts)}). Skipping: {line.strip()}")
                    continue # Skip this malformed line

                # Extract the first three columns
                chrom = parts[0]
                chrom_start = parts[1]
                chrom_end = parts[2]

                # Generate the new columns
                peak_id = f"peak{line_num}" # Unique ID for each peak
                field_zero_1 = "0"          # The first "0" column
                field_zero_2 = "0"          # The second "0" column
                field_plus = "+"            # The "+" column

                # Construct the new line for the output file
                # The BED format typically has: chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb
                # The request maps to: chrom, chromStart, chromEnd, peak_id, score="0", strand="0", and an additional field "+".
                # For simplicity, we'll append directly in the order requested by the user.
                # If this needs to strictly conform to BED6+ with specific fields, adjustments would be needed.
                output_parts = [
                    chrom,
                    chrom_start,
                    chrom_end,
                    peak_id,
                    field_zero_1,
                    field_zero_2,
                    field_plus
                ]

                outfile.write('\t'.join(output_parts) + '\n')
                lines_written += 1

        print(f"\n--- Processing Summary ---")
        print(f"Total lines processed from input: {lines_processed}")
        print(f"Total lines written to output: {lines_written}")
        print(f"Expanded data saved to: {output_filename}")

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python expand_bed.py <input_bed_file>")
        print("Example: python expand_bed.py my_3_column_file.bed")
        sys.exit(1)
    
    input_file = sys.argv[1]
    expand_bed_file(input_file)
