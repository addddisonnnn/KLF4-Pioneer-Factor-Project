import os
import sys

def remove_lines_from_bed(input_filename, start_line_to_remove, end_line_to_remove):
    """
    Removes lines within a specified range (inclusive) from a BED file.

    Args:
        input_filename (str): The path to the input BED file.
        start_line_to_remove (int): The 1-based starting line number to remove.
        end_line_to_remove (int): The 1-based ending line number to remove.
    """

    # Generate an output filename by appending "_filtered"
    base_name = os.path.splitext(input_filename)[0]
    output_filename = base_name + "_filtered.bed"

    lines_read = 0
    lines_kept = 0
    lines_removed = 0

    print(f"Processing input file: {input_filename}")
    print(f"Lines to remove: from {start_line_to_remove} to {end_line_to_remove} (inclusive)")
    print(f"Output will be saved to: {output_filename}")

    try:
        with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
            for line_num, line in enumerate(infile, 1): # enumerate starts from 0, so add 1 for 1-based line numbers
                lines_read += 1
                
                # Check if the current line number is within the removal range
                if start_line_to_remove <= line_num <= end_line_to_remove:
                    lines_removed += 1
                    # print(f"Removing line {line_num}: {line.strip()}") # Optional: uncomment to see removed lines
                else:
                    outfile.write(line)
                    lines_kept += 1

        print(f"\n--- Processing Summary ---")
        print(f"Total lines read: {lines_read}")
        print(f"Lines removed: {lines_removed}")
        print(f"Lines kept and written to output: {lines_kept}")
        print(f"Filtered data saved to: {output_filename}")

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python remove_bed_lines.py <input_bed_file>")
        print("Example: python remove_bed_lines.py Group1_expanded.bed")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # Define the specific line numbers to remove as requested by the user
    LINES_START = 5861
    LINES_END = 5985

    remove_lines_from_bed(input_file, LINES_START, LINES_END)
