import pandas as pd

def combine_csv_files(file1_path, file2_path, output_file_path):
    """
    Combines two CSV files with different headers.

    Retains all columns from file1 (Human_HN_Histone_Overlaps.csv),
    renames 'H3K4ME3' to 'H3K4me3' and moves it to the end.
    Adds specific histone modification columns from file2 (Human_Histone_Overlaps.csv).

    Args:
        file1_path (str): Path to the first CSV file
                          (e.g., "Human_HN_Histone_Overlaps.csv").
        file2_path (str): Path to the second CSV file
                          (e.g., "Human_Histone_Overlaps.csv").
        output_file_path (str): Path for the combined output CSV file.
    """
    print(f"Loading first CSV file: '{file1_path}'...")
    df1 = pd.read_csv(file1_path)
    print(f"Loading second CSV file: '{file2_path}'...")
    df2 = pd.read_csv(file2_path)

    # Define common columns for merging. These are assumed to uniquely identify rows.
    merge_cols = ['chrom', 'start', 'end']

    # --- Step 1: Prepare df1 (Human_HN_Histone_Overlaps.csv) ---
    # Rename H3K4ME3 column in df1. We use a temporary name to avoid conflicts
    # if 'H3K4me3' already exists before the final reordering.
    if 'H3K4ME3' in df1.columns:
        df1 = df1.rename(columns={'H3K4ME3': 'H3K4me3_temp_final_pos'})
    else:
        print(f"Warning: 'H3K4ME3' column not found in '{file1_path}'. It will not be renamed or moved.")

    # --- Step 2: Select specific columns from df2 (Human_Histone_Overlaps.csv) ---
    # These are the columns to be added from the second file.
    cols_to_add_from_df2 = [
        'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me2',
        'H3K79me2', 'H3K9ac', 'H3K9me2', 'H3K9me3', 'H4K20me1'
    ]

    # Filter for columns that actually exist in df2 to avoid KeyError
    existing_cols_df2 = [col for col in cols_to_add_from_df2 if col in df2.columns]
    missing_cols_df2 = [col for col in cols_to_add_from_df2 if col not in df2.columns]

    if missing_cols_df2:
        print(f"Warning: The following requested columns were not found in '{file2_path}': {missing_cols_df2}. They will not be added.")

    # Create a subset of df2 containing only the merge_cols and the desired additional columns
    df2_subset = df2[merge_cols + existing_cols_df2]

    # --- Step 3: Merge the two dataframes ---
    # An inner merge ensures that only rows present in both files (based on merge_cols)
    # are included in the final output. This assumes 'chrom', 'start', 'end'
    # uniquely identify the rows and that the files have the "same rows" as stated.
    print("Merging dataframes based on 'chrom', 'start', 'end'...")
    combined_df = pd.merge(df1, df2_subset, on=merge_cols, how='inner')

    # --- Step 4: Reorder columns and finalize 'H3K4me3' ---
    # Get the current list of columns in the merged DataFrame
    current_cols = combined_df.columns.tolist()

    final_column_order = []
    h3k4me3_col_present = False

    # Build the new column order
    for col in current_cols:
        if col == 'H3K4me3_temp_final_pos':
            h3k4me3_col_present = True
            # Don't add it yet, we'll add it at the very end
            continue
        final_column_order.append(col)

    # Add the renamed H3K4me3 column to the very end if it was present
    if h3k4me3_col_present:
        final_column_order.append('H3K4me3')
        # Assign the values from the temporary column to the final named column
        combined_df['H3K4me3'] = combined_df['H3K4me3_temp_final_pos']
        # Drop the temporary column
        combined_df = combined_df.drop(columns=['H3K4me3_temp_final_pos'])

    # Reindex the DataFrame with the new, desired column order
    combined_df = combined_df[final_column_order]

    # --- Step 5: Save the combined DataFrame to CSV ---
    print(f"Saving combined data to '{output_file_path}'...")
    combined_df.to_csv(output_file_path, index=False)
    print("Combination complete!")
    print(f"Final combined CSV saved to: '{output_file_path}'")


# --- Main execution block ---
if __name__ == "__main__":
    # Define your input and output file paths here
    file1_path = "Human_HN_Histone_Overlaps.csv"
    file2_path = "Human_Histone_Overlaps.csv"
    output_file_path = "Combined_Human_Histone_Data.csv"

    combine_csv_files(file1_path, file2_path, output_file_path)
