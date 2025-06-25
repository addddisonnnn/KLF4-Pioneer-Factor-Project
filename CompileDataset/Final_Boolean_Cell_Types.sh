#!/bin/bash
#SBATCH --partition=general-compute
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --qos=general-compute
#SBATCH --job-name="Bool_Cell_Lines" 


# This script generates the final bed file with all the boolean values for whether or not this peak is found in each cell line specified here
# Author: Addison Yam

module load gcc/11.2.0
module load bedtools

echo "Generating final boolean table"

# Keeps track of all the unique cell type names, these correspond to the names of the original BED files and are going to the column headers
unique_cell_types=(BJ H9_ESC iPSC LIS49_hESC MCF7 U87)

# Array that maps the cell type name to their original filenames
declare -A cell_type_files
cell_type_files[BJ]="BJ_1.bed BJ_2.bed BJ_3.bed BJ_4.bed BJ_5.bed"
cell_type_files[H9_ESC]="H9_ESC.bed"
cell_type_files[iPSC]="iPSC.bed"
cell_type_files[LIS49_hESC]="LIS49_hESC.bed"
cell_type_files[MCF7]="MCF7_1.bed MCF7_2.bed MCF7_3.bed"
cell_type_files[U87]="U87.bed"

# Write the header row
header="chrom\tstart\tend"
for cell_type_name in "${unique_cell_types[@]}"; do
	header="{$header}\t${cell_type_name}"
done
echo -e "header" > final_boolean_peaks.bed

#Process each peak in the master bed file
while IFS=$'\t' read -r chrom start end; do
	output_line="$chrom\t$start\t$end" # Start building the output line with BED coordinates

	# Create a temporary file containing only the current master peak for bedtools.
	echo -e "$chrom\t$start\t$end" > current_master_peak.bed

	# For each unique cell type, check for overlap with its *original* replicate files
	for cell_type_name in "${unique_cell_types[@]}"; do
		original_files="${cell_type_files[$cell_type_name]}"

		is_present=0 # Assume 0 (not present) initially
		# Iterate through each original file for the current cell type
		for original_file in $original_files; do
			# Check if the current master peak overlaps with any peak in this original file
			if bedtools intersect -a current_master_peak.bed -b "${original_file}" -u | grep -q .; then
				is_present=1 # Found an overlap in at least one replicate
				break      # No need to check other replicates for this cell type; it's present
			fi
		done
		# Append the boolean result for this cell type to the output line
		output_line="${output_line}\t${is_present}"
	done

	# Write the complete output line for the current master peak to the .bed file
	echo -e "$output_line" >> final_boolean_peaks.bed

done < master_non_overlapping.bed

rm -f current_master_peak.bed

echo "Done with generating final boolean table"
