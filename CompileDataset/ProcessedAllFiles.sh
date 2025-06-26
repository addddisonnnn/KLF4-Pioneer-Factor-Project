#!/bin/bash
#SBATCH --partition=general-compute
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --qos=general-compute
#SBATCH --job-name="Process BEDs" 

# Script Purpose: This script has three major parts. First, it organizes the bed file to be in proper format. Then, it removes any overlapping peaks by taking the first. Lastly, it checks and documents if each peak is overlapping with which cell type
# Author: Addison Yam

module load gcc/11.2.0
module load bedtools

############################################
################# PHASE ONE ################
############################################

# This phase prepares the master bed file to be in the file format of interest by concatenating all of the bed files of interest, sorting the peaks, removing unnecessary columns, and then removing unnecessary chromosomes (so it can be inputted for the GREAT tool

echo "Starting with Phase One"
date

echo "Creating a master list of all distinct peaks from all original files"
# Every file is being concatenated one after the other
cat bone_marrow.bed Mbd3f_6.bed MEF_5.bed SCC_1.bed V6_5_ESC_5.bed CD19.bed Mbd3f_7.bed MEF_6.bed SCC_2.bed V6_5_ESC_6.bed ESC.bed Mbd3f_8.bed MEF_7.bed SCC_3.bed V6_5_ESC_7.bed Mbd3f_1.bed Mbd3f_9.bed MEF_8.bed SCC_4.bed V6_5_ESC_8.bed Mbd3f_2.bed MEF_1.bed pre_iPSC_1.bed V6_5_ESC_1.bed V6_5_ESC_9.bed Mbd3f_3.bed MEF_2.bed pre_iPSC_2.bed V6_5_ESC_2.bed Mbd3f_4.bed MEF_3.bed V6_5_ESC_3.bed Mbd3f_5.bed MEF_4.bed Prostate_Stem.bed V6_5_ESC_4.bed > all_original_peaks.bed
echo "Done creating master bed file of all the original peaks"

echo "Sorting original peaks"
# Sorts first by the first column (chromosome), then the start, and finally end coordinates, also removes any duplicates 
sort -k1,1V -k2,2n -k3,3n -u all_original_peaks.bed > master_all_distinct_peaks.bed
echo "Done with sorting"

echo "Remove unnecessary columns" 
# Keeps the first, second, and third columns
cut -f1-3 master_all_distinct_peaks.bed > master_trimmed.bed
echo "Done trimming"

echo "Filtering out unnecessary chromosomes"
# List of standard mm10 chromosomes ===
VALID_CHROMS=$(echo -e "chr1\nchr2\nchr3\nchr4\nchr5\nchr6\nchr7\nchr8\nchr9\nchr10\nchr11\nchr12\nchr13\nchr14\nchr15\nchr16\nchr17\nchr18\nchr19\nchrX\nchrY\nchrM")
# Converts list into grep pattern ===
PATTERN=$(echo "$VALID_CHROMS" | paste -sd'|' -)
# Filters BED file
grep -E "^($PATTERN)[[:space:]]" master_trimmed.bed > concat_sorted_trimmed.bed
echo "Filtered BED written to concat_sorted_trimmed.bed"
echo "Done with Phase One"
date

############################################
################# PHASE TWO ################
############################################

# This phase removes any peaks that overlap with each other and only takes the first peak.
echo "Starting with Phase Two"
date
echo "Generating master bed file of non-overlapping regions"
# Create empty output file
> master_non_overlapping.bed

# Uses awk so that we keep the first peak by using the first peak's start and end and not recording what is in between
awk 'BEGIN {OFS="\t"; prev_chrom=""; prev_end=0}
        {
        # If this peak is on a new chromosome OR it starts after the previous peak ends,
        # then it is a new, non-overlapping peak to add to our master list.
        if ($1 != prev_chrom || $2 >= prev_end) {
                print $1, $2, $3; # Print chr, start, end of the current peak
                prev_chrom = $1;
                prev_end = $3;
        } else {
                # If this peak overlaps with the previous one, we have already included the first
                # peak (the one that came before it in the sorted file) that covers this region.
                # We update prev_end if the current peak extends further, to ensure future check
                # are against the furthest extent of the currently tracked "first" region.
                if ($3 > prev_end) {
                        prev_end = $3;
                }
        }
        }' concat_sorted_trimmed.bed \
    > master_non_overlapping.bed

echo "Done with creating list of non-overlapping peaks" 
echo "Done with Phase Two"
date

############################################
############### PHASE THREE ################
############################################

# This phase generates the final bed file with all the boolean values for whether or not this peak is found in each cell line specified here

echo "Starting Phase Three"
date
echo "Generating final boolean table"

# Keeps track of all the unique cell type names, these correspond to the names of the original BED files and are going to the column headers
unique_cell_types=(bone_marrow Mbd3f MEF SCC V6_5_ESC CD19 ESC pre_iPSC Prostate_Stem)

# Array that maps the cell type name to their original filenames
declare -A cell_type_files
cell_type_files[bone_marrow]="bone_marrow.bed"
cell_type_files[Mbd3f]="Mbd3f_1.bed Mbd3f_2.bed Mbd3f_3.bed Mbd3f_4.bed Mbd3f_5.bed Mbd3f_6.bed Mbd3f_7.bed Mbd3f_8.bed Mbd3f_9.bed"
cell_type_files[MEF]="MEF_1.bed MEF_2.bed MEF_3.bed MEF_4.bed MEF_5.bed MEF_6.bed MEF_7.bed MEF_8.bed"
cell_type_files[SCC]="SCC_1.bed SCC_2.bed SCC_3.bed SCC_4.bed"
cell_type_files[V6_5_ESC]="V6_5_ESC_1.bed V6_5_ESC_2.bed V6_5_ESC_3.bed V6_5_ESC_4.bed V6_5_ESC_5.bed V6_5_ESC_6.bed V6_5_ESC_7.bed V6_5_ESC_8.bed V6_5_ESC_9.bed"
cell_type_files[CD19]="CD19.bed"
cell_type_files[ESC]="ESC.bed"
cell_type_files[pre_iPSC]="pre_iPSC_1.bed pre_iPSC_2.bed"
cell_type_files[Prostate_Stem]="Prostate_Stem.bed"

# Create the header row
header="chrom\tstart\tend"
for cell_type_name in "${unique_cell_types[@]}"; do
        header="{$header}\t${cell_type_name}"
done
# Writes the header to and creates the output file
echo -e "$header" > master.bed

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
	echo -e "$output_line" >> master.bed

done < master_non_overlapping.bed

rm -f current_master_peak.bed all_original_peaks.bed master_all_distinct_peaks.bed master_trimmed.bed concat_sorted_trimmed.bed master_non_overlapping.bed

echo "Done with generating final boolean table"
echo "Done with Phase Three"
date
