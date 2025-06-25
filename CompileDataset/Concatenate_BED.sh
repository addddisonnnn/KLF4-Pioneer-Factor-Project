# Script Purpose: Concatenate all of the bed files of interest. Then sort the peaks.
# Author: Addison Yam

echo "Creating a master list of all distinct peaks from all original files"
# Every file is being concatenated one after the other
cat BJ_1.bed BJ_2.bed BJ_3.bed BJ_4.bed BJ_5.bed H9_ESC.bed iPSC.bed LIS49_hESC.bed MCF-7_1.bed MCF-7_2.bed MCF-7_3.bed U87.bed > all_original_peaks.bed
echo "Done creating master bed file of all the original peaks"

echo "Sorting original peaks"
# Sorts first by the first column (chromosome), then the start, and finally end coordinates, also removes any duplicates 
sort -k1,1 -k2,2n -k3,3n -u all_original_peaks.bed > master_all_distinct_peaks.bed
echo "Done with sorting"

echo "Remove unnecessary columns" 
# Keeps the first, second, and third columns
cut -f1-3 master_all_distinct_peaks.bed > master_trimmed.bed
echo "Done trimming"

echo "Filtering out unnecessary chromosomes"
# List of standard hg38 chromosomes ===
VALID_CHROMS=$(echo -e "chr1\nchr2\nchr3\nchr4\nchr5\nchr6\nchr7\nchr8\nchr9\nchr10\nchr11\nchr12\nchr13\nchr14\nchr
15\nchr16\nchr17\nchr18\nchr19\nchr20\nchr21\nchr22\nchrX\nchrY\nchrM")
# Converts list into grep pattern ===
PATTERN=$(echo "$VALID_CHROMS" | paste -sd'|' -)
# Filters BED file
grep -E "^($PATTERN)[[:space:]]" master_trimmed.bed > master.bed
echo "Filtered BED written to master.bed"
