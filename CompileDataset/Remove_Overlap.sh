# This script removes any peaks that overlap with each other and only takes the first peak 
# Author: Addison Yam

echo "Generating master bed file of non-overlapping regions"
# Create empty output file
> master_non_overlapping.bed

# Uses awk so that we keep the first peak by using the first peak's start and 
end and not recording what is in between
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
        }' master.bed \
    > master_non_overlapping.bed

echo "Done with creating list of non-overlapping peaks"
