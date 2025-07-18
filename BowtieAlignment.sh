# Script: Intended to align fasta files with reference genome via bowtie alignment

echo "Loading Modules . . ."

module load gcc/11.2.0
module load bowtie2/2.4.4
module load samtools/1.16.1

echo "Running bowtie2 . . . "
bowtie2 -p 16 -f -x /projects/academic/mjbuck/References/igenome-041714/Bowtie2_Indexes/mm9 
-1 SRR18970853_1.fasta -2 SRR18970853_2.fasta 2> MEFs_MNase4U_alignmentstats | samtools sort
 -@ 12 -O bam -o MEFs__MNase4U_sorted.bam

# -x : location of where the reference genomes are
# -1 : first fasta file (usually from SRR/NCBI)
# -2 : second fasta file (usually from SRR/NCBI)
# samtools converts sam to bam file

# the following is if there is one fasta input, not two
echo "Running bowtie2 . . . "
bowtie2 -p 16 -f -x /projects/academic/mjbuck/References/genomes-01202020/human/gencode/R36.GRCh38.p13/Sequence/Bowtie2/GRCh38.p13.genome 
-U H3K27me3.fasta 2> H3K27me3_alignmentstats2 | samtools sort -@ 12 -O bam -o H3K27me3_sorted.bam
