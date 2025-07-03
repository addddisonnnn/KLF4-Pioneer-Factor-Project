# Script to run fimo and the steps necessary to prepare the inputs for using MEME-Suite's FIMO tool

# Commands to run fimo:
module load gcc
module load openmpi
module load meme/5.5.7

# Get the meme files from searching KLF4 on Jasper. And used the MA0039.2  for mouse and ​​MA0039.5 for human. 
# human.meme and mouse.meme

# Need to download reference genome files to convert BED to FASTA, genomes from https://hgdownload.soe.ucsc.edu/downloads.html
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
mv hg38.fa humansequence.fa
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz
mv mm10.fa mousesequence.fa

# Convert the bed files to fasta. 
module load bedtools
bedtools getfasta -fi hg38.fa -bed humanmaster.bed -fo humansequence.fa
bedtools getfasta -fi mm10.fa -bed mousemaster.bed -fo mousesequence.fa 

# Run FIMO
# fimo -oc humanmotif human.meme humansequence.fa (more simple version if 100,000 is enough motif scores but I needed more)
fimo --thresh 1e-3 --max-stored-scores 1000000 --parse-genomic-coord -oc humanmotif human.meme humansequence.fa
# fimo -oc mousemotif mouse.meme mousesequence.fa
fimo --thresh 1e-3 --max-stored-scores 1000000 --parse-genomic-coord -oc mousemotif mouse.meme mousesequence.fa

# Output saved to humanmotif and mousemotif directories, relevant files are named fimo.tsv within these directories
