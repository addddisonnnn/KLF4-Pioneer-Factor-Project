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
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz

# Convert the bed files to fasta. The bed files contain the chrom, start, and end which is all needed as well as the reference genomes to get the sequences of interest
module load bedtools
bedtools getfasta -fi hg38.fa -bed humanmaster.bed -fo humansequence.fa
bedtools getfasta -fi mm10.fa -bed mousemaster.bed -fo mousesequence.fa 

# Run FIMO
# fimo -oc humanmotif human.meme humansequence.fa (more simple version if 100,000 is enough motif scores but I needed more)
fimo --thresh 1e-3 --max-stored-scores 1000000 --parse-genomic-coord -oc humanmotif human.meme humansequence.fa
# fimo -oc mousemotif mouse.meme mousesequence.fa
fimo --thresh 1e-3 --max-stored-scores 1000000 --parse-genomic-coord -oc mousemotif mouse.meme mousesequence.fa

# Output saved to humanmotif and mousemotif directories, relevant files are named fimo.tsv within these directories

# More running of FIMO with more potential motifs
fimo --thresh 1e-3 --max-stored-scores 1000000 --parse-genomic-coord -oc ZNF384 ZNF384.meme humansequence.fa
fimo --thresh 1e-3 --max-stored-scores 1000000 --parse-genomic-coord -oc OTX2 OTX2.meme humansequence.fa
fimo --thresh 1e-3 --max-stored-scores 1000000 --parse-genomic-coord -oc TP53 TP53.meme humansequence.fa
fimo --thresh 1e-3 --max-stored-scores 1000000 --parse-genomic-coord -oc ZBED4 ZBED4.meme humansequence.fa
fimo --thresh 1e-3 --max-stored-scores 1000000 --parse-genomic-coord -oc ZNF460 ZNF460.meme humansequence.fa
fimo --thresh 1e-3 --max-stored-scores 1000000 --parse-genomic-coord -oc ZBTB6 ZBTB6.meme humansequence.fa
fimo --thresh 1e-3 --max-stored-scores 1000000 --parse-genomic-coord -oc FOSL2 FOSL2.meme humansequence.fa
fimo --thresh 1e-3 --max-stored-scores 1000000 --parse-genomic-coord -oc ZNF213 ZNF213.meme humansequence.fa
fimo --thresh 1e-3 --max-stored-scores 1000000 --parse-genomic-coord -oc POU6F2 POU6F2.meme humansequence.fa
fimo --thresh 1e-3 --max-stored-scores 1000000 --parse-genomic-coord -oc FOXC2 FOXC2.meme humansequence.fa
