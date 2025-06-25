# Script designed to download SRA samples
# Author: Addison Yam

module load gcc/11.2.0 
module load openmpi/4.1.1 
module load sra-toolkit

nohup fastq-dump --split-files --fasta SRR18970927 # (Sample GSM6077127 Klf4-ChIP-seq, mESCs, biol rep 1)
nohup fastq-dump --split-files --fasta SRR18970928 # (Sample GSM6077128 Klf4-ChIP-seq, mESCs, biol rep 2)
