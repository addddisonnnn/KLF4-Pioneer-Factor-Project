# Analysis Performed on several genomic files

## ArchTEx
Architectural Tag Extender - identifies the average length of the DNA fragments that were sequenced using cross-correlation of single-read sequencing. 

Ran ArchTEx on all MNase-seq experiment results to determine which is the most appropriate for analysis. The following ArchTEx runs were run with a reference genome of mm9 build. These MNase-seq datasets were found in a paper titled [Nucleosome fibre topology guides transcription factor binding to enhancers](https://www.nature.com/articles/s41586-024-08333-9). Downloading these MNase files gave us fasta files which we needed to run through BowtieAlignment to align these peaks with a reference genome, which gave us SAM files, which were converted with samtools to BAM files. 

### MNase-seq concentration of 1U ![MEFs_MNase1U](https://github.com/user-attachments/assets/910ac2cd-8476-4555-84ca-6d3bbee69f80)

### MNase-seq concentration of 4U ![MEFs_MNase4U](https://github.com/user-attachments/assets/4a4a8abb-d1bd-4bdf-89c5-f451fd8a0458)

### MNase-seq concentration of 16U ![MEFs_MNase16U](https://github.com/user-attachments/assets/a3421028-1bc5-4aef-9ba2-191215b850fe)


### MNase-seq concentration of 64U ![MEFs_MNase64U](https://github.com/user-attachments/assets/00333699-65b6-448f-983f-1dfcad2c81c0)

MNase-seq 64U proved to be the best at getting the exact positions of nucleosomes. 

### Ran ArchTEx with the MNase-seq 64U data as the bam file and KLF4 ChIP-seq data as the reference file. ![MEF_MNase64U_KLF4_ChIP](https://github.com/user-attachments/assets/6683f4b9-a71f-490d-88ff-04b7da48f0e0)
Since it's hard to see where the nucleosomes are in the picture, here is the output from ArchTEx visualized in Excel. ![KLF4 ChIP-seq with MNase64U](https://github.com/user-attachments/assets/3c38a851-5086-47ae-8a29-92d597cd9bfb)
