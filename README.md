# KLF4-Pioneer-Factor-Project
Documenting my progress towards discovering relevant determinants and bettering our understanding the mechanisms and limitations of KLF4 pioneer factor ability


# Introduction
Transcription factors have an essential role in regulating gene expression by binding to particular genomic sites. A subclass of TFs called pioneer transcription factors possesses the characteristic of binding to and opening inaccessible condensed chromatin, which is crucial for a myriad of developmental processes such as cell differentiation and reprogramming. Chromatin, composed of DNA wrapped around core histone proteins forming nucleosomes, is organized into complex 3D structures. PFs bind to regions of high amounts of motifs at the ends and loop junctions of nucleosomes. Due to nucleosome structure, PFs are directed toward enhancers, thereby activating transcription.

![image](https://github.com/user-attachments/assets/8a5ae1d3-1b3f-4741-9dd3-78b433b19b2b)

(Clarke, Thomas & Mostoslavsky, Raul. (2022). DNA repair as a shared hallmark in cancer and ageing. Molecular Oncology. 16. 10.1002/1878-0261.13285.)

# Objective
This project aims to mechanistically understand how transcription factors selectively bind and regulate gene expression, with a focus on Krüppel-like factor 4 (KLF4). The objective is to highlight key predictive variables, such as motifs, co-factor motifs, chromatin accessibility, histone modifications, and other genomic features.


# Main Steps (folders with more information available)
## Step 1: Compile Datasets
I got the majority of my data from the Cistrome database, allowing for simple search and filtering for KLF4. All of the datasets on the database are ChIP-seq data, making our lives so much easier. They are also converted to the most up-to-date build, and all of their data has been converted into BED files. We chose datasets that met our standards (>1K peaks). 

## Step 2: Process & Clean our Data
This is where we figure out potential variables, factors, or determinants of interest. Histone methylation/modifications, other transcription factors, and more are examples. I ran MEME-Suite's FIMO tool to get information about whether the KLF4 motif is bound at the region, KLF4 motif score, strand (+,-), and total KLF4 motif score (because sometimes one region can have several KLF4 motifs bound potentially).

## Step 3: ML Model Training
Utilized machine learning techniques to study KLF4 binding in different cell types based on a myriad of determinants. Random forests and decision trees were the first techniques used. Also, two cell types we focused on were Head-Neck Sarcoma cells and MCF7 cells.

## Step 4: Employ Predictive Models

# Abstract / Overall Summary
Transcription factors play an essential role in regulating gene expression by binding to particular genomic sites, and a subclass called pioneer transcription factors possesses the characteristic of binding to and opening inaccessible condensed chromatin within nucleosomes, crucial for various developmental processes such as cell differentiation and reprogramming. This project aims to mechanistically understand how TFs selectively bind to genomic sites and regulate gene expression, with a specific focus on Krüppel-like factor 4 (KLF4), a pioneer factor and one of four Yamanaka factors capable of reprogramming somatic cells into stem cells. To explore the determinants of KLF4 binding specificity, we apply supervised learning techniques to highlight key predictive variables, including motifs, co-factor motifs, histone modifications, and other genomic features impacting pioneer factor ability. The project is currently curating a compendium of KLF4 CHiP-seq data to train subsequent models. Visualizations of nucleosome positioning and establishments of motifs of KLF4 binding in the cell line MCF-7 associated with breast cancer have occurred. Main goals include curation of a KLF4 CHiP-seq dataset compendium across cell lines in human and mice, training models, and developing a computational pipeline. A better understanding of KLF4 may lead to more experiments on stem cells, TFs, and gene expression.
