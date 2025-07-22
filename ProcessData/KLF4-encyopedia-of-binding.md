## Encyclopedia info KLF4BindingMasterFile.csv
This keeps track of all the columns in the Master file, what it should store, what it means, and where it came from

|Column Name|Data Type|Description|Reference|
|---|---|---|--|
|chromosome name|str|---|--|
|start|int|--|--|
|stop|int|---|--|
|KLF4 motif|bool|---|--|
|KLF4 motif score|float|---|--|
|Number of KLF4 motif|int|---|--|

### Progression of columns of the human master file: 
chrom    start   end    BJ     H9_ESC iPSC   LIS49_hESC     MCF7   U87    HN_SCC
chrom,start,end,BJ_epi,BJ_iPSC,H9_ESC,iPSC,LIS49_hESC,MCF7_endo,MCF7_induc,U87,HN_SCC,Klf4 motif,Klf4 motif score,Number of Klf4 motifs,Strand,Total motif score,H3K4ME3
chrom,start,end,BJ_epi,BJ_iPSC,H9_ESC,iPSC,LIS49_hESC,MCF7_endo,MCF7_induc,U87,HN_SCC,Klf4 motif,Klf4 motif score,Number of Klf4 motifs,Strand,Total motif score,H3K4ME3,FOSL2,FOXC2,OTX2,POU6F2,TP53,ZBED4,ZBTB6,ZNF213,ZNF384,ZNF460,H2AFZ,H3F3A,H3K27ac,H3K27me3,H3K36me3,H3K4me1,H3K4me2,H3K79me2,H3K9ac,H3K9me2,H3K9me3,H4K20me1
chrom,start,end,BJ_epi,BJ_iPSC,H9_ESC,iPSC,LIS49_hESC,MCF7_endo,MCF7_induc,U87,HN_SCC,Klf4 motif,Klf4 motif score,Number of Klf4 motifs,Strand,Total motif score,H3K4ME3,FOSL2 motif,FOXC2 motif,OTX2 motif,POU6F2 motif,TP53 motif,ZBED4 motif,ZBTB6 motif,ZNF213 motif,ZNF384 motif,ZNF460 motif,FOSL2 motif score,FOXC2 motif score,OTX2 motif score,POU6F2 motif score,TP53 motif score,ZBED4 motif score,ZBTB6 motif score,ZNF213 motif score,ZNF384 motif score,ZNF460 motif score,H3K27ac HN,H3K27me3 HN,H3K36me3 HN,H3K4me1 HN,H3K4me3 HN

Final: chrom,start,end,BJ_epi,BJ_iPSC,H9_ESC,iPSC,LIS49_hESC,MCF7_endo,MCF7_induc,U87,HN_SCC,Klf4 motif,Klf4 motif score,Number of Klf4 motifs,Strand,Total motif score,FOSL2 motif,FOXC2 motif,OTX2 motif,POU6F2 motif,TP53 motif,ZBED4 motif,ZBTB6 motif,ZNF213 motif,ZNF384 motif,ZNF460 motif,FOSL2 motif score,FOXC2 motif score,OTX2 motif score,POU6F2 motif score,TP53 motif score,ZBED4 motif score,ZBTB6 motif score,ZNF213 motif score,ZNF384 motif score,ZNF460 motif score,H3K27ac HN,H3K27me3 HN,H3K36me3 HN,H3K4me1 HN,H3K4me3 HN,H2AFZ,H3F3A,H3K27ac,H3K27me3,H3K36me3,H3K4me1,H3K4me2,H3K79me2,H3K9ac,H3K9me2,H3K9me3,H4K20me1,H3K4me3

Grouping of the variables:
- chrom, start, end
- Cell lines: BJ_epi, BJ_iPSC, H9_ESC, iPSC, LIS49_hESC, MCF7_endo, MCF7_induc, U87, HN_SCC
- KLF4 info: Klf4 motif, Klf4 motif score, Number of Klf4 motifs, Strand, Total motif score
- Comotifs present: FOSL2 motif, FOXC2 motif, OTX2 motif, POU6F2 motif, TP53 motif, ZBED4 motif, ZBTB6 motif, ZNF213 motif, ZNF384 motif, ZNF460 motif
- Comotifs scores: FOSL2 motif score, FOXC2 motif score, OTX2 motif score, POU6F2 motif score, TP53 motif score, ZBED4 motif score, ZBTB6 motif score, ZNF213 motif score, ZNF384 motif score, ZNF460 motif score
- MCF7 histone modifications: H2AFZ, H3F3A, H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K4me2, H3K4me3, H3K79me2, H3K9ac, H3K9me2, H3K9me3, H4K20me1
- HN SCC histone modifications: H3K27ac HN, H3K27me3 HN, H3K36me3 HN, H3K4me1 HN, H3K4me3 HN


### Progression of columns of the mouse master file: 
chrom  start   end    bone_marrow    Mbd3f  MEF    SCC    V6_5_ESC       CD19   ESC    pre_iPSC       Prostate_Stem
chrom,start,end,bone_marrow,Mbd3f_MEF,Mbd3f_iPSC,MEF,SCC,V6_5_ESC,CD19,ESC,pre_iPSC,Prostate_Stem,OSKM_MEF,Klf4 motif,Klf4 motif score,Number of Klf4 motifs,Strand,Total motif score
