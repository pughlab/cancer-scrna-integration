## A comparison of data integration methods for single-cell RNA sequencing of cancer samples

*This repository contains custom code used for analysis and plotting in Richards et al., bioRxiv (2021).*

Tumours are routinely profiled with single-cell RNA sequencing (scRNA-seq) to characterize their diverse cellular ecosystems of malignant, immune, and stromal cell types.  When combining data from multiple samples or studies, batch-specific technical variation can confound biological signals. However, scRNA-seq batch integration methods are often not  designed for, or benchmarked, on datasets containing cancer cells. Here, we compare 5 data integration tools applied to 171,206 cells from 5 tumour scRNA-seq datasets. Based on our results, STACAS and fastMNN are the most suitable methods for integrating tumour datasets, demonstrating robust  batch effect correction while preserving relevant biological variability in the malignant compartment. This comparison provides a framework for evaluating how well single-cell integration methods correct for technical variability while preserving biological heterogeneity of malignant and non-malignant cell populations.

##
### Data Availability 
The majority of data used during this study was obtained from publicly available sources. All datasets have been made available in their re-processed form through CReSCENT (https://crescent.cloud/; Study IDs CRES-P24, CRES-P25, CRES-P26, CRES-P27, CRES-P28).
##
### Contact Information
Laura Richards (lauram.richards@mail.utoronto.ca)  
Dr. Trevor Pugh (trevor.pugh@utoronto.ca)  
