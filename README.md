## A benchmark of batch-effect correction methods for single-cell RNA sequencing of cancer samples

This repository contains custom code used for analysis and plotting in Richards et al. (2021).

> Tumours are routinely profiled with single-cell RNA sequencing to characterize the diverse cellular ecosystems of malignant, immune, and stromal cell types.  When combining data from multiple samples or studies, it is crucial to correct for batch-specific variations. However, single cell batch integration methods are almost never designed for or benchmarked on datasets containing malignant cells. Here, we compare 5 data integration tools on 5 cancer single-cell RNA-seq datasets. Based on our results, STACAS is the most suitable method for integration of tumour datasets. However, due to its significantly shorter run time and ability to run without the need for parameter tuning, fastMNN is a competitive alternative. This benchmark provides a framework for evaluating how well single-cell integration methods correct for technical variability while preserving biological heterogeneity between cancer patients, in both malignant and microenvironmental cell compartments. 


###
###
###
##
## 1) Datasets

| DatasetID | Cancer Type  | Sequencing Technology | Input Material | Num. Samples  | Num. Patients | Num. Cells | Cell Types | Data Download |
|-----------|-------------|----------------------|---------------|---------------|---------------|------------|-----------|------|
| Richards-GBM-LGG | Glioblastoma,<br />Oligodendroglioma | 10x Genomics, 3' (v2) | Nuclei | 8 | 3 | 35,549 | Malignant, Astrocytes, Oligodendrocytes, Neurons, Tcells, Myeloid, Vascular cells | NA |
| Yost-BCC [1] | Basal cell carcinoma | 10x Genomics, 5' | Cells | 22 | 11 | 53,030 | B_cells, CAFs, T_cells, DCs, Endothelial, Macrophages, Malignant, Melanocytes, Myofibroblasts, NK_cells, pDCs, Plasma_cells | [GSE123813](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123813) |
| Ma-LIHC [2] | Hepatocellular carcinoma,<br />Intrahepatic cholangiocarcinoma | 10x Genomics, 3' (v2) | Cells | 19 | 19 | 9,752 | B_cells, CAFs, HPCs, Malignant, T_cells, Macrophages, Endothelial | [GSE125449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125449) |
| Caron-ALL [3] | Childhood acute lymphoblastic leukemia | 10x Genomics, 3' (v2) | Cells | 11 | 11 | 38,827 | Malignant, T_cells, Erythrocytes, B_cells, NK_cells, Macrophages | [GSE132509](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132509) |
| Bi-RCC [4] | Clear cell renal cell carcinoma | 10x Genomics, 3' (v2) | Cells | 8 | 8 | 34,048 | T_cells, Macrophages, NK_cells, Malignant, Endothelial, B_cells, DCs, Mast_cells, Plasma_cells, Fibroblast | [Broad Single Cell Portal: SCP1288](https://singlecell.broadinstitute.org/single_cell/study/SCP1288/tumor-and-immune-reprogramming-during-immunotherapy-in-advanced-renal-cell-carcinoma) |

###
###
###
##
## 2) Data Integration Tools

| Method | Reference | Description | Rank | RunTime | MemoryUsage |
|--------|-----------|-------------|------|---------|-------------|
|        |           |             |      |         |             |

###
###
###
##
## 3) References

1. Yost, K.E., Satpathy, A.T., Wells, D.K. et al. Clonal replacement of tumor-specific T cells following PD-1 blockade. Nat Med 25, 1251–1259 (2019). https://doi.org/10.1038/s41591-019-0522-3
2. Ma, L., Hernandez, M.O., Zhao, Y. et al. Tumor Cell Biodiversity Drives Microenvironmental Reprogramming in Liver Cancer. Cancer Cell 36, 418-430.e6 (2019). https://doi.org/10.1016/j.ccell.2019.08.007
3. Caron, M., St-Onge, P., Sontag, T. et al. Single-cell analysis of childhood leukemia reveals a link between developmental states and ribosomal protein expression as a source of intra-individual heterogeneity. Sci Rep 10, 8079 (2020). https://doi.org/10.1038/s41598-020-64929-x
4. Bi, K., He, M.X., Bakouny, Z. et al. Tumor and immune reprogramming during immunotherapy in advanced renal cell carcinoma. Cancer Cell 39, 649-661.e5 (2021). https://doi.org/10.1016/j.ccell.2021.02.015

###     
###      
###     
##     
## 4) Contact Information
Laura Richards (lauram.richards@mail.utoronto.ca)  
Dr. Trevor Pugh (trevor.pugh@utoronto.ca)  
