---
title: "hCBO merge data"
author: "Qiang (Danny) Chen"
date: "2025-09-19"
output:
  html_document:
    keep_md: true
    highlight: tango
    code_folding: hide
    number_sections: yes
    theme: united
    toc: yes
    toc_float:
      collapsed: yes
      smooth_scroll: no
  word_document:
    toc: yes
  pdf_document:
    toc: yes
  bibliography: biblio.bib
vignette: |
  %\VignetteIndexEntry{RNASeq} %\VignetteEngine{knitr::rmarkdown} %\VignetteEncoding{UTF-8} %\usepackage[utf8]{inputenc}
---

We integrated *GM25256, NCRM5, GM23913* at day60 in our 2 month-old hCBO data with the following two datasets and identified DEGs among cerebellar markers across protocols. 

For the quadrato lab, we focused on three samples below (they generated three replicates at day 60):  
[GSE247974](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247974)
GSM7904010: organoid D1, scRNA-seq  
GSM7904011: organoid D2, scRNA-seq  
GSM7904012: organoid D3, scRNA-seq  

For the pasca lab, we only focused on one sample (at day 72) below:  
[GSE233574](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233574)
GSM7430706: OrganoidScreen; sample11; FGF2-50  

There is only one sample (GSM7430706) in GSE233574. We did sub-sampling in cells of GSM7430706 and calculated psudo-bulk of sub-samples for comparisons. 



``` r
# knitr::opts_chunk$set(tidy=FALSE, cache=FALSE, echo=TRUE, dev="png", message=FALSE, error=FALSE, warning=FALSE)
knitr::opts_chunk$set(echo=TRUE, warning=F, message=F)

library("tidyr")
library("dplyr")
library("ggplot2")
library("knitr")
library("pheatmap")
library("DT")

library("Seurat")
library("DESeq2")

rm(list=ls())

d_geo  = "/data/qiangchen/Projects/NCATS/SCTL/SeungMiRyu/hCBO/GEO/"
d_sctl = "/data/qiangchen/Projects/NCATS/SCTL/SeungMiRyu/hCBO/ISB038/htseq_count/"

d_proj = "/data/qiangchen/Projects/NCATS/SCTL/SeungMiRyu/hCBO/SCTLday60_GSE247974_GSE233574/"
d_dat  = paste0(d_proj, "Data/")
d_res  = paste0(d_proj, "Results/")
```

# Load SCTL data


``` r
f_sctl = paste0(d_sctl, "ISB038.counts.Rdata")
load(f_sctl)

countData_SCTL = as.data.frame(df)

colnames = colnames(countData_SCTL)
colnames = gsub("_htseq_counts", "", colnames)
colnames = gsub("-", "_", colnames)
colnames(countData_SCTL) = colnames

cols_selected = c(1, grep("day60", colnames))

countData_SCTL = countData_SCTL[, cols_selected]

countTable_SCTL = countData_SCTL[, 2:ncol(countData_SCTL)]
print(nrow(countData_SCTL))
```

# Load GSE247974 data


``` r
samples = c("GSM7904010_D1", "GSM7904011_D2", "GSM7904012_D3")

bulkList = lapply(samples, function(cur_sample) {

    d_sample = paste0(d_geo, "GSE247974/Raw/", cur_sample, "/")
    cur_data = Read10X(data.dir = d_sample)
    cur_so  = CreateSeuratObject(counts=cur_data, project=cur_sample)

    cur_metadata = cur_so@meta.data
    cur_group = as.numeric(as.factor(cur_metadata$orig.ident))

    count_matrix = GetAssayData(cur_so, assay = "RNA", slot = "counts")

    pseudobulk = data.frame(t(rowsum(t(as.matrix(count_matrix)), group=cur_group)))
    colnames(pseudobulk) = cur_sample 
    
    return(pseudobulk)

})

countTable_GSE247974 = as.data.frame(do.call(cbind, bulkList))
colnames(countTable_GSE247974) = c("GSE247974_day60_1", "GSE247974_day60_2", "GSE247974_day60_3")
countData_GSE247974 = cbind(Gene=rownames(countTable_GSE247974), countTable_GSE247974)

print(nrow(countData_GSE247974))
```

# Load GSE233574 data

We used only sample11 in this study. Sample11 includes 29910 features and 1028 cells. To conduct DEG analysis, we generated 3 bulk-RNA data using subsmples of cells. Each subsample includes randomly selected 80% of total cells. 



``` r
f_GSE233574 = paste0(d_geo, "GSE233574/Raw/GSE233574_OrganoidScreen_processed_SeuratObject.rds")

so_GSE233574 = readRDS(f_GSE233574)

so_GSE233574 = subset(so_GSE233574, subset=sample=="sample11")
metadata_GSE233574 = so_GSE233574@meta.data

count_matrix = GetAssayData(so_GSE233574, assay = "RNA", slot = "counts")

n_total = ncol(count_matrix)

# Calculate 80% of indices
n_cells = floor(0.8 * n_total)

group = 1:n_cells

# Generate random subsample of indices
subsample_ind1 = sample(1:n_total, size = n_cells, replace = FALSE)
subsample_ind2 = sample(1:n_total, size = n_cells, replace = FALSE)
subsample_ind3 = sample(1:n_total, size = n_cells, replace = FALSE)

count1 = rowSums(count_matrix[, subsample_ind1])
count2 = rowSums(count_matrix[, subsample_ind2])
count3 = rowSums(count_matrix[, subsample_ind3])

countData_GSE233574 = data.frame(Sub1=count1, Sub2=count2, Sub3=count3)
countData_GSE233574 = cbind(rownames(countData_GSE233574), countData_GSE233574)
colnames(countData_GSE233574) = c("Gene", "GSE233574_day72_1", "GSE233574_day72_2", "GSE233574_day72_3")

print(nrow(countData_GSE233574))
```

# Merge counts


``` r
merged = merge(countData_SCTL, countData_GSE247974, by="Gene")
print(nrow(merged))

merged = merge(merged, countData_GSE233574, by="Gene")
print(nrow(merged))

count_matrix = merged[, 2:ncol(merged)]
rownames(count_matrix) = merged$Gene

f_counts = paste0(d_dat, "raw_counts.SCTL_GSE247974_GSE233574.tsv")
write.table(count_matrix, f_counts, sep="\t", row.names=F, quote=F)

pca_result = prcomp(t(count_matrix)) 
pca_scores = as.data.frame(pca_result$x)
pca_scores$Sample = rownames(pca_scores)
pca_scores$Group = pca_scores$Sample
pca_scores$Group = gsub("_1", "", pca_scores$Group)
pca_scores$Group = gsub("_2", "", pca_scores$Group)
pca_scores$Group = gsub("_3", "", pca_scores$Group)
pca_scores$Group = as.factor(pca_scores$Group)

p = ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group, label=Group)) +
    geom_point(size = 5) +
    geom_text() + 
    theme_minimal() +
    labs(title = "PCA plot for merged counts",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) + 
    theme(legend.position = "none")

print(p)
```

![](step01.merge_samples_files/figure-html/merge_count-1.png)<!-- -->

``` r
f_fig = paste0(d_res, "PCA_plot.merged.D60_D72.before_normalization.pdf")
ggsave(f_fig, p)
```

# Boxplot of raw counts.


``` r

boxplot(count_matrix,main="Raw counts", ylab = "read counts", las=2)
```

![](step01.merge_samples_files/figure-html/box_raw-1.png)<!-- -->

# QC

QC was done following the same approach in previous analysis for this project


``` r

sampleInfo = data.frame(ID = colnames(count_matrix))
sampleInfo = sampleInfo %>% separate(ID, into=c("Sample", "Day", "Replicate"), remove=F)
sampleInfo$Protocol = "SCTLday60"
sampleInfo$Protocol[sampleInfo$Sample=="GSE247974"] = "GSE247974"
sampleInfo$Protocol[sampleInfo$Sample=="GSE233574"] = "GSE233574"

sampleInfo$Protocol = factor(sampleInfo$Protocol, levels=c("SCTLday60", "GSE247974", "GSE233574"))

table(sampleInfo$Protocol, sampleInfo$Sample)
##            
##             GM23913 GM25256 GSE233574 GSE247974 NCRM
##   SCTLday60       3       3         0         0    3
##   GSE247974       0       0         0         3    0
##   GSE233574       0       0         3         0    0

f_sample = paste0(d_dat, "sample_info.SCTLday60_GSE247974_GSE233574.tsv")
write.table(sampleInfo, f_sample, sep="\t", row.names=F, quote=F)

q = which(rowSums(count_matrix) > 10)
count_matrix = count_matrix[q, sampleInfo$ID]

f_count = paste0(d_dat, "count_matrix.SCTLday60_GSE247974_GSE233574.qced.rda")
save(count_matrix, sampleInfo, file=f_count)

```

# Session Information


``` r

sessionInfo()
## R version 4.4.3 (2025-02-28)
## Platform: x86_64-pc-linux-gnu
## Running under: Rocky Linux 8.7 (Green Obsidian)
## 
## Matrix products: default
## BLAS/LAPACK: /usr/local/intel/2022.1.2.146/mkl/2022.0.2/lib/intel64/libmkl_rt.so.2;  LAPACK version 3.9.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: America/New_York
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] DESeq2_1.46.0               SummarizedExperiment_1.36.0
##  [3] Biobase_2.66.0              MatrixGenerics_1.18.1      
##  [5] matrixStats_1.5.0           GenomicRanges_1.58.0       
##  [7] GenomeInfoDb_1.42.3         IRanges_2.40.1             
##  [9] S4Vectors_0.44.0            BiocGenerics_0.52.0        
## [11] Seurat_5.2.1                SeuratObject_5.0.2         
## [13] sp_2.2-0                    DT_0.33                    
## [15] pheatmap_1.0.12             knitr_1.50                 
## [17] ggplot2_3.5.1               dplyr_1.1.4                
## [19] tidyr_1.3.1                
## 
## loaded via a namespace (and not attached):
##   [1] RColorBrewer_1.1-3      rstudioapi_0.17.1       jsonlite_2.0.0         
##   [4] magrittr_2.0.3          spatstat.utils_3.1-3    farver_2.1.2           
##   [7] rmarkdown_2.29          ragg_1.3.3              zlibbioc_1.52.0        
##  [10] vctrs_0.6.5             ROCR_1.0-11             spatstat.explore_3.4-2 
##  [13] S4Arrays_1.6.0          htmltools_0.5.8.1       SparseArray_1.6.2      
##  [16] sass_0.4.9              sctransform_0.4.1       parallelly_1.43.0      
##  [19] KernSmooth_2.23-26      bslib_0.9.0             htmlwidgets_1.6.4      
##  [22] ica_1.0-3               plyr_1.8.9              plotly_4.10.4          
##  [25] zoo_1.8-13              cachem_1.1.0            igraph_2.1.4           
##  [28] mime_0.13               lifecycle_1.0.4         pkgconfig_2.0.3        
##  [31] Matrix_1.7-2            R6_2.6.1                fastmap_1.2.0          
##  [34] GenomeInfoDbData_1.2.13 fitdistrplus_1.2-2      future_1.34.0          
##  [37] shiny_1.10.0            digest_0.6.37           colorspace_2.1-1       
##  [40] patchwork_1.3.0         tensor_1.5              RSpectra_0.16-2        
##  [43] irlba_2.3.5.1           textshaping_1.0.0       labeling_0.4.3         
##  [46] progressr_0.15.1        spatstat.sparse_3.1-0   httr_1.4.7             
##  [49] polyclip_1.10-7         abind_1.4-8             compiler_4.4.3         
##  [52] withr_3.0.2             BiocParallel_1.40.1     fastDummies_1.7.5      
##  [55] R.utils_2.13.0          MASS_7.3-65             DelayedArray_0.32.0    
##  [58] tools_4.4.3             lmtest_0.9-40           httpuv_1.6.15          
##  [61] future.apply_1.11.3     goftest_1.2-3           R.oo_1.27.0            
##  [64] glue_1.8.0              nlme_3.1-167            promises_1.3.2         
##  [67] grid_4.4.3              Rtsne_0.17              cluster_2.1.8          
##  [70] reshape2_1.4.4          generics_0.1.3          gtable_0.3.6           
##  [73] spatstat.data_3.1-6     R.methodsS3_1.8.2       data.table_1.17.0      
##  [76] XVector_0.46.0          spatstat.geom_3.3-6     RcppAnnoy_0.0.22       
##  [79] ggrepel_0.9.6           RANN_2.6.2              pillar_1.10.1          
##  [82] stringr_1.5.1           spam_2.11-1             RcppHNSW_0.6.0         
##  [85] later_1.4.1             splines_4.4.3           lattice_0.22-6         
##  [88] survival_3.8-3          deldir_2.0-4            tidyselect_1.2.1       
##  [91] locfit_1.5-9.12         miniUI_0.1.1.1          pbapply_1.7-2          
##  [94] gridExtra_2.3           scattermore_1.2         xfun_0.52              
##  [97] stringi_1.8.7           UCSC.utils_1.2.0        lazyeval_0.2.2         
## [100] yaml_2.3.10             evaluate_1.0.3          codetools_0.2-20       
## [103] tibble_3.2.1            cli_3.6.4               uwot_0.2.3             
## [106] systemfonts_1.2.1       xtable_1.8-4            reticulate_1.40.0      
## [109] munsell_0.5.1           jquerylib_0.1.4         Rcpp_1.0.14            
## [112] globals_0.16.3          spatstat.random_3.3-3   png_0.1-8              
## [115] spatstat.univar_3.1-2   parallel_4.4.3          dotCall64_1.2          
## [118] listenv_0.9.1           viridisLite_0.4.2       scales_1.3.0           
## [121] ggridges_0.5.6          crayon_1.5.3            purrr_1.0.4            
## [124] rlang_1.1.5             cowplot_1.1.3

rm(list=ls())

```

