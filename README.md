# Data and scripts to analyze eukaryotic plankton communities differences from eight reefs in Bocas del Toro

# Data files
1) licor_data.csv == Light logger data
2) minmeanmax_tempdata.csv == Temperature logger data
3) orig_raw_data_total.mapping.txt == ASV sample information

Get the sequencing files at: https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA507270

# Scripts
1) licorData_code.R == Light data analysis
2) tempData_code.R == Temperature analysis
3) Complete_analysis.R == Plankton community genetic analysis, including analysis using dada2 (analysis through taxonomy assignment), MCMC.OTU (ASV filtering), vegan (PCoA), and DESeq2 (differential abundance testing)

# sessionInfo()
R version 3.5.3 (2019-03-11)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.5

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] DESeq2_1.22.2               pheatmap_1.0.12             cluster_2.1.0               MCMC.OTU_1.0.10            
 [5] MCMCglmm_2.29               ape_5.3                     coda_0.19-3                 Matrix_1.2-18              
 [9] vegan_2.5-6                 lattice_0.20-38             permute_0.9-5               forcats_0.4.0              
[13] stringr_1.4.0               dplyr_0.8.3                 purrr_0.3.3                 readr_1.3.1                
[17] tidyr_1.0.0                 tibble_2.1.3                ggplot2_3.2.1               tidyverse_1.3.0            
[21] phyloseq_1.24.2             ShortRead_1.40.0            GenomicAlignments_1.18.1    SummarizedExperiment_1.12.0
[25] DelayedArray_0.8.0          matrixStats_0.55.0          Biobase_2.42.0              Rsamtools_1.34.1           
[29] GenomicRanges_1.34.0        GenomeInfoDb_1.18.2         Biostrings_2.50.2           XVector_0.22.0             
[33] IRanges_2.16.0              S4Vectors_0.20.1            BiocParallel_1.16.6         BiocGenerics_0.28.0        
[37] dada2_1.10.1                Rcpp_1.0.3                  labdsv_2.0-1                mgcv_1.8-31                
[41] nlme_3.1-143               

loaded via a namespace (and not attached):
 [1] cubature_2.0.4         Rtsne_0.15             colorspace_1.4-1       hwriter_1.3.2         
 [5] htmlTable_1.13.3       corpcor_1.6.9          base64enc_0.1-3        fs_1.3.1              
 [9] rstudioapi_0.10        bit64_0.9-7            AnnotationDbi_1.44.0   fansi_0.4.1           
[13] lubridate_1.7.4        xml2_1.2.2             codetools_0.2-16       splines_3.5.3         
[17] geneplotter_1.60.0     knitr_1.27             zeallot_0.1.0          ade4_1.7-13           
[21] Formula_1.2-3          jsonlite_1.6           annotate_1.60.1        broom_0.5.3           
[25] dbplyr_1.4.2           compiler_3.5.3         httr_1.4.1             backports_1.1.5       
[29] assertthat_0.2.1       lazyeval_0.2.2         cli_2.0.1              htmltools_0.4.0       
[33] acepack_1.4.1          tools_3.5.3            igraph_1.2.4.2         gtable_0.3.0          
[37] glue_1.3.1             GenomeInfoDbData_1.2.0 reshape2_1.4.3         cellranger_1.1.0      
[41] vctrs_0.2.1            multtest_2.36.0        iterators_1.0.12       tensorA_0.36.1        
[45] xfun_0.12              rvest_0.3.5            lifecycle_0.1.0        XML_3.99-0.1          
[49] zlibbioc_1.28.0        MASS_7.3-51.5          scales_1.1.0           hms_0.5.3             
[53] biomformat_1.8.0       rhdf5_2.24.0           RColorBrewer_1.1-2     memoise_1.1.0         
[57] gridExtra_2.3          rpart_4.1-15           RSQLite_2.2.0          latticeExtra_0.6-28   
[61] stringi_1.4.5          genefilter_1.64.0      foreach_1.4.7          checkmate_1.9.4       
[65] rlang_0.4.2            pkgconfig_2.0.3        bitops_1.0-6           Rhdf5lib_1.2.1        
[69] htmlwidgets_1.5.1      bit_1.1-15.1           tidyselect_0.2.5       plyr_1.8.5            
[73] magrittr_1.5           R6_2.4.1               generics_0.0.2         Hmisc_4.3-0           
[77] DBI_1.1.0              foreign_0.8-74         pillar_1.4.3           haven_2.2.0           
[81] withr_2.1.2            nnet_7.3-12            survival_3.1-8         RCurl_1.95-4.13       
[85] modelr_0.1.5           crayon_1.3.4           locfit_1.5-9.1         grid_3.5.3            
[89] readxl_1.3.1           data.table_1.12.8      blob_1.2.0             digest_0.6.23         
[93] reprex_0.3.0           xtable_1.8-4           RcppParallel_4.4.4     munsell_0.5.0   
