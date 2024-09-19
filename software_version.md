
SUPPA version 2.3, installed via conda, full environment `conda list`:
```
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                       2_gnu    conda-forge
brotli                    1.0.9                h166bdaf_8    conda-forge
brotli-bin                1.0.9                h166bdaf_8    conda-forge
bzip2                     1.0.8                h7f98852_4    conda-forge
ca-certificates           2023.5.7             hbcca054_0    conda-forge
certifi                   2023.5.7           pyhd8ed1ab_0    conda-forge
charset-normalizer        3.1.0              pyhd8ed1ab_0    conda-forge
idna                      3.4                pyhd8ed1ab_0    conda-forge
joblib                    1.2.0              pyhd8ed1ab_0    conda-forge
ld_impl_linux-64          2.40                 h41732ed_0    conda-forge
libblas                   3.9.0           16_linux64_openblas    conda-forge
libbrotlicommon           1.0.9                h166bdaf_8    conda-forge
libbrotlidec              1.0.9                h166bdaf_8    conda-forge
libbrotlienc              1.0.9                h166bdaf_8    conda-forge
libcblas                  3.9.0           16_linux64_openblas    conda-forge
libexpat                  2.5.0                hcb278e6_1    conda-forge
libffi                    3.4.2                h7f98852_5    conda-forge
libgcc-ng                 12.2.0              h65d4601_19    conda-forge
libgfortran-ng            12.2.0              h69a702a_19    conda-forge
libgfortran5              12.2.0              h337968e_19    conda-forge
libgomp                   12.2.0              h65d4601_19    conda-forge
liblapack                 3.9.0           16_linux64_openblas    conda-forge
libnsl                    2.0.0                h7f98852_0    conda-forge
libopenblas               0.3.21          pthreads_h78a6416_3    conda-forge
libsqlite                 3.42.0               h2797004_0    conda-forge
libstdcxx-ng              12.2.0              h46fd767_19    conda-forge
libuuid                   2.38.1               h0b41bf4_0    conda-forge
libzlib                   1.2.13               h166bdaf_4    conda-forge
ncurses                   6.3                  h27087fc_1    conda-forge
numpy                     1.24.3          py311h64a7726_0    conda-forge
openssl                   3.1.0                hd590300_3    conda-forge
packaging                 23.1               pyhd8ed1ab_0    conda-forge
pandas                    2.0.1           py311h320fe9a_1    conda-forge
patsy                     0.5.3              pyhd8ed1ab_0    conda-forge
pip                       23.1.2             pyhd8ed1ab_0    conda-forge
platformdirs              3.5.1              pyhd8ed1ab_0    conda-forge
pooch                     1.7.0              pyha770c72_3    conda-forge
pysocks                   1.7.1              pyha2e5f31_6    conda-forge
python                    3.11.3          h2755cc3_0_cpython    conda-forge
python-dateutil           2.8.2              pyhd8ed1ab_0    conda-forge
python-tzdata             2023.3             pyhd8ed1ab_0    conda-forge
python_abi                3.11                    3_cp311    conda-forge
pytz                      2023.3             pyhd8ed1ab_0    conda-forge
readline                  8.2                  h8228510_1    conda-forge
requests                  2.31.0             pyhd8ed1ab_0    conda-forge
scikit-learn              1.2.2           py311h103fc68_1    conda-forge
scipy                     1.10.1          py311h64a7726_3    conda-forge
setuptools                67.7.2             pyhd8ed1ab_0    conda-forge
six                       1.16.0             pyh6c4a22f_0    conda-forge
statsmodels               0.14.0          py311h1f0f07a_1    conda-forge
suppa                     2.3                        py_2    bioconda
threadpoolctl             3.1.0              pyh8a188c0_0    conda-forge
tk                        8.6.12               h27826a3_0    conda-forge
typing-extensions         4.6.1                hd8ed1ab_0    conda-forge
typing_extensions         4.6.1              pyha770c72_0    conda-forge
tzdata                    2023c                h71feb2d_0    conda-forge
urllib3                   2.0.2              pyhd8ed1ab_0    conda-forge
wheel                     0.40.0             pyhd8ed1ab_0    conda-forge
xz                        5.2.6                h166bdaf_0    conda-forge
```


R used for `src/split_events.R` and `R/str_q_remove_outliers.R` version 4.2.0, sessionInfo:
```
R version 4.2.0 (2022-04-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux 8.8 (Ootpa)

Matrix products: default
BLAS/LAPACK: /vast/palmer/apps/avx2/software/OpenBLAS/0.3.12-GCC-10.2.0/lib/libopenblas_haswellp-r0.3.12.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] purrr_1.0.2   stringr_1.4.0 getopt_1.20.3

loaded via a namespace (and not attached):
[1] compiler_4.2.0  magrittr_2.0.3  cli_3.6.2       tools_4.2.0
[5] vctrs_0.6.5     stringi_1.7.6   lifecycle_1.0.4 rlang_1.1.3
```

R used for subsequent analyses in `R/psi.R` and `R/dpsi.R`:
```
R version 4.4.0 (2024-04-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 22631)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] wbData_0.9.2    lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4    
 [6] purrr_1.0.2     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1  
[11] tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] gld_2.6.6           gtable_0.3.5        beeswarm_0.4.0      ggrepel_0.9.5      
 [5] RApiSerialize_0.1.2 agricolae_1.3-7     lattice_0.22-6      tzdb_0.4.0         
 [9] vctrs_0.6.5         tools_4.4.0         generics_0.1.3      parallel_4.4.0     
[13] proxy_0.4-27        fansi_1.0.6         cluster_2.1.6       AlgDesign_1.2.1    
[17] pkgconfig_2.0.3     pheatmap_1.0.12     Matrix_1.7-0        data.table_1.15.4  
[21] polylabelr_0.2.0    RColorBrewer_1.1-3  ggridges_0.5.6      RcppParallel_5.1.7 
[25] readxl_1.4.3        lifecycle_1.0.4     rootSolve_1.8.2.4   farver_2.1.1       
[29] compiler_4.4.0      egg_0.4.5           Exact_3.3           munsell_0.5.1      
[33] qs_0.26.1           vipor_0.4.7         DescTools_0.99.55   class_7.3-22       
[37] crayon_1.5.2        pillar_1.9.0        MASS_7.3-60.2       boot_1.3-30        
[41] nlme_3.1-164        tidyselect_1.2.1    mvtnorm_1.2-4       stringi_1.8.3      
[45] labeling_0.4.3      polyclip_1.10-6     conover.test_1.1.6  grid_4.4.0         
[49] colorspace_2.1-0    lmom_3.0            expm_0.999-9        cli_3.6.2          
[53] magrittr_2.0.3      eulerr_7.0.2        utf8_1.2.4          e1071_1.7-14       
[57] withr_3.0.0         scales_1.3.0        rappdirs_0.3.3      bit64_4.0.5        
[61] ggbeeswarm_0.7.2    timechange_0.3.0    httr_1.4.7          matrixStats_1.3.0  
[65] bit_4.0.5           gridExtra_2.3       cellranger_1.1.0    hms_1.1.3          
[69] stringfish_0.16.0   rlang_1.1.3         Rcpp_1.0.12         printMat_0.1.0     
[73] glue_1.7.0          vroom_1.6.5         rstudioapi_0.16.0   R6_2.5.1 
```





