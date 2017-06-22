This is the Shiny app created in the Hackathon to demo new features to potentially merge with the old Shiny app (`Proteomic_Correlation_Shiny/DepLab/inst/shiny/`)

## Files:
- `app.R`: The main Shiny app, run this from RStudio or with `runApp()` from R in this dir
- `demo-data.R`: Script that sets up the data and plots to be used in the app
- `user_settings_io.R`: Shiny module for the user's RDS file settings import/export

## Software

### Dependencies

Some of the most crucial dependencies at the moment (there are likely to be more):

```
library(DepLabData)
library(DepLab)
library(ggplot2)
library(shiny)
library(plotly)
```

Currently tested under the following configuration:

```
R version 3.3.2 (2016-10-31)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: macOS Sierra 10.12.4

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] bindrcpp_0.2         NMF_0.20.6           Biobase_2.34.0       BiocGenerics_0.20.0  cluster_2.0.6       
 [6] rngtools_1.2.4       pkgmaker_0.22        registry_0.3         DepLabData_0.1.0     DepLab_0.1.1.9000   
[11] RSQLite_2.0          dplyr_0.7.0          DBI_0.7              data.table_1.10.4    shinyHeatmaply_0.1.0
[16] heatmaply_0.10.1     viridis_0.4.0        viridisLite_0.2.0    plotly_4.6.0         ggplot2_2.2.1       
[21] shiny_1.0.3         

loaded via a namespace (and not attached):
 [1] httr_1.2.1         tidyr_0.6.3        bit64_0.9-7        jsonlite_1.5       foreach_1.4.3     
 [6] gtools_3.5.0       assertthat_0.2.0   stats4_3.3.2       blob_1.1.0         yaml_2.1.14       
[11] robustbase_0.92-7  lattice_0.20-35    glue_1.1.0         digest_0.6.12      RColorBrewer_1.1-2
[16] colorspace_1.3-2   htmltools_0.3.6    httpuv_1.3.3       plyr_1.8.4         purrr_0.2.2.2     
[21] xtable_1.8-2       mvtnorm_1.0-6      scales_0.4.1       gdata_2.18.0       whisker_0.3-2     
[26] tibble_1.3.3       nnet_7.3-12        lazyeval_0.2.0     magrittr_1.5       mime_0.5          
[31] mclust_5.3         memoise_1.1.0      doParallel_1.0.10  MASS_7.3-47        gplots_3.0.1      
[36] class_7.3-14       tools_3.3.2        gridBase_0.4-7     trimcluster_0.1-2  stringr_1.2.0     
[41] kernlab_0.9-25     munsell_0.4.3      fpc_2.1-10         caTools_1.17.1     rlang_0.1.1       
[46] grid_3.3.2         iterators_1.0.8    htmlwidgets_0.8    crosstalk_1.0.0    labeling_0.3      
[51] bitops_1.0-6       gtable_0.2.0       codetools_0.2-15   flexmix_2.3-14     reshape2_1.4.2    
[56] TSP_1.1-5          R6_2.2.2           seriation_1.2-2    gridExtra_2.2.1    prabclus_2.2-6    
[61] bit_1.1-12         bindr_0.1          KernSmooth_2.23-15 dendextend_1.5.2   stringi_1.1.5     
[66] modeltools_0.2-21  Rcpp_0.12.11       gclus_1.3.1        DEoptimR_1.0-8     diptest_0.75-7    
```
