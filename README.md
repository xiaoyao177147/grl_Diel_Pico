# Code and data for: Antiphase Diel Patterns of Abundance and Biomass of Picophytoplankton in the Oligotrophic Ocean

## Root

* Field_SCS.R: R script for Fig. 1 (the field study in the northern SCS) and Table S3 of the paper;

* model.R: R script for Fig. 3 (A simple model simulating the daily patterns of cell numbers, cell size, and biomass of Prochlorococcus) of the paper;

* Table_S1.R: R script for Table S1 of the paper;

* Calibration_FSC.R: R script for Fig. S2 (Picophytoplankton culture-based calibration between cell diameter and FSC) of the paper;

* bac_virus.R: R script for Fig. S4 (Diel variation of biomass of heterotrophic bacteria and viruses) of the paper;

* rVolume.R: R script for Fig. S11 (Comparison of daily percent increase in biovolume in culture and SCS) of the paper;

* Four_carbon_conversion.R: R script for Fig. S12 (Diel variation of picophytoplankton in biomass derived from different volume-to-carbon conversion methods) of the paper;

## Data

* Picoauto.csv: Time-series measurements of picophytoplakton data;

* Sampletime.csv: Sampling time;

* FSC_size.csv: Size-FSC calibration data;

* Bac_virus.csv: Time-series measurements of bacteria and virus data;

* Chla.csv: Chla data;

* Dilution.csv: Modified dilution experiments data; 
 
## Figure

* Plots of apparent growth rate (9h^-1^) versus dilution factor in modified dilution experiments;

## SeaFlow

R scripts (4 steps) and data (SeaFlow data v1.3 http://doi.org/10.5281/zenodo.2678021 and SCS data) for Fig. 2 and Fig. S8--S10 of the paper;

### Data 

Raw data and data generated by a series of intermediate steps;

### Fig 

Partial plots generated in meta-analysis.

```
> sessionInfo()

R version 4.1.2 (2021-11-01)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.936  LC_CTYPE=Chinese (Simplified)_China.936   
[3] LC_MONETARY=Chinese (Simplified)_China.936 LC_NUMERIC=C                              
[5] LC_TIME=Chinese (Simplified)_China.936    

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggbeeswarm_0.6.0     ggh4x_0.2.1          car_3.0-12           carData_3.0-4       
 [5] forcats_0.5.1        stringr_1.4.0        dplyr_1.0.7          purrr_0.3.4         
 [9] readr_2.1.1          tidyr_1.1.4          tibble_3.1.6         tidyverse_1.3.1     
[13] lubridate_1.8.0      bpDir_0.1.2          circular_0.4-93      suncalc_0.5.0       
[17] data.table_1.14.2    ggpubr_0.4.0         scales_1.1.1         egg_0.4.5           
[21] ggplot2_3.3.5        gridExtra_2.3        lemon_0.4.5          patchwork_1.1.1     
[25] cowplot_1.1.1        gtable_0.3.0         plyr_1.8.6           investr_1.4.0       
[29] basicTrendline_2.0.5

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.7       mvtnorm_1.1-2    lattice_0.20-45  assertthat_0.2.1 utf8_1.2.2      
 [6] R6_2.5.1         cellranger_1.1.0 backports_1.2.1  reprex_2.0.1     httr_1.4.2      
[11] pillar_1.6.3     rlang_0.4.11     readxl_1.3.1     rstudioapi_0.13  munsell_0.5.0   
[16] broom_0.7.9      vipor_0.4.5      compiler_4.1.2   modelr_0.1.8     xfun_0.29       
[21] pkgconfig_2.0.3  tidyselect_1.1.1 fansi_0.5.0      crayon_1.4.1     tzdb_0.1.2      
[26] dbplyr_2.1.1     withr_2.4.3      MASS_7.3-54      nlme_3.1-153     jsonlite_1.7.2  
[31] lifecycle_1.0.1  DBI_1.1.1        magrittr_2.0.1   cli_3.1.0        stringi_1.7.6   
[36] ggsignif_0.6.3   fs_1.5.0         xml2_1.3.2       ellipsis_0.3.2   generics_0.1.0  
[41] vctrs_0.3.8      boot_1.3-28      tools_4.1.2      beeswarm_0.4.0   glue_1.4.2      
[46] hms_1.1.1        abind_1.4-5      plotrix_3.8-2    colorspace_2.0-2 rstatix_0.7.0   
[51] rvest_1.0.2      knitr_1.37       haven_2.4.3 
```
