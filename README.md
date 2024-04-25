# ag-soil-anaerobe

This repository houses the data and code for Lacroix et al. 2024: Microbial proxies for anoxic microsites vary with management and partially explain soil carbon concentration. 

## DATA

The following files represent the raw ddPCR data:

* `intact_core_grav.csv`
* `dna_extraction.csv`
* `mcra_data.xlsx`
* `glta_SHI_28oct_manual.csv`
* `lacroix_glta_02nov_manual.csv`
* `EML_nirK_SHI_samples_16jan_auto.csv` 
* `EML_nirS_SHI_samples_21jan_MANUAL_6000.csv` 
* `EML_dsrAB_SHI_samples_15mar2022_auto.csv`

`anaerobe_pca.csv` contains the first principal component values of anaerobe abundances, as calculated in `ddPCR.Rmd`. 

`AllData_FINAL.xlsx` contains all of the data from the study, with different tabs representing different datatypes. More information can be found in the "README" tab of this spreadsheet.

## CODE

For each `.Rmd` file there is an associated `.md` file that shows the code (with a clickable table of contents) and output formatted for web.

`ddPCR.Rmd` shows the raw data processing fo ddPCR data and contains all of the univariate statistical tests applied to the ddPCR data.

`PredictSoilC.Rmd` shows the tests that predicted Soil C concentration. 

`PredictSoilC-withPCAresult.Rmd`shows all of the tests performed in `PredictSoilC.Rmd` but uses PC1 of anaerobe abundance in place of absolute anaerobe abundance.






