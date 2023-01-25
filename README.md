# ag-soil-anaerobe

This repository houses the data and code for Lacroix et al. 2023: Anoxic microsites enhance soil carbon storage and respond to management. 

## DATA

The following files represent the raw ddPCR data:

* `intact_core_grav.csv`
* `dna_extraction.csv`
* `mcra_data.xlsx`
* `EML_16S_SHI_redo_08MAR2022_auto.csv`
* `EML_16S_SHI_redo_08MAR2022_auto.csv`
* `glta_SHI_28oct_manual.csv`
* `lacroix_glta_02nov_manual.csv`
* `EML_nirK_SHI_samples_16jan_auto.csv` 
* `EML_nirS_SHI_samples_21jan_MANUAL_6000.csv` 
* `EML_dsrAB_SHI_samples_15mar2022_auto.csv`

`AllData_03AUG.xlsx` contains all of the data from the study, with different tabs representing different datatypes. More information can be found in the "README" tab of this spreadsheet.

## CODE

For each `.Rmd` file there is an associated `.md` file that shows the code (with a clickable table of contents) and output formatted for web.

`ddPCR.Rmd` shows the raw data processing fo ddPCR data and contains all of the univariate statistical tests applied to the ddPCR data.

`PredictAnaerobeAbundance.Rmd` contains all of the code for the exercises which predicted anaerobe abundance.

`PredictSoilC.Rmd` shows the tests that predicted Soil C concentration. 






