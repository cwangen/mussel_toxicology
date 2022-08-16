# mussel_toxics  (as of 8/16/22)
A repo containing data and analysis of mussel toxicology data from the Puget Sound. The main repo contains a .gitignore file that has exceptions for the outputs and reports folders, as well as a license, and .Rproj file, and this README.md.

# Summary by Folder

## data
* Clean: 
  * Contains totals_all.csv, which is created by 01_data-wrangling.R in the R folder. The R script should be updated in order for this file to reflect the most recently received data.
* Raw: 
  * Contains raw Excel files as received from WDFW. 
## outputs
* Contains plots and maps created by code contained in the R folder.
## R
* 01_data-wrangling.R: Takes raw Excel file and performs data cleaning.
* 02_initial-summary-stats.R: Performs simple data analysis.
* 02a_maps.R: Creates multi-panel maps for all analytes. 
* 03a, 03b, 03c: Individual analyte analysis.
## reports
* Contains report .RMD files and supporting files.
