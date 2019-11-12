# restoration-planting-regrowth-comparison
Code used to generate results, figures and tables for DOI: 10.1111/rec.13077 [![DOI](https://zenodo.org/badge/221116714.svg)](https://zenodo.org/badge/latestdoi/221116714)

This project uses **forest inventory data** and **plant functional trait data**, processing them in **R scripts** to generate results and figures.

**Forest inventory data**
The forest inventory data used to generate this publication were obtained from Greening Australia and Greenfleet. These data were obtained via data agreements and cannot be reproduced here. I have provided anonymised plot-level aggregate data, the minimum required to reproduce all results and figures used in this project. I have not provided any plant-level data, site coordinates or identifying information.

**Functional trait data**
Functional trait data was obtained from the TRY Plant Trait Database and Australian herbarium records (for maximum height only). Citations to trait databases and publications for TRY functional trait data are made in the Methods section of the paper. Raw traits from TRY are not reproduced here.

**R scripts** Two scripts are provided in this repository, along with summarised data used for for the statistical models and analyses in the manuscript. These scripts makes extensive use of the code folding functionality in R-studio (Alt + O is the default shortcut on Window and Linux machines to collapse all folds).

**reproduce-results.R** is a minimum working script that contains all of the data and analyses to produce the figures and tables found in the manuscript. It does not, in most cases, show how the raw data was transformed from the raw functional trait and forest inventory data. It makes use of files in the "Data" and "Functions" sub-folder, so please point your working directory to the location of this script and include all other directories in the repository.

**full-processing-and-analysis.R** shows the complete data analysis pathway, but will require the raw forest inventory and functional trait data acquired and aggregated by the authors of this publication. These are not provided here. This script is provided for transparency and reproducibility, that data processing steps are fully documented.
