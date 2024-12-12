# **Weak meta-analytic evidence for the evolution of  body size, fecundity, and survival in changing thermal environments**

This repository contains the data and code to reproduce the results from Pottier et al. (2025). Weak meta-analytic evidence for the evolution of   body size, fecundity, and survival in changing thermal environments. *In prep* 

A rendered version of the code is available in R/data_processing_and_analysis.html, or this webpage: https://p-pottier.github.io/Exp_evol_temp/. This document contains headers and tabs to help navigate through the code. On the left hand side, you also have a toggle to change from dark to light theme.

Importantly, note that this repository does not contain model outputs (RData/ folder). This is because the model outputs were too large to be shared in Github, but please contact *p.pottier@unsw.edu.au* if you would like access to some of these files. Note also that his code requires large computational power and all statistical models ran on the computational cluster *Katana* supported by Research Technology Services at UNSW Sydney (https://research.unsw.edu.au/katana). The R code needed to run these models are provided in the R/ folder, and the resouces required to run each R file are provided in the pbs/ folder. These resources can be adapted for different supercomputers. 

If you need any assistance with the data and code, or identify a mistake, please feel free to contact Patrice Pottier at *p.pottier@unsw.edu.au* 

------------

Below is an overview of the different folders in this repository and their content. 

# **bibliographic_searches/**
*This folder contains all the bibliographic files used to conduct the systematic review* 

## **files_for_screening/**
*This folder contains bibliographic files that were assigned to each contributor for screening* 

## **final_searches/** 
*This folder contains the final bibliographic files compiled from the different databases* 

### **all_results_combined/**
* `all_results_evol_repro_temp.ris`: Bibliographic file containing all bibliographic records, before deduplication

### **all_results_combined_deduplicated/**
* `all_records_deduplicated_evol_repro_temp.ris`: Bibliographic file containing all bibliographic records, after deduplication in R 
* `all_records_deduplicated_in_Rayyan_evol_repro_temp.ris`: Bibliographic file containing all bibliographic records, after deduplication in R and in Rayyan

### **proquest/**
*This folder contains the bibliographic files compiled from ProQuest (Dissertations and Theses)* 

### **scopus/**
*This folder contains the bibliographic files compiled from Scopus* 

### **web_of_science/**
*This folder contains the bibliographic files compiled from Web of Science (core collection)* 

## **pilot/** 
* `all_records_deduplicated_evol_repro_temp.ris`: Bibliographic file from a naive pilot search, which was re-used to identify relevant key terms for the systematic review* 

---

# **data/**
* `processed_data.csv`: Processed data used for the analyses. To find the data before data processing, see *data_extraction/all_extracted_data.csv*. Steps of data curation and processing are provided in the code.

---

# **data_extraction/**
* `all_extracted_data.xlsx`: Data extracted from all collaborators. Note that this data was cleaned and processed in Excel, as well as in R to ensure consistency in extractions among contributors. Note that this file contains the **metadata** in one of the tabs. 
* `all_extracted_data.csv`: Same as above, but in .csv format

## **by_authors/** 
*This folder contains the data extracted by each of the contributors (see initials).*
*Note that all these data sheets were checked by the lead author (see "_checked" suffix), and then combined into a single file "all_extracted_data.xlsx".* 
*Additional checks and corrections of typos were done after the data was combined.*  

---

# **fig/**
*This folder contains output file figures from the code. The figures presented in the manuscript were edited in Powerpoint a posteriori* 

-- 

# **pbs/**
*This folder contains all resource files needed to run the statistical models in a high-performance computing environment*
*The name of each individual file corresponds to an individual R file in the **R/models/** folder*

---

# **R/**
*This folder contains all R files*

## **data_processing_and_analysis_files/**
*File data for the rendered html file*

## **models/**
*All R files used to run statistical models. The name of each file corresponds to the effect size used as the response variable (lnRR, lnVR, or lnCVR), along with the moderator variables tested (e.g., assay_temp_diff)*

* `data_processing_and_analysis.Rmd`: Rmd file summarising all the code used to process the data, calculate the effect sizes, run statistical models, analyse model outputs, and produce the figures.
* `data_processing_and_analysis.html`: Knitted version of *data_processing_and_analysis.Rmd*  
* `deduplication.Rmd`: Rmd file used to deduplicate bibliographic records
* `keyword_search_design.Rmd`: Rmd file used to find additional key terms for the systematic review search
* `splitting_bibliographic_files_for_screening.Rmd`: Rmd file used to divide bibliographic records and assign them to each contributor for screening. 

# **RData/**
*Folder containing all model outputs. This folder is empty in Github because the model outputs are too large to be share through Github. However, please feel free to contact Patrice Pottier to request access to these files. The output from each model is also displayed in the webpage and "data_processing_and_analysis.html"*
