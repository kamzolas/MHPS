# MHPS (MASH Human Proximity Score)
Ranking of commonly used rodent MASLD models, based on their proximity to human MASLD according to their phenotypic, molecular, and histological profiles. 

This pipeline supports the analyses published in our Nature Metabolism paper (https://doi.org/10.1038/s42255-024-01043-6)

Running the MHPS pipeline
--------------------------------------------
Copy the MHPS pipeline folder to a local directory. 
To execute the MHPS pipeline use the script "RUN_MHPS.R" which will source and load all needed functions and input data from the MHPS pipeline folders.


Input files for MHPS pipeline (already included in MHPS pipeline folder)
--------------------------------------------
TableS1B_HHPS.xlsx
Copied columns with model names and HHPS scores of "Metabolic relevance" and "Ability to induce MASH fibrosis" from Supplementary Table S1B

BioStudies: https://www.ebi.ac.uk/biostudies/studies/S-BSST1361
Degs_EPoS.csv, Degs_UCAM_VCU.csv, FGSEA_humanEPOSvsNAFL_merged.csv, FGSEA_humanUCAM_VCUvsNAFL_merged.csv

Microarrays_D.E.Analysis_results.csv, Mouse_DEGs.csv, Rat_DEGs_with_Log2CPMs_AfterFilteringLowLog2CPMs.csv

Table S9 - METADATA per animal.xlsx

Utils:
MouseHuman_orthologyMapping.txt, rat_to_human.txt, Final Liver-related pathways.csv

Output
--------------------------------------------
List object "MHPS_ranking.RDS" including the two types of ranking "Metabolic relevance" and "Ability to induce MASH fibrosis"
