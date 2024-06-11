# load libraries ---------------------------------------------------------------------------

library(openxlsx)
library(dplyr)


# Define working directory -----------------------------------------------------------------

path_dir = ""


# Prepare HHPS -----------------------------------------------------------------------------

source(paste0(path_dir,"/MHPS_Pipeline/R_scripts/cleanTables_histology.R"))

tbls1b_hhps = read.xlsx(paste0(path_dir,"/MHPS_Pipeline/Data/Rodent_histology/TableS1B_HHPS.xlsx"))

prep_HHPS(path_dir,tbls1b_hhps)


# groupGENES_AB_AC.Rs -----------------------------------------------------------------------

source(paste0(path_dir,"/MHPS_Pipeline/R_scripts/groupGENES_AB_AC.R"))

degs_epos = read.csv(paste0(path_dir,"MHPS_Pipeline/Data/Human_DEGs_KEGG/Degs_EPoS.csv"))
degs_ucam_vcu = read.csv(paste0(path_dir,"MHPS_Pipeline/Data/Human_DEGs_KEGG/Degs_UCAM_VCU.csv"))

degs <- mergeTables_degs(degs_epos,degs_ucam_vcu)
divideDEGs_AB_AC(degs,effectSizeCol="log2FoldChange_",sigCol="pvalue_",log2_effectSizeCutoff=log2(1.5),path_dir=path_dir)


# groupKEGG_AB_AC.R -------------------------------------------------------------------------

source(paste0(path_dir,"/MHPS_Pipeline/R_scripts/groupKEGG_AB_AC.R"))

kegg_epos <- read.csv2(paste0(path_dir,"/MHPS_Pipeline/Data/Human_DEGs_KEGG/FGSEA_humanEPOSvsNAFL_merged.csv"),header=T,sep=",",stringsAsFactors = F)
kegg_ucam_vcu <- read.csv2(paste0(path_dir,"/MHPS_Pipeline/Data/Human_DEGs_KEGG/FGSEA_humanUCAM_VCUvsNAFL_merged.csv"),header=T,sep=",",stringsAsFactors = F)

keggs <- mergeTables_pathways(kegg_epos,kegg_ucam_vcu)
divideKEGGs_AB_AC(keggs,pvalcutoff=0.05,path_dir=path_dir)


# Create_PathwayRepository_KEGG_GENES.R ------------------------------------------------------

source(paste0(path_dir,"/MHPS_Pipeline/R_scripts/Create_Repository_KEGG_GENES.R"))

# This steps takes some minuttes
get_kegg_pathways(path_dir,kegg_filename="keggpathways.RDS")

keggpathways = readRDS(paste0(path_dir,"/MHPS_Pipeline/Data/DHPS_inputs/keggpathways.RDS"))

create_KEGG_repos(keggpathways,path_dir=path_dir)

mk_DEGs_Repo(modSetFiles=paste0(path_dir,"/MHPS_Pipeline/Data/DHPS_inputs/HumanGENES_AB_AC.RDS"),repoDir=paste0(path_dir,"/MHPS_Pipeline/Data/DHPS_inputs/Repos/"))


# Compute_RodentRegulationTable.R -------------------------------------------------------------

source(paste0(path_dir,"/MHPS_Pipeline/R_scripts/Compute_RodentRegulationTable.R"))

mmu_mic_de <- read.csv(paste0(path_dir,"/MHPS_Pipeline/Data/Rodent_DEGs/Microarrays_D.E.Analysis_results.csv"))
mmu_de <- read.csv(paste0(path_dir,"/MHPS_Pipeline/Data/Rodent_DEGs/Mouse_DEGs.csv"))
rt <- read.csv(paste0(path_dir,"/MHPS_Pipeline/Data/Rodent_DEGs/Rat_DEGs.csv"))

# Filter intervention models
mmu_de <- mmu_de %>% dplyr::select(-contains("6J.AMLN.C2.28W.AMLN.C2.CR.8W"),-contains("6J.CDAHFD.F45.4W.REV.8W"))

fc_rodents <- cleanTables_rodents(mmu_mic_de,mmu_de,rt,path_dir=path_dir)
fc_to_ranks(fc_rodents,path_dir=path_dir)


# buildPEPs_Repository_KEGG_GENES.R -----------------------------------------------------------

source(paste0(path_dir,"/MHPS_Pipeline/R_scripts/buildPEPs_Repository_KEGG_GENES.R"))

build_repoKEGG(path_dir=path_dir)
build_repoDEGs(path_dir=path_dir)


# ComputeDHPS.R -------------------------------------------------------------------------------

source(paste0(path_dir,"/MHPS_Pipeline/R_scripts/ComputeDHPS.R"))

compute_DHPS(path_dir=path_dir)

paste0(path_dir,"/MHPS_Pipeline/Data/MHPS_inputs/DHPS_AB_Metabolic")
dhps_ab = readRDS("/novo/users/lmhd/Projects/LITMUS/Analyses/D6_1_paper/Rscripts/litmus-t6.1-manuscript/Pipeline_published_version//MHPS_Pipeline/Data/MHPS_inputs/DHPS_AB_Metabolic.RDS")


# ComputeROC_ALT_AST.R -------------------------------------------------------------------------

source(paste0(path_dir,"/MHPS_Pipeline/R_scripts/ComputeROC_ALT_AST.R"))

phdb <- read.xlsx(paste0(path_dir,"/MHPS_Pipeline/Data/Rodent_phenotype/Table S9 - METADATA per animal.xlsx"))
# Filter intervention models
phdb <- phdb %>% filter(!phdb$PublicationModelName %in% c("6J-AMLN-C2-28W-AMLN-C2-CR-8W","6J-CDAHFD-F45-4W-REV-8W"))

phdb_upd <- cleanTables_createDiseaseGroups(phdb)

tresholds_alt_ast(path_dir=path_dir,phdb_upd)

# ComputePHPS.R -------------------------------------------------------------------------------

source(paste0(path_dir,"/MHPS_Pipeline/R_scripts/ComputePHPS.R"))

phdb <- read.xlsx(paste0(path_dir,"/MHPS_Pipeline/Data/Rodent_phenotype/Table S9 - METADATA per animal.xlsx"))
# Filter intervention models
phdb <- phdb %>% filter(!phdb$PublicationModelName %in% c("6J-AMLN-C2-28W-AMLN-C2-CR-8W","6J-CDAHFD-F45-4W-REV-8W"))

compute_PHPS(path_dir=path_dir,phdb)


# ComputeMHPS.R -------------------------------------------------------------------------------

source(paste0(path_dir,"/MHPS_Pipeline/R_scripts/ComputeMHPS.R"))

compute_MHPS(path_dir=path_dir)





