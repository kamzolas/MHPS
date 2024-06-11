library(repo)
library(gep2pep)

prepare_ranks <- function(path_dir){
  alldat <- readRDS(paste0(path_dir,"/MHPS_Pipeline/Data/DHPS_inputs/RodentRankedTable.RDS"))
  
  # Remove the dots in the model names
  colnames(alldat) <- gsub("\\.","-",colnames(alldat))
  colnames(alldat) <- gsub("6J-WD-C0-3-STZ-9W-WD-5W","6J-WD-C0.3-STZ-9W-WD-5W",colnames(alldat))
  colnames(alldat) <- gsub("6J-WD-C1-3-FG-CCL4-24W","6J-WD-C1.3-FG-CCL4-24W",colnames(alldat))
  colnames(alldat) <- gsub("6J-MCD-5-5W","6J-MCD-5.5W",colnames(alldat))
  colnames(alldat) <- gsub("C0-2","C0.2",colnames(alldat))
  
  # Separate the models into microarray, RNA-seq rat and RNA-seq mouse
  mmu_mic <- alldat[,grep("6N.WD.C0.2",colnames(alldat))]
  rat <- alldat[,c(grep("RZ",colnames(alldat)),grep("R.CDA",colnames(alldat)))]
  mmu <- alldat[,!colnames(alldat) %in% c(colnames(mmu_mic),colnames(rat))]
  
  mmu_mic <- mmu_mic[rowSums(t(apply(mmu_mic,1,is.na)))==0,];nrow(mmu_mic)
  mmu_mic <- apply(mmu_mic,2,rank,ties.method="random")
  rat <- rat[rowSums(t(apply(rat,1,is.na)))==0,];nrow(rat)
  rat <- apply(rat,2,rank,ties.method="random")
  mmu <- mmu[rowSums(t(apply(mmu,1,is.na)))==0,];nrow(mmu)
  mmu <- apply(mmu,2,rank,ties.method="random")
  
  return(list(mmu_mic=mmu_mic,rat=rat,mmu=mmu))
}

# Open the repositories storing the complete set of KEGG pathways and calculate the ES for all the animal models
build_repoKEGG <- function(path_dir){
  rodents_list <- prepare_ranks(path_dir)
  mmu_mic <- rodents_list$mmu_mic
  rat <- rodents_list$rat
  mmu <- rodents_list$mmu
  
  rp_path_mmu_mic <- file.path(paste0(path_dir,"/MHPS_Pipeline/Data/DHPS_inputs/Repos/gep2pep_KEGG_AB_AC_Microarray_mmu/"))
  rp_mmu_mic <- repo_open(rp_path_mmu_mic)
  
  buildPEPs(rp_mmu_mic,mmu_mic,replace_existing = T,min_size = 1)
  
  rp_path_rat <- file.path(paste0(path_dir,"/MHPS_Pipeline/Data/DHPS_inputs/Repos/gep2pep_KEGG_AB_AC_RNAseq_rat/"))
  rp_rat <- repo_open(rp_path_rat)
  
  buildPEPs(rp_rat,rat,replace_existing = T,min_size = 1)
  
  rp_path_mmu <- file.path(paste0(path_dir,"/MHPS_Pipeline/Data/DHPS_inputs/Repos/gep2pep_KEGG_AB_AC_RNAseq_mmu/"))
  rp_mmu <- repo_open(rp_path_mmu)
  
  buildPEPs(rp_mmu,mmu,replace_existing = T,min_size = 1)  
}

# Open the repositories storing the gene set of the human DEGs and calculate the ES for all the animal models
run_build_DEG <- function(dbDir,dbName,grp_name){
  rodents_list <- prepare_ranks(path_dir)
  mmu_mic <- rodents_list$mmu_mic
  rat <- rodents_list$rat
  mmu <- rodents_list$mmu
  
  dbName <- paste0(dbDir,paste0(dbName,"_"))
  for(i in 1:length(grp_name)){
    rp_path_mmu_mic <- file.path(paste0(dbName[1],grp_name[i]))
    rp_mmu_mic <- repo_open(rp_path_mmu_mic)
    
    buildPEPs(rp_mmu_mic,mmu_mic,replace_existing = T,min_size = 5, max_size = 3000)
    
    rp_path_rat <- file.path(paste0(dbName[2],grp_name[i]))
    rp_rat <- repo_open(rp_path_rat)
    
    buildPEPs(rp_rat,rat,replace_existing = T,min_size = 5, max_size = 3000)
    
    rp_path_mmu <- file.path(paste0(dbName[3],grp_name[i]))
    rp_mmu <- repo_open(rp_path_mmu)
    
    buildPEPs(rp_mmu,mmu,replace_existing = T,min_size = 5, max_size = 3000) 
  }  
}

build_repoDEGs <- function(path_dir){
  grp_name_list <- list(
    c("AB","AC"))
  
  db_nm_list <- list(
    c("gep2pep_DEGs_Microarray_mmu","gep2pep_DEGs_RNAseq_rat","gep2pep_DEGs_RNAseq_mmu"))
  
  for(i in 1:length(grp_name_list)){
    run_build_DEG(dbDir=paste0(path_dir,"/MHPS_Pipeline/Data/DHPS_inputs/Repos/"),dbName=db_nm_list[[i]],grp_name=grp_name_list[[i]])  
  } 
}

