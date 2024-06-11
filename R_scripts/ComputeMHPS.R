library(purrr)

merge_MHPS <- function(phps,dhps,hhps){
  phps = cbind.data.frame(as.data.frame(phps),"models"=rownames(phps))

  list_df = list(phps,dhps,hhps)
  all_layers <- Reduce(function(x,y) merge(x,y,all=F,by="models"), list_df)
  
  all_layers$MHPS <- rowSums(all_layers[,c("DHPS","PHPS","HHPS")]*c(1/3,1/3,1/3),na.rm = T)
  
  all_layers <- all_layers[order(all_layers$MHPS,decreasing = T),]
  
  return(all_layers)
}

compute_MHPS <- function(path_dir){
  # The PHPS (normalised) is found under phenotype_heatmap and the column "PHPS"
  rank_res_phenotype = readRDS(paste0(path_dir,"/MHPS_Pipeline/Data/MHPS_inputs/PHPS_ResultTable.RDS"))
  
  # HHPS (normalised) is found under the column "HHPS"
  rank_res_histology = readRDS(paste0(path_dir,"/MHPS_Pipeline/Data/MHPS_inputs/HHPS_Results.RDS"))
  
  # DHPS (normalised) is found under the column "DHPS"
  rank_res_transcriptomics <- list("AB_Metabolic"=readRDS(paste0(path_dir,"/MHPS_Pipeline/Data/MHPS_inputs/DHPS_AB_Metabolic.RDS")),"AC_Fibrosis"=readRDS(paste0(path_dir,"/MHPS_Pipeline/Data/MHPS_inputs/DHPS_AC_Fibrosis.RDS")))
  
  rank_type=c("AB_Metabolic",rkt="AC_Fibrosis")
  
  nhps_res_list = lapply(rank_type,function(rkt){merge_MHPS(phps=rank_res_phenotype[[rkt]]$phenotype_heatmap,
                                                            dhps=rank_res_transcriptomics[[rkt]],
                                                            hhps=rank_res_histology[[rkt]])
  })
  names(nhps_res_list) = rank_type
  saveRDS(nhps_res_list,paste0(path_dir,"/MHPS_Pipeline/Data/MHPS_ranking.RDS"))  
}






