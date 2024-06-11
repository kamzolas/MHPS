prep_HHPS <- function(path_dir,tbl_hhps){
  tbl_hhps = tbl_hhps %>% rename("Model.Name"="models")
  
  HHPS_list = list("AB_Metabolic"=tbl_hhps[,c("models","HHPS.Metabolic.relevance")] %>% rename("HHPS.Metabolic.relevance"="HHPS"),
                   "AC_Fibrosis"=tbl_hhps[,c("models","HHPS.Ability.to.induce.MASH.fibrosis")] %>% rename("HHPS.Ability.to.induce.MASH.fibrosis"="HHPS"))
  
  saveRDS(HHPS_list,paste0(path_dir,"/MHPS_Pipeline/Data/MHPS_inputs/HHPS_Results.RDS"))
}