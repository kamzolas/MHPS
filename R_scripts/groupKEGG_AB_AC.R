# Clean the data tables with human KEGG pathways from the two data sets
cleanTables_pathways <- function(kegg_epos,kegg_ucam_vcu){
  kegg_epos <- kegg_epos[,!grepl("vsCTRL",colnames(kegg_epos))]
  colnames(kegg_ucam_vcu) <- sub("_EPoS","-UCMSAN",colnames(kegg_ucam_vcu))
  
  ep_nafl <- Reduce(function(x,y) merge(x=x,y=y,by="X",all=T),list(kegg_ucam_vcu,kegg_epos))
  
  colnames(ep_nafl) <- sub("NAFL.F0","MildNAFLD",colnames(ep_nafl))
  colnames(ep_nafl) <- sub("NAFL.NASHF0","MildNAFLD",colnames(ep_nafl))
  colnames(ep_nafl) <- sub("NASHF12","ModNAFLD",colnames(ep_nafl))
  colnames(ep_nafl) <- sub("NASHF34","SevNAFLD",colnames(ep_nafl))
  colnames(ep_nafl) <- sub("_EPoS","-EPoS",colnames(ep_nafl))  
  return(ep_nafl)
}

mergeTables_pathways <- function(kegg_epos,kegg_ucam_vcu){
  kegg_epos <- kegg_epos %>% rename_all(~ifelse(. == "X", ., paste0(.,";epos")))
  kegg_ucam_vcu <- kegg_ucam_vcu %>% rename_all(~ifelse(. == "X", ., paste0(.,";ucamvcu")))
  
  keggs_human <-full_join(kegg_epos,kegg_ucam_vcu,by="X")

  colnames(keggs_human) <- sub("NAFL.F0","MildNAFLD",colnames(keggs_human))
  colnames(keggs_human) <- sub("NAFL.NASHF0","MildNAFLD",colnames(keggs_human))
  colnames(keggs_human) <- sub("NASHF12","ModNAFLD",colnames(keggs_human))
  colnames(keggs_human) <- sub("NASHF34","SevNAFLD",colnames(keggs_human))
  return(keggs_human)
}

# Divide into A, B and C.
getSignifPaths5 <- function(greshm,sigm){
  
  # Identify the position of the relevant columns
  idc1 <- grepl("vsCTRL",colnames(sigm))
  idc2 <- grepl(";ucamvcu",colnames(sigm)) & !idc1
  idc3 <- grepl(";epos",colnames(sigm)) & !idc1
  
  # A
  all.up.idc <- rowSums(sign(greshm)>0) == ncol(greshm) & sigm[,idc1] == "*" & rowSums(sigm[,idc2] == "*") > 0 & rowSums(sigm[,idc3] == "*") > 0
  all.dwn.idc <- rowSums(sign(greshm)<0) == ncol(greshm) & sigm[,idc1] == "*" & rowSums(sigm[,idc2] == "*") > 0 & rowSums(sigm[,idc3] == "*") > 0
  
  # B
  naflctrl.up.idc <- sigm[,idc1]=="*" & greshm[,idc1]>0 & rowSums(sigm[,idc2] == "*") == 0 & rowSums(sigm[,idc3] == "*") == 0
  naflctrl.dwn.idc <- sigm[,idc1]=="*" & greshm[,idc1]<0 & rowSums(sigm[,idc2] == "*") == 0 & rowSums(sigm[,idc3] == "*") == 0
  
  # B
  naflnash.up.idc <- rowSums(sign(greshm[,!idc1])>0) == (ncol(greshm)-1) & sigm[,idc1] !="*" & rowSums(sigm[,idc2] == "*") > 0 & rowSums(sigm[,idc3] == "*") > 0
  naflnash.dwn.idc <- rowSums(sign(greshm[,!idc1])<0) == (ncol(greshm)-1) & sigm[,idc1] !="*" & rowSums(sigm[,idc2] == "*") > 0 & rowSums(sigm[,idc3] == "*") > 0
  
  # Wrap up
  idcList <- list(list(all.up.idc,all.dwn.idc),list(naflctrl.up.idc,naflctrl.dwn.idc),list(naflnash.up.idc,naflnash.dwn.idc))
  
  pathList <- lapply(1:length(idcList), function(x) {  
    idc.up <- idcList[[x]][[1]]
    idc.dwn <- idcList[[x]][[2]]
    if(sum(idc.up)>0){
      df.up <- data.frame("ID"=paste0("KEGG_",rownames(greshm[idc.up,])),"NES"=1)  
    } else {df.up <- data.frame("ID"=NA,"NES"=1)}
    if(sum(idc.dwn)>0){
      df.dwn <- data.frame("ID"=paste0("KEGG_",rownames(greshm[idc.dwn,])),"NES"=-1)  
    } else {df.dwn <- data.frame("ID"=NA,"NES"=-1)}
    df.all <- rbind(df.up,df.dwn)
    df.all <- df.all[!is.na(df.all$ID),]
    df.all$pvalue <- 0
    return(df.all)})
  names(pathList) <- c("a.All","b.NAFLvsCTRL","c_strict.NASHvsNAFL")
  
  return(pathList)
}

# Function to extract the relevant comparisons
KEGGres_FilteredTables <- function(cohortList,colkeep=NULL,pvalcutoff,path_dir){
  livrel <- read.csv2(paste0(path_dir,"/MHPS_Pipeline/Data/Utils/Final Liver-related pathways.csv"),header=T,sep=",",stringsAsFactors =F)
  mergdf <- Reduce(function(x,y) merge(x=x,y=y,by="X",all=T),cohortList)
  mergdf <- as.data.frame(mergdf,stringsAsFactors=T)
  rownames(mergdf) <- mergdf$X
  mergdf$X <- NULL
  mergdf <- mergdf[rownames(mergdf) %in% livrel$Pathways,]
  
  rownames(mergdf) <- sub("KEGG_","",rownames(mergdf))
  nrow(mergdf);mergdf <- mergdf[rowSums(t(apply(mergdf,1,is.na)))==0,];nrow(mergdf)
  
  # Keep only selected comparisons
  if(!is.null(colkeep)){mergdf <- mergdf[,colkeep]}
  
  # Keep only NES and p-value
  greshm <- mergdf[,grep("NES",colnames(mergdf))]
  sigm <- mergdf[,grep("pval",colnames(mergdf))]
  
  # Filter pathways with no significance
  nrow(greshm);greshm <- greshm[rowSums(sigm <= pvalcutoff)>0,];nrow(greshm)
  sigm <- sigm[rowSums(sigm <= pvalcutoff)>0,]
  sigm[sigm<=0.05] <- "*"
  sigm[sigm>0.05] <- ""
  sigm[is.na(sigm)] <- ""
  
  # Prepare output matrix
  for(i in 1:ncol(greshm)){greshm[,i] <- as.numeric(greshm[,i])}
  colnames(greshm) <- sub("NES_","",colnames(greshm))
  colnames(greshm) <- sub("_Human","",colnames(greshm))
  return(list(greshm,sigm))
}

# Function to create the human KEGG module sets divided into AB and AC
divideKEGGs_AB_AC <- function(sanucm,pvalcutoff,path_dir){
  selcols=c("vsCTRL","vsMildNAFLD")
  
  # GSEA result from the selected cohort(s)
  if(class(sanucm) != "list"){cohortList <- list(sanucm)}
  
  # Define the comparisons to include and create the module sets
  colkeep <- c()
  for(i in 1:length(selcols)){colkeep <- c(colkeep,grep(selcols[i],colnames(sanucm)[-1]))}
  
  greshm <- KEGGres_FilteredTables(cohortList,colkeep,pvalcutoff,path_dir)[[1]]
  sigm <- KEGGres_FilteredTables(cohortList,colkeep,pvalcutoff,path_dir)[[2]]
  
  pathList <- getSignifPaths5(greshm,sigm)
  
  pathList_upd <- lapply(1:length(pathList),function(x) {list("PW"=pathList[[x]])})
  names(pathList_upd) <- names(pathList)
  
  hs_modSet <- pathList_upd
  hs_modSet$b.NAFLvsCTRL$PW <- rbind.data.frame(pathList_upd$a.All$PW,pathList_upd$b.NAFLvsCTRL$PW)
  hs_modSet$c_strict.NASHvsNAFL$PW <- rbind.data.frame(pathList_upd$a.All$PW,pathList_upd$c_strict.NASHvsNAFL$PW)
  hs_modSet$a.All <- NULL
  
  names(hs_modSet) <- c("AB","AC")
  
  # Save the module set as input for DSEA
  saveRDS(hs_modSet,paste0(path_dir,"MHPS_Pipeline/Data/DHPS_inputs/HumanKEGG_AB_AC.RDS"))
}