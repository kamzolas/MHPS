# Functions to clean data tables and divide DEGs into group A, B and C

DEG_FilteredTables <- function(de_genes,effectSizeCol,sigCol,log2_effectSizeCutoff=NULL){
  mergdf <- de_genes

  # Keep only the values needed for the grouping
  greshm <- mergdf[,grep(effectSizeCol,colnames(mergdf))]
  sigm <- mergdf[,grep(sigCol,colnames(mergdf))]

  # Filter pathways with no significant chanve
  if(!is.null(log2_effectSizeCutoff)){
      filt_idc <- sigm <= 0.05 & abs(greshm) >= log2_effectSizeCutoff
  } else { filt_idc <- sigm <= 0.05}
  
  # Divide into one table with fold-changes and one with presence of significance indicated with an asterisk
  greshm <- greshm[rowSums(filt_idc,na.rm = T)>0,]
  sigm <- sigm[rowSums(filt_idc,na.rm=T)>0,]
  sigm[sigm<=0.05] <- "*"
  sigm[sigm>0.05] <- ""
  sigm[is.na(sigm)] <- ""

  # Prepare output matrix
  for(i in 1:ncol(greshm)){greshm[,i] <- as.numeric(greshm[,i])}
  colnames(greshm) <- sub(effectSizeCol,"",colnames(greshm))
  return(list(greshm,sigm))
}

sortDEG_ABC <- function(greshm,sigm){
  
  idc1 <- grepl("vsCTRL",colnames(sigm))
  idc2 <- grepl(";ucamvcu",colnames(sigm)) & !idc1
  idc3 <- grepl(";epos",colnames(sigm)) & !idc1
  
  # A
  all.up.idc <- rowSums(sign(greshm)>0,na.rm=T) == ncol(greshm) & sigm[,idc1] == "*" & rowSums(sigm[,idc2] == "*",na.rm=T) > 0 & rowSums(sigm[,idc3] == "*",na.rm=T) > 0
  all.dwn.idc <- rowSums(sign(greshm)<0,na.rm=T) == ncol(greshm) & sigm[,idc1] == "*" & rowSums(sigm[,idc2] == "*",na.rm=T) > 0 & rowSums(sigm[,idc3] == "*",na.rm=T) > 0

  # B
  naflctrl.up.idc <- sigm[,idc1]=="*" & greshm[,idc1]>0 & rowSums(sigm[,idc2] == "*") == 0 & rowSums(sigm[,idc3] == "*",na.rm=T) == 0
  naflctrl.dwn.idc <- sigm[,idc1]=="*" & greshm[,idc1]<0 & rowSums(sigm[,idc2] == "*") == 0 & rowSums(sigm[,idc3] == "*",na.rm=T) == 0

  # C
  naflnash.up.idc <- rowSums(sign(greshm[,!idc1])>0,na.rm=T) == (ncol(greshm)-1) & sigm[,idc1] !="*" & rowSums(sigm[,idc2] == "*",na.rm=T) > 0 & rowSums(sigm[,idc3] == "*",na.rm=T) > 0
  naflnash.dwn.idc <- rowSums(sign(greshm[,!idc1])<0,na.rm=T) == (ncol(greshm)-1) & sigm[,idc1] !="*" & rowSums(sigm[,idc2] == "*",na.rm=T) > 0 & rowSums(sigm[,idc3] == "*",na.rm=T) > 0

  # Wrap up
  idcList <- list(list(all.up.idc,all.dwn.idc),list(naflctrl.up.idc,naflctrl.dwn.idc),list(naflnash.up.idc,naflnash.dwn.idc))
  pathList <- lapply(1:length(idcList), function(x) {  
    idc.up <- idcList[[x]][[1]]
    idc.dwn <- idcList[[x]][[2]]
    if(sum(idc.up)>0){
      df.up <- data.frame(dlayer=sub(".*_","",rownames(greshm[idc.up,])),"direction"=1)
    } else {df.up <- data.frame(dlayer=NA,"direction"=1)}
    if(sum(idc.dwn)>0){
      df.dwn <- data.frame(dlayer=sub(".*_","",rownames(greshm[idc.dwn,])),"direction"=-1)
    } else {df.dwn <- data.frame(dlayer=NA,"direction"=-1)}
    df.all <- rbind(df.up,df.dwn)
    df.all <- df.all[!is.na(df.all$dlayer),]
    return(list(df.all))})
  names(pathList) <- c("A","B","C")

  return(pathList)
}

mergeTables_degs <- function(gn_ep,gn_uv){
  gn_ep <- gn_ep %>% rename_all(~ifelse(. == "ENSEMBL_GENE_ID", ., paste0(.,";epos")))
  gn_uv <- gn_uv %>% rename_all(~ifelse(. == "ENSEMBL_GENE_ID", ., paste0(.,";ucamvcu")))
  
  degs_human <-full_join(gn_ep,gn_uv,by="ENSEMBL_GENE_ID")
  degs_human$`external_gene_name;epos` <- NULL
  colnames(degs_human) <- gsub("external_gene_name;ucamvcu","external_gene_name",colnames(degs_human))
  rownames(degs_human) <- paste(1:nrow(degs_human),degs_human$external_gene_name,sep="_")
  degs_human$external_gene_name <- NULL
  degs_human$ENSEMBL_GENE_ID <- NULL
  return(degs_human)
}

divideDEGs_AB_AC <- function(degs,effectSizeCol,sigCol,log2_effectSizeCutoff,path_dir){

  deg_filt_list = DEG_FilteredTables(degs,effectSizeCol,sigCol,log2_effectSizeCutoff)

  greshm <- deg_filt_list[[1]]
  sigm <- deg_filt_list[[2]]

  genes_abc_strict <- sortDEG_ABC(greshm,sigm)

  # Create gene sets merging into AB and AC
  genes_ab_ac <- genes_abc_strict
  genes_ab_ac$B[[1]] <- rbind.data.frame(genes_ab_ac$A[[1]],genes_ab_ac$B[[1]])
  genes_ab_ac$C[[1]] <- rbind.data.frame(genes_ab_ac$A[[1]],genes_ab_ac$C[[1]])
  genes_ab_ac$A <- NULL
  names(genes_ab_ac) <- c("AB","AC")
  genes_ab_ac$AB[[1]] <- genes_ab_ac$AB[[1]][!duplicated(genes_ab_ac$AB[[1]]$dlayer),]
  genes_ab_ac$AC[[1]] <- genes_ab_ac$AC[[1]][!duplicated(genes_ab_ac$AC[[1]]$dlayer),]

  saveRDS(genes_ab_ac,paste0(path_dir,"/MHPS_Pipeline/Data/DHPS_inputs/HumanGENES_AB_AC.RDS"))
}

