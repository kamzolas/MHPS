library(readxl)
library(dplyr)

cleanTables_rodents <- function(lly_de,mmu_de,rt,path_dir){
  mmu_de$SYMBOL <- mmu_de$external_gene_name
  mc <- plyr::join_all(list(lly_de,mmu_de), by="SYMBOL", type="full")
  
  # Convert mouse gene symbols into their human orthologs and create matrix of logfc ranks
  mc$humanOrth <- NA
  rt$humanOrth <- NA
  
  mtoh <- read.table(paste0(path_dir,"/MHPS_Pipeline/Data/Utils/MouseHuman_orthologyMapping.txt"),sep="\t",header=T)
  rtoh <- read.table(paste0(path_dir,"/MHPS_Pipeline/Data/Utils/rat_to_human.txt"),sep="\t",header=T)

  mtohsub <- mtoh[mtoh$Mouse.Marker.Symbol %in% mc$external_gene_name,]
  rtohsub <- rtoh[rtoh$ensembl_gene_id %in% rt$ENSGENEID,]
  
  mc$humanOrth[match(mtohsub$Mouse.Marker.Symbol,mc$external_gene_name)] <- as.character(mtohsub$Human.Marker.Symbol)
  rt$humanOrth[match(rtohsub$ensembl_gene_id,rt$ENSGENEID)] <- as.character(rtohsub$hsapiens_homolog_associated_gene_name)
  
  mc.medexpr <- apply(mc[,grep("Log2CPM",colnames(mc))],1,median)
  rt.medexpr <- apply(rt[,grep("Log2CPM",colnames(rt))],1,median)
  
  mc2 <- mc[order(mc.medexpr,decreasing=T),]
  rt2 <- rt[order(rt.medexpr,decreasing=T),]
  
  nrow(mc);mc2 <- mc2[!duplicated(mc2$humanOrth),];nrow(mc2)
  mc2 <- mc2[!is.na(mc2$humanOrth),];nrow(mc2)
  
  nrow(rt);rt2 <- rt2[!duplicated(rt2$humanOrth),];nrow(rt2)
  rt2 <- rt2[!is.na(rt2$humanOrth),];nrow(rt2)
  
  # Combine the mouse and rat datasets by the human ortholog genes
  anm <- plyr::join_all(list(mc2,rt2),by="humanOrth",type="full")
  saveRDS(anm,paste0(path_dir,"/MHPS_Pipeline/Data/Rodent_DEGs/DE_allRodents.RDS"))
  
  anm.filt <- anm[,grepl("log2FoldChange",colnames(anm)) | grepl("pvalue",colnames(anm)) | grepl("padj",colnames(anm))]
  rownames(anm.filt) <- anm$humanOrth
  
  saveRDS(anm.filt,paste0(path_dir,"/MHPS_Pipeline/Data/Rodent_DEGs/DE_allRodents_filtered.RDS"))
  return(anm.filt)
}

fc_to_ranks <- function(anm.filt,path_dir){
  # Convert all the fold changes into ranks
  anm.filt[1:2,grep("log2FoldChange",colnames(anm.filt))[1:2]]
  
  anm.ranked <- anm.filt
  for(i in 1:(ncol(anm.filt)/3)){
    idc <- grep("log2FoldChange",colnames(anm.ranked))[i]
    anm.ranked[,idc] <- rank(-1*anm.ranked[,idc],na.last = "keep",ties.method="random")
  }
  
  anm.ranked <- anm.ranked[,grep("log2FoldC",colnames(anm.ranked))]
  colnames(anm.ranked) <- sub("log2FoldChange_","",colnames(anm.ranked))
  
  saveRDS(anm.ranked,paste0(path_dir,"/MHPS_Pipeline/Data/DHPS_inputs/RodentRankedTable.RDS"))
}

