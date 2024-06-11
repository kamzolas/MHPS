library(repo)
library(gep2pep)
library(GSEABase)
library(readxl)
library(stringr)
library(dplyr)

# Functions to extract ES matrices from the repo
extrESM <- function(fpv,repo_db){
  rp_path <- file.path(fpv[1])
  rpnew <- repo_open(rp_path)
  esmall <- as.data.frame(loadESmatrix(rpnew, repo_db),stringsAsFactors=F);esmall$ids <- rownames(esmall)
  pvmall <- as.data.frame(loadPVmatrix(rpnew, repo_db),stringsAsFactors=F);pvmall$ids <- rownames(pvmall)
  for(i in 2:length(fpv)){
    rp_path <- file.path(fpv[i])
    rpnew <- repo_open(rp_path)
    esm <- as.data.frame(loadESmatrix(rpnew, repo_db),stringsAsFactors=F)
    esm$ids <- rownames(esm)
    esmall <- full_join(esmall,esm,by="ids")
    pvm <- as.data.frame(loadPVmatrix(rpnew, repo_db),stringsAsFactors=F)
    pvm$ids <- rownames(pvm)
    pvmall <- full_join(pvmall,pvm,by="ids")
  }
  psig <- pvmall < 0.05
  psig[psig==T] <- "*";psig[psig=="FALSE"] <- ""
  return(list(esmall,pvmall,psig))
}

# Function to prepare the gene sets for DSEA
finalmod <- function(module1res,repo_db,setnm,filterRows=F,fpv,nrowToFilter=NULL,resall=NULL,n=NULL){
  m1res <- module1res[[1]]
  m1res$ID <- paste0(m1res$ID,setnm)
  m1res <- list(m1res)
  names(m1res) <- repo_db
  
  if(filterRows==T){
    print("Find number of rows to filter...")
    pruneRes <- pruneHPS(fpv,fmod,resall,repo_db,setnm,n)
    nrowToFilter <- pruneRes[[1]]
    esmheat.res <- pruneRes[[2]]
    plm <- esmheat.res[[n]]
    
    print("Starting to filter the pathwayset")
    
    upreg.tokeep <- rownames(plm)[1:nrowToFilter]
    print(m1res[[repo_db]]$ID)
    filt.idc <- m1res[[repo_db]]$ID %in% upreg.tokeep

    m1res.up <- m1res[[repo_db]][m1res[[repo_db]]$NES == 1,]
    m1res.dwn <- m1res[[repo_db]][m1res[[repo_db]]$NES == -1,]
    if(n==2){
      m1res.up <- m1res[[repo_db]][m1res[[repo_db]]$NES == 1 & filt.idc,]
    } else {
      m1res.dwn <- m1res[[repo_db]][m1res[[repo_db]]$NES == -1 & filt.idc,]
    }
    m1res <- list(rbind(m1res.up,m1res.dwn))
    names(m1res) <- repo_db
  }
  return(m1res)
}

# Function to merge the DSEA results to obtain only one list with DSEA results across all rodent groups
combResTables <- function(mod1all,excl=NULL){
  mod1all.up <- lapply(1:length(mod1all), function(x) mod1all[[x]][[1]]$Upregulated)
  tbl.up <- Reduce(rbind,mod1all.up)
  tbl.up <- try(tbl.up[order(sign(tbl.up$ES)*(-log10(tbl.up$PV)),decreasing=T),])
  mod1all.dwn <- lapply(1:length(mod1all), function(x) mod1all[[x]][[1]]$Downregulated)
  tbl.dwn <- Reduce(rbind,mod1all.dwn)
  tbl.dwn <- try(tbl.dwn[order(sign(tbl.dwn$ES)*(-log10(tbl.dwn$PV)),decreasing=F),])
  resb <- list(tbl.dwn,tbl.up);names(resb) <- c("Downregulated","Upregulated")
  res <- c(resb,mod1all[[1]][[2]])
  return(res)
}

# Function to aggregate result from up- and downregulated pathwayset analyses into one table based on summed ES
combES2 <- function(resall){
  resmf.dwn <- resall$Downregulated
  resmf.dwn <- try(resmf.dwn[order(rownames(resmf.dwn)),])
  resmf.up <- resall$Upregulated
  resmf.up <- resmf.up[order(rownames(resmf.up)),]
  if(class(resmf.dwn)=="try-error"){
    allmf <- cbind(matrix(NA,ncol=2,nrow=nrow(resmf.up)),resmf.up)
  } else {allmf <- cbind(resmf.dwn,resmf.up)}
  
  names(allmf) <- c("ES.dwn","PV.dwn","ES.up","PV.up")
  allmf$ES.dwn[allmf$PV.dwn > 0.05] <- 0
  allmf$ES.up[allmf$PV.up > 0.05] <- 0
  allmf$ES.combined <- rowSums(cbind(-1*allmf[,"ES.dwn"],allmf[,"ES.up"]),na.rm=T)
  allmf <- allmf[order(allmf$ES.combined,decreasing = T),]
  return(allmf)
}

# Function to transform scores to the interval [0,1]
scoreTransf <- function(x){
  x1 <- x-min(x,na.rm=T)
  return(x1/max(x1,na.rm=T))
}

# Functions to rank rodents according to KEGG
pthsea <- function(input_pathways,rp,collection_name,pv.cut){
  cat(sprintf("Running drug2gene using %s \n", collection_name))
  GNsets <- loadCollection(rp, collection = collection_name)
  ksn <- sapply(GNsets,setName)
  
  input_pathways <- input_pathways[!is.na(input_pathways$ID),]
  ftry.dwn <-input_pathways$ID[input_pathways$pvalue < pv.cut & input_pathways$NES < 0]
  ftry.up <- input_pathways$ID[input_pathways$pvalue < pv.cut & input_pathways$NES > 0]
  cat(sprintf("Number of downregulated pathways surviving the cutoff: %d \n",length(ftry.dwn)))
  cat(sprintf("Number of upregulated pathways surviving the cutoff: %d \n",length(ftry.up)))
  allpaths <- list(ftry.dwn,ftry.up)
  
  shared_pathways <- list(c(),c())
  allres <- list(c(),c())
  names(allres) <- names(shared_pathways) <- c("Downregulated","Upregulated")
  for(i in 1:length(allpaths)){
    if(length(allpaths[[i]]) > 4){
      pathways <- sapply(GNsets,setIdentifier)[which(ksn %in% allpaths[[i]])]
      cat(sprintf("Shared %s pathways: %d \n", names(allres)[i],length(pathways)))
      subdb <- GNsets[sapply(GNsets, setIdentifier) %in% pathways]
      psea <- try(PathSEA(rp, subdb))
      res <- try(getResults(psea, collection_name))
      print(class(res))
      if(class(res)!="try-error"){
        print("ok")
        rl1p <- res[res$ES > 0,];rl1p <- rl1p[order(rl1p$ES,decreasing=T),]
        rl1n <- res[res$ES < 0,];rl1n <- rl1n[order(rl1n$ES,decreasing=F),]
        if(i == 1){res <- as.data.frame(rbind(rl1n,rl1p[order(rl1p$ES),]))} 
        if(i == 2){res <- as.data.frame(rbind(rl1p,rl1n[order(rl1n$ES),]))}
        allres[[i]] <- res
      } else {allres[[i]] <- res}
      shared_pathways[[i]] <- ksn[ksn %in% allpaths[[i]]]
    }     
  }
  return(list(allres,shared_pathways))
}

lay1_kegg <- function(human_module_sets,fpr,fpv,grp1=NULL,grp2=NULL,combineGrps=F){
  pathway_dbList <- c("KEGG_PATHWAYS")
  txtadd <- ""
  
  resList <- list()
  for(i in 1:length(human_module_sets[[1]])){
    finres <- lapply(1:length(human_module_sets), function(x) {
      pathw_set <- human_module_sets[[x]]
      fmod <- finalmod(pathw_set,repo_db = pathway_dbList, setnm = txtadd,filterRows=F,fpv)
      mod1all <- lapply(1:length(fpv), function(x) {
        rp_path <- file.path(fpv[x])
        rpnew <- repo_open(rp_path)
        d2gres <- pthsea(fmod[[1]],rpnew,names(fmod),0.05)
      })
      resall <-combResTables(mod1all,excl)
      resall <- combES2(resall)
      return(resall)
    })
    names(finres) <- names(human_module_sets)
    resList[[i]] <- finres
    for(j in 1:length(resList[[i]])){
      resList[[i]][[j]] <- resList[[i]][[j]][order(rownames(resList[[i]][[j]])),]
    }
  }
  
  # Create the table with average ES.combined obtained from A + C or A + B
  phw <- resList[[1]][[grp1]]
  
  phw$ES.combined <- scoreTransf(phw$ES.combined)
  
  phw <- phw[order(rownames(phw)),]
  colnames(phw) <- paste0("KEGG.",colnames(phw))
  
  allres <- phw
  allres$models <- rownames(allres)  
  return(allres)
}

# Function to rank rodents according to DEGs
lay3_deg <- function(fpr,fpv,grp1=NULL,grp2=NULL,layer_nms){
  fpv_list <- lapply(1:length(layer_nms),function(x) fpv[grep(paste0("_",layer_nms)[x],str_sub(fpv,-3))])
  names(fpv_list) <- layer_nms
  rp_path <- file.path(fpv_list[[1]][2])
  rp.deg <- repo_open(rp_path)
  repo_db = as.character(rp.deg$print()[4,1])
  
  deg_list <- lapply(1:length(fpv_list),function(x){
    fpv_x <- fpv_list[[x]]
    epos.degres <- extrESM(fpv=fpv_x,repo_db=repo_db)
    
    degres <- t(epos.degres[[1]])
    colnames(degres) <- paste0("ES.",sub(".*_","",degres["ids",]))
    degres <- as.data.frame(degres[rownames(degres) != "ids",],stringsAsFactors=F)
    
    degresSig <- t(epos.degres[[2]])
    colnames(degresSig) <- paste0("PV.",sub(".*_","",degresSig["ids",]))
    degresSig <- as.data.frame(degresSig[rownames(degresSig) != "ids",],stringsAsFactors=F)

    epdeg <- as.data.frame(cbind(degres,degresSig),stringsAsFactors=F)
    for(i in 1:ncol(epdeg)){epdeg[,i] <- as.numeric(epdeg[,i])}
    epdeg <- epdeg[,c(2,4,1,3)]
    names(epdeg) <- c("DEGS.ES.dwn","DEGS.PV.dwn","DEGS.ES.up","DEGS.PV.up")
    epdeg$DEGS.ES.dwn[epdeg$DEGS.PV.dwn > 0.05] <- 0
    epdeg$DEGS.ES.up[epdeg$DEGS.PV.up > 0.05] <- 0
    epdeg$DEGS.ES.combined <- rowSums(cbind(-1*epdeg[,"DEGS.ES.dwn"],epdeg[,"DEGS.ES.up"]),na.rm=T)
    epdeg <- epdeg[order(epdeg$DEGS.ES.combined,decreasing = T),]
    
    epdeg$models <- rownames(epdeg)
    return(epdeg)
  })
  names(deg_list) <- names(fpv_list)
  
  # Ensure that the order of models is the same across A, B and C
  for(i in 1:length(deg_list)){
    deg_list[[i]] <- deg_list[[i]][order(rownames(deg_list[[i]])),]
  }
  
  epdeg <- deg_list[[grp1]]
  epdeg$DEGS.ES.combined <- scoreTransf(epdeg$DEGS.ES.combined)
  
  return(epdeg)
}

# Combine the rankings from KEGG and DEGs
combineLayers <- function(human_module_sets,fpr,fpv,grp1,grp2,layer_nms){
  allres <- lay1_kegg(human_module_sets=human_module_sets,fpr[1],fpv=fpv$KEGG,grp1=grp1,grp2=grp2)
  epdeg <- lay3_deg(fpr=fpr,fpv=fpv$DEG,grp1=grp1,grp2=grp2,layer_nms=layer_nms)
  
  # Merge all the DSEA layers into one table
  allres2 <- Reduce(function(x,y) merge(x=x,y=y,by="models",all=T),list(allres,epdeg))
  
  # Create the merged score of KEGG pathways and DE genes
  allres2$DSEA.combined <- rowSums(allres2[,grepl("KEGG.ES.combined",colnames(allres2)) | grepl("DEGS.ES.combined",colnames(allres2))],na.rm=T)/2
  allres2$DHPS <- scoreTransf(allres2$DSEA.combined)
  
  return(allres2)
}

# Compute DHPS and save the result
compute_DHPS <- function(path_dir){
  hs_modSet_list <- readRDS(paste0(path_dir,"/MHPS_Pipeline/Data/DHPS_inputs/HumanKEGG_AB_AC.RDS"))
  layer_names <- c("AB","AC")
  combGrps <- c(F)
  fpr_vec <- paste0(path_dir,"/MHPS_Pipeline/Data/DHPS_inputs/Repos/")
  fpvType_list <- list(
    "KEGG"=paste0(fpr_vec,c("gep2pep_KEGG_AB_AC_RNAseq_mmu/","gep2pep_KEGG_AB_AC_RNAseq_rat/","gep2pep_KEGG_AB_AC_Microarray_mmu/")),
    "DEG"=paste0(fpr_vec,list.files(fpr_vec)[grepl("DEGs_",list.files(fpr_vec))])
  )
  rank_type <- list("AB_Metabolic"=list("grp1"="AB","grp2"="AB"),
                    "AC_Fibrosis"=list("grp1"="AC","grp2"="AC"))
  ranking_type = c("AB_Metabolic","AC_Fibrosis")
  path_dir = paste0(path_dir,"/MHPS_Pipeline/Data/MHPS_inputs/")
  
  # Calculate DHPS
  resdf_list = lapply(ranking_type,function(i){
    resdf <- combineLayers(human_module_sets = hs_modSet_list,
                           fpr=fpr_vec,
                           fpv=fpvType_list,
                           grp1=rank_type[[i]]$grp1,
                           grp2=rank_type[[i]]$grp2,
                           layer_nms=layer_names)
    output_dir = paste0(path_dir,"DHPS_",i,".RDS")
    saveRDS(resdf,output_dir)
  })  
}


