library(KEGGREST)
library(GSEABase)
library(gep2pep)
library(repo)
library(parallel)

# Extract human KEGG pathways through KEGG rest API
filtKGGnms <- function(qname){
  qname <- sub(" - Homo sapiens \\(human\\)","",qname)
  qname <- paste("KEGG",toupper(qname),sep="_")
  qname <- gsub(" ","_",qname)
  qname <- gsub("-","_",qname)
  qname <- gsub("\\(","",qname)
  qname <- gsub("\\)","",qname)
  qname <- gsub("_/","",qname)
  qname <- gsub(",","",qname)
  qname <- gsub("___","_",qname)
  qname <- gsub("__","_",qname)
  return(qname)
}

get_kegg_pathways <- function(path_dir,kegg_filename){
  kpathway_hsa <- keggList("pathway","hsa")
  
  # This takes some time to run
  testkgl <- lapply(1:length(kpathway_hsa),function(x){
    print(x)
    query <- keggGet(names(kpathway_hsa)[x])  
    if(length(grep("GENE",names(query[[1]])))!=0){
      qname <- query[[1]]$NAME
      qname <- filtKGGnms(qname)
      gns <- sub(";.*","",query[[1]]$GENE[seq(2,length(query[[1]]$GENE),2)])
      gns <- gns[!duplicated(gns)]
      testset <- GeneSet(gns,geneIdType=SymbolIdentifier(),setName=qname,collectionType=CategorizedCollection(category="KEGG",subCategory = "PATHWAYS"),setIdentifier = qname)  
      return(testset)
    }
  })
  testkgl <- testkgl[-which(sapply(testkgl,is.null))]
  
  keggpathway.list <- sapply(testkgl,geneIds)
  names(keggpathway.list) <- sapply(testkgl,setName)
  
  saveRDS(object=testkgl,file=paste0(path_dir,"/MHPS_Pipeline/Data/DHPS_inputs/",kegg_filename))
}

# Create the repos storing the human KEGG pathways
create_KEGG_repos <- function(keggpathways,path_dir){
  # Convert to a GeneSetCollection object
  db <- GeneSetCollection(keggpathways)
  
  # Define the path for the repo repository and create the rp object. Fill in kegg data.
  repoRoot_mmu <- file.path(paste0(path_dir,"/MHPS_Pipeline/Data/DHPS_inputs/Repos/gep2pep_KEGG_AB_AC_RNAseq_mmu"))
  rp_mmu <- createRepository(repoRoot_mmu, db)
  
  repoRoot_rat <- file.path(paste0(path_dir,"/MHPS_Pipeline/Data/DHPS_inputs/Repos/gep2pep_KEGG_AB_AC_RNAseq_rat"))
  rp_rat <- createRepository(repoRoot_rat, db)
  
  repoRoot_lly <- file.path(paste0(path_dir,"/MHPS_Pipeline/Data/DHPS_inputs/Repos/gep2pep_KEGG_AB_AC_Microarray_mmu"))
  rp_lly <- createRepository(repoRoot_lly, db)
}

# Create the repos storing the human DE genes

# Function to create the gene set in the right format
createEpoSDEGdb <- function(eposDEGlist,dbname){
  db.eposDEG <- c()
  for(i in 1:length(eposDEGlist)){
    eposDEG <- eposDEGlist[[i]]
    # Divide gene sets into up- and downregulation
    eposDEG.dwn <- eposDEG[eposDEG$direction == -1,]
    epos_gs_dwn <- GeneSet(as.character(eposDEG.dwn[,1]),geneIdType=SymbolIdentifier(),setName=paste0(dbname[i],"_Dwn"),collectionType=CategorizedCollection(category=dbname[i],subCategory = "DEG"),setIdentifier = paste0(dbname[i],"_Dwn")) 
    
    eposDEG.up <- eposDEG[eposDEG$direction == 1,]
    epos_gs_up <- GeneSet(as.character(eposDEG.up[,1]),geneIdType=SymbolIdentifier(),setName=paste0(dbname[i],"_Up"),collectionType=CategorizedCollection(category=dbname[i],subCategory = "DEG"),setIdentifier = paste0(dbname[i],"_Up")) 
    db.eposDEG <- c(db.eposDEG,GeneSetCollection(epos_gs_up),GeneSetCollection(epos_gs_dwn))
  }
  ### Gather all types of gene sets as one GeneSetCollection
  db <- GeneSetCollection(db.eposDEG)  
  return(db)
}

# Save the DEGs gene sets as repositories divided by technology and species
mk_DEGs_Repo <- function(modSetFiles,repoDir,path_dir){
  degList <- readRDS(modSetFiles)
  repoNames <- c("gep2pep_DEGs_RNAseq_mmu","gep2pep_DEGs_RNAseq_rat","gep2pep_DEGs_Microarray_mmu")
  repoNames <- paste0(repoDir,paste0(repoNames,"_"))
  for(i in 1:length(degList)){
    db <- createEpoSDEGdb(degList[[i]],names(degList))
    repoRoot_mmu <- file.path(paste0(repoNames[1],names(degList)[i]))
    rp_mmu <- createRepository(repoRoot_mmu, db)
    
    repoRoot_rat <- file.path(paste0(repoNames[2],names(degList)[i]))
    rp_rat <- createRepository(repoRoot_rat, db)
    
    repoRoot_lly <- file.path(paste0(repoNames[3],names(degList)[i]))
    rp_lly <- createRepository(repoRoot_lly, db)
  } 
}


