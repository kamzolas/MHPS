library(pheatmap)
library(readxl)
library(tibble)
library(dplyr)

# Functions
set_scores <- function(feat,score_list){
  scoretable <- data.frame(matrix(NA,ncol=3,nrow=length(feat)))
  colnames(scoretable) <- c("Up","NoReg","Dwn")
  rownames(scoretable) <- feat
  for(i in 1:length(score_list)){scoretable[i,] <- score_list[[i]]}
  return(scoretable)
}

# Calculate the statistics as a wilcoxon two group comparison or for ALT and AST using a global threshold.
calc_diff <- function(scoretable,mdls,thlist,phenotype_data){
  
  # Define the output tables
  model_scores <- data.frame(matrix(NA,nrow=nrow(scoretable),ncol=length(mdls)))
  rownames(model_scores) <- rownames(scoretable)
  colnames(model_scores) <- mdls
  model_signif <- model_scores
  model_effectsize <- model_scores
  model_pvals <- model_scores
  model_dietMeans <- model_scores
  model_ctrlMeans <- model_scores
  
  # Iterate over all the different phenotype variables (if phenotype variable is missing, NA will be inserted)
  for(j in 1:nrow(scoretable)){
    biochem <- rownames(scoretable)[j]
    
    # Phenotype variables evaluated as ctrl vs treatment
    if(biochem %in% colnames(phenotype_data)){
      for(i in 1:length(mdls)){
        mdl <- mdls[i]
        idc <- which(phenotype_data$PublicationModelName %in% mdl)
        valm <- as.numeric(as.data.frame(phenotype_data)[,biochem][idc])
        valc <- as.numeric(as.data.frame(phenotype_data)[,biochem][grepl(mdl,phenotype_data$PublicationModelName_ComparisonModel)])

        sgn <- sign(median(valm,na.rm=T)-median(valc,na.rm = T))
        med <- median(valm,na.rm=T)-median(valc,na.rm = T)
        fold_change <- mean(valm,na.rm=T)/mean(valc,na.rm = T)
        pct_change <- med/median(valc,na.rm = T)
        pvl <- try(wilcox.test(valm,valc)$p.value,silent=TRUE)
        rowidc <- which(rownames(model_scores) %in% biochem)
        model_scores[rowidc,mdl] <- med
        if(class(pvl)=="try-error"){
          model_signif[rowidc,mdl] <- NA
          model_effectsize[rowidc,mdl] <- NA
          model_pvals[rowidc,mdl] <- NA
          model_dietMeans[rowidc,mdl] <- NA
          model_ctrlMeans[rowidc,mdl] <- NA}
        if(class(pvl)!="try-error"){
          model_effectsize[rowidc,mdl] <- fold_change
          model_pvals[rowidc,mdl] <- pvl
          model_dietMeans[rowidc,mdl] <- mean(valm,na.rm=T)
          model_ctrlMeans[rowidc,mdl] <- mean(valc,na.rm = T)
          
          b1 <- sgn*(-log10(pvl)) >= -log10(0.05)
          b2 <- sgn*(-log10(pvl)) <= log10(0.05)
          if(b1){model_signif[rowidc,mdl] <- 1}
          if(b2){model_signif[rowidc,mdl] <- -1}
          if(!b1 & !b2){model_signif[rowidc,mdl] <- 0}
        }
      }
    }
    
    # Evaluate ALT and AST through a global cutoff
    if(biochem == "ALT&AST NAFLD" | biochem == "ALT&AST NASH"){
      for(i in 1:length(mdls)){
        mdl <- mdls[i]
        idc <- which(phenotype_data$PublicationModelName %in% mdl)
        val_ast <- as.numeric(as.data.frame(phenotype_data)[idc,"AST.(U/L)"])
        val_alt <- as.numeric(as.data.frame(phenotype_data)[idc,"ALT.(U/L)"])
        ast_alt_thresh <- median(val_ast,na.rm=T) >= thlist[["ALT&AST"]]["ast"] & median(val_alt,na.rm=T) >= thlist[["ALT&AST"]]["alt"]
        rowidc <- which(rownames(model_scores) %in% biochem)
        model_signif[rowidc,mdl] <- ast_alt_thresh
      }
    }
  }
  
  return(list("model_signif"=model_signif,"model_effectsize"=model_effectsize,"model_pvals"=model_pvals,"model_dietMeans"=model_dietMeans,"model_ctrlMeans"=model_ctrlMeans))
}

# For inspection, plot the resulting phenotype ranking as a heatmap
plot_pheatmap <- function(scoretable,model_signif,mdls_incl){
  models_tpl <- model_signif[which(rowSums(!is.na(scoretable))!=0),]
  scoretable <- scoretable[rowSums(!is.na(scoretable))!=0,]
  feat_idc <- rownames(scoretable)
  for(i in 1:length(feat_idc)){
    idc_up <- which(models_tpl[feat_idc[i],] == 1)
    idc_no <- which(models_tpl[feat_idc[i],] == 0)
    idc_dwn <- which(models_tpl[feat_idc[i],] == -1)
    models_tpl[feat_idc[i],idc_up]  <- scoretable[feat_idc[i],"Up"]
    models_tpl[feat_idc[i],idc_no]  <- scoretable[feat_idc[i],"NoReg"]
    models_tpl[feat_idc[i],idc_dwn]  <- scoretable[feat_idc[i],"Dwn"]
  }
  
  idx <- which(apply(apply(models_tpl,2,is.na),2,sum) == nrow(models_tpl))
  models_tpl <- rbind(models_tpl,apply(models_tpl,2,sum,na.rm=T))
  models_tpl[nrow(scoretable)+1,idx] <- NA
  rownames(models_tpl)[nrow(scoretable)+1] <- "SUM"
  
  toplot <- models_tpl[,order(models_tpl[nrow(scoretable)+1,],decreasing = F,na.last = F)]
  toplot <- toplot[,colnames(toplot) %in% mdls_incl]
  scorerange <- as.character(seq(min(toplot,na.rm=T),max(toplot,na.rm=T),by=0.5))
  
  clrs_set <- c("dodgerblue4","dodgerblue3","lightblue","white","lightgoldenrod1","gold2","darkorange","brown1","brown4","black")
  brks <- seq(-3,6,by=0.5)
  clrs <- colorRampPalette(clrs_set,space="rgb")(length(brks))
  names(clrs) <- brks
  
  toplot <- toplot[,ncol(toplot):1]
  colgaps <- cumsum(table(match(unique(toplot[rownames(toplot) == "SUM",]),toplot[rownames(toplot) == "SUM",])))
  toplot <- t(toplot)
  phm <- pheatmap(toplot,cluster_cols = F,cluster_rows = F,fontsize_col = 9,fontsize_row = 9,color = clrs[names(clrs) %in% scorerange],na_col="grey80",cellwidth = 11.5, cellheight = 9,gaps_col = (ncol(toplot)-1),gaps_row = colgaps,angle_col = 45,silent=T)
  
  toplot <- cbind(toplot,"PHPS"=scoreTransf(toplot[,"SUM"]))
  return(toplot)
}

# Score and rank the models
score_and_rank <- function(scoring_scheme,score_list,feat,phenotype_data,thlist,mdls_incl){
  names(score_list[[scoring_scheme]]) <- feat
  scoretable <- set_scores(feat,score_list[[scoring_scheme]])
  calc_diff_res <- calc_diff(scoretable,mdls=mdls_incl,thlist,phenotype_data)
  model_signif <- calc_diff_res$model_signif
  model_effectsize <- calc_diff_res$model_effectsize
  model_pvals <- calc_diff_res$model_pvals
  model_dietMeans <- calc_diff_res$model_dietMeans
  model_ctrlMeans <- calc_diff_res$model_ctrlMeans
  ph <- plot_pheatmap(scoretable,model_signif,mdls_incl)
  return(list(scoretable=scoretable,mdls=mdls_incl,model_signif=model_signif,model_effectsize=model_effectsize,phenotype_heatmap=ph,"model_pvals"=model_pvals,"model_dietMeans"=model_dietMeans,"model_ctrlMeans"=model_ctrlMeans))
}

# Function to transform scores to the interval [0,1]
scoreTransf <- function(x){
  x1 <- x-min(x,na.rm=T)
  return(x1/max(x1,na.rm=T))
}

compute_PHPS <- function(path_dir,phdb){
  # Define alt and ast thresholds, scoring schemes for PHPS and phenotype variables to be analysed
  
  # Define ALT and AST thresholds calculated by ROC curve analysis
  alt_thresh = readRDS(paste0(path_dir,"/MHPS_Pipeline/Data/PHPS_inputs/ALT_threshold.RDS"))
  ast_thresh = readRDS(paste0(path_dir,"/MHPS_Pipeline/Data/PHPS_inputs/AST_threshold.RDS"))
  
  thlist <- list("ALT&AST"=c(ast=ast_thresh,alt=alt_thresh))
  
  # Scoring schemes:
  score_list <- list(
    "AC_Fibrosis"=list("BW.(gram)"=c(0.5,0,-0.5),
                       "TGs.(mmol/L)"=c(0.5,0,-0.5),
                       "Cholesterol.(mmol/L)"=c(0.5,0,-0.5),
                       "LW/BW%"=c(1,0,-1),
                       "ALT.(U/L)"=c(NA,NA,NA),
                       "AST.(U/L)"=c(NA,NA,NA),
                       "ALT&AST NAFLD"=c(NA,NA,NA),
                       "ALT&AST NASH"=c(2,0,0)),
    "AB_Metabolic"=list("BW"=c(2,0,-2),
                        "TGs.(mmol/L)"=c(1,0,-1),
                        "Cholesterol.(mmol/L)"=c(1,0,-1),
                        "LW/BW%"=c(1,0,-1),
                        "ALT.(U/L)"=c(NA,NA,NA),
                        "AST.(U/L)"=c(NA,NA,NA),
                        "ALT&AST NAFLD"=c(0.5,0,0),
                        "ALT&AST NASH"=c(NA,NA,NA)))
  
  # Phenotype variables included in the two group tests
  feat <- c("BW.(gram)","TGs.(mmol/L)","Cholesterol.(mmol/L)","LW/BW%","ALT.(U/L)","AST.(U/L)","ALT&AST NAFLD","ALT&AST NASH")
  
  # Score and calculate PHPS for all disease models
  models_include = unique(phdb$PublicationModelName[phdb$PublicationDietGroup != "CTRL"])
  
  # Human Proximity scoring types
  hps_types = c("AB_Metabolic","AC_Fibrosis")
  
  # Run the 'Metabolic' or 'Fibrotic' phenotype scoring and ranking
  rank_res_list <- lapply(hps_types,function(x){
    rank_res <- score_and_rank(scoring_scheme=x,
                               score_list=score_list,
                               feat=feat,
                               phenotype_data=phdb,
                               thlist=thlist,
                               mdls_incl=models_include)  
  })
  names(rank_res_list) = hps_types
  
  # Save the result tables used for NHPS
  saveRDS(rank_res_list,paste0(path_dir,"/MHPS_Pipeline/Data/MHPS_inputs/PHPS_ResultTable.RDS"))  
}

