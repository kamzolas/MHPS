library(readxl)
library(caret)
library(PRROC)
library(dplyr)

# Divide rodents into disease categories CTRL, NAFL, NASHF2+

ds_ctrl <- function(ds_tbl,stea="S",infl="I(CRN)",bal="B(CRN)",fib="F(CRN)"){
  b1 <- ds_tbl[,stea] == 0
  b2 <- ds_tbl[,infl] == 0
  b3 <- ds_tbl[,bal] == 0
  b4 <- ds_tbl[,fib] == 0
  return(which(b1 & b2 & b3 & b4))
}
ds_nafl <- function(ds_tbl,stea="S",infl="I(CRN)",bal="B(CRN)",fib="F(CRN)"){
  b1 <- ds_tbl[,stea] > 0
  b2 <- ds_tbl[,infl] < 2
  b3 <- ds_tbl[,bal] == 0
  b4 <- ds_tbl[,fib] == 0
  return(which(b1 & b2 & b3 & b4))
}
ds_nash <- function(ds_tbl,stea="S",infl="I(CRN)",bal="B(CRN)",fib="F(CRN)"){
  b1 <- ds_tbl[,stea] > 0
  b2 <- (ds_tbl[,infl] >= 2 | ds_tbl[,bal] > 0) & ds_tbl[,fib] == 0
  b3 <- (ds_tbl[,infl] >= 0 | ds_tbl[,bal] >= 0) & ds_tbl[,fib] > 0
  return(which(b1 & (b2 | b3)))
}

cleanTables_createDiseaseGroups <- function(phd){
  phd$Disease_status <- NA
  
  phd$Disease_status[ds_nash(phd)] <- "NASH"
  phd$Disease_status[ds_ctrl(phd)] <- "CTRL"
  phd$Disease_status[ds_nafl(phd)] <- "NAFL"
  
  phd$Disease_status_fibrosis <- phd$Disease_status
  
  fib_cut <- 2
  phd$Disease_status_fibrosis[which(phd$Disease_status_fibrosis == "NASH" & as.numeric(phd$`F(CRN)`) >= fib_cut)] <- "NASHF2+"
  phd$Disease_status_fibrosis[which(phd$Disease_status_fibrosis == "NAFL")] <- "NAFL"  
  
  # Define the models/animals that should be used for evaluating the ALT and AST threshold
  # Some models are not included in the evaluation of the optimal ALT and AST thresholds due to technical limitations or disease relevance 
  
  not_eval <- c("OB-CHOW-8W","6J-CCl4-4W","6J-CCl4-6W","6J-CTRL-4W","6J-CTRL-6W")
  ds_tbl <- phd[phd$Sex=="M" & !phd$PublicationModelName %in% not_eval,]
  
  return(ds_tbl)
}

# Optimal threshold of ALT and AST ROC curves are found using Youden's index:
get_roc_gg <- function(phd,ds_ref,ds_comp,predcol,mdl_rm=NULL,ds_col,plot_title){
  phf1_ab <- phd[!is.na(phd$`B(CRN)`) & !is.na(phd$`I(CRN)`) & !is.na(phd[,ds_col]) & !is.na(phd[,predcol]) & (as.data.frame(phd)[,ds_col] %in% ds_ref | as.data.frame(phd)[,ds_col] %in% ds_comp),]
  phf1_ab <- phf1_ab[!phf1_ab$PublicationModelName %in% mdl_rm,]
  
  wc0 <- 1*(as.data.frame(phf1_ab)[,ds_col] %in% ds_comp)
  sc0 <- pull(phf1_ab,predcol)
  
  PRROC_obj <- roc.curve(scores.class0 = sc0, weights.class0=wc0,curve=TRUE)
  auc <- signif(PRROC_obj$auc,2)
  roc_curve_dat <- as.data.frame(PRROC_obj$curve)
  
  colnames(roc_curve_dat) <- c("1-Specificity","Sensitivity","-")
  
  best_thresh <- PRROC_obj$curve[which.max(PRROC_obj$curve[,2] - PRROC_obj$curve[,1]),3]
  val_at_best <- PRROC_obj$curve[which.max(PRROC_obj$curve[,2] - PRROC_obj$curve[,1]),]  
  roc_curve_dat$Best_thresh <- as.character(round(val_at_best[3],2))
  roc_curve_dat$AUC <- as.character(auc)
  roc_curve_dat$FPR <- val_at_best[1]
  roc_curve_dat$TPR <- val_at_best[2]
  
  gg_roc <- ggplot(data=roc_curve_dat,aes(x=`1-Specificity`,y=Sensitivity)) + 
    geom_point(aes(x=FPR,y=TPR,fill=AUC,color=Best_thresh),size=2) +
    scale_color_discrete(name=paste0(predcol," threshold"),labels=paste0(roc_curve_dat$Best_thresh," U/L")) +
    scale_fill_manual(values="white",name=paste0("AUC: ",roc_curve_dat$AUC),labels="") +
    geom_line(color="grey10",alpha=0.8) + 
    theme_bw() + 
    ggtitle(plot_title) +
    theme(legend.spacing.y=unit(0,"cm"),
          legend.text = element_text(size=9),
          legend.title = element_text(size=9),
          legend.box.margin = margin(0,0,0,0),
          legend.margin=margin(-3,-5,-3,-5),
          legend.justification = "right",
          aspect.ratio=1,
          plot.title =element_text(size=11, face='bold')) +
    geom_abline(slope=1,intercept=0,color=alpha("grey40",0.5)) +
    guides(fill = guide_legend(override.aes = list(color="white"),order = 1)) 
  gg_roc
  
  return(list("gg_roc_curve"=gg_roc,"PRROC_obj"=PRROC_obj,"Best_threshold"=best_thresh,"Val_at_bestTresh"=val_at_best))
}

get_sens_spec <- function(path_dir,phd,gg_alt,gg_ast){
  # Calculating sensitivity and specificity from combining AST and ALT thresholds
  phf2 <- phd[!is.na(phd$Disease_status) & !is.na(phd$`AST.(U/L)`),]
  alt_th2 <- gg_alt$Best_threshold
  ast_th2 <- gg_ast$Best_threshold
  
  # NASHF2+ vs CTRL
  d0 <- 1*(as.numeric(phf2$`ALT.(U/L)`) > alt_th2 & as.numeric(phf2$`AST.(U/L)`) > ast_th2)
  d0 <- factor(d0,levels=c("0","1"))
  ref0 <- phf2$Disease_status_fibrosis
  ref0[ref0 %in% c("CTRL")] <- 0
  ref0[ref0=="NASHF2+"] <- 1
  ref0 <- factor(ref0,levels=c("0","1"))
  
  res_list = list("alt_and_ast"=c("sensitivity"=sensitivity(d0,ref0,positive="1"),"specificity"=specificity(d0,ref0,negative="0")))
  
  # NASHF2+ vs CTRL: ALT and AST separately
  # ALT:
  d0 <- 1*(as.numeric(phf2$`ALT.(U/L)`) > alt_th2)
  d0 <- factor(d0,levels=c("0","1"))
  ref0 <- phf2$Disease_status_fibrosis
  ref0[ref0 %in% c("CTRL")] <- 0
  ref0[ref0=="NASHF2+"] <- 1
  ref0 <- factor(ref0,levels=c("0","1"))
  
  res_list = append(res_list,list("alt"=c("sensitivity"=sensitivity(d0,ref0,positive="1"),"specificity"=specificity(d0,ref0,negative="0"))))
  
  # AST:
  d0 <- 1*(as.numeric(phf2$`AST.(U/L)`) > ast_th2)
  d0 <- factor(d0,levels=c("0","1"))
  ref0 <- phf2$Disease_status_fibrosis
  ref0[ref0 %in% c("CTRL")] <- 0
  ref0[ref0=="NASHF2+"] <- 1
  ref0 <- factor(ref0,levels=c("0","1"))
  
  res_list = append(res_list,list("ast"=c("sensitivity"=sensitivity(d0,ref0,positive="1"),"specificity"=specificity(d0,ref0,negative="0"))))
  saveRDS(res_list,paste0(path_dir,"/MHPS_Pipeline/Data/PHPS_inputs/Sensitivity_specificity.RDS"))
}

tresholds_alt_ast <- function(path_dir,phd){
  gg_alt <- get_roc_gg(phd,ds_ref="CTRL",ds_comp="NASHF2+",predcol="ALT.(U/L)",mdl_rm=NULL,ds_col="Disease_status_fibrosis",plot_title="ALT: NASHF2+ vs CTRL")
  gg_ast <- get_roc_gg(phd,ds_ref="CTRL",ds_comp="NASHF2+",predcol="AST.(U/L)",mdl_rm=NULL,ds_col="Disease_status_fibrosis",plot_title="AST: NASHF2+ vs CTRL")

  saveRDS(gg_alt$Best_threshold,paste0(path_dir,"/MHPS_Pipeline/Data/PHPS_inputs/ALT_threshold.RDS"))
  saveRDS(gg_ast$Best_threshold,paste0(path_dir,"/MHPS_Pipeline/Data/PHPS_inputs/AST_threshold.RDS")) 
  
  get_sens_spec(path_dir,phd,gg_alt,gg_ast)
}




















