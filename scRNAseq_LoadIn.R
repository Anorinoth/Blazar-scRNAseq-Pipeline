#Clear current environment------------------------------------------------------------------------------

rm(list=ls())
graphics.off()

#Define and generate directories------------------------------------------------------------------------------
identifier<-"Blazar_Project_017"
str<-strsplit(Sys.info()[4], "\\.")[[1]]
ROOT_DIR="/mnt/v/Transcriptomics/Blazar_Project_017/Working Directory"
assignment<-"ATT6"

if (dir.exists(paste(ROOT_DIR,"/analysis/",sep=""))==FALSE) {
  dir.create(paste(ROOT_DIR,"/analysis/",sep=""))
}
if (dir.exists(paste(ROOT_DIR,"/analysis/", assignment,sep=""))==FALSE) {
  dir.create(paste(ROOT_DIR,"/analysis/", assignment,sep=""))
}

DATA_DIR <- file.path("/mnt/v/Transcriptomics/Blazar_Project_017/Data/")      # SPECIFY HERE
PROG_DIR <- file.path(ROOT_DIR, "analysis", assignment,"prog")     # SPECIFY HERE
if (dir.exists(paste(PROG_DIR))==FALSE) {
  dir.create(paste(PROG_DIR))
}
RES_DIR  <- file.path(ROOT_DIR, "analysis", assignment,"res")     # SPECIFY HERE
if (dir.exists(paste(RES_DIR))==FALSE) {
  dir.create(paste(RES_DIR))
}
RMD_DIR  <- file.path(ROOT_DIR, "analysis", assignment,"rmd")     # SPECIFY HERE
if (dir.exists(paste(RMD_DIR))==FALSE) {
  dir.create(paste(RMD_DIR))
}

#Define and generate Sub directories------------------------------------------------------------------------------
if (dir.exists(paste(RES_DIR,"/QC_Plots",sep=""))==FALSE) {
  dir.create(paste(RES_DIR,"/QC_Plots",sep=""))}
if (dir.exists(paste(RES_DIR,"/Dim_scaling",sep=""))==FALSE) {
  dir.create(paste(RES_DIR,"/Dim_scaling",sep=""))}
if (dir.exists(paste(RES_DIR,"/UMAP",sep=""))==FALSE) {
  dir.create(paste(RES_DIR,"/UMAP",sep=""))}
if (dir.exists(paste(RES_DIR,"/ClusterIDs",sep=""))==FALSE) {
  dir.create(paste(RES_DIR,"/ClusterIDs",sep=""))}
if (dir.exists(paste(RES_DIR,"/SampMatch",sep=""))==FALSE) {
  dir.create(paste(RES_DIR,"/SampMatch",sep=""))}
if (dir.exists(paste(RES_DIR,"/DeMulti",sep=""))==FALSE) {
  dir.create(paste(RES_DIR,"/DeMulti",sep=""))}

if (dir.exists(paste(RES_DIR,"/DeepClustAnalysis",sep=""))==FALSE) {
  dir.create(paste(RES_DIR,"/DeepClustAnalysis",sep=""))}

if (dir.exists(paste(RES_DIR,"/BioMarkers",sep=""))==FALSE) {
  dir.create(paste(RES_DIR,"/BioMarkers",sep=""))}
if (dir.exists(paste(RES_DIR,"/BioMarkers/TopMarkers",sep=""))==FALSE) {
  dir.create(paste(RES_DIR,"/BioMarkers/TopMarkers",sep=""))}
if (dir.exists(paste(RES_DIR,"/BioMarkers/ClusterMarkers",sep=""))==FALSE) {
  dir.create(paste(RES_DIR,"/BioMarkers/ClusterMarkers",sep=""))}

#Custom Functions------------------------------------------------------------------------------
using=function(pack) {
  libs<-unlist(list(pack))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  n<-length(need)
  if (n>0) {
    for (i in 1:n) {
      needing=need[i]
      deps <- tools::package_dependencies(packages = needing,
                                          recursive = TRUE)
      deps=unlist(list(deps))
      reqdep<-unlist(lapply(deps,require,character.only=TRUE))
      install.packages(need[i], INSTALL_opts = '--no-lock')
      if (system.file(package = need[i])=="") {
        BiocManager::install(need[i])
      }
      if (system.file(package = need[i])=="") {
        print(paste("Error auto-installing package ",need[i],". Please install manually",sep=""))
      }
    }
  }
  }

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

iDEA.louis.edit <- function(object){
  num_core <- object@num_core
  
  ## using parallel to correct
  library(doSNOW)
  if(num_core > 1){
    if(num_core > detectCores()){warning("LOUIS:: the number of cores you're setting is larger than detected cores!");num_core = detectCores()}
  }#end fi
  ## modified here if need parallel, by default num_core=1
  if(.Platform$OS.type == "windows"){num_core <- 1}
  cl <- makeCluster(num_core)
  registerDoSNOW(cl)
  num_annot <- length(object@de)
  
  pb <- txtProgressBar(max = num_annot, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  #parallel
  
  res_all <- pbmclapply(1:num_annot, FUN = function(i){
    res <- object@de[[i]]
    if (is.null(object@de[[i]])) {
      NULL
    } else {
      Annot <- rep(0, object@num_gene)
      Annotind = which(names(object@annotation) == names(object@de)[i])
      Annot[object@annotation[[Annotind]]] <- 1
      Annot <- Annot - mean(Annot)
      Annot <- as.matrix(data.frame(rep(1, object@num_gene), Annot) )
      ##################
      ## Louis function
      LouisMethod <- function(res, A){
        numer <- (1-res$pip)*res$pip*pnorm(res$beta, mean=0, sd=sqrt(res$sigma2_beta*res$sigma2_e))*pnorm(res$beta, mean=0, sd=sqrt(res$sigma2_beta2*res$sigma2_e))
        denom <- res$pip*pnorm(res$beta, mean=0, sd=sqrt(res$sigma2_beta*res$sigma2_e)) +
          (1-res$pip)*pnorm(res$beta, mean=0, sd=sqrt(res$sigma2_beta2*res$sigma2_e))
        sel_index <- which(!is.na(numer/denom^2)) ### delet the na genes
        numer <- numer[sel_index]
        denom <- denom[sel_index]
        A <- A[sel_index,]
        # define missing information matrix
        T <- matrix(0, ncol=2, nrow=2)
        T[1,1] <- sum( as.numeric((numer/denom^2) * A[,1]*A[,1] ) )
        T[1,2] <- sum( as.numeric((numer/denom^2) * A[,1]*A[,2] ) )
        T[2,1] <- sum( as.numeric((numer/denom^2) * A[,1]*A[,2] ) )
        T[2,2] <- sum( as.numeric((numer/denom^2) * A[,2]*A[,2] ) )
        
        annot_var <- diag( solve(res$info_mat - T) )
        # return the results
        return(annot_var)
      }# end function
      #################
      
      if(!is.na(res$annot_coef[2])){
        try_test <- try( louis_var <- LouisMethod(res, Annot), silent=T)
        if(class(try_test)=="try-error"){
          print("try-error!")
          resTemp = NULL
        }else{
          resTemp = data.frame(annot_id=object@annot_id[Annotind], annot_coef=res$annot_coef[2], annot_var=res$annot_var[2], annot_var_louis=louis_var[2], sigma2_b=res$sigma2_beta)
        }
      }else{
        resTemp = NULL
      }
    }}) #### end fi
  res_all <- do.call(rbind, res_all)    
  ## remove the negative variance genes
  pos_index <- which(res_all$annot_var_louis>0)
  zscore <- res_all$annot_coef[pos_index]/sqrt(res_all$annot_var[pos_index])
  zscore_louis <- res_all$annot_coef[pos_index]/sqrt( res_all$annot_var_louis[pos_index] )
  ## only keep positive variance corrected by louis method, modified by sun, 2019-10-8 16:28:30
  res_all <- res_all[pos_index, ]
  res_all$pvalue_louis <- 2*pnorm(-abs(zscore_louis))
  res_all$pvalue <- 2*pnorm(-abs(zscore))
  object@gsea <- res_all
  rm(res_all)
  # return the results
  return(object)
} #end funcs

# Package Load In------------------------------------------------------------------------------
packages = c("Seurat","RColorBrewer","ggplot2","Matrix","gridExtra","dplyr","slickR","cowplot",
             "tictoc","patchwork","HGNChelper","ggraph","igraph","data.tree","scater","clustree",
             "future","prettydoc","limma","openxlsx","imager","devtools","tidyverse","cowplot",
             "Matrix.utils","edgeR","magrittr","purrr","reshape2","S4Vectors","tibble",
             "SingleCellExperiment","pheatmap","apeglm","png","DESeq2","piano","pbmcapply","ggrepel",
             "readxl","knitr","fgsea","data.table","tools","randomcoloR","ggpubr","BiocParallel",
             "pracma","biomaRt","furrr","grid","vsn","PCAtools","ggalt","ggforce","kableExtra",
             "msigdbr","GSEABase","qusage","rlist","htmlwidgets","DelayedMatrixStats","ggnewscale",
             "scales","gprofiler2","snowfall","raster","cate","doSNOW","doParallel")

using(packages)

if(!require(iDEA)){
  devtools::install_github('xzhoulab/iDEA')
  library(iDEA)
} else {
  library(iDEA)
}

#Load Databases------------------------------------------------------------------------------
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";

gsc = getBroadSets(paste(DATA_DIR,"/Enrichment/msigdb_v7.4.xml",sep=""))

#Set options for rendering------------------------------------------------------------------------------
options(stringsAsFactors=FALSE)
knitr::opts_chunk$set(fig.width=7, fig.height=5,
                      echo=FALSE, warning=FALSE, message=FALSE)
knitr::opts_chunk$set(dev.args=list(bg="transparent"))

tc<-theme_bw()+
  theme(
    panel.background = element_rect(fill = "transparent") # bg of the panel
    , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
    , legend.background = element_rect(fill = "transparent") # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )
theme_set(tc)

plan("multisession", workers = 4)

