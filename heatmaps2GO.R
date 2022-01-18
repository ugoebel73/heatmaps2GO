library(AnnotationDbi)
library(org.Mm.eg.db)

library(limma)
library(edgeR)

library("pheatmap")
library(cowplot)

## have these last, to not have to qualify each tidyverse verb (?better to only attach them?)
library(dplyr)
library(tidyr)
library(stringr)

library("SummarizedExperiment") ## how can I avoid this? It masks a lot of basic functions.


source("functions_data.R")
source("functions_transformations.R")
source("functions_visualization.R")


## Load the Carien Niessen & Oana Persa tximport results (from kallisto), transformed into a SummarizedExperiment:
## Here, the organism name and AnnotationDbi database are carried with the object as  metadata(se)$org_name,metadata(se)$org_db_name.
se <- make_niessen_se(use_all_counts=TRUE, model="MODEL2")


## The GO terms for which we want to plot heatmaps of the count data
GO_terms <- read.csv("GO_terms.csv",comment.char="#")[,"GO_term"] ## ASSUME that the file contains a column "GO_term" 

## The direct "is a" children of the GO terms of interest -- possibly of use later
GO_isa_children <- sapply(list(BP=GO.db::GOBPCHILDREN,
                               MF=GO.db::GOMFCHILDREN,
                               CC=GO.db::GOCCCHILDREN),
                          function(ch) {
                              l <- sapply(as.list(ch)[GO_terms],
                                          function(x) {y <- x[names(x)=="isa"]; names(y) <- NULL;y},
                                          simplify=FALSE)
                              l[sapply(l,length)>0]
                          },simplify=FALSE)

## Associate GO terms with gene names
genes2GO <- enrichR2GO(GO_terms)

## If the type of genes2GO$gene_id differs from metadata(se)$gene_id_type, then we need a mapping between the tewo types:
if(genes2GO$type != metadata(se)$gene_id_type) { 
    ##tmp <- AnnotationDbi::mapIds(x = eval(as.symbol(metadata(se)$org_db_name)), 
    ##                             keys = rownames(se),
    ##                             column = genes2GO$type,
    ##                             keytype = metadata(se)$gene_id_type,
    ##                             multiVals = "list") ## We don't expect a 1:1 mapping !
    ##gene_id_map <- tibble(old=names(tmp),new=tmp) %>% unnest(new); rm(tmp)
    
    gene_id_map <- tibble(old=names(metadata(se)$gene_map),new=metadata(se)$gene_map)
    ## more natural, but !!se should caryy along the type of the mapped IDs, which should be cross-checked with genes2GO$type!!!
    
} else {
    gene_id_map <- NULL ## means "don't map in remap_rownames()
}

## Retain in genes2GO only those genes which have a mapping type1<->type2 in the gene_id_map (that is, which can be mapped):
genes2GO$tbl <- genes2GO$tbl %>% filter(if(is.null(gene_id_map)) TRUE else (gene_id  %in% toupper(gene_id_map$new)))

## NOTE that the symbols in the Enrichr DBs are all uppercase (like human gene symbols, but at the time 
##      of the initial analysis of Carien's data the Enricher website used these DBs for {human, mouse, rat}).
##      Mouse symbols (only first character uppercase) must be made all-uppercase before using the gene_id_map !
##      (One could make the filtering into a function and set the conversion function as a parameter.)


############## RECOMPUTE FROM HERE ON IF SOMETHING HAS CHANGED:
## Rename and **collapse rows with ambiguous name mapping before scaling**:
tmp <- remap_rownames(assays(se)$counts,  ## if rows with amiguous name mapping must be collapsed, do it before scaling!
                      gene_id_map,
                      merge_function=colSums)
mtrx  <-        tmp$tbl
mtrx_mapping <- tmp$grp


## Read out the GSEA (here really:  Enrichr GO ORA) results of the original analysis, group genes by significant GO term
tmp <- sapply(names(metadata(se)$Enrichr_GO_minFC1.5_maxQ0.05_MODEL2),
              function(nm) {
                  v        <- as.data.frame(metadata(se)$Enrichr_GO_minFC1.5_maxQ0.05_MODEL2[[nm]])[,"DGEgenes_in_GSEAset"]
                  names(v) <- as.data.frame(metadata(se)$Enrichr_GO_minFC1.5_maxQ0.05_MODEL2[[nm]])[,"GSEAid"]

                  sapply(v[intersect(names(v),GO_terms)],strsplit,",")
              })
GO_groups <- sapply(GO_terms,
                    function(this_term) {
                        Reduce(union,
                               sapply(names(which(sapply(tmp,function(x)this_term %in% names(x)))),
                                      function(nm)tmp[[nm]][[this_term]],
                                      simplify=FALSE)
                               )
                    },simplify=FALSE)

## define the final set of "genes of interest", now in the namespace of the new gene ID type, in case there was a mapping
if(!is.null(gene_id_map)) {
    ## gene_ids <- unique((tibble(new=names(mtrx_mapping),
    ##                            old=      mtrx_mapping) %>% unnest(old) %>% filter(old %in% metadata(se)$gene_ids))$new
    ##                   )
    
    gene_ids <- Reduce(union,GO_groups) ## use the mapped names used in GSEA -- they do have the correct type because they
                                        ## came from the same analysis as metadata(se)$gene_map,
                                        ## but this is not obvious in the current logical flow of the script! CLEAN UP!
} else {
    gene_ids <- metadata(se)$gene_ids
}

## =======================================================================================================================
## -----------------------------------------------------------------------------------------------------------------------
library(DESeq2)
library(ggplot2)

my_plotPCA <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE, returnScree=FALSE,
                        PCs=c(x=1,y=2),
                        pointSize=3,
                        hjust=-0.75,vjust=0.5,
                        sampleMap=NULL,main="",
                        ###xlim=c(-35,35),ylim=c(-35,35))  {
                        abs_xmax=NA,abs_ymax=NA)  {
    ## from /home/ugoebel/CECAD/Reznick/Project_Berlin_RNA-seq/Analysis/DGE/deseq2_plotPCA2.R
    ## Not sure where I found the code of DESeq2::plotPCA() originally, but one version is here:
    ## https://github.com/mikelove/DESeq2/blob/master/R/plots.R
    ## Or see here: https://rdrr.io/github/mikelove/DESeq2/man/plotPCA.html
    ## (-> type "DESeq2:::plotPCA.DESeqTransform") [this is how I likely got it initially]
    
    library(ggplot2)
    ##source("~/CECAD/Programming/Sharable_code/General/general_functions.R")

    
    rv <- rowVars(assay(object)) ## originally expected to be the output of rlog(dds)
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
        drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }
    else {
        colData(object)[[intgroup]]
    }
    d <- data.frame(pca$x[, PCs[1]], pca$x[, PCs[2]], group = group, 
                    intgroup.df, name = colnames(object))
    colnames(d)[1:2] <- paste("PC",PCs,sep="")
    if (returnData) {
        attr(d, "percentVar") <- percentVar[PCs]
        ##return(d)
        return(list(d=d,
                    raw=pca))
    } else if (returnScree) {
        ## copied from https://support.bioconductor.org/p/83626/
        scree_plot <- data.frame(percentVar*100)

        scree_plot[,2]<- c(1:nrow(scree_plot))
        
        colnames(scree_plot)<-c("percent_variance_explained","principle_component_number")
        return(ggplot(scree_plot, mapping=aes(x=principle_component_number, y=percent_variance_explained))+geom_bar(stat="identity"))
    }
    if(!is.null(sampleMap)) {
        rownames(d) <- sampleMap[rownames(d)]
    }
    g <- ggplot(data = d, aes_string(x = paste("PC",PCs["x"],sep=""), y = paste("PC",PCs["y"],sep=""), color = "group")) + 
                geom_point(size = pointSize) + geom_text(aes(label=rownames(d),hjust=hjust, vjust=vjust)) +
                xlab(paste0("PC",PCs["x"],": ", round(percentVar[PCs["x"]] * 100), "% variance",sep="")) + 
                ylab(paste0("PC",PCs["y"],": ", round(percentVar[PCs["y"]] * 100), "% variance",sep="")) + 
                coord_fixed()
    if(!is.na(abs_xmax)) g <- g + xlim(-abs_xmax,abs_xmax)
    if(!is.na(abs_ymax)) g <- g + ylim(-abs_ymax,abs_ymax)

    g + ggtitle(main)
}


dds <- estimateSizeFactors(DESeqDataSetFromMatrix(countData=round(mtrx[rownames(mtrx)[rownames(mtrx) %in% gene_ids],]),
                                                  colData = colData(se), design=metadata(se)$design))
rld <- rlog(dds, blind=TRUE)

## see /home/ugoebel/CECAD/Reznick/Project_Berlin_RNA-seq/Analysis/DGE/DGE_v3.R

## These are PCAs on a subset of a pre-normalized full dataset

p <- list()
l <- c()
k <- 1
for(i in 1:3) {
    for(j in (i+1):4) {
        cat(i,j,"\n")
        p[[k]] <- my_plotPCA(rld, intgroup=c("Condition"), sampleMap=as.matrix(colData(se))[,"SampleName"],
                             PCs=c(x=i,y=j),
                             hjust=0.2,vjust=-0.5)
        ##main=paste0("PC_",i," vs PC_",j)))
        l[k] <- paste0("PC_",i," vs PC_",j)
        k <- k+1
    }
}
cowplot::plot_grid(p[[1]],p[[2]],p[[3]],
                   p[[4]],p[[5]],p[[6]],
                   ncol = 3, nrow = 2,
                   labels=l)

dev.copy2pdf(file="mtrxPCA_PC1-4.pdf")

v <- colSums(assays(rld)[[1]]); names(v) <- as.matrix(colData(se))[,"SampleName"][names(v)]
##v
##  ALU_584   ALU_713   ALU_714   ALU_715   ALU_720      90-4      90-5     274-5  ANK_2492      95-1      95-4     174-4 
## 2725.667  2731.006  2731.304  2729.972  2764.972  2753.604  2745.221  2744.958  2759.432  2752.254  2739.851  2741.030 
## ANK_2488      95-2     14_06     176_1 ANK_275-1 
## 2757.476  2759.260  2768.935  2761.223  2766.697

tmp <-  my_plotPCA(rld, intgroup=c("Condition"), returnData=TRUE)
PCs_gene_ids <- tmp$d
PCs_gene_ids_raw <- tmp$raw



## -----------------------------------------------------------------------------------------------------------------------
dds_full <- estimateSizeFactors(DESeqDataSetFromMatrix(countData=round(mtrx),
                                                  colData = colData(se), design=metadata(se)$design))
rld_full <- rlog(dds_full, blind=TRUE)

layout(matrix(1:16,nrow=4,byrow=TRUE))
for(i in 1:(ncol(mtrx)-1)) {
    hist(assays(rld)[[1]][,i],n=10,main=colnames(mtrx)[i],ylim=c(0,70))
}
layout(1)
dev.copy2pdf(file="rlog_distribution_gene_ids.pdf") ## not very normal, but also not negative binomial

layout(matrix(1:16,nrow=4,byrow=TRUE))
for(i in 1:(ncol(mtrx)-1)) {
    hist(assays(rld_full)[[1]][,i],n=10,main=colnames(mtrx)[i])
}
layout(1)
dev.copy2pdf(file="rlog_distribution_full_mtrx.pdf") ## this now looks very similar to https://www.biostars.org/p/429952/ (negative binomial)


## -----------------------------------------------------------------------------------------------------------------------
p <- list()
l <- c()
k <- 1
for(i in 1:3) {
    for(j in (i+1):4) {
        cat(i,j,"\n")
        p[[k]] <- my_plotPCA(rld_full, intgroup=c("Condition"), sampleMap=as.matrix(colData(se))[,"SampleName"],
                             PCs=c(x=i,y=j),
                             hjust=0.2,vjust=-0.5,
                             abs_xmax=ifelse(k<=3,35,ifelse(k==4,25,20)),
                             abs_ymax=ifelse(k<=4,20,12))
                             
        ##main=paste0("PC_",i," vs PC_",j)))
        l[k] <- paste0("PC_",i," vs PC_",j)
        k <- k+1
    }
}
cowplot::plot_grid(p[[1]],p[[2]],p[[3]],
                   p[[4]],p[[5]],p[[6]],
                   ncol = 3, nrow = 2,
                   labels=l)

dev.copy2pdf(file="mtrxPCA_full500_PC1-4.pdf")
## system("mv mtrxPCA_full_PC1-4.pdf mtrxPCA_full500_PC1-4.pdf") ## from original filename mtrxPCA_full_PC1-4.pdf
## -----------------------------------------------------------------------------------------------------------------------
p <- list()
l <- c()
k <- 1
for(i in 1:3) {
    for(j in (i+1):4) {
        cat(i,j,"\n")
        p[[k]] <- my_plotPCA(rld_full, intgroup=c("Condition"), sampleMap=as.matrix(colData(se))[,"SampleName"],
                             PCs=c(x=i,y=j),
                             ntop=Inf,
                             hjust=0.2,vjust=-0.5)
                             
        ##main=paste0("PC_",i," vs PC_",j)))
        l[k] <- paste0("PC_",i," vs PC_",j)
        k <- k+1
    }
}
cowplot::plot_grid(p[[1]],p[[2]],p[[3]],
                   p[[4]],p[[5]],p[[6]],
                   ncol = 3, nrow = 2,
                   labels=l)

dev.copy2pdf(file="mtrxPCA_full_PC1-4.pdf")

## -----------------------------------------------------------------------------------------------------------------------
mtrx_filterByExpr <- mtrx_full.save[filterByExpr(DGEList(mtrx_full.save), metadata(se)$design),]
dim(mtrx_filterByExpr)
##[1] 17946    17
dds_full_filterByExpr <- estimateSizeFactors(DESeqDataSetFromMatrix(countData=round(mtrx_filterByExpr),
                                                                    colData = colData(se), design=metadata(se)$design))
rld_full_filterByExpr <- rlog(dds_full_filterByExpr, blind=TRUE)

## -----------------------------------------------------------------------------------------------------------------------
p <- list()
l <- c()
k <- 1
for(i in 1:3) {
    for(j in (i+1):4) {
        cat(i,j,"\n")
        p[[k]] <- my_plotPCA(rld_full_filterByExpr, intgroup=c("Condition"), sampleMap=as.matrix(colData(se))[,"SampleName"],
                             PCs=c(x=i,y=j),
                             ntop=Inf,
                             hjust=0.2,vjust=-0.5)
                             
        ##main=paste0("PC_",i," vs PC_",j)))
        l[k] <- paste0("PC_",i," vs PC_",j)
        k <- k+1
    }
}
cowplot::plot_grid(p[[1]],p[[2]],p[[3]],
                   p[[4]],p[[5]],p[[6]],
                   ncol = 3, nrow = 2,
                   labels=l)
dev.copy2pdf(file="mtrxPCA_full_filterByExpr_PC1-4.pdf")

my_plotPCA(rld_full_filterByExpr, intgroup=c("Condition"), sampleMap=as.matrix(colData(se))[,"SampleName"],ntop=Inf,returnScree=TRUE)
dev.copy2pdf(file="mtrxScree_full_filterByExpr_PC1-4.pdf")

# .......................................................................................................................
cnd <- tapply(rownames(colData(se)),colData(se)$Condition,c)

p <- list()
for(cond in names(cnd)) {
    p[[cond]] <- my_plotPCA(rld_full_filterByExpr[,cnd[[cond]]], intgroup=c("Condition"), sampleMap=as.matrix(colData(se))[,"SampleName"],
                             PCs=c(x=1,y=2),
                             ntop=Inf,
                             hjust=0.2,vjust=-0.5)
}
    
cowplot::plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],
                   ncol = 2, nrow = 2,
                   labels=names(cond))
dev.copy2pdf(file="mtrxPCA_full_filterByExpr_PC1-2_byCondition.pdf")

# -----------------------------------------------------------------------------------------------------------------------
p <- list()
l <- c()
k <- 1
for(i in 1:3) {
    for(j in (i+1):4) {
        cat(i,j,"\n")
        p[[k]] <- my_plotPCA(rld_full[gene_ids,], intgroup=c("Condition"), sampleMap=as.matrix(colData(se))[,"SampleName"],
                             PCs=c(x=i,y=j),
                             hjust=0.2,vjust=-0.5,
                             abs_xmax=ifelse(k<=3,15,ifelse(k<6,8,5)),
                             abs_ymax=ifelse(k<=3,10,5))
                             
        ##main=paste0("PC_",i," vs PC_",j)))
        l[k] <- paste0("PC_",i," vs PC_",j)
        k <- k+1
    }
}
cowplot::plot_grid(p[[1]],p[[2]],p[[3]],
                   p[[4]],p[[5]],p[[6]],
                   ncol = 3, nrow = 2,
                   labels=l)

dev.copy2pdf(file="mtrxPCA_full_gene_ids_PC1-4.pdf")
## -----------------------------------------------------------------------------------------------------------------------
p <- list()
l <- c()
k <- 1
for(i in 1:3) {
    for(j in (i+1):4) {
        cat(i,j,"\n")
        p[[k]] <- my_plotPCA(rld_full_filterByExpr[gene_ids,], intgroup=c("Condition"), sampleMap=as.matrix(colData(se))[,"SampleName"],
                             PCs=c(x=i,y=j),
                             hjust=0.2,vjust=-0.5,
                             abs_xmax=ifelse(k<=3,15,ifelse(k<6,8,5)),
                             abs_ymax=ifelse(k<=3,10,5))
                             
        ##main=paste0("PC_",i," vs PC_",j)))
        l[k] <- paste0("PC_",i," vs PC_",j)
        k <- k+1
    }
}
cowplot::plot_grid(p[[1]],p[[2]],p[[3]],
                   p[[4]],p[[5]],p[[6]],
                   ncol = 3, nrow = 2,
                   labels=l)

dev.copy2pdf(file="mtrxPCA_full_filterByExpr_gene_ids_PC1-4.pdf")

my_plotPCA(rld_full_filterByExpr[gene_ids,], intgroup=c("Condition"), sampleMap=as.matrix(colData(se))[,"SampleName"],returnScree=TRUE)
dev.copy2pdf(file="mtrxScree_full_filterByExpr_gene_ids_PC1-4.pdf")

# .......................................................................................................................
cnd <- tapply(rownames(colData(se)),colData(se)$Condition,c)

p <- list()
for(cond in names(cnd)) {
    p[[cond]] <- my_plotPCA(rld_full_filterByExpr[gene_ids,cnd[[cond]]], intgroup=c("Condition"), sampleMap=as.matrix(colData(se))[,"SampleName"],
                             PCs=c(x=1,y=2),
                             ntop=Inf,
                             hjust=0.2,vjust=-0.5)
}
    
cowplot::plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],
                   ncol = 2, nrow = 2,
                   labels=names(cond))
dev.copy2pdf(file="mtrxPCA_full_filterByExpr_gene_ids_PC1-2_byCondition.pdf")

l <- list()
for(cond in names(cnd)) {
    l[[cond]] <- my_plotPCA(rld_full_filterByExpr[gene_ids,cnd[[cond]]], intgroup=c("Condition"), sampleMap=as.matrix(colData(se))[,"SampleName"],returnData=TRUE)
}

## -----------------------------------------------------------------------------------------------------------------------
##library("rrcov")
m <- t(assays(rld_full_filterByExpr)[[1]][rownames(assays(rld_full_filterByExpr)[[1]]) %in% gene_ids,])
rownames(m) <-  as.matrix(colData(se))[,"SampleName"][rownames(m)]
m_full <- t(assays(rld_full_filterByExpr)[[1]])
rownames(m_full) <-  as.matrix(colData(se))[,"SampleName"][rownames(m_full)]

## Using the filtered, rlog-normalized dataset, restricted to the GO term genes
pc <- rrcov::PcaGrid(m[1:4,],k=2) ## cre
names(which(!pc$flag))
##character(0)

pc <- rrcov::PcaGrid(m[5:8,],k=2) ## iko
names(which(!pc$flag))
##[1] "ALU_720"

pc <- rrcov::PcaGrid(m[9:12,],k=2) ## p53
names(which(!pc$flag))
##[1] "ANK_2492"
pc <- rrcov::PcaGrid(m[13:17,],k=2) ## dko
names(which(!pc$flag))
##character(0)

## Using the complete filtered, rlog-normalized dataset
pc <- rrcov::PcaGrid(m_full[1:4,],k=2) ## cre
names(which(!pc$flag))
##character(0)

pc <- rrcov::PcaGrid(m_full[5:8,],k=2) ## iko
names(which(!pc$flag))
##[1] "ALU_720"

pc <- rrcov::PcaGrid(m_full[9:12,],k=2) ## p53
names(which(!pc$flag))
##[1] "ANK_2492"
pc <- rrcov::PcaGrid(m_full[13:17,],k=2) ## dko
names(which(!pc$flag))
##[1] "ANK_275-1"








##
pc <- rrcov::PcaGrid(m)
table(pc$flag)
##TRUE 
##  17

pc_4 <- rrcov::PcaGrid(m,k=4)
table(pc_4$flag)
##TRUE 
##  17 

pc_full <- rrcov::PcaGrid(m_full)
##FALSE  TRUE 
##    1    16 
names(which(!pc_full$flag))
##[1] "ALU_720"

pc_full_4 <- rrcov::PcaGrid(m_full,k=4)
table(pc_full_4$flag)
##FALSE  TRUE 
##    2    15 
names(which(!pc_full_4$flag))
##[1] "ALU_715" "ALU_720"

pcH <- rrcov::PcaHubert(m)
table(pcH$flag)
##TRUE 
##  17 

pcH_4 <- rrcov::PcaHubert(m,k=4)
table(pcH_4$flag)
##FALSE  TRUE 
##    2    15 
names(which(!pcH_4$flag))
##[1] "ALU_720"   "ANK_275-1"

pcH_full <- rrcov::PcaHubert(m_full)
table(pcH_full$flag)
##FALSE  TRUE 
##    4    13 
names(which(!pcH_full$flag))
##[1] "ALU_714"  "ALU_715"  "ALU_720"  "ANK_2492"

pcH_full_4 <- rrcov::PcaHubert(m_full,k=4)
table(pcH_full_4$flag)
##FALSE  TRUE 
##    3    14 
names(which(!pcH_full_4$flag))
##"ALU_714" "ALU_715" "ALU_720"

## -----------------------------------------------------------------------------------------------------------------------
## =======================================================================================================================
## scale and subset only after the PCAs!
mtrx_full.save <- mtrx

## Scale the collapsed matrix:
scalefunc_name   <- "CPM"
##scalefunc_name   <- "CPM_voom"
mtrx <-scale_functions[[scalefunc_name]](mtrx, ## first scale as a whole ...
                                         do_log=TRUE,pseudocount=0.1) [rownames(mtrx)[rownames(mtrx) %in% gene_ids],] ## .. then subset


## -----------------------------------------------------------------------------------------------------------------------
load("groups_941.RData") ## was saved "Thu Sep 30 15:45:36 2021"
## -----------------------------------------------------------------------------------------------------------------------

## Explore ways to make a scaled matrix without duplicate genes:
intersections <- intersection_groups(GO_groups)
d <- as.data.frame(intersections[-which(colnames(intersections)=="set")]) ; d
##   GO.0070498 GO.0038061 GO.0033209 GO.0019221 GO.0071357 GO.0050727 GO.0002479 GO.0043312 GO.0071456 GO.0036294 setSize
##1           0          0          0          0          0          0          0          1          0          0      75
##2           0          0          0          1          0          0          0          0          0          0      71
##3           0          0          0          1          1          0          0          0          0          0      26
##4           0          0          0          0          0          1          0          0          0          0      25
##5           1          1          1          1          0          0          1          0          1          0      15
##6           0          0          1          1          0          0          0          0          0          0      12
##7           0          0          0          0          0          0          0          0          1          1      10
##8           0          0          0          1          0          1          0          0          0          0       9
##9           0          0          0          0          0          0          0          0          1          0       8
##10          1          1          1          1          0          0          1          1          1          0       7
##11          1          0          0          1          0          0          0          0          0          0       6
##12          0          0          0          1          0          0          0          1          0          0       6
##13          0          0          0          0          0          1          0          1          0          0       4
##14          0          1          0          0          0          0          0          0          0          0       3
##15          0          0          0          0          0          0          1          0          0          0       3
##16          0          1          1          1          0          0          0          0          0          0       2
##17          0          0          0          0          0          0          1          1          0          0       2
##18          0          0          0          1          0          0          1          1          0          0       2
##19          0          0          0          1          0          0          0          0          1          1       2
##20          0          0          1          1          1          0          0          0          0          0       1
##21          1          0          0          1          0          1          0          0          0          0       1
##22          0          0          1          1          0          1          0          0          0          0       1
##23          1          1          1          1          0          1          0          0          0          0       1
##24          1          0          0          1          1          1          0          0          0          0       1
##25          0          0          1          1          0          0          0          1          0          0       1
##26          0          0          0          1          1          0          0          1          0          0       1
##27          0          0          0          1          1          0          1          1          0          0       1
##28          0          0          0          1          0          0          0          0          1          0       1
##29          1          1          1          1          0          0          0          0          1          0       1
##30          1          1          1          1          1          0          1          0          1          0       1
##31          1          1          1          1          0          1          1          0          1          0       1

## -----------------------------------------------------------------------------------------------------------------------
library("ComplexHeatmap")
source("~/CECAD/Pipeline/Git/Heatmaps/Ali_GO.Heatmap.Code.And.Example/My_adaptations/functions_UG.R")
load("groups_941.RData")

### Alternative samples labels.
sample.labels <-  colData(se)$SampleName

conditions <- colData(se)$Condition
names(conditions) <- colnames(mtrx)
condition_colors <- c(cre="grey", iko="yellow",p53="darkblue",dko="darkgreen")  

go_terms_colors_941        <- colors(length(groups_941))
names(go_terms_colors_941) <- names(groups_941)

go_terms_colors_by_GO        <- colors(length(GO_groups))
names(go_terms_colors_by_GO) <- names(GO_groups)

## -----------------------------------------------------------------------------------------------------------------------
## pheatmaps code, from https://www.biostars.org/p/223532/
scale_rows  <-  function(x){
    m  <-  apply(x, 1, mean, na.rm = TRUE)
    s  <-  apply(x, 1, sd, na.rm = TRUE)
    return((x - m) / s)
}

scale_mat  <-  function(mat, scale){
    if(!(scale %in% c("none", "row"))){
        stop("scale argument shoud take values: 'none' or 'row'")
    }
    mat  <-  switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
    return(mat)
}

#####################colnames(mtrx) <-  colData(se)[,"nice_colnames"]

## Pre-cluster mtrx by group
mtrx_by_groups_941 <- sapply(c("none","row"),
                             function(scale) {
                                 Reduce(rbind,
                                        sapply(groups_941,
                                               function(set) {
                                                   m <- scale_mat(mtrx[set,],scale)
                                                   m[hclust(dist(m))$order,]
                                               },simplify=FALSE)
                                        )
                             },simplify=FALSE)

chm <- list()
for(scale in c("none","row")) {
    chm[[scale]] <- heatmap_UG(mtrx_by_groups_941[[scale]], groups_941,
                               conditions,
                               go_terms_colors = go_terms_colors_941, condition_colors = condition_colors,sample.labels = sample.labels,
                               how_to_scale="none", ### scale, !! scaling is already done!!
                               border=NA,
                               fontsize_row=2,
                               matrix_is_ordered=TRUE,
                               returned.object="ComplexHeatmap") 
}
chm[["none"]]
dev.copy2pdf(file="complexHeatmap_no_scaling.pdf")
##dev.copy2pdf(file="complexHeatmap_no_scaling_CPM_voom.pdf")
chm[["row"]]
dev.copy2pdf(file="complexHeatmap_row_scaling.pdf")
##dev.copy2pdf(file="complexHeatmap_row_scaling_CPM_voom.pdf")


## -----------------------------------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------
## Try how a simple stacking of "GO_mtrx"s looks like if they are pre-scaled and pre-clustered:

mtrx_by_GO <- sapply(c("none","row"),
                     function(scale) {
                         Reduce(rbind,
                                sapply(GO_groups,
                                       function(set) {
                                           m <- scale_mat(mtrx[set,],scale)
                                           m[hclust(dist(m))$order,]
                                       },simplify=FALSE)
                                )
                     },simplify=FALSE)

go_terms_colors_by_GO        <- colors(length(GO_groups))
names(go_terms_colors_by_GO) <- names(GO_groups)

chm_by_GO <- list()
for(scale in c("none","row")) {
    chm_by_GO[[scale]] <- heatmap_UG(mtrx_by_GO[[scale]], GO_groups,
                                     conditions,
                                     go_terms_colors = go_terms_colors_by_GO,
                                     condition_colors = condition_colors,sample.labels = sample.labels,
                                     how_to_scale="none", ## scale, !! pre-scaled !! (added Sep 16, 2021)
                                     border=NA,
                                     fontsize_row=2,
                                     pre_blown=TRUE,
                                     matrix_is_ordered=TRUE,
                                     returned.object="ComplexHeatmap") 
}
chm_by_GO[["none"]]
dev.copy2pdf(file="complexHeatmap_byGO_no_scaling.pdf")
##dev.copy2pdf(file="complexHeatmap_byGO_no_scaling_CPM_voom.pdf")
chm_by_GO[["row"]]
dev.copy2pdf(file="complexHeatmap_byGO_row_scaling.pdf")
##dev.copy2pdf(file="complexHeatmap_byGO_row_scaling_CPM_voom.pdf")



## -----------------------------------------------------------------------------------------------------------------------
## Sep 30, 2021: filter and scale the input matrix as in the limma/voom workflow:

matrix_for_limma <- function(cnt,design,use.cpm="edgeR",prior=5,do.filter=TRUE,
                             normalize.method="none") ## "none" is the default in voom() -> *no* further normalization
{
    y <- DGEList(cnt) 

    if(do.filter) { ## this is how it was done in the limma workflow -- but here I finally only use the genes
                    ## which have been selected by that workflow, so I could (for comparison with my workflow above)
                    ## take cpm() without filtering at this stage
        
        keep <- filterByExpr(y, design)
        y <- y[keep,]
    }
    y <- calcNormFactors(y)

    if        (use.cpm=="edgeR") {
        edgeR::cpm(y, log=TRUE, prior.count=prior) ## https://support.bioconductor.org/p/100731/ (recommends 5)
                                               ## see on the use of prior.count in cpm(): https://support.bioconductor.org/p/107719/
                                               ## https://rdrr.io/bioc/edgeR/src/R/cpm.R
                                               ## https://rdrr.io/bioc/edgeR/src/R/makeCompressedMatrix.R
    } else if (use.cpm=="voom") {
        ## as in function "voom" of packageVersion("limma") = 3.46.0
        lib.size <- y$samples$lib.size * y$samples$norm.factors

        normalizeBetweenArrays(t(log2(t(y$counts + 0.5)/(lib.size + 1) * 1e+06)),
                               method=normalize.method) ## NO additional scaling, unless a different "normalize.method" is specified
    } else {
        stop("Unknown scaling method in matrix_for_limma() !\n")
    }
    
}
prior_count <- 0.1 ##0.5
do_filter <- FALSE# TRUE ##FALSE
##use_cpm <- "edgeR"
use_cpm <- "voom"

###### limma_m <- matrix_for_limma(assays(se)$counts,metadata(se)$design,
######                             do.filter= do_filter,
######                             use.cpm= use_cpm,
######                             prior= prior_count) 
###### 
###### ## Rename and collapse rows with ambiguous name mapping --
###### ## NOTE: log cpm scaling has already been done here. logCPMs are not additive (adding would amount to multiplying on the original scale),
###### ##       so take the *average* of rows mapping to the same symbol, instead of the colsum. 
###### tmp <- remap_rownames(limma_m,  ## if rows with amiguous name mapping must be collapsed, do it before scaling!
######                       gene_id_map,
######                       merge_function=colMeans)
###### limma_mtrx  <-        tmp$tbl
###### limma_mtrx_mapping <- tmp$grp
###### 
###### limma_mtrx <- limma_mtrx[rownames(limma_mtrx) %in% gene_ids,]
###### 
###### colnames(limma_mtrx) <-  colData(se)[,"nice_colnames"]

## paramters for use in PCA
prior_count <- 0.5; do_filter <- TRUE; use_cpm <- "voom"

## alternative: run the limma workflow on the pre-mapped/collapsed counts
limma_mtrx <- matrix_for_limma(remap_rownames(assays(se)$counts,
                                              gene_id_map,
                                              merge_function=colSums)$tbl,
                               metadata(se)$design,
                               do.filter=do_filter,
                               prior= prior_count)
limma_mtrx_full.save <- limma_mtrx

limma_mtrx <- limma_mtrx[rownames(limma_mtrx) %in% gene_ids,]

## Use a subset of a pre-normalized full dataset, like for the PCAs


## -----------------------------------------------------------------------------------------------------------------------
layout(matrix(1:16,nrow=4,byrow=TRUE))
for(i in 1:(ncol(limma_mtrx)-1)) {
    hist(limma_mtrx[,i],n=10,main=colnames(limma_mtrx)[i],ylim=c(0,70))
}
layout(1)
dev.copy2pdf(file="voom_distribution_gene_ids.pdf") ## not very normal, even more skewed than the rlog distribution

layout(matrix(1:16,nrow=4,byrow=TRUE))
for(i in 1:(ncol(limma_mtrx)-1)) {
    hist(limma_mtrx_full.save[,i],n=10,main=colnames(limma_mtrx_full.save)[i])
}
layout(1)
dev.copy2pdf(file="voom_distribution_full_mtrx.pdf") ## sort of normal (but skewed towards higher counts)


## -----------------------------------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------
## cleaned code 13.10.2021
## .......................................................................................................................
prior_count <- 0.5; do_filter <- TRUE; use_cpm <- "voom"

limma_mtrx <- matrix_for_limma(remap_rownames(assays(se)$counts,
                                              gene_id_map,
                                              merge_function=colSums)$tbl,
                               metadata(se)$design,
                               do.filter=do_filter,
                               prior= prior_count)
limma_mtrx_full.save <- limma_mtrx

limma_mtrx <- limma_mtrx[rownames(limma_mtrx) %in% gene_ids,]
## .......................................................................................................................
rlog_mtrx <- assays(rld_full_filterByExpr)[[1]][rownames(assays(rld_full_filterByExpr)[[1]]) %in% gene_ids,]
## .......................................................................................................................
## -----------------------------------------------------------------------------------------------------------------------
GOs_to_leave_out <- c("GO.0071456","GO.0036294")
## "groups_941" which contain at least one of the GO terms that should be left out:
groups_to_leave_out <- names(which(sapply(names(groups_941),function(n) grepl(paste(GOs_to_leave_out,collapse="|"),n))))

## -----------------------------------------------------------------------------------------------------------------------
source("~/CECAD/Pipeline/Git/Heatmaps/Ali_GO.Heatmap.Code.And.Example/My_adaptations/functions_UG.R")
load("groups_941.RData")

## for whatever reason, I had a "." here instead of ":" (necessary only for the non-redundant groups, which have "." as a separator)
GOs_to_leave_out <- c("GO:0071456","GO:0036294") 

### Alternative samples labels.
sample.labels <-  colData(se)$SampleName

conditions <- colData(se)$Condition
names(conditions) <- colnames(mtrx)
condition_colors <- c(cre="grey", iko="yellow",p53="darkblue",dko="darkgreen")  

m_vars <- c(log2CPM="mtrx", voom="limma_mtrx",rlog="rlog_mtrx")
group_vars <- c(redundant="GO_groups",nonredundant="groups_941")
leaveout_vars <- c(redundant="GOs_to_leave_out", nonredundant="groups_to_leave_out")
col_vars <- c(redundant="go_terms_colors_by_GO", nonredundant="go_terms_colors_941")

for(grouping in names(group_vars)) {
    g <- eval(as.symbol(group_vars[grouping]))
    l <- eval(as.symbol(leaveout_vars[grouping]))
    col <- eval(as.symbol(col_vars[grouping]))
              
              
    for(normalization in names(m_vars)) {
        M <- eval(as.symbol(m_vars[normalization]))
        
        M_by_group <- sapply(c("full","main","rest"),
                             function(how) {
                                 if(how=="full") {
                                     go <- g
                                 } else if(how=="partial") {
                                     go <- g[-which(names(g) %in% l)]
                                 } else if(how=="rest") {
                                     go <- g[ which(names(g) %in% l)]
                                 }
                                 go_col <- colors(length(go))
                                 names(go_col) <-  names(go)
                                 
                                 sapply(c("none","row"),
                                        function(scale) {
                                            cat(grouping,normalization,how,scale,"\n")
                                            ##if(how=="rest") browser()
                                            
                                            m <- Reduce(rbind,
                                                        sapply(go,
                                                               function(set) {
                                                                   m <- scale_mat(M[set,],scale)
                                                                   m[hclust(dist(m))$order,] ## this clusters *rows*
                                                               },simplify=FALSE))
                                            p <- heatmap_UG(m, go,
                                                            conditions,
                                                            go_terms_colors = go_col,
                                                            condition_colors = condition_colors,
                                                            sample.labels = sample.labels,
                                                            how_to_scale="none", ### !! scaling is already done!!
                                                            pre_blown=TRUE,
                                                            border=NA,
                                                            put_legend=TRUE,
                                                            fontsize_row=ifelse(how=="rest",8,2),
                                                            matrix_is_ordered=TRUE,
                                                            returned.object="ComplexHeatmap")  ## clustering of *columns* is done here
                                            print(p)
                                            
                                            ##browser()
                                            dev.copy2pdf(file=paste0(normalization,"_",
                                                                     grouping,"_how.",how,"_scale.",scale,".pdf"))
                                            
                                        },simplify=FALSE)
                             },simplify=FALSE)
    }
}

## -----------------------------------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------
