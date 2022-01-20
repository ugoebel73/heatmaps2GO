library(org.Mm.eg.db)
##library(tidyverse)
library(dplyr)

source("functions.R")

org_name <- "Mus musculus"
ensembl_version <- 102

## --------------------------------------------------------------------------
## Richard Acton's solution for downloading an EnsDb:
extractor_function <- make_ensdb_extractor(org_name,
                                           ensembl_version)
## note that the extractor_function() remembers the parameter values
## of the calling function that were active upon its creation 
## (in the function environment? as promises?)
## (although they are not shown in the code when it is printed!)
EnsDb_v102  <- get_tx2gene_from_extractor(extractor_function())

## --------------------------------------------------------------------------
## all v102 Ensembl IDs via the EnsDb:
ensdb_v102_genes   <- unique(EnsDb_v102[,c("GENEID", "SYMBOL")])

## --------------------------------------------------------------------------
## the DGE result on which the heatmaps will be based:
path <- "DATA/DGE_plusVOOM_limma_CONTRASTS_FROM_MEANS_MODEL_all_Dec21_2021.xlsx"
DGE_results <- sapply(readxl::excel_sheets(path),
                      function(n) readxl::read_xlsx(path,sheet=n),
                      simplify=FALSE)
## --------------------------------------------------------------------------
## the v102 Ensembl IDs actually used in the DGE analysis:
DGE_genes <- DGE_results[[1]]$ensembl_geneid ## same for all comparisons 
                                             ## -> use [[1]]
## --------------------------------------------------------------------------
## get an Ensembl ID <-> Gene Ontology relation from the Mouse OrgDb package
## (note that the OrgDb genome builds are from NCBI, not Ensembl,
## and that moreover ENSOURCEDATE: 2021-Apr13 suggests that the Ensembl data
## in the object are > v102 and likely from GRCm39, not GRCm38 as v102!)
org_genes2GO <- AnnotationDbi::select(org.Mm.eg.db,
                                      keys=ensdb_v102_genes$GENEID,
                                      columns=c("SYMBOL",
                                                "GOALL",
                                                "EVIDENCEALL",
                                                "ONTOLOGYALL"),
                                      keytype="ENSEMBL")
table(is.na(unique(org_genes2GO$GOALL)))
##FALSE  TRUE 
##22763     1 ## but nearly all v102 Ensembl Ids do have a GO annotation 

## add an indicator column for whether or not a gene is in DGE_genes:
org_genes2GO$in_DGE <- org_genes2GO$ENSEMBL %in% DGE_genes
## --------------------------------------------------------------------------
## For comparison: Ensembl v102 data from biomaRt
mart102 <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                           dataset="mmusculus_gene_ensembl",
                           host=biomaRt::listEnsemblArchives() %>% 
                             dplyr::filter(version=="102") %>% 
                             pull(url))
biomaRt_v102_genes <-    biomaRt::getBM(attributes = c("ensembl_gene_id",
                                                       "external_gene_name",
                                                       "go_id"),
                                        mart=mart102)                        
## --------------------------------------------------------------------------
## GO evidence code definitions:
path <- "DATA/GO_evidence_codes.tsv"
GO_evidence <- read.csv(path,sep="\t",comment.char="#")

## --------------------------------------------------------------------------

## the GO terms of interest:
GO_terms <- read.csv("GO_terms.csv",comment.char="#")[,"GO_term"]
## .. and the associated gene lists:
GO_genesets <- sapply(list(all=TRUE,
                           dge=org_genes2GO$in_DGE),
                      function(cnd) {
                          org_genes2GO[cnd,] %>% 
                                 filter(GOALL %in% GO_terms) %>% 
                                 rename(GO_term=GOALL) %>% 
                                 rename(Ontology=ONTOLOGYALL) %>%
                                 rename(Evidence=EVIDENCEALL) %>%
                                 select(-in_DGE) %>% 
                                 nest_by(GO_term,
                                         Ontology,.key="gene_set")
                      },simplify=FALSE)
##sapply(GO_genesets[["all"]]$gene_set,function(x)length(unique(x$ENSEMBL)))
##[1]   2 395  73 130 116  15 336  23  50 105
##sapply(GO_genesets[["dge"]]$gene_set,function(x)length(unique(x$ENSEMBL)))
##[1]   2 320  66  96 108  11 272  22  44  92
## too similar??
##i <- which(org_genes_plusGOALL$GOALL %in% GO_terms)
##length(unique(org_genes_plusGOALL[i,"ENSEMBL"]))
##[1] 867 ## this is really not much

## accessor function for the columns of the nested gene info tibbles:
get_geneset_column <- function(tag="all",field="ENSEMBL") {
  l <- eval(parse(text=paste0('purrr::map(GO_genesets[["',
                              tag,
                              '"]]$gene_set,function(x) x %>% dplyr::pull(',
                              field,
                              '))')))
  names(l) <- GO_genesets[[tag]]$GO_term
  l
}  

## --------------------------------------------------------------------------



