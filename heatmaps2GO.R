library(tidyverse)
library(org.Mm.eg.db)

source("functions.R")

org_name <- "Mus musculus"
ensembl_version <- 102

## --------------------------------------------------------------------------
## Richard Acton's solution:
extractor_function <- make_ensdb_extractor(org_name,
                                           ensembl_version)
## note that the extractor_function() remembers the parameter values
## of the calling function that were active upon its creation as promises
## (although they are not shown in the code when it is printed!)
EnsDb_v102  <- get_tx2gene_from_extractor(extractor_function())
ens_genes   <- unique(EnsDb_v102[,c("GENEID", "SYMBOL")])
## --------------------------------------------------------------------------

org_genes_plusGOALL <- AnnotationDbi::select(org.Mm.eg.db,
                                             keys=ens_genes$GENEID,
                                             columns=c("SYMBOL",
                                                       "GOALL",
                                                       "EVIDENCEALL",
                                                       "ONTOLOGYALL"),
                                             keytype="ENSEMBL")
table(is.na(unique(org_genes_plusGOALL$GOALL)))
##FALSE  TRUE 
##22763     1

## --------------------------------------------------------------------------

GO_terms <- read.csv("GO_terms.csv",comment.char="#")[,"GO_term"]



