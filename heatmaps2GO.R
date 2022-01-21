library(org.Mm.eg.db)
##library(tidyverse)
library(dplyr)

source("functions.R")

org_name <- "Mus musculus"
ensembl_version <- 102
## --------------------------------------------------------------------------
## biomaRt as a source for Ensembl gene and GO data 
## (alternative to the NCBI based OrgDb package org.Mm.eg.db):
mart102 <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                            dataset="mmusculus_gene_ensembl",
                            host=biomaRt::listEnsemblArchives() %>% 
                              
                              filter(version=="102") %>% 
                              pull(url))

## --------------------------------------------------------------------------

## currently valid GO terms from GO.db:
GO_terms <- DBI::dbGetQuery(GO.db::GO_dbconn(), 
                            "SELECT go_id as GO_term, 
                             ontology as Ontology FROM go_term") 
print(GO.db::GO.db) ## later also used to get parent/child relations
##GODb object:
#  | GOSOURCENAME: Gene Ontology
#| GOSOURCEURL: http://current.geneontology.org/ontology/go-basic.obo
#| GOSOURCEDATE: 2021-09-01
#| Db type: GODb
#| package: AnnotationDbi
#| DBSCHEMA: GO_DB
#| GOEGSOURCEDATE: 2021-Sep13
#| GOEGSOURCENAME: Entrez Gene
#| GOEGSOURCEURL: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
#| DBSCHEMAVERSION: 2.1

## --------------------------------------------------------------------------
## GO evidence code definitions:
path <- "DATA/GO_evidence_codes.tsv"
GO_evidence <- read.csv(path,sep="\t",comment.char="#")

## --------------------------------------------------------------------------
## GO terms for which we want to draw heatmaps:
our_GO_terms <- read.csv("GO_terms.csv",comment.char="#")[,"GO_term"]

## --------------------------------------------------------------------------
## the DGE result on which the heatmaps will be based:
path <- "DATA/DGE_plusVOOM_limma_CONTRASTS_FROM_MEANS_MODEL_all_Dec21_2021.xlsx"
DGE_results <- sapply(readxl::excel_sheets(path),
                      function(n) readxl::read_xlsx(path,sheet=n),
                      simplify=FALSE)

## --------------------------------------------------------------------------
## the v100/102 Ensembl IDs actually used in the DGE analysis:
DGE_genes <- DGE_results[[1]]$ensembl_geneid ## same for all comparisons 
                                             ## -> use [[1]]
## --------------------------------------------------------------------------
## get an Ensembl ID (from Ens102 through biomaRt) <-> Gene Ontology relation 
## from the Mouse OrgDb package
## (note that the OrgDb genome builds are from NCBI, not Ensembl,
## and that moreover ENSOURCEDATE: 2021-Apr13 suggests that the Ensembl data
## in the object are > v102 and likely from GRCm39, not GRCm38 as v102!)
org_genes2GO <- AnnotationDbi::select(org.Mm.eg.db,
                                      keys=biomaRt::getBM(
                                        attributes = "ensembl_gene_id",  
                                        mart=mart102)[,1],
                                      columns=c("SYMBOL",
                                                "GOALL",
                                                "EVIDENCEALL",
                                                "ONTOLOGYALL"),
                                      keytype="ENSEMBL")                   %>%
  
                     filter(!is.na(GOALL) & (GOALL %in% GO_terms$GO_term)) %>%
                     rename(Evidence=EVIDENCEALL)                          %>%
                     rename(GO_term=GOALL)                                 %>%
                     rename(Ontology=ONTOLOGYALL)             

## add an indicator column for whether or not a gene is in DGE_genes:
org_genes2GO$in_DGE <- org_genes2GO$ENSEMBL %in% DGE_genes
## --------------------------------------------------------------------------
## get an Ensembl ID <-> Gene Ontology relation from biomRt:
biomaRt_v102_genes2GO <-    biomaRt::getBM(attributes = c("ensembl_gene_id",
                                                          "external_gene_name",
                                                          "go_id",
                                                          "go_linkage_type"),
                                        mart=mart102)                 %>%
                            
                                filter(go_id %in% GO_terms$GO_term)   %>% 
                                rename(ENSEMBL=ensembl_gene_id)       %>%
                                rename(SYMBOL=external_gene_name)     %>%
                                rename(Evidence=go_linkage_type)      %>%
                                rename(GO_term=go_id)                 %>%
                            left_join(GO_terms,by="GO_term") ## join ontology

biomaRt_v102_genes2GO$in_DGE <- biomaRt_v102_genes2GO$ENSEMBL %in% DGE_genes

## Why biomaRt?
## pro: - version can be selected (while org.Mm.eg.de is current version only)
##      - contains GO data (missing in the equivalent EnsDb package
##        from AnnotationHub, which would be native Ensembl like biommaRt) 
## con: - GO data seem to be a GO slim (at least the Ensembl REST API 
##        returns goslim data when queried for GO)
##      - not entirely clear how the gene <-> GO term relation is defined:
##        In compliance with the Ensembl model 
##        (see https://support.bioconductor.org/p/122315/) the terms seem to
##        directly associated with genes (<- transcripts), the relation not
##        being propagated to children of a term. I infer this from the fact
##        that when the genes associated with the children are manually 
##        added, then the number of genes annotated to  "our_GO_terms"
##        becomes similar to the number from org.Mm.eg.db, 
##        while otherwise it is much smaller for most terms.
## --------------------------------------------------------------------------
## "genes per GO term" lists for our_GO_terms by org.Mm.eq.db and by biomaRt:

GO_genesets <- sapply(c("org", "biomaRt_v102"),
                      function(var) {
                        var <- paste0(var, "_genes2GO")
                        lapply(list(all=TRUE,
                                    dge=eval(as.symbol(var))$in_DGE),
                                function(cnd) {
                                    eval(as.symbol(var))[cnd,]             %>% 
                                        filter(GO_term %in% our_GO_terms)  %>% 
                                        select(-in_DGE)                    %>% 
                                        nest_by(GO_term,
                                                Ontology,.key="gene_set")
                       })
               },simplify = FALSE) ## sapply sets the names -> simpler

## note that nest_by sorts by GO term, 
## so order is different "our_GO_terms"!

## --------------------------------------------------------------------------
## accessor function to extract a column of the nested tibbles 
## in GO_genesets[[]]:

get_geneset_column <- function(tag="all",
                               obj="GO_genesets$org",
                               field="ENSEMBL") {
  l <- eval(parse(text=paste0('purrr::map(',obj,'[["',
                              tag,
                              '"]]$gene_set,function(x) x %>% dplyr::pull(',
                              field,
                              '))')))
  names(l) <- GO_genesets[[tag]]$GO_term
  l
}  
## --------------------------------------------------------------------------
## Explore how gene sets defined by or.Mm.eg.de and biomaRt relate:

setcounts <- tibble(GO_term=sort(our_GO_terms),
                    org_all=lapply(get_geneset_column(), unique),
                    org_dge=lapply(get_geneset_column(tag="dge"), unique),
                    biomaRt_all=lapply(get_geneset_column(
                      obj = "GO_genesets$biomaRt_v102"), unique),
                    biomaRt_dge=lapply(get_geneset_column(
                      tag="dge",
                      obj = "GO_genesets$biomaRt_v102"),unique)
                    )
                      
                    

##...........................................................................
## augment Ensembl v102 <-> GO associations  by the Ensembl v102 genes 
## associated with the children of each of "our_GO_terms" (as defined by GO.db):

child_terms <-  sapply(our_GO_terms,
                       function(x) {
                         cat(x,"\n")
                         more_specific_GOs(x,GO2entrez=org.Mm.egGO2EG)
                       },simplify=FALSE)

setcounts$biomaRt_all_inclusive <- 
       lapply(sort(our_GO_terms),
              function(x) {
                cnd <- biomaRt_v102_genes2GO$GO_term %in% child_terms[[x]]
                unique(biomaRt_v102_genes2GO$ENSEMBL[cnd])
              })

setcounts$biomaRt_dge_inclusive <- 
       lapply(sort(our_GO_terms),
              function(x) {
                cnd <- biomaRt_v102_genes2GO$in_DGE &
                       (biomaRt_v102_genes2GO$GO_term %in% child_terms[[x]])
                unique(biomaRt_v102_genes2GO$ENSEMBL[cnd])
              })
## > setcounts
## # A tibble: 10 × 7
## GO_term    org_all     org_dge     biomaRt_all biomaRt_dge biomaRt_all_inclusive biomaRt_dge_inclusive
## <chr>      <list>      <list>      <list>      <list>      <list>                <list>               
##   1 GO:0002479 <chr [2]>   <chr [2]>   <chr [2]>   <chr [2]>   <chr [2]>             <chr [2]>            
##   2 GO:0019221 <chr [395]> <chr [320]> <chr [144]> <chr [105]> <chr [423]>           <chr [334]>          
##   3 GO:0033209 <chr [73]>  <chr [66]>  <chr [40]>  <chr [36]>  <chr [74]>            <chr [66]>           
##   4 GO:0036294 <chr [130]> <chr [96]>  <chr [3]>   <chr [3]>   <chr [152]>           <chr [132]>          
##   5 GO:0038061 <chr [116]> <chr [108]> <chr [10]>  <chr [10]>  <chr [119]>           <chr [107]>          
##   6 GO:0043312 <chr [15]>  <chr [11]>  <chr [4]>   <chr [4]>   <chr [12]>            <chr [9]>            
##   7 GO:0050727 <chr [336]> <chr [272]> <chr [85]>  <chr [76]>  <chr [347]>           <chr [261]>          
##   8 GO:0070498 <chr [23]>  <chr [22]>  <chr [18]>  <chr [17]>  <chr [24]>            <chr [23]>           
##   9 GO:0071357 <chr [50]>  <chr [44]>  <chr [3]>   <chr [3]>   <chr [45]>            <chr [38]>           
##  10 GO:0071456 <chr [105]> <chr [92]>  <chr [119]> <chr [101]> <chr [142]>           <chr [122]>          

## Note that lists with ~identical counts need not be identical (but similar):
## > length(intersect(setcounts[[3,"org_all"]][[1]], setcounts[[3,"biomaRt_all_inclusive"]][[1]]))
## [1] 62
## > length(setcounts[[3,"org_all"]][[1]])
## [1] 73
## > length(setcounts[[3,"biomaRt_all_inclusive"]][[1]])
## [1] 74


## > length(intersect(setcounts[[4,"org_all"]][[1]], setcounts[[4,"biomaRt_all_inclusive"]][[1]]))
## [1] 99
## > length(setcounts[[4,"org_all"]][[1]])
## [1] 130
## > length(setcounts[[4,"biomaRt_all_inclusive"]][[1]])
## [1] 152
##---------------------------------------------------------------------------