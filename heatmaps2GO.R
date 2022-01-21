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

GO_synonyms <- as.list(GO.db::GOSYNONYM)
names(GO_synonyms) <- sapply(GO_synonyms,function(x)x@GOID)

## can get definitions from
## as.list(GO.db::GOTERM)  <<<--- explore! 

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

GO2genes <- sapply(c("org", "biomaRt_v102"),
                  function(var) {
                    var <- paste0(var, "_genes2GO")
                    lapply(list(all=TRUE,
                                dge=eval(as.symbol(var))$in_DGE),
                            function(cnd) {
                                eval(as.symbol(var))[cnd,]              %>% 
                                     filter(GO_term %in% our_GO_terms)  %>% 
                                     select(-in_DGE)                    %>% 
                                     nest_by(GO_term,
                                             Ontology,.key="gene_set")
                       })
            },simplify = FALSE) ## sapply sets the names -> simpler

## note that nest_by sorts by GO term, 
## so the order is different from "our_GO_terms"!

## --------------------------------------------------------------------------
## accessor function to extract a column of the nested tibbles 
## in GO2genes[[]]:

get_geneset_column <- function(tag="all",
                               obj="GO2genes$org",
                               field="ENSEMBL") {
  l <- eval(parse(text=paste0('purrr::map(',obj,'[["',
                              tag,
                              '"]]$gene_set,function(x) x %>% dplyr::pull(',
                              field,
                              '))')))
  names(l) <- GO2genes[[tag]]$GO_term
  l
}  
## --------------------------------------------------------------------------
## Explore how gene sets defined by or.Mm.eg.de and biomaRt relate:

GO_genesets <- tibble(GO_term=sort(our_GO_terms),
                    org_all=lapply(get_geneset_column(), unique),
                    org_dge=lapply(get_geneset_column(tag="dge"), unique),
                    biomaRt_all=lapply(get_geneset_column(
                      obj = "GO2genes$biomaRt_v102"), unique),
                    biomaRt_dge=lapply(get_geneset_column(
                      tag="dge",
                      obj = "GO2genes$biomaRt_v102"),unique)
                    )
##...........................................................................
## augment Ensembl v102 <-> GO associations  by the Ensembl v102 genes 
## associated with the children of each of "our_GO_terms" (as defined by GO.db):

child_terms <-  sapply(our_GO_terms,
                       function(x) {
                         cat(x,"\n")
                         more_specific_GOs(x,GO2entrez=org.Mm.egGO2EG)
                       },simplify=FALSE)

for(tag in c("all","dge")) {
  GO_genesets[[paste0("biomaRt_",tag,"_inclusive")]] <-
        lapply(sort(our_GO_terms),
               function(x) {
                 ##browser()
                 cnd <- biomaRt_v102_genes2GO$GO_term %in% child_terms[[x]]
                 if(tag=="dge") {
                   cnd <- cnd &  biomaRt_v102_genes2GO$in_DGE
                 } 
                 unique(biomaRt_v102_genes2GO$ENSEMBL[cnd])
               })
}

## > GO_genesets
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

## Note that lists with ~identical counts need not be identical 
## (but they do share a non-trivial core):
## > length(intersect(GO_genesets[[3,"org_all"]][[1]], 
##                    GO_genesets[[3,"biomaRt_all_inclusive"]][[1]]))
## [1] 62
## > length(GO_genesets[[3,"org_all"]][[1]])
## [1] 73
## > length(GO_genesets[[3,"biomaRt_all_inclusive"]][[1]])
## [1] 74


## > length(intersect(GO_genesets[[4,"org_all"]][[1]], 
##                    GO_genesets[[4,"biomaRt_all_inclusive"]][[1]]))
## [1] 99
## > length(GO_genesets[[4,"org_all"]][[1]])
## [1] 130
## > length(GO_genesets[[4,"biomaRt_all_inclusive"]][[1]])
## [1] 152
## ..........................................................................
## Where are the proteasome-related genes that are in the Enrichr heatmaps,
## but no longer in the org.Mm.eg.db-based ones?
load("DATA/groups_464.RData")
proteasome_symbols <- groups_464[[3]]

biomaRt_v102_genes2GO %>% 
  filter((SYMBOL %in% proteasome_symbols) &
         (ENSEMBL %in% Reduce(union,
                              GO_genesets$biomaRt_all_inclusive)))
## a single gene (same for "all" and "dge")
## ENSMUSG00000005779  Psmb4

org_genes2GO %>% 
  filter((SYMBOL %in% proteasome_symbols) &
           (ENSEMBL %in% Reduce(union,
                                GO_genesets$org_all)))
## the same single gene

## ..........................................................................
## Where is the Enrichr annotation coming from?

load("/home/ugoebel/CECAD/Pipeline/Git/Heatmaps/heatmaps2GO_test/Reproduce/enrichr_proteasome_mapping.RData")
enrichr_proteasome_mapping <- unlist(enrichr_proteasome_mapping) ## all length 1
length(enrichr_proteasome_mapping)
##[1] 25 

enrichr_proteasome_GOs <- 
all(enrichr_proteasome_mapping %in% DGE_genes)
##[1] TRUE
load("/home/ugoebel/CECAD/Pipeline/Git/Heatmaps/heatmaps2GO_test/Reproduce/enrichr_proteasome_GOs.RData")
enrichr_proteasome_GOs <- sub("\\.",":",sub("\\s+","", enrichr_proteasome_GOs))
##1] "GO:0070498:100%" "GO:0038061:100%" "GO:0033209:100%" "GO:0019221:100%"
##[5] "GO:0071456:100%" "GO:0002479:96%"  "GO:0043312:28%"  "GO:0071357:4%"  
##[9] "GO:0050727:4%"  

all(sub(":\\d+%","",enrichr_proteasome_GOs) %in% GO_terms$GO_term)
##[1] TRUE  ## should be TRUE, because this is a subset of our_GO_terms

nrow(biomaRt_v102_genes2GO %>% filter(GO_term %in% 
                                      sub(":\\d+%","",
                                          enrichr_proteasome_GOs)) %>%
                               distinct(ENSEMBL))
##[1] 392
g <- biomaRt_v102_genes2GO %>% filter(ENSEMBL %in% 
                                      enrichr_proteasome_mapping)  %>%
                                      distinct(ENSEMBL) %>% pull(ENSEMBL)
length(g)
##[1] 25 ## all 25 genes are actually represented

table(g %in% DGE_genes) ## and they are all in our analysis set
##TRUE 
##  25 

biomaRt_v102_genes2GO %>% filter(ENSEMBL %in% g) %>%
                          distinct(SYMBOL) %>% pull(SYMBOL)
##[1] "Psmc4"  "Psme1"  "Psme2"  "Psmb3"  "Psmb4"  "Psmb5"  "Psmd12" "Psma7" 
##[9] "Psmd4"  "Psma5"  "Psmb2"  "Psmb6"  "Psmb7"  "Psma3"  "Psmb9"  "Psmf1" 
##[17] "Psmb8"  "Psma6"  "Uba52"  "Psmd3"  "Psmc1"  "Psmc2"  "Psma2"  "Psmb10"
##[25] "Psmb1" 

t <- biomaRt_v102_genes2GO %>% filter(ENSEMBL %in% g) %>%
                               distinct(GO_term) %>% pull(GO_term)
length(unique(t))
##[1] 97
length(intersect(our_GO_terms,t))
##[1] 0
## the 25 genes are associated with 97 different GO terms, 
## but not with any term of our_GO_terms

length(which(t %in% names(GO_synonyms)))
##[1] 22
sapply(GO_synonyms[intersect(t, names(GO_synonyms))],
       function(x) any(our_GO_terms %in% x@Synonym))
#GO:0016887 GO:0030163 GO:0005515 GO:0008233 GO:0005839 GO:0004175 GO:0036064 
#FALSE      FALSE      FALSE      FALSE      FALSE      FALSE      FALSE 
#GO:0003674 GO:0006511 GO:0019882 GO:0003723 GO:0007519 GO:0003735 GO:0006412 
#FALSE      FALSE      FALSE      FALSE      FALSE      FALSE      FALSE 
#GO:0005840 GO:0022625 GO:0022627 GO:0030234 GO:0050790 GO:0016020 GO:0000902 
#FALSE      FALSE      FALSE      FALSE      FALSE      FALSE      FALSE 
#GO:0042098 
#FALSE 
sapply(GO_synonyms[intersect(t, names(GO_synonyms))],
       function(x) any(our_GO_terms %in% x@Secondary))
#GO:0016887 GO:0030163 GO:0005515 GO:0008233 GO:0005839 GO:0004175 GO:0036064 
#FALSE      FALSE      FALSE      FALSE      FALSE      FALSE      FALSE 
#GO:0003674 GO:0006511 GO:0019882 GO:0003723 GO:0007519 GO:0003735 GO:0006412 
#FALSE      FALSE      FALSE      FALSE      FALSE      FALSE      FALSE 
#GO:0005840 GO:0022625 GO:0022627 GO:0030234 GO:0050790 GO:0016020 GO:0000902 
#FALSE      FALSE      FALSE      FALSE      FALSE      FALSE      FALSE 
#GO:0042098 
#FALSE 

## like the information transfer was via
## ENSEMBL from DGE (1)-> symbol (2)-> mapped to (?human) Enrichr symbol 
##                  (3)-> associated with GO term in Enrichr database
##                  (4)-> back-associated with original Ensembl ID
##                        (implicitly by matrix row of origin)
## after (1)->, the information on the original Ensembl ID is lost;
## GO terms from other genes which happen to map to an (Enrichr flavor!) symbol 
## will be transferred to the Ensembl ID.


## see also GO.db.pdf under GOTERM:
##All the obsolete GO terms are under the nodes "obsolete molecular function" (GO:0008369), "ob-
##solete cellular component" (GO id GO:0008370), and "obsolete biological process" (GO:0008371).
##Each of these GO identifiers has a group of GO identifiers as their direct children with GO terms
##that were defined by GO but are deprecated in the current build. These deprecated GO terms were
##appended by "(obsolete)" when the data package was built.
##Mappings were based on data provided by: Gene Ontology http://current.geneontology.org/ontology/go-
##  basic.obo With a date stamp from the source of: 2021-09-01
## EXPLORE THIS <<<------------


## copied from ws02:/data/public/ugoebel/Analysis/Test_Ensembl_versions
## (should be 102 through biomaRt):
load("DATA/biomart_descriptors_102.RData") ## varname: descriptors102

length(unique(descriptors_102$ensembl_gene_id))
##[1] 56305
length(unique(biomaRt_v102_genes2GO$ENSEMBL))our_G
##[1] 22302 ## much less .. like due to removed empty GO terms:

tmp <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                     "external_gene_name",
                                     "go_id",
                                     "go_linkage_type"),
                      mart=mart102)
length(unique(tmp$ensembl_gene_id))
##[1] 56305
length(unique(tmp$ensembl_gene_id[tmp$go_id!=""]))
##[1] 22303


## --------------------------------------------------------------------------
