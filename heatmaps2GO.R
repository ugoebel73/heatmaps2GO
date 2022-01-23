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
GO_info <- as.list(GO.db::GOTERM)
names(GO_info) <- sapply(GO_info,function(x)x@GOID)

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
## (this is the group 3 of groups_464, consisting of proteasome genes only)
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

g <- biomaRt_v102_genes2GO %>% filter(ENSEMBL %in% 
                                      enrichr_proteasome_mapping)  %>%
                                      distinct(ENSEMBL) %>% pull(ENSEMBL)
length(g)
##[1] 25 ## all 25 genes are actually represented

table(g %in% DGE_genes) ## and they are all in our analysis set
##TRUE 
##  25 

g_symbols <- biomaRt_v102_genes2GO %>% filter(ENSEMBL %in% g) %>%
                                       distinct(SYMBOL) %>% pull(SYMBOL)
##[1] "Psmc4"  "Psme1"  "Psme2"  "Psmb3"  "Psmb4"  "Psmb5"  "Psmd12" "Psma7" 
##[9] "Psmd4"  "Psma5"  "Psmb2"  "Psmb6"  "Psmb7"  "Psma3"  "Psmb9"  "Psmf1" 
##[17] "Psmb8"  "Psma6"  "Uba52"  "Psmd3"  "Psmc1"  "Psmc2"  "Psma2"  "Psmb10"
##[25] "Psmb1" 

## But they are not associated with our_GO_terms in biomaRt or in org.Mm.eg.db,
## except for one gene:
tmp <- setNames(sapply(GO_genesets[["biomaRt_all_inclusive"]],
                       function(x)intersect(x,g)),
                GO_genesets$GO_term)
tmp[sapply(tmp,length)>0]
##$`GO:0050727`
##[1] "ENSMUSG00000005779"

tmp <- setNames(sapply(GO_genesets[["org_all"]],
                       function(x)intersect(x,g)),
                GO_genesets$GO_term)
tmp[sapply(tmp,length)>0]
##$`GO:0050727`
##[1] "ENSMUSG00000005779"


## With which terms were the 25 genes actually associated in Enrichr?
source("functions_enrichr.R")
Enrichr_dbs <- sapply(read_Enrichr_dbs(Enrichr_base="DATA"),
                      function(x) {
                        y <- x$genes
                        names(y) <- extractCapturedSubstrings(
                                       pattern="\\((GO:\\d+)\\)$",
                                       string=names(y))    
                        y 
                      })
sapply(toupper(g_symbols),function(this_g)  
  names(which(
    sapply(Enrichr_dbs[["GO_Biological_Process_2018"]][our_GO_terms],
           function(x)this_g %in% x)
    ))
)
## $PSMC4
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0071456"
## $PSME1
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0071456"
## $PSME2
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0071456"
## $PSMB3
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0071456"
## $PSMB4
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0071456"
## $PSMB5
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0071456"
## $PSMD12
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0043312" "GO:0071456"
## $PSMA7
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0071456"
## $PSMD4
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0071456"
## $PSMA5
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0043312" "GO:0071456"
## $PSMB2
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0071456"
## $PSMB6
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0071456"
## $PSMB7
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0043312" "GO:0071456"
## $PSMA3
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0071456"
## $PSMB9
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0071456"
## $PSMF1
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0071456"
## $PSMB8
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0071357" "GO:0002479" "GO:0071456"
## $PSMA6
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0050727" "GO:0002479" "GO:0071456"
## $UBA52
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0071456"
## $PSMD3
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0043312" "GO:0071456"
## $PSMC1
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0071456"
## $PSMC2
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0043312" "GO:0071456"
## $PSMA2
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0043312" "GO:0071456"
## $PSMB10
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0071456"
## $PSMB1
## [1] "GO:0070498" "GO:0038061" "GO:0033209" "GO:0019221" "GO:0002479" "GO:0043312" "GO:0071456"                        


t <- biomaRt_v102_genes2GO %>% filter(ENSEMBL %in% g) %>%
                               distinct(GO_term) %>% pull(GO_term)
length(unique(t))
##[1] 97
length(intersect(our_GO_terms,t))
##[1] 0
## the 25 genes *are* associated with 97 different GO terms, 
## but non is even a synonym of our_GO_terms:

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

## these are the biological themes 
## with which the 25 genes are associated in GO.db:

##sapply(t,function(x)GO_info[[x]]@Term)
tmp <- sapply(t,
              function(x) {
                if (GO_info[[x]]@Ontology=="BP") GO_info[[x]]@Term
                else                             NULL
             })
tmp[sapply(tmp,length)>0]
## $`GO:0030163`
## [1] "protein catabolic process"
## $`GO:1901800`
## [1] "positive regulation of proteasomal protein catabolic process"
## $`GO:0001824`
## [1] "blastocyst development"
## $`GO:0045899`
## [1] "positive regulation of RNA polymerase II transcription preinitiation complex assembly"
## $`GO:0019884`
## [1] "antigen processing and presentation of exogenous antigen"
## $`GO:0010950`
## [1] "positive regulation of endopeptidase activity"
## $`GO:0061136`
## [1] "regulation of proteasomal protein catabolic process"
## $`GO:2000045`
## [1] "regulation of G1/S transition of mitotic cell cycle"
## $`GO:0043161`
## [1] "proteasome-mediated ubiquitin-dependent protein catabolic process"
## $`GO:0006508`
## [1] "proteolysis"
## $`GO:0051603`
## [1] "proteolysis involved in cellular protein catabolic process"
## $`GO:0010498`
## [1] "proteasomal protein catabolic process"
## $`GO:0010499`
## [1] "proteasomal ubiquitin-independent protein catabolic process"
## $`GO:0002862`
## [1] "negative regulation of inflammatory response to antigenic stimulus"
## $`GO:0006979`
## [1] "response to oxidative stress"
## $`GO:0006511`
## [1] "ubiquitin-dependent protein catabolic process"
## $`GO:0043248`
## [1] "proteasome assembly"
## $`GO:0010243`
## [1] "response to organonitrogen compound"
## $`GO:0014070`
## [1] "response to organic cyclic compound"
## $`GO:0052548`
## [1] "regulation of endopeptidase activity"
## $`GO:0002376`
## [1] "immune system process"
## $`GO:2000116`
## [1] "regulation of cysteine-type endopeptidase activity"
## $`GO:1901799`
## [1] "negative regulation of proteasomal protein catabolic process"
## $`GO:0030154`
## [1] "cell differentiation"
## $`GO:0045444`
## [1] "fat cell differentiation"
## $`GO:0019882`
## [1] "antigen processing and presentation"
## $`GO:0051092`
## [1] "positive regulation of NF-kappaB transcription factor activity"
## $`GO:0007519`
## [1] "skeletal muscle tissue development"
## $`GO:0006412`
## [1] "translation"
## $`GO:0016567`
## [1] "protein ubiquitination"
## $`GO:0019941`
## [1] "modification-dependent protein catabolic process"
## $`GO:0042176`
## [1] "regulation of protein catabolic process"
## $`GO:0050790`
## [1] "regulation of catalytic activity"
## $`GO:1901215`
## [1] "negative regulation of neuron death"
## $`GO:0000902`
## [1] "cell morphogenesis"
## $`GO:0042098`
## [1] "T cell proliferation"



## likely the information transfer was via
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
length(unique(biomaRt_v102_genes2GO$ENSEMBL))
##[1] 22302 ## much less .. like due to removed empty GO terms:

tmp <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                     "external_gene_name",
                                     "go_id",
                                     "go_linkage_type"),
                      mart=mart102)
length(unique(tmp$ensembl_gene_id))
##[1] 56305
length(unique(tmp$ensembl_gene_id[tmp$go_id!=""]))
##[1] 22303 ## OK

## --------------------------------------------------------------------------
## --------------------------------------------------------------------------
## --------------------------------------------------------------------------
## prepare input for heatmaps:
indicators <- sapply(colnames(GO_genesets)[-1],
                     function(n) make_setMembership_vectors(
                                            setNames(GO_genesets[[n]],
                                                     GO_genesets[["GO_term"]])
                                            ),
                    simplify=FALSE)

geneset <- "biomaRt_dge_inclusive"
m <- as.matrix(indicators[[geneset]][1:(ncol(indicators[[geneset]])-2)])
nclus <- 8
clusters <-  cluster_setMembership_vectors(by="hamming", 
                                           v=m,nclus=nclus)
  
sapply(clusters,function(i) indicators[[geneset]][i,], simplify=FALSE)

## not very nice ...


