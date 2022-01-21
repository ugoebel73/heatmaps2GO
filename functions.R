## -----------------------------------------------------------------------------
## functions by Richard Acton (slightly modified and renamed by UG)
#' make_ensdb_extractor
#' 
#' @param species Scientific name of a species
#' @param ensembl_version the version of ensemble to use
#' @param cache the cache directory to be used by AnnotationHub
#' 
make_ensdb_extractor <- function(species, ensembl_version,
                                         cache="~/.cache/R/AnnotationHub") {
    txdb <- function() {
      AH <- AnnotationHub::AnnotationHub(cache = cache)
      AnnotationHub::query(
        AH, pattern = c(species, "EnsDb", ensembl_version)
    )[[1]]
  }
  return(txdb)
}
#' get_tx2gene_from_extractor
#' 
#' @param ensDb An EnsDb object 
#' 
get_tx2gene_from_extractor <- function(ensDb) {
  ## the return value of the function returned by make_ensdb_extractor
  ## can be used as parameter 'ensDb'
  
    k <- AnnotationDbi::keys(ensDb, 
                             keytype = "TXNAME")
    tx2gene <- AnnotationDbi::select(ensDb, 
                                     k, 
                                     columns=c("GENEID", "SYMBOL"), 
                                     keytype="TXNAME"
                              ) %>% 
                              dplyr::mutate(TXNAME = paste0("transcript:", 
                                                            TXNAME))
}
## -----------------------------------------------------------------------------
## Functions for defining sets of genes or terms based on the GO DAG structure:

more_specific_genes <- function(go_term,GO2entrez=org.Mm.egGO2EG) {
  ## get the genes directly associated with the query term
  genes <- unique(GO2entrez[[go_term]]) 
  
  ## get the direct children of the query term
  ch <- as.list(GO.db::GOBPCHILDREN)[[go_term]]
  
  ## recurse for each child term, merge results
  for(g in ch[!is.na(ch)]) {
    genes <- union(genes,
                   more_specific_genes(g,GO2entrez=GO2entrez))
  }
  genes
}

less_specific_GOs <- function(go_term,GO2entrez=org.Hs.egGO2EG) {
  terms <- go_term
  ch <- setdiff(as.list(GO.db::GOBPPARENTS)[[go_term]], "all") 
  
  for(g in ch[!is.na(ch)]) {
    terms <- union(terms,
                   less_specific_GOs(g,GO2entrez=GO2entrez)) 
    
  }
  terms
}
more_specific_GOs <- function(go_term,GO2entrez=org.Hs.egGO2EG) {
  terms <- go_term
  ch <- as.list(GO.db::GOBPCHILDREN)[[go_term]]
  
  for(g in ch[!is.na(ch)]) {
    terms <- union(terms,
                   more_specific_GOs(g,GO2entrez=GO2entrez)) 
    
  }
  terms
}



## -----------------------------------------------------------------------------
## Functions for partitioning genes sets annotated to GO terms, such that
## (1) the gene sets of the partition do not overlap
## (2) the partition is such that 
##   each set has a high (if possible: maximal) number of GO terms 
##   shared by all set members,
##   and also a high (if possible: maximal) number of GO terms 
##   which are NOT assigned to any set member
##   (where shared presence should have a higher weight than shared absence)

#' make_setMembership_partition 
#' @param setList a named lists of items (of any atomic type)
#' 
make_setMembership_partition <- function(setList) {
  
  ## get all possible combinatorial combinations of the input list elements
  combinations <- sapply(1:(2^length(setList)-1),
                         function(i) as.logical(intToBits(i))[1:length(setList)])
  rownames(combinations) <- names(setList)
  
  ## for each combinatorial combination of list elements,
  ## find the items that occur 
  ## (1) in all elements of the combination
  ## (2) in no other list element
  sets <- apply(combinations,2,
                function(x) {
                  setdiff(Reduce(intersect,setList[names(x)[which( x)]]), 
                          Reduce(union    ,setList[names(x)[which(!x)]])  
                  ) 
                })
  
  storage.mode(combinations) <- "integer" ## for readability in the output
  
  tibble(data.frame(t(combinations)), ## columns = list elements;
         ## rows = 0/1 indicator vectors representing a column combination
         set=sets,                    ## the items in each combination
         setSize=sapply(sets,length)  ## the number of items in each combination
  ) %>% filter(setSize>0) %>% arrange(desc(setSize)) 
  
}



#' cluster_setMembership_vectors
#' @param v a matrix or data.frame with columns=list elements, rows = 0/1 indicator vectors
#' 
cluster_setMembership_vectors <- function(v,by="GOSemSim",...) {
  if(by=="GOSemSim") {
    cluster_GO_vectors(v,semData) 
  } else {
    stop('Currently "GOSemSim" is the only clustering method supported by cluster_setMembership_vectors()!')
  }
  
}

#' cluster_GO_vectors
#' @param v a matrix or data.frame with columns=list elements, rows = 0/1 indicator vectors
#' @param semData (parameter of GOSemSim::mgoSim) GOSemSimDATA object
#' @param measure (parameter of GOSemSim::mgoSim) One of "Resnik", "Lin", "Rel", "Jiang", "TCSS", "Wang" 
#' @param combine (parameter of GOSemSim::mgoSim) One of "max", "avg", "rcmax", "BMA" 
#' @param nclus cut the tree of GO similarities such that nclus clusters are returned
#'
cluster_GO_vectors <- function(v,semData, measure = "Wang", combine = "BMA",nclus=10) {
  v <- as.matrix(v)
  row2vec <- function(i,v) colnames(v)[as.logical(v[i,])]
  
  m <- matrix(NA,nrow=nrow(v),ncol=nrow(v))
  for(i in 1:nrow(v)) {
    for(j in 1:nrow(v)) m[i,j] <- mgoSim(row2vec(i,v), row2vec(j,v), semData)
  }
  if(max(m)<1) {
    stop("Similarity value >1 found in cluster_GO_vectors!\n")
  }
  
  tmp <- cutree(hclust(as.dist(1-m)),k=nclus)
  
  tapply(1:length(tmp),tmp,c)
}

