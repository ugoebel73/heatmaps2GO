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
make_setMembership_vectors <- function(setList) {
  
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
  if(by == "GOSemSim") {
    cluster_setMembership_byGOSemSim(v,semData, ...) 
  } else if (by =="hamming") {
    cluster_setMembership_by_Hamming(v, ...)
  } else {
    stop('Unknown clustering function in cluster_setMembership_vectors()!')
  }
  
}

#' cluster_GO_vectors
#' @param v a matrix or data.frame with columns=list elements, rows = 0/1 indicator vectors
#' @param semData (parameter of GOSemSim::mgoSim) GOSemSimDATA object
#' @param measure (parameter of GOSemSim::mgoSim) One of "Resnik", "Lin", "Rel", "Jiang", "TCSS", "Wang" 
#' @param combine (parameter of GOSemSim::mgoSim) One of "max", "avg", "rcmax", "BMA" 
#' @param nclus cut the tree of GO similarities such that nclus clusters are returned
#'
cluster_setMembership_byGOSemSim <- function(v,semData, measure = "Wang", combine = "BMA",nclus=10) {
  v <- as.matrix(v)

  ## return the colnames(assumed here: GO IDs) which have a 1 (presence) in row i 
  row2vec <- function(i,v) colnames(v)[as.logical(v[i,])]

  ## compute a similarity matrix of the GO terms associated with pairs of rows
  m <- matrix(NA,nrow=nrow(v),ncol=nrow(v)) 
  for(i in 1:nrow(v)) {
    for(j in 1:nrow(v)) m[i,j] <- mgoSim(row2vec(i,v), row2vec(j,v), semData)
  }
  if(max(m)<1) {
    stop("Similarity value >1 found in cluster_GO_vectors!\n")
  }

  ## convert to distance matrix and cut the hclust tree
  ## at a level yielding "nclus" clusters 
  tmp <- cutree(hclust(as.dist(1-m)),k=nclus)
  
  tapply(1:length(tmp),tmp,c)
}

hammingLike_matrix_score <- function(m,exponent=4, zero_penalty=0.01,
                                     combineScores=sum)  {
    ## "m" is a matrix with values in {0,1}
    ##     with rows to be read as presence(1)/absence(0) indicators
    ## "zero_penality" is added to all-zero columns of "m"
    ## "exponent" will be used as exponent y in score_colMean().

    ## For each column c_i of "m", the function calls score_colMean() to
    ## score the similarity of presence/absence values in c_i.
    
    score_colMean <- function(x,y=exponent) {
        if((y %% 2) > 0) {
            stop("exponent must be even in score_colMean()")
        }
        (2*x-1)^y
    }

    if(nrow(m) > 1) {
      collapse_cols <- colMeans
    } else {
      collapse_cols <- function(x) x
    }
    
    combineScores(sapply(collapse_cols(apply(m, 
                                       2, ## for each column of input matrix m:
                                       function(x) x+ (all(x==0) * zero_penalty)) ## if it is a column of shared zeroes,
                                                                                  ## add the zero_penalty
                                       ), score_colMean) ## Compute (2*[mean of (penalized) column] -1)^(exponent) :
                                                   ## The score is the closer to 1 the closer the colMean
                                                   ## is either to 1 or to 0;
                                                   ## at colMean=0.5 the score is 0;
                                                   ## multiplying all-zero columns by "zero_penalty" before scoring
                                                   ## drives their score away from 1 (= weight shared absences
                                                   ## less than shared presences)
    ) ## return a combination of the per-column scores of m (the sum by default) 

}
## to consider: is there a way to use this function as a (similarity->)distance function in a kmeans-like algorithm? 

cluster_setMembership_by_Hamming <- function(v, exponent=4, zero_penalty=0.01, nclus=8) {
    ## UPGMA clustering of the indicator vectors in v (combined as matrix rows),
    ## where pairwise mergers of clusters are scored by function hammingLike_matrix_score()
    
    v <- as.matrix(v)
    score_colMean <- function(x,y=exponent) (2*x-1)^y
    
    row_clusters <- sapply(1:nrow(v),function(i)list(i))
    
    while(length(row_clusters)>nclus) {

        this_max <- -Inf
        
        for(i in 2:length(row_clusters)) { ## get similarity scores of all pairwise mergers of the current clusters
            for(j in 1:(i-1)) {
                s <- hammingLike_matrix_score(v[unlist(row_clusters[c(j,i)]),])
                
                
                if(s > this_max) {
                    this_max <- s
                    i_max <- i
                    j_max <- j
                }
            }
        }
        row_clusters <- c(row_clusters[-c(j_max,i_max)], ## join the two clusters with highest similarity
                          list(c(row_clusters[[j_max]],
                                 row_clusters[[i_max]])
                               ))
    }
    row_clusters
}
cluster_setMembership_by_kmeans <- function(v,  
                                            nstart=1, nruns=500, 
                                            iter_max=100, 
                                            nclus=8, 
                                            min_shared_frac=0.999,
                                            max_fail_shared=1,
                                            verbose=FALSE) {
  v <- as.data.frame(v) ## this will retain the rownames 
                        ## of the extracted sub-matrices (not really needed) 
  
  sol <- sapply(1:nruns,
                function(n) {
                  ##browser()
                  if(verbose) cat(n,"\n")
                  
                  i <- kmeans(as.matrix(v),
                              nclus, nstart=nstart,
                              iter.max=iter_max)$cluster
                  clusters  <- tapply(1:length(i),i, 
                                      function(j) v[j,,drop=FALSE])
                  
                  freq <- sapply(clusters,
                                 function(this_clus) apply(this_clus,
                                                           2,
                                                           function(x)sum(x)/
                                                                   length(x)))
                  ##browser()
                  if(sum(!apply(freq,2,function(x)any(x>=min_shared_frac))) 
                     > max_fail_shared)  { 
                    ## maximally "max_fail_shared" clusters are allowed to 
                    ## *not* have at least one fixed term
                    NULL
                    
                  } else {
                    ##browser()
                    tapply(1:length(i),i, c)
                  }
                },simplify=FALSE)
  
  sol[!sapply(sol,is.null)]
}

