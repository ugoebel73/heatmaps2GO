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

## -----------------------------------------------------------------------------
## Functions for transforming a matrix to be displayed as a heatmap:
## .............................................................................

collapse_matrix_rows <- function(m,rowsets,merge_function=colMeans) {
  ## rowsets is a list, 
  ##   with element value = a vector of rownames or row indices of matrix m
  ##        element name  = rowname of a new row in m 
  ##                        formed by collapsing the corresponding rows
  ## merge_function() is used to collapse a set of rows 
  ##                  contained in a given "rowsets" entry
  
  if(length(rowsets)==0) {
    return(list(collapsed=NULL,remaining=m))
    
  } else if (any(duplicated(names(rowsets)))) {
    cat('Found duplicated  "rowsets" names-- returning NULL!\n' )
    return(NULL)
  } else {
    if        (all(sapply(rowsets,class)=="character")) {
      by_name <- TRUE
    } else if (all(sapply(rowsets,class)=="numeric")) {
      by_name <- FALSE
    } else {
      cat('"rowsets" must be either all "character" or all "numeric"',
          ' -- returning NULL!\n')
    }
    
    OK  <- ( by_name && 
              (all(sapply(rowsets,function(s) all(s %in% rownames(m)))))) ||
           (!by_name && 
              (all(sapply(rowsets,function(s) all(s %in% (1:nrow(s)))))))
    if(!OK) {
      cat('All vector entries of each list element of "rowsets" must either map',
          ' to rownames or to row indices of matrix m -- returning NULL!\n')
    }
  } 
  collapsed <- t(sapply(rowsets,
                        function(s) merge_function(m[s,,drop=FALSE])))
  colnames(collapsed) <- colnames(m)
  
  ## handle any remaining, not collapsed rows of m 
  if(by_name) {
    i_remains <- which(!(rownames(m) %in% Reduce(union,rowsets)))
  } else {
    i_remains <- setdiff(1:nrow(n),Reduce(union,rowsets))
  }
  if(length(i_remains)>0) { 
    r <- m[i_remains,,drop=FALSE]
    rownames(r) <- rownames(m)[i_remains] 
    colnames(r) <- colnames(m)
  } else {
    r <- NULL
  }
  
  ## return the collapsed and the "remaining", uncollapsed part of the input matrix as separate pieces:
  list(collapsed=collapsed, remaining=r)
  
}

remap_rownames <- function(m, new2old, merge_function=colMeans,
                           remove_NA_rownames=TRUE ) {
  ## This function maps the original rownames of matrix m to new rownames.
  ## Parameter new2old is a tibble with columns "new" and "old", 
  ## describing the mapping (which is NOT necessarily 1:1).
  ## 
  ## The function handles the case where >1 old name maps to the same new name:
  ## Each such group of "clashing" matrix rows is merged to yield 
  ## a single new row (using function merge_function()),
  ## and the original rows are removed.
  
  ## --- is.null(new2old) is taken to mean "don't map":
  if(is.null(new2old)) {
    return(m)
  }
  
  ## --- First make sure that all rownames of m do appear as "old" names:
  if(!all(rownames(m) %in% new2old$old)) {
    cat("** The new2old mapping must cover all rownames of matrix m ! **")
    return(m)
  }
  ## remove any "old" names which are not represented in m
  new2old <- new2old %>% 
               filter(old %in% rownames(m)) 
  
  
  ## --- Then do the actual mapping:
  
  ## group old names by new name:
  grp <- tapply(new2old$old, new2old$new, c, 
                simplify=FALSE)
  
  ## collapse groups of rows mapping to a shared new rowname:
  res <- collapse_matrix_rows(m,
                              grp[sapply(grp, length) > 1],
                              merge_function)
  
  ## deal with the remaining rows 
  if(!is.null(res$remaining)) {
    ## apply the mapping 
    ## (which is now necessarily 1:1 **except for possible NAs in "new"**) 
    ## to this part of the matrix:
    new2old <- new2old %>% filter(old %in% rownames(res$remaining))
    
    mp        <- new2old$new
    names(mp) <- new2old$old
    
    if(any(is.na(mp))) {
      mp[which(is.na(mp))] <- "NA" ## this ASSUMEs that "NA" is not a valid new name!
    }
    rownames(res$remaining) <- mp[rownames(res$remaining)] 
    ## here, rownames with missing "new" names are still NA (not "NA")!
    
    if(remove_NA_rownames) {
      res$remaining <- res$remaining[!is.na(rownames(res$remaining)),]
    }
    
    ## append the 1:1 mapping (with collapsed NAs as group "NA") 
    ## to the collapsed groups, and return grp, 
    ## allowing to trace the mapping back in the calling environment:
    grp <- c(grp[sapply(grp, length) > 1],
             tapply(names(mp),mp,c))
    
  }
  ## return the collapsed and the original part of the input matrix separately,
  ## allowing the calling function to handle them separately:
  return(list(tbl=rbind(res$collapsed,res$remaining),grp=grp))
  
}

