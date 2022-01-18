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