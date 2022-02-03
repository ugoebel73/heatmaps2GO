
map_Ensembl_versions_from_IDMapper_output <- function(IDMapper_output,
                                                      target_release="104",
                                                      query_IDs) {

    ## parse an output file created by and downloaded from
    ## https://www.ensembl.org/Mus_musculus/Tools/IDMapper
    ## (or the equivalent for another species)

    ## The IDMapper result is expected to derive from a query 
    ## with the Ensembl IDs in parameter "queryIDs".

    ## NOTE that the downloaded results file is in some respects different
    ## from the results on the web page displayed on completion of the search:
    ## - The results file contains genes with "New stable ID"=="<retired>";
    ##   there is no entry for these genes on the page
    ##
    ## - The results file contains some empty sections 
    ##   ("missing" in the code below).
    ##   They are not attributable to specific query Ensembl IDs,
    ##   but the number of empty sections equals the number of query genes 
    ##   with no non-empty section in the results file,
    ##   so one very likely generates the other.
    ##   On the web page, genes with no non-empty section are not shown.
    ##
    ## - Very few entries in the results file contain a line with release="105"
    ##   (for results file KLjMWDn3CtXEFwcB-7995102.idmapper.txt,
    ##    where "105" is the current release).
    ##   On the web page, all displayed entries contain such a line.
    ##   ** So obviously what is displayed on the web page is also 
    ##      in the current release. 
    ##      But the ID may have been changed from the query ID.
    ##   **
    ##  --> distinguish in output of this function:
    ##  (A1) <retired>
    ##  (A2) missing in results file
    ##  (B) ID changed in most recent line of the result file entry
    ##  (C) the rest (should be OK)
    
    l <- readLines(IDMapper_output)
    header <- "Old stable ID, New stable ID, Release, Mapping score"

    from <-   which(l==header)    +1
    to   <- c(which(l==header)[-1]-1, length(l))
    
    cols <- strsplit(header,",\\s*")[[1]]
    mappings <- sapply(1:length(from),
                       function(i) {
                         ##cat(i,"\n")
                         this_l <-  strsplit(setdiff(l[from[i]:to[i]],""),
                                             ",\\s*")
                         if        (length(this_l)==0) {
                             return(NA) ## was: NULL
                         } else if (length(this_l)==1) {
                             d <- matrix(this_l[[1]],nrow=1)
                         } else {
                             d <- Reduce(rbind,
                                         strsplit(setdiff(l[from[i]:to[i]],""),
                                                  ",\\s*")
                                         )
                         }
                           
                         d <- as.data.frame(d)
                         colnames(d) <- cols
                           
                         d
                     },simplify=FALSE)
    missing <- is.na(mappings)

    ## get the Ensembl ID (without version) 
    ## corresponding to each "mappings" entry
    n <- sapply(mappings[!missing], 
                function(m) {
                    id <- unique(sub("\\.\\d+$","",m[,"Old stable ID"]))
                    if(length(id)>1) browser() 
                    ## >1 old ID per entry -- should not happen

                    id
                })
    names(mappings)[!missing] <- n

    if(any(missing)) {
        ## missing IDs should be those query IDs which had no hit:
        missed_query_IDs <- setdiff(query_IDs,n)
        if(sum(missing)!=length(missed_query_IDs)) {
            stop(paste("Inconsistent number of missing IDs in",
                       "map_Ensembl_versions_from_IDMapper_output: please check!"))
        }
        ## assign dummy result entries for the missing IDs
        mappings[missing] <-
            sapply(missed_query_IDs,
                   function(i) {
                       d <- data.frame(i,"<missing>","-1","-1")
                       colnames(d) <- c("Old stable ID",
                                        "New stable ID",
                                        "Release", "Mapping score")
                       d
                   }, simplify=FALSE)
        ## the names are not set by the assignment, why?
        names(mappings)[missing] <- missed_query_IDs 
    }



    ## for each Ensembl ID, get the last assigned "New stable ID",
    ## and the release in which it was assigned:
    res <- data.frame(t(sapply(mappings,
                               function(m) {
                                   ## get the last entry referring to the
                                   ## most recent release which is equal to or
                                   ## before target_release
                                   ## (ignore newer releases):
                                   l <- which.max(m$Release<=target_release)

                                   ## extract the valid ID and the release number:
                                   n <- m[l,"New stable ID"]
                                   if(grepl("^ENS.+\\.\\d+$",n)) {
                                       n <- sub("\\.\\d+$","",n)
                                   }
                                   c(newID=n,release=m[l,"Release"])
                               })))
    res$retired <- res$newID == "<retired>"
    res$missing <- res$newID == "<missing>"
    res$changed <- !res$retired & !res$missing & (res$newID != rownames(res))

    res

}







