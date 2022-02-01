## parsing
## sais:~/CECAD/Niessen/Project_Persa/Analysis/NewData_opersa/SignatureSearch/MSigDB/GdqrLddoe51HQaM8-7836507.idmapper.txt
## ==
## ugoebel-X570-AORUS-ELITE:~/Projects/Niessen/Project_Persa/SignatureSearch/MSigDB/November2021/GdqrLddoe51HQaM8-7836507.idmapper.txt
## (R-4.1.0)


## produced by https://www.ensembl.org/Mus_musculus/Tools/IDMapper/Results?tl=GdqrLddoe51HQaM8-7836507
## ID History Converter results
## Job details
## Job summary
## ID mapping of ids100.txt
## Species Mouse (Mus musculus)
## Assembly GRCm39

source("~/CECAD/Programming/Sharable_code/General/general_functions.R")
refseq_IDs_from_EnsemblGTF <- function(ensembl_GTF_file="~/CECAD/Pipeline/Bioc_annotation/Ensembl103/Mus_musculus.GRCm39.103.gtf") {
    v <- extractCapturedSubstrings("(?:^|/)([^_/]+_[^_/]+)\\.(GRCm\\d+)\\.(\\d+)\\.",
                                   ensembl_GTF_file)
    if(any(sapply(v,length)==0)){
        stop("Expecting an Ensembl GTF file name!")
    }
    organism <- v[1]
    genome_build <- v[2]
    ensembl_version <- v[3]

    gtf <- rtracklayer::readGFF(ensembl_GTF_file,version=2)

    tibble(organism=organism,
           genome_build=genome_build,
           ensembl_version=ensembl_version,
           refseqIDs=unique(gtf[gtf$source=="RefSeq","gene_id"]) )
                                
        
}

problematic_IDs  <- function(mapping=mappings[!is.na(mappings)],
                             target_release="103",
                             keep=mitoIDs) {
    
    ## discard releases > target_release
    mapping <- mapping[sapply(mapping,
                              function(m)as.numeric(m[nrow(m),"Release"])<=as.numeric(target_release))]
    
    ## Did the ID switch away from the initial ID (which is the name of its list entry)
    ## in any release (previous or target)?
    switched <-  sapply(names(mapping),
                        function(n)any(n != sub("\\.\\d+$","",
                                                mapping[[n]][,"New stable ID"]))
                        )
    
    ## new ID and release at the most current entry for each Ensembl ID
    last_newID   <- sapply(mapping,function(m)sub("\\.\\d+$","",m[nrow(m),"New stable ID"]))
    last_release <- sapply(mapping,function(m)as.numeric(m[nrow(m),"Release"]))
    
    
    ## the Ensembl IDs which will be returned as "potentially bad"
    ## (to be discarded, or potentially renamed to a more current ID): 
    potentially_bad <- setdiff(names(mapping)[(last_release < as.numeric(target_release) | switched)],
                               keep)
    
    
    
    ## break the "potentially bad" IDs up into subsets:
    
    retired <- names(which(last_newID[potentially_bad]=="<retired>"))
    current_switch <- names(which( (last_release[potentially_bad]==target_release)
                                  &
                                  ((last_newID[potentially_bad] != names(last_newID[potentially_bad]))  ## an actual switch
                                      & 
                                      (last_newID[potentially_bad] != "<retired>"))                     ## not a "retirement"
                                  )
                            )
    list(retired=retired, ## clearly out of the game
         current_switch=last_newID[current_switch], ## these may be kept and renamed
         rest=setdiff(potentially_bad,
                      union(retired,current_switch)))
    
}



map_Ensembl_versions_from_IDMapper_output <- function(IDMapper_output,
                                                      target_release="105") {
    ## parse an output file created by and downloaded from
    ## https://www.ensembl.org/Mus_musculus/Tools/IDMapper
    ## (or the equivalent for another species)

    l <- readLines(IDMapper_output)
    header <- "Old stable ID, New stable ID, Release, Mapping score"

    from <-   which(l==header)    +1
    to   <- c(which(l==header)[-1]-1, length(l))
    
    cols <- strsplit(header,",\\s*")[[1]]
    mappings <- sapply(1:length(from),
                       function(i) {
                           ##cat(i,"\n")
                           this_l <-  strsplit(setdiff(l[from[i]:to[i]],""),",\\s*")
                           if        (length(this_l)==0) {
                               return(NA) ## was: NULL
                           } else if (length(this_l)==1) {
                               d <- matrix(this_l[[1]],nrow=1)
                           } else {
                               d <- Reduce(rbind,
                                           strsplit(setdiff(l[from[i]:to[i]],""),",\\s*")
                                           )
                           }
                           
                           d <- as.data.frame(d)
                           colnames(d) <- cols
                           
                           d
                       },simplify=FALSE)
    missing <- which(is.na(mappings))


    ## get the Ensembl ID (without version) corresponding to each "mappings" entry
    n <- sapply(mappings[-missing], 
                function(m) {
                    id <- unique(sub("\\.\\d+$","",m[,"Old stable ID"]))
                    if(length(id)>1) browser() ## >1 old ID per entry -- should not happen

                    id
                })
    nm <- rep(NA,length(mappings))
    nm[-missing] <- n
    names(mappings) <- nm

    ## The mitochondrial IDs of the mouse Ensembl annotations come from the early release v64,
    ## but they are seemingly not updated in newer releases -> they should be kept when mapping to a newer release.
    ## Luckily they are (at least in v103 and v105, others not checked) the only genes with source=="RefSeq",
    ## so they are easy to identify.
    ## Extract them here ...
    fn <-"~/CECAD/Pipeline/Bioc_annotation/Ensembl103/Mus_musculus.GRCm39.103.gtf"
    mitoIDs <- refseq_IDs_from_EnsemblGTF(ensembl_GTF_file=fn)$refseqIDs

   
    ## .. and then pass them to problematic_IDs, which returns IDs which have been retired, have no entries in newer
    ## releases, or have been assigned a new ID in the target release:
    problematic_IDs(mapping=mappings[!is.na(mappings)],
                    target_release="103",
                    keep=mitoIDs)






}







