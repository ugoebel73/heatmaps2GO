##-------------------------------------------------------------------------------------------------------------------------

extractCapturedSubstrings <- function(pattern, string, global_search=FALSE) {
  if(global_search) {
    tmp <- attributes(gregexpr(pattern,string,perl=TRUE)[[1]])[c("capture.start","capture.length")] ## only 1 input string -->> [[1]]
  } else {
    tmp <- attributes( regexpr(pattern,string,perl=TRUE)     )[c("capture.start","capture.length")]
  }
  substring(string,first=as.vector(tmp$capture.start),last=as.vector(tmp$capture.start+tmp$capture.length-1))
}

##----------------- parsers for the individual Enrichr databases ----------------------------------------------------------
##----------------- NOTE: the database names are implicitly specified here, as the names of the parsers!! -----------------
##-----------------       For each parser name, a file with this name is expected in directory Enrichr_base ---------------
Enrichr_parsers <- list()
Enrichr_parsers[["ChEA_2016"]]                      <- function(x) {
                                                       y <- toupper(strsplit(x[1],"\\s+")[[1]])
                                                       ## split the leading part on whitespece
                                                       has_organism <- y[length(y)] %in% c("HUMAN","MOUSE","RAT")
                                                       c(term=y[1],
                                                         id=y[2],
                                                         method=y[3],
                                                         cell_type=paste(y[4:(length(y)-as.numeric(has_organism))],collapse="_"),
                                                         genome=ifelse(has_organism,y[length(y)],NA))
                                                    }
##-------------------------------------------------------------------------------------------------------------------------
Enrichr_parsers[["ENCODE_TF_ChIP-seq_2015"]]        <- function(x) {
                                                       y <- toupper(strsplit(x[1],"\\s+")[[1]]) 
                                                       c(term=y[1],
                                                         cell_type=paste(y[-c(1,length(y))],collapse="_"), ## pure guess ..
                                                         method="ChIP-seq",
                                                         genome=y[length(y)])
                                                    }
##-------------------------------------------------------------------------------------------------------------------------
Enrichr_parsers[["GO_Biological_Process_2018"]] <-
    Enrichr_parsers[["GO_Molecular_Function_2018"]] <-
    Enrichr_parsers[["GO_Cellular_Component_2018"]] <- function(x) {
                                                       y <- strsplit(x[1],"\\s+")[[1]]
                                                       c(term=paste(y[1:(length(y)-1)],collapse="_"),
                                                         id=extractCapturedSubstrings("\\(([^()]+)\\)$",x[1]))
                                                    }
##-------------------------------------------------------------------------------------------------------------------------
Enrichr_parsers[["KEGG_2019_Mouse"]]                <- function(x) {
                                                       y <- strsplit(x[1],"\\s+")[[1]]
                                                       c(term=paste(y,collapse="_"))
                                                    }
##-------------------------------------------------------------------------------------------------------------------------
Enrichr_parsers[["Reactome_2016"]]                  <- function(x) {
                                                       y <- strsplit(x[1],"\\s+")[[1]]
                                                       c(term=paste(y[1:(length(y)-3)],collapse="_"),
                                                         genome=paste(y[(length(y)-2):(length(y)-1)],collapse="_"), ## always?
                                                         id=y[length(y)])
                                                    }
##-------------------------------------------------------------------------------------------------------------------------
Enrichr_parsers[["WikiPathways_2019_Mouse"]]       <- function(x) {
                                                      y <- strsplit(x[1],"\\s+")[[1]]
                                                      c(term=paste(y[1:(length(y)-1)],collapse="_"),
                                                        id=extractCapturedSubstrings("(WP\\d+)$",x[1]))
                                                    }
##-------------------------------------------------------------------------------------------------------------------------

##--------------------------------------------------------------------------------------------------------------------------------------
##----------------- use the parsers to read the gene sets and additional info from the database files ----------------------------------
read_Enrichr_dbs <- function(Enrichr_base= ".",
                             databases=c("GO_Biological_Process_2018",
                                         "GO_Cellular_Component_2018",
                                         "GO_Molecular_Function_2018"),
                             parsers=Enrichr_parsers ## for each element of
                                                     ## vector "databases"
                                                     ## a file with this name 
                                                     ## is expected in 
                                                     ## directory Enrichr_base
                             ) {
    Enrichr_dbs <- list()
    
    for(n in databases) {
        cat(n,"\n")
        
        v <- readLines(paste(Enrichr_base, n,sep="/"))
        if(!all(sapply(v,function(x)grepl("\t\t",x)))) {
            cat("*** Unexpected database format for", db, "-- skipping! ***\n")
            next
        }
        Enrichr_dbs[[n]] <- list()
        
        tmp <- strsplit(v,'\t\t')
        tag <- sapply(tmp,function(x)gsub("\\s+","_",x[1]))              ## x[1], "as is", "_"-separated
        
        Enrichr_dbs[[n]]$info <- sapply(tmp, parsers[[n]])               ## x[1], parsed
        if(!is.null(dim(Enrichr_dbs[[n]]$info))) {
            Enrichr_dbs[[n]]$info <- t(Enrichr_dbs[[n]]$info)
            
        } else {
            cn <- unique(names(Enrichr_dbs[[n]]$info))
            cn <- ifelse(length(cn)==1,cn,"term")
            Enrichr_dbs[[n]]$info <- matrix(Enrichr_dbs[[n]]$info,ncol=1)
            colnames(Enrichr_dbs[[n]]$info) <- cn
        }
        rownames(Enrichr_dbs[[n]]$info) <- tag
        
        Enrichr_dbs[[n]]$genes <- sapply(tmp,function(x)Reduce(union,sapply(x[2:length(x)],function(y)strsplit(y,"\\t"))))
        names(Enrichr_dbs[[n]]$genes) <- tag                             ## x[2:(length(x)],each split into genes on "\t",
                                                                         ## then merged.
                                                                         ## The GO dbs have >1 "\t\t" in some rows
                                                                         ## (?? indicating hierarchical levels??).
                                                                         ## Simply ignoring the additional "\t\t"s
                                                                         ## yields correct set sizes as reported by Enrichr.
                                                                         ##--------->>   But: if possible, clarify the meaning of the format!!
        
    }
    
    Enrichr_dbs
}
##--------------------------------------------------------------------------------------------------------------------------------------
## The following function was written in the course of the very initial analysis of the Niessen data.
## It was copied to this file on Dec 2, 2021 from ugoebel-X570-AORUS-ELITE:/home/ugoebel/Programming/Sharable_code/DGE/dge_functions.R
DGEoutput2Enrichr <- function(d, rule="QVAL", qval_col="Adjusted.P.value", t_col="t", logFC_col="logFC", gene_col="__ROWNAMES",
                              gene_map=NULL, N_t=500, maxQ=0.05, cmp=">",geneset_dbs=biojupies_dbs) {
    ## d: the result data.frame of a differential gene expression program
    ## qval_col:  the column holding the adjusted p-value (or any other measure to be handled by the QVAL   rule) 
    ## t_col:     the column holding the t-value          (or any other measure to be handled by the TRULE  rule)
    ## logFC_col: the column holding the effect strength AND EFFECT DIRECTION (logFC, or the estimated model coefficient)
    ## rule==QVAL:   return genes (with d[,qval_col]<=maxQ) & get(cmp)(d[,logFC_col],0)]
    ##                                                      ## The get() function takes a string as input and returns a reference
    ##                                                      ## to the function corresponding to that string, if it exists.
    ##                                                      ## Used here to pass the comparison operator (to 0) as a parameter.
    ## rule==TRULE:  return the N_t genes with the highest (cmp==">") or lowest (cmp=="<") values in d[,t_col]
    ##               (used by BioJupies on the limma result tables, to generate input for Enrichr)
    
    ## gene_col:  the column holding the gene names to be used in the output
    ##            If gene_col=="__ROWNAMES" --> use rownames(d)
    ##
    ## gene_map:  a vector v with v[gene_col]=external which maps gene names in d to mouse MGI names. If NULL, use as is.

    if(gene_col=="__ROWNAMES") {
        genes <- rownames(d)
    } else {
        genes <- d[,gene_col]
    }
    if(!is.null(gene_map)) {
        genes <- gene_map[genes]
    }
    if("THIS_GENE_NAMES" %in% colnames(d)) {
        cat("*** Column name conflict in bioJupiesList2enrichR -- returning NULL!")
    } else {
        d$THIS_GENE_NAMES <- genes  ## append a column with the output gene names
    }
    
    if(rule=="QVAL") {
        l <- unique(d$THIS_GENE_NAMES[(d[,qval_col]<=maxQ) & get(cmp)(d[,logFC_col],0)])
    } else if (rule=="TRULE") { 
        l <- head(unique(d[order(d[,t_col],decreasing=(cmp==">")),]$THIS_GENE_NAMES),
                  n=500)
    }

    list(e=enrichr(l,geneset_dbs),l=l)
}

        
query_Enrichr <- function(DGE_result, gene_map, ## from gene names in DGE_result to gene symbols
                          dbs, ## the parsed database files
                          use_dbs, ## names of the databases to use (Enrichr output from these dbs will be combined in the function output)
                          additional_columns, ## additional parsed information from the database file
                          DGE_maxq=0.05,Enrichr_maxq=0.05,
                          rule="QVAL",
                          qval_col="adj.P.Val",
                          logFC_col= "logFC",
                          gene_col ="__ROWNAMES",
                          t_col="t",
                          effect="up") {
                          ##logFC_thrsh=0) {

    ## NOTE here: when function DGEoutput2Enrichr() queries the Enrichr server with the database names,
    ##            then the *server-side copies* of the databases are used for the actual search!
    ##
    ##            I use the downloaded database flatfiles only to extract additional information.
    ##            ! This ASSUMES that a given database name is associate with stable content during
    ##              the entire evolution of the Enrichr web server !
    ##
    ##            The gene sets from the downloaded database files are only used to compute an overlap 
    ##            of the uppercase-d input gene list and the individual gene sets.
    ##            Currently I stop with an error if the uppercase-d overlap differs from the overlap reported by Enrichr.
    ##            With Carien Niessen's "limma" and "edgeR" regulated mouse genes I do not run into this error,
    ##            so Enrichr indeed seems to simply match the uppercase-d input against the databases
    ##            (which misses mouse genes with any small deviation from the human gene symbol).  

    raw <- DGEoutput2Enrichr(DGE_result, 
                             rule=rule,
                             maxQ=DGE_maxq,
                             qval_col =qval_col,
                             logFC_col=logFC_col,
                             gene_col =gene_col,
                             gene_map =gene_map,
                             cmp=ifelse(effect=="up",">","<"),
                             geneset_dbs=use_dbs) 
                             ##min_absLFC=logFC_thrsh)
    genes <- raw$l ## the differentially expressed query genes

    result <- list(inlist=genes,
                   Enrichr=sapply(use_dbs, function(x)c()))
    
    for(db in use_dbs) {
        ## raw Enrichr result for database db
        d <- raw$e[[db]]; if(is.null(d)) next ## skip if no result

        ## are all requested additional output columns indeed present in the passed "dbs" ?
        additional_columns <- intersect(additional_columns,colnames(dbs[[db]]$info))
        
        ##---------------------------------------------------------------------------------------------------------------
        ## filter by p-value
        d <- d[d[,"Adjusted.P.value"]<=Enrichr_maxq,]; if(nrow(d)==0) next ## skip if no significant result

        ##---------------------------------------------------------------------------------------------------------------
        ## check whether all "Term" entries (column 1) of the Enrichr result do appear as rownames of the "info" data.frame
        ## of the current database -- if not, there must have been an upstream database parsing error of some sort
        dbtag <- gsub("\\s+","_",d[,"Term"]) ## whitespace has been replaced by "_" in the rownames -- apply to "Term" column, too

        for(colmn in additional_columns) {
            if(length(setdiff(dbtag,rownames(dbs[[db]]$info))) > 0) {
                stop("Some Terms in the Enrichr result are missing in the parsed databases -- check for parsing errors!")
            }
        }
        
        this_G2g <- data.frame(Enrichr_db=db,
                               Enrichr_geneset=dbtag, ## (not necessary at least for GO: replicated in "additional columns"!)
                               DGE_effect=effect,
                               Enrichr_pval=d[,"P.value"],
                               Enrichr_qval=d[,"Adjusted.P.value"],
                               Enrichr_odds=d[,"Odds.Ratio"]) 
        
        ##---------------------------------------------------------------------------------------------------------------
        ## append additional columns (from the parsing of the database file) to the raw Enrichr result
        for(colmn in additional_columns) {
            if(colmn %in% colnames(dbs[[db]]$info)) {
                this_G2g[,paste("Enrichr_",colmn,sep="")] <- dbs[[db]]$info[dbtag,colmn]
            } else {
                this_G2g[,paste("Enrichr_",colmn,sep="")] <- NA
            }
        }
        ##---------------------------------------------------------------------------------------------------------------
        ## append the numbers of genes "in input & in term"(=_list) and "in term"(=_term) explicitly,
        ## by parsing the "Overlap" column
        n <- sapply(d[,"Overlap"],extractCapturedSubstrings,pattern="^(\\d+)/(\\d+)$")
        this_G2g[,"Enrichr_ngenes_list"] <- as.numeric(n[1,])
        this_G2g[,"Enrichr_ngenes_term"] <- as.numeric(n[2,])
        
      
        ## NOTE that the gene IDs in the vertebrate Enrichr GO databases are all uppercase (likely human)!
        ## What I do here is to simply check the input gene list converted to all-uppercase against a term's gene set,
        ## and return the genes which have a match.
        ## NOTE that it may not(??) be guaranteed that all gene IDs returned by Enrichr do have a mouse counterpart --
        ##      but only if Enrichr does something more sophisticated than a simple string match!
        ids <- sapply(dbtag,
                      function(tag) {
                          g1 <- sort(genes[which(toupper(genes) %in% dbs[[db]]$genes[[tag]])]) ## uppercase-d mouse genes contained in gene set
                          g2 <- sort(strsplit(d[gsub("\\s+","_",d$Term)==tag,"Genes"],";")[[1]]) ## overlapping genes reported by Enrichr
                          if(!all(toupper(g1)==g2)) {
                               stop(paste0("Term ", tag, ": ",
                                           "Gene overlaps (a) returned by Enrichr and ",
                                           "(b) from the intersection of the gene set and the input genes are not the same!\n"))
                          }
                          g1 ## return the input version (can differ from the reported version only by case)
                      }, 
                      simplify=FALSE)
        ## have a separate column for the number of input genes actually present in each gene set:
        ## (with the "stop" solution above, Enrichr_ngenes_list==DGE_ngenes_list  by construction)
        this_G2g[,"DGE_ngenes_list"] <- sapply(ids,length)
       
        ## append the lists of input genes actually present in the gene sets:
        this_G2g[,"DGE_genes_in_geneset"] <- sapply(ids,paste,collapse=",")

        ##---------------------------------------------------------------------------------------------------------------
        ## finally, to allow an easy conversion to class "enrichResult" (-> use with clusterProfiler):
        ## (copied from sais:~/CECAD/Ladwig/Project/Analysis/DGE/analysis.R)
        ## On the meaning of the columns in the ORA output tables:
        ## See https://www.biostars.org/p/220465/ :
        ## My interpretation:
        ## BgRatio,   M/N = ratio of (genes in my "universe" AND in GO ontology AND in term) /
        ##                           (genes in my "universe" AND in GO ontology) 
        ## GeneRatio, k/n = ratio of (genes in my "genes" AND in GO ontology AND in term) / ##_list
        ##                           (genes in my "genes" AND in GO ontology)
        
        ## Because the Enrichr result includes no concept of a "universe", assume that the universe is just the entire database:
        ## BgRatio,   M/N = ratio of (genes in GO ontology AND in term) / ##_term
        ##                           (genes in GO ontology)

        
        flat_db <- Reduce(union,dbs[[db]]$genes)
        this_G2g[,"Enrichr_dbsize"] <- length(flat_db) ## genes in any of the gene sets of the current database
        this_G2g[,"DGE_list_in_db"] <- length(intersect(toupper(genes),flat_db)) ## note that the intersection can be small!!

        this_G2g[,"DGE_listsize"] <- length(genes) ## for comparison with "DGE_list_in_db" -- ideally the two numbers should agree

        this_G2g[,"GeneRatio"] <- paste(this_G2g[,"DGE_ngenes_list"],
                                        this_G2g[,"DGE_list_in_db"],sep="/")
        this_G2g[,"BgRatio"] <- paste(this_G2g[,"Enrichr_ngenes_term"],
                                      this_G2g[,"Enrichr_dbsize"],sep="/")
        
        
        ##---------------------------------------------------------------------------------------------------------------
        result$Enrichr[[db]] <- this_G2g 
    }
    
    result
    
}

##-------------------------------------------------------------------------------------------------------------------------
##----------------- Visualization using clusterProfiler functions ---------------------------------------------------------
enrichr2enrichResult <- function(x,
                                 colnames_map =c(Enrichr_id="ID",
                                                 Enrichr_term="Description",
                                                 GeneRatio="GeneRatio",
                                                 BgRatio="BgRatio",
                                                 Enrichr_pval="pvalue",
                                                 Enrichr_qval="p.adjust",
                                                 Enrichr_ngenes_list="Count",
                                                 DGE_genes_in_geneset="geneID"),
                                 use_ont="GO_Biological_Process_2018",
                                 organism="Mus musculus",
                                 keytype="SYMBOL",
                                 delimiter=","
                                 ) {
    ## "x" is an output object of query_Enrichr()
    ## The function converts x$Enrichr to an object of class enrichResult, thereby also using the x$inlist element.
    ## The output object is usable by clusterProfiler visualization functions expecting an enrichResult
    ## (tested so far: dotplot(), emapplot(pairwise_termsim()).
    ## For emapplot, use the default method="JC" in pairwise_termsim()).

    ## (other databases will be added)
    ont_map <- c(GO_Biological_Process_2018="BP",
                 GO_Molecular_Function_2018="MF",
                 GO_Cellular_Component_2018="CC")

    if(!(use_ont %in% names(ont_map))) {
        stop(paste("enrichr2enrichResult currently only supports the following Enrichr databases:",
                   paste(names(ont_map),collapse=","), "!"))
    }
    if(!(use_ont %in% names(x$Enrichr))) {
        stop(paste("The supplied x$Enrichr does has no element",use_ont,"!"))
    }

    if(!all(c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","Count","geneID") %in% colnames_map)) {
        stop("The supplied colnames_map does not create all required column names!")
    }
     if(!all(names(colnames_map) %in% colnames(x$Enrichr[[use_ont]]))) {
        stop("Some columns listed in the supplied colnames_map are not present in the supplied x$Enrichr[[use_ont]]!")
    }
       
    df <- x$Enrichr[[use_ont]][names(colnames_map)]
    colnames(df) <- colnames_map[colnames(df)]
    rownames(df) <- df$ID ## this is important for pairwise_termsim()!

    if(delimiter !="/") {
        df$geneID <- gsub(delimiter,"/",df$geneID) ## geneInCategory() splits on "/" only
    }
    
    new("enrichResult", result = df,
        gene=x$inlist,
        ontology=ont_map[use_ont],
        organism=organism,
        keytype=keytype)
    
}
