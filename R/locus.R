#' @title Genomic locus for CNV calls
#' 
#' @description 
#' \code{locus} Extract genomic locus information to each 
#' row of a dataframe consisting of CNV-like entries
#' 
#' @details 
#' This function takes as input a dataframe containg CNV calls
#' or any similar entries and returns the same dataframe with 
#' three additional columns: "loc.str", "loc.end", and "locus".
#' This can be useful per se and it is a required step for
#' comparing two datasets with the function \code{\link{inter_comp}}.
#' Input must possess the following columns: 
#' \item{"chr", chromosome of the call in GRCh format (i.e. "1", 
#' not "chr1")}
#' \item{"start" start of the call}
#' \item{"end" end of the call}
#' \newline
#' By default the function will attempt to download the required 
#' cytobands file of the selected assembly (default is "hg19"), 
#' it is possible to pass a local file as \code{bands} setting
#' the \code{whichCyto} parameter to \code{"local"} instead. 
#' The function uses a for loop and this is its major bottleneck. 
#' In order to speed up the process the input dataset is splitted 
#' according to the \code{n.cores} parameter and the splits are 
#' processed in parallel. As an example, processing a dataframe with 
#' ~10500 entries takes about 10 seconds using 4 cores, and about 4 
#' seconds using 16 cores on our system, while the same work using 
#' only one core takes around 31 seconds.
#' Default number of cores is 4, in this way it should 
#' work with default parameters even on a laptop.
#' 
#' @param cnv.in a dataframe consisting of CNVcalls
#' @param whichCyto "remote" or "local"
#' @param bands dataframe containg genomic locus reference for the 
#' selected assembly
#' @param assembly genomic assobly, either "hg18", "hg19" or "hg38"
#' @param n.cores number of usable CPU cores
#' 
#' @export
#' 
#' @return cnv.out
#' 
#' @author Simone Montalbano simone.montalbano@protonmail.com


locus <- 
    function(cnv.in, whichCyto = "remote" , bands, assembly="hg19", n.cores = 4) {
        require(dplyr)
        require(magrittr)
        require(data.table)
        require(foreach)
        require(doParallel)
        
        registerDoParallel(cores = n.cores)
        if  (whichCyto == "remote") bands <- fread(paste0("https://hgdownload.cse.ucsc.edu/goldenPath/", 
                                                          assembly, "/database/cytoBand.txt.gz"))
        if (whichCyto == "local") cat("Using local cytoBands file! \n")
        colnames(bands) <- c("chr", "start", "end", "locus", "buh")
        bands %<>% dplyr::select(chr, start, end, locus) %>% 
            arrange(start) %>% arrange(chr) %>% 
            mutate(chr = substr(chr,4,5))
        cnv.in %>% mutate(loc.str = NA, loc.end = NA, locus = NA)
        
        len <- round(nrow(cnv.in)/n.cores)
        splits <- split.data.frame(cnv.in, rep(1:ceiling(n.cores), each=len, len.out=nrow(cnv.in)))
        
        loop <- function(n) {
            for (i in 1:nrow(splits[[n]])) {
                tmp <- bands %>% filter(chr == splits[[n]]$chr[i])
                if (nrow(tmp > 0)) {
                    tmp.str <- tmp %>% filter(start <= splits[[n]]$start[i] &
                                                  end >= splits[[n]]$start[i])
                    tmp.end <- tmp %>% filter(start <= splits[[n]]$end[i] &
                                                  end >= splits[[n]]$end[i])
                    splits[[n]]$loc.str[i] <- tmp.str$locus
                    splits[[n]]$loc.end[i] <- tmp.end$locus
                    if (splits[[n]]$loc.str[i] == splits[[n]]$loc.end[i]) {
                        splits[[n]]$locus[i] <- paste0(splits[[n]]$chr[i], splits[[n]]$loc.str[i]) }
                    else splits[[n]]$locus[i] <- paste0(splits[[n]]$chr[i], splits[[n]]$loc.str[i], 
                                                        "-", splits[[n]]$loc.end[i])
                }
                else {
                    splits[[n]]$loc.str[i] <- NA
                    splits[[n]]$loc.end[i] <- NA
                    splits[[n]]$locus[i] <- NA
                    
                }
                
            }
            return(splits[[n]])
        }
        
        cnv.out <- foreach(n=1:n.cores, .combine=rbind) %dopar% loop(n)
        
        cnv.out %<>% filter(!is.na(locus))
        
        return(cnv.out)
    }