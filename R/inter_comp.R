#' @title Inter results CNV calls comparison
#' 
#' @description 
#' \code{inter_comp} compare results from different 
#' methods to find the common calls, or find the unique calls 
#' between groups in cases/controls scenario. 
#' 
#' @details
#' This function compare the results of CNV calling methods. It can
#' be useful when merging results from different pipelines on the same 
#' data, in order to highlight the common calls (of higher 
#' confidence) and avoid duplicate calls.
#' It can also handle case/control situations (selecting the
#' calls present in cases only), and family-based 
#' studies.
#' \newline
#' The function is specifically designed to work on SNP array data
#' and require SNPs position information to process calls. In particular,
#' the actual comparison is made on the markers rather than on the raw
#' genomic coordinates. As an example, in default settings, two calls 
#' will be treated as the same if they share 50% or more of the 
#' markers. 
#' Inside the function there is a filter based on the number of markers,
#' default behavior is to eliminate the calls with less than 10
#' markers. If this is undesired simply set \code{min.markers} to 0.
#' \newline
#' There are three possible modes that can be set with the parameter
#' \code{comp.type}. If set to "inter" the function will scan \code{in.a} 
#' and \code{in.b} for replicate calls, it will add a new column, 
#' "uniq" (if 0 call is replicated) and then it will eliminate the 
#' replicated calls from \code{in.b} before merging the two datasets.
#' In contrast, "matched" and "case/control" assume that \code{in.a}
#' contains calls for the case(s) and \code{in.b} for the control(s). 
#' It will then attempt to select the  calls of the case that are not
#' replicated in control. The difference between the two is that 
#' "case/control" assume one-sample-one-object while "matched" account
#' for family ID and can handle more than one sample per object.
#' \newline
#' Required input files are: 
#' \code{markers} is a data-frame containing information about the SNP 
#' markers of the array used (required columns: "chr" "position"
#' "snp");
#' \code{in.a in.b} are two data-frame containing the actual CNV calls
#' (required columns: "chr" "start" "end" "CN" "loc.start" "loc.end").
#' \newline
#' The function uses a for loop and this is its major bottleneck. 
#' In order to speed up the process the input dataset is splitted 
#' according to the \code{n.cores} parameter and the splits are 
#' processed in parallel.
#' Default number of cores is 4, in this way it should 
#' work with default parameters even on a laptop.
#' 
#' @param in.b First dataset or Case(s)
#' @param in.b Second dataset or Control(s)
#' @param markers Data-frame containing the SNP informations:
#' "chr", "pos" (chromosomal coordinate), "snp" (name of the marker)
#' @param threshold Desired threshold for the comparison, e.g. 
#' if 0.5 two calls will be treated as the same if sharing 50% or 
#' more of the markers
#' @param met.a Name of the algorithm/pipeline of the first dataset, 
#' useful to keep the information after the merge
#' @param met.b Name for the second dataset
#' @param comp.type type of comparison, can be either "inter", "matched"
#' or "case/control"
#' @param min.markers minimum number of marker a calls need to contains 
#' @param n.cores number of CPU core to use 
#' @param keepCols set it to F in order to discard the intermediate 
#' colums from \code{inter_comp} and \code{locus}.
#'  
#' @export
#' 
#' @return res
#' 
#' @author Simone Montalbano simone.montalbano@protonmail.com

inter_comp <- 
    function(in.a, in.b, markers, threshold = 0.5, met.a = "methodA", 
             met.b = "methodB", comp.type = "inter", min.markers = 10,
             n.cores = 4, keepCols = T) {
        require(dplyr)
        require(magrittr)
        require(foreach)
        require(doParallel)
        
        registerDoParallel(cores = n.cores)
        if (comp.type == "matched") {
            dt1 <- in.a %>% arrange(start) %>% arrange(chr) %>%
                mutate(srt.ix = NA, end.ix = NA, num.snp = NA,
                       uniq = NA)
            dt2 <- in.b %>% arrange(start) %>% arrange(chr)%>%
                mutate(srt.ix = NA, end.ix = NA, num.snp = NA,
                       uniq = NA)
        }
        else {
            dt1 <- in.a %>% arrange(start) %>% arrange(chr) %>%
                mutate(srt.ix = NA, end.ix = NA, num.snp = NA,
                       uniq = NA, method = NA)
            dt2 <- in.b %>% arrange(start) %>% arrange(chr)%>%
                mutate(srt.ix = NA, end.ix = NA, num.snp = NA,
                       uniq = NA, method = NA)
        }
        
        
        snps <- markers %>% arrange(pos) %>% arrange(chr)
        tr <- threshold
        
        ## Functions
        
        start.end <- function(input) {
            snps %<>% mutate(chrposstart=paste(chr, pos, sep=":"))
            tmp <- input %>% mutate(chrposstart=paste(chr, start, sep=":"),
                                    chrposend=paste(chr, end, sep=":"))
            
            ixst <- match(tmp$chrposstart, snps$chrposstart)
            ixend <- match(tmp$chrposend, snps$chrposstart)
            
            input$srt.ix[!is.na(ixst)] <- ixst[!is.na(ixst)]               # first snp index
            input$end.ix[!is.na(ixend)] <- ixend[!is.na(ixend)]            # last snp index
            input$num.snp <- input$end.ix - input$srt.ix                   # total snp
            
            return(input)
        }
        
        comp <- function(a, b, n) {
            
            a$uniq <- NA
            
            if (comp.type == "matched") {
                for (i in 1:nrow(a)) {
                    
                    tmp <- b %>% filter(family == a$family[i] & CN == a$CN[i] & chr == a$chr[i] &
                                            (loc.str == a$loc.str[i] | loc.end == a$loc.end[i]))
                    
                    if (nrow(tmp) == 0) a$uniq[i] <- 1
                    else {
                        comm <- c()
                        markA <- a$srt.ix[i]:a$end.ix[i]
                        for (ix in 1:nrow(tmp)) {
                            markB <- tmp$srt.ix[ix]:tmp$end.ix[ix]
                            comm.tmp <- intersect(markA, markB)
                            comm <- c(comm, as.integer(length(comm.tmp)))
                        }
                        
                        len <- tr*length(markA)
                        ind <- which.max(comm)
                        if (comm[ind] >= len) a$uniq[i] <- 0
                        else                  a$uniq[i] <- 1
                    }
                }
            }
            
            if (comp.type == "case/control") {
                for (i in 1:nrow(a)) {
                    
                    tmp <- b %>% filter(CN == a$CN[i] & chr == a$chr[i] &
                                            (loc.str == a$loc.str[i] | loc.end == a$loc.end[i]))
                    
                    if (nrow(tmp) == 0) a$uniq[i] <- 1
                    else {
                        comm <- c()
                        markA <- a$srt.ix[i]:a$end.ix[i]
                        for (ix in 1:nrow(tmp)) {
                            markB <- tmp$srt.ix[ix]:tmp$end.ix[ix]
                            comm.tmp <- intersect(markA, markB)
                            comm <- c(comm, as.integer(length(comm.tmp)))
                        }
                        
                        len <- tr*length(markA)
                        ind <- which.max(comm)
                        if (comm[ind] >= len) a$uniq[i] <- 0
                        else                  a$uniq[i] <- 1
                    }
                }
            }
            
            if (comp.type == "inter") {
                chroms <- unique(a$chr)
                a %<>% filter(chr == chroms[n])
                for (i in 1:nrow(a)) {
                    tmp <- b %>% filter(sample.name == a$sample.name[i] & CN == a$CN[i] & chr == a$chr[i] &
                                            (loc.str == a$loc.str[i] | loc.end == a$loc.end[i]))
                    if (nrow(tmp) == 0) a$uniq[i] <- 1
                    else {
                        comm <- c()
                        markA <- a$srt.ix[i]:a$end.ix[i]
                        for (ix in 1:nrow(tmp)) {
                            markB <- tmp$srt.ix[ix]:tmp$end.ix[ix]
                            comm.tmp <- intersect(markA, markB)
                            comm <- c(comm, as.integer(length(comm.tmp)))
                        }
                        len <- tr*length(markA)
                        ind <- which.max(comm)
                        if (comm[ind] >= len) a$uniq[i] <- 0
                        else                  a$uniq[i] <- 1
                    }
                }
            }
            
            return(a)
        }
        
        ## Application
        
        dt1.1 <- start.end(dt1)
        dt2.1 <- start.end(dt2)
        dt1.1 %<>% filter(num.snp >= min.markers)
        dt2.1 %<>% filter(num.snp >= min.markers)
        
        if (comp.type == "matched") dt1.2 <- comp(dt1.1, dt2.1)
        
        if (comp.type == "inter") {
            chroms.a <- as.character(unique((in.a %>% arrange(chr))$chr))
            chroms.b <- as.character(unique((in.b %>% arrange(chr))$chr))
            # if (!all.equal(chroms.a, chroms.b)) cat("Different number of chromosome between input,
            #                                         maybe one method does not support chr X?\n")
            # str(chroms.a)
            # str(chroms.b)
            dt1.2 <- foreach(n=1:length(chroms.a), .combine = rbind) %dopar% comp(dt1.1, dt2.1, n)
            dt2.2 <- foreach(n=1:length(chroms.b), .combine = rbind) %dopar% comp(dt2.1, dt1.1, n)
        }
        
        # select
        if (comp.type == "matched") res <- dt1.2 %>% filter(uniq == 1)
        
        if (comp.type == "inter") {
            dt1.2$method <- met.a
            dt2.3 <- dt2.2 %>% filter(uniq == 1)
            dt2.3$method <- met.b
            res <- rbind(dt1.2, dt2.3)
        }
        
        if (keepCols == T) res %<>% arrange(start) %>% arrange(chr) %>% arrange(sample.name)
        if (keepCols == F) res %<>% select(everything(), -loc.str, -loc.end, -end.ix, -srt.ix) %>%
            arrange(start) %>% arrange(chr) %>% arrange(sample.name)
        
        
        return(res)
    }