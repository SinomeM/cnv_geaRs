#' @title CNV raw data plot
#' 
#' @description 
#' \code{lrr_plot} create scatter plots of LRR and MAF 
#' in a selected region
#' 
#' @details 
#' This function can be used to inspect the raw data (LRR and MAF)
#' of putative CNV calls. It also calculate the mean LRR in a 
#' sliding window with default size of 100 SNPs and default step of
#' 50 SNPs. A bar in each plot indicate the position of the call.
#' \newline
#' For convenience consider to select the data from the entire FinalREport
#' outside R, e.g. using awk. As an example, \code{awk -F"\t" '{if ($1 == "sampleID" && $3 == 1) 
#' print $0}' FinalReport.txt > chr1_sampleID.tsv} will select the SNPs 
#' in chromosome 1 of a desired sample and save it as a separate file 
#' that can be more easily loaded in R.
#' \newline 
#' In any case \code{final.rep} must be a data-frame containing the data 
#' of one single sample and the following columns: "Chr", "Position",
#' "Log.R.Ratio", "B.allele.Freq". Order is not relevant.
#'
#' @param chr Chromosome where the putative CNV is present (e.g "1", not "chr1")
#' @param cnv.start Start position of the CNV
#' @param cnv.end End position of the CNV
#' @param region Size of the region to be plotted on each side of the CNV
#' @param final.rep Data frame in Illumina Final report format
#' @param window Sliding window size (in SNPs)
#' @param step Sliding window step (in SNPs)
#'  
#' @export
#' 
#' @return pl
#' 
#' @author Simone Montalbano simone.montalbano@protonmail.com

plot_lrr_maf <- function(chr, cnv.start, cnv.end, region = 5, final.rep,
                         window = 100, step = 50 ) {
    
    require(dplyr)
    require(magrittr)
    require(data.table)
    require(ggplot2)
    require(RcppRoll)
    require(stats)
    
    # Select only the SNPs in the correct chr. A region 1 Mbp larger than the one that will
    # be plotted
    tmp <- final.rep %>% filter(Chr == chr & (Position >= (cnv.start - (region+1)*1000000) & 
                                                  Position <= (cnv.end + (region+1)*1000000) )) %>% arrange(Position)
    
    # Sliding Window Mean on LRR 
    pos <- tmp$Position[seq(1, nrow(tmp), by=step)]
    pos2 <- tmp$Position[seq(1, nrow(tmp), by=(step/5))]
    x <- tmp$Log.R.Ratio
    val <- roll_mean(x, n = window, by = step)
    appr <- approx(x = pos[1:(length(pos)-2)], y = val, method = "constant", xout = pos2)
    appr <- data.frame(appr$x, appr$y)
    
    # Region to be plotted
    dat <- tmp %>% filter(Position >= (cnv.start - (region)*1000000) & 
                              Position <= (cnv.end + (region)*1000000))
    appr %<>% filter(appr.x >= (cnv.start - (region)*1000000) & 
                         appr.x <= (cnv.end + (region)*1000000))
    
    # Plot LRR 
    pl1 <- ggplot(dat, aes(x = Position/1000000, y = Log.R.Ratio)) +
        geom_point(size = 1, colour = "#FF9999") +
        geom_line(data = appr, aes(x = appr.x/1000000, y = appr.y), colour = "black") + 
        geom_segment(aes(x = cnv.start/1000000, xend = cnv.end/1000000, y = 0, yend = 0),
                     size = 0.5, colour = "#FF0000") +
        xlab("Position (Mbp)") +
        ylab("LRR") + 
        theme_bw() + theme(panel.border=element_blank()) +
        scale_x_continuous(labels = scales::comma)
    
    # Plot MAF 
    pl2 <- ggplot(dat, aes(x = Position/1000000, y = B.allele.Freq)) +
        geom_point(size = 1, colour = "#33CCCC") +
        geom_segment(aes(x = cnv.start/1000000, xend = cnv.end/1000000, y = 0.5, yend = 0.5),
                     size = 0.5, colour = "#33CCFF") +
        xlab("Position (Mbp)") +
        ylab("BAF") +
        theme_bw() + 
        theme(panel.border=element_blank()) +
        scale_x_continuous(labels = scales::comma)
    
    # Combine 
    pl <- cowplot::plot_grid(pl1, pl2, nrow = 2)
    
    return(pl)
}
