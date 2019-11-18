#' @title LRR plot
#' 
#' @description \code{lrr_plot} crate scatter plots of LRR values
#' 
#' @details This function can be used to inspect the raw data (LRR)
#' of putative CNV calls 
#' \newline
#' DA FINIRE!!
#' 
#' @author Simone Montalbano simone.montalbano@protonmail.com

## fare documentazione 
##     ## AGGIUNGI UNITA' DI MISURA su x
##     aggiungi le assi 

plot_lrr_maf <- function(sam, chr, cnv.start, cnv.end, region = 5, final.rep,
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
    
    # Plot LRR 
    pl1 <- ggplot(dat, aes(x = Position, y = Log.R.Ratio)) +
        geom_point(size = 1, colour = "#FF9999") +
        geom_line(data = appr, aes(x = appr.x, y = appr.y), colour = "black") + 
        geom_segment(aes(x = cnv.start, xend = cnv.end, y = 0, yend = 0),
                     size = 0.2, colour = "#FF0000")
    
    # Plot MAF 
    pl2 <- ggplot(dat, aes(x = Position, y = B.allele.Freq)) +
        geom_point(size = 1, colour = "#33CCCC") +
        geom_segment(aes(x = cnv.start, xend = cnv.end, y = 0.5, yend = 0.5),
                     size = 0.2, colour = "#33CCFF")
    
    # Combine 
    pl <- cowplot::plot_grid(pl1, pl2, nrow = 2)
    
    return(pl)
}
