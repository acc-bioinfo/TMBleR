#' .fetch_filtered_SNV
#' 
#' This is just a helper function that can be used to retrieved the SNV
#' that have been filtered out
#'
#'
#' @param original_vcf Granges object
#' @param filtered_vcf Granges object
#'
#' @return vcfcollabpsed object
#' 
.fetch_filtered_SNV <- function(original_vcf, filtered_vcf){
    
    out <- IRanges::subsetByOverlaps(original_vcf, filtered_vcf, invert = TRUE)

    # sanity check
    if((length(out) + length(unique(filtered_vcf))) != length(original_vcf)){
        stop("Unexpected problem while trying to retrieve SNV in non exonic regions")
    }
    
    # make it into a dataframe
    out <- as.data.frame(out@rowRanges) %>%
        dplyr::rename("CHR" = .data$seqnames) %>%
        dplyr::rename("START" = .data$start) %>%
        dplyr::rename("END" = .data$end) %>%
        dplyr::mutate("REF" = ".") %>%
        dplyr::mutate("ALT" = ".") %>%
        dplyr::select(-.data$width, -.data$strand, -.data$paramRangeID)
    return(out)
}
