#' Provide a suitable input to the callTMB function for TMB quantification
#'
#' This function reads in input an object containing mutations and generates in
#' output a list with the following elements: filter description,
#' an object with sequencing design and an object containing variants.
#'
#' @param vcf  a \code{CollapsedVCF} object containing variants
#' @param design a \code{GRanges} object containing WES or panel sequencing
#' design
#' @param filter a \code{character string} describing the type of filter
#' previously applied to input variants
#'
#' @return Returns a list with filter description, a design object
#' and an object containing variants, as required in input by the
#' \code{\link{callTMB}} function
#'
#' @author Laura Fancello
#'
inputToTMB <- function(vcf, design, filter){

    # Sanity Checks  -----------------------------------------------------------
    ## Check input arguments
    if (is.null(vcf)) {
        stop("argument 'vcf' is missing, with no default")
    }
    if (!(methods::is(vcf)[1] == "CollapsedVCF")) {
                        stop("Input does not contain valid object for variants:
                             please provide valid CollapsedVCF object")
    }
    if (is.null(design)) {
        stop("argument 'design' is missing, with no default: please provide a GRanges object containing
             the sequencing design")
    }
    if (is.null(filter)) {
        stop("argument 'filter' is missing, with no default")
    }

    # Preprocess input  ---------------------------------------------------------
    ## Remove off target regions by overlapping with GRanges object with design
    # Transform output object of filtered variants in GRanges, to be able to give it as input to subsetByOverlaps
    vcf_GR <- SummarizedExperiment::rowRanges(vcf)
    # Remove off targets
    vcf_ok <- IRanges::subsetByOverlaps(vcf_GR, design)
    # Go back to data.frame
    vcf_ok <- as.data.frame(vcf_ok, row.names = NULL)
    vcf_ok=vcf_ok[,c(-4,-5,-6,-9,-10)]
    colnames(vcf_ok) <- c("CHR", "START", "END", "REF", "ALT")
    
    # Generate suitable input to applyTMB function  ----------------------------
    ## Generate a list with the following elements: filter description, design
    ## and object containing variants
    input_TMB <- list(variants=vcf_ok, design=design, filter=filter)
    
    return(input_TMB)
}
