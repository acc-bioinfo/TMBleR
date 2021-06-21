#' Filter variants according to variant allele frequency
#'
#' This function removes from the input variants those with variant allele
#' frequency (VAF) inferior to the cutoff value provided. The AF field is
#' required in the input vcf file to apply this function.
#'
#' @param vcf a \code{CollapsedVCF} object containing somatic variants
#' @param vaf.cutoff minimum value of variant allele frequency accepted
#'
#' @return Returns a \code{CollapsedVCF} object containing only variants with
#' variant allele frequency above the cutoff
#'
#' @author Laura Fancello
#'
filterByVAF <- function(vcf, vaf.cutoff){

    # Sanity Checks  -----------------------------------------------------------
    ## Check the input arguments
    if (!(methods::is(vcf)[1] == "CollapsedVCF")) {
        stop("No valid vcf provided.")
    }
    if (is.null(vaf.cutoff)) {
        stop("argument 'vaf.cutoff' is missing, with no default")
    }
    if(!(methods::is(vaf.cutoff)[1]=="numeric")){
        stop("No valid vaf.cutoff provided: please indicate a numeric value")
    }

    # Filter by VAF  -----------------------------------------------------------
    ## Remove variants not passing VAF filter
    if (!(is.null(VariantAnnotation::geno(vcf)$AF))) { # Check if allele frequency field (AF) is present
        # make sure that the metafield "number" in the input vcf  is set to "A"
        #       ##FORMAT=<ID=DP,Number=A
        vaf <- unname(rapply(VariantAnnotation::geno(vcf)$AF[,1], function(x) x[1], how="unlist"))
        vcf_filtered <- vcf[vaf > vaf.cutoff]
    } else {stop("No AF (Allele Frequency) field found in vcf file")}
    
    return(vcf_filtered)
}

