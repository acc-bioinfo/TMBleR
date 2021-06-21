#' Remove cancer variants
#'
#' This function removes from the vcf dataset those variants defined as pathogenic
#'  (cancer) in COSMIC or truncating variants in tumor suppressor genes.
#'
#' @param vcf a \code{CollapsedVCF} object containing variants
#' @param assembly human genome assembly: hg19 or hg38
#' @param tsList path to file containg list of tumor suppressors. If not
#' provided a list of 1217 tumor suppressors from the TSgene2 database
#'  is used.
#'
#' @return Returns a \code{data.frame} object containing only those variants
#' which passed the filter
#'
#' @author Laura Fancello
#'
removeCancerVariants=function(vcf, assembly, tsList){

    # Sanity Checks  -----------------------------------------------------------
    ## Check input arguments
    if (is.null(vcf)) {
        stop("argument 'vcf' is missing, with no default")
    }
    if (!(methods::is(vcf)[1] == "CollapsedVCF")) {
        stop("No valid vcf provided: please provide a CollapsedVCF object")
    }
    if (is.null(assembly)) {
        stop("argument 'assembly' is missing, with no default: please specify 'hg19' or 'hg38'")
    }
    if ((assembly != "hg19") && (assembly != "hg38")) {
        stop("No valid genome assembly: please specify 'hg19' or 'hg38'")
    }
    if (is.null(tsList)) {
        stop("No tsList provided.")
    }

    # Remove cancer variants  --------------------------------------------------
    ## Apply function to remove truncating mutations in tumor suppressor genes
    vcf_NoTruncatingInTS <- filterTruncatingTumorSuppr(vcf, assembly, tsList)

    ## Apply function to remove variants described in COSMIC
    if(dim(vcf_NoTruncatingInTS)[1]==0){  # empty output
      vcf_NoCosmic <- vcf_NoTruncatingInTS  
    }
    else{
      vcf_NoCosmic <- filterCOSMIC(vcf_NoTruncatingInTS, assembly)
    }
    
    return(vcf_NoCosmic)
}
