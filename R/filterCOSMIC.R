#' Filter variants described in COSMIC
#'
#' This function removes from the input variants those described in COSMIC
#'
#' @param vcf a \code{data.frame} object containing variants
#' @param assembly human genome assembly to use: "hg19" or "hg38"
#'
#' @return Returns a \code{CollapsedVCF} object with variants not described in COSMIC
#'
#' @author Laura Fancello
#'
filterCOSMIC <- function(vcf, assembly){

    .filter_cosmic_helper <- function(original_vcf, COSMIC_file){
        # Convert cosmic into VRanges
        COSMIC_vrange = VariantAnnotation::VRanges(seqnames = COSMIC_file$CHR,
                                                   ranges  = IRanges::IRanges(COSMIC_file$START, COSMIC_file$END),
                                                   ref = COSMIC_file$REF,
                                                   alt = COSMIC_file$ALT
        )
        # remove cosmic
        vcf_out <- IRanges::subsetByOverlaps(original_vcf, COSMIC_vrange, invert = TRUE)
        return(vcf_out)
    }
    
    # Remove COSMIC variants  --------------------------------------------------
    if (assembly == "hg19") {
        
        # Is the cosmic dataset already loaded?
        if(exists("COSMIC_hg19")){
            ## Remove from input variants those present in COSMIC
            vcf_NoCosmic <- .filter_cosmic_helper(vcf, COSMIC_hg19)
        } else {
            warning("data object COSMIC_hg19 was not found. Loading a COSMIC demo dataset. Please do not use this for TMB analysis but only for demoing TMBleR")
            # Load demo data COSMIC_hg19_demo. This dataset is imported in the 
            if(!exists("COSMIC_hg19_demo")){utils::data(COSMIC_hg19_demo)}
            ## Remove from input variants those present in COSMIC
            vcf_NoCosmic <- .filter_cosmic_helper(vcf, COSMIC_hg19_demo)
        }
        
    }
    if (assembly == "hg38") {
        
        # Is the cosmic dataset already loaded?
        if(exists("COSMIC_hg38")){
            ## Remove from input variants those present in COSMIC
            vcf_NoCosmic <- .filter_cosmic_helper(vcf, COSMIC_hg38)
        } else {
            warning("data object COSMIC_hg38 was not found. Loading a COSMIC demo dataset. Please do not use this for TMB analysis but only for demoing TMBleR")
            # Load internal demo dataset COSMIC_hg38_demo. This dataset is imported
            # in the environment as COSMIC_hg38
            if(!exists("COSMIC_hg38_demo")){utils::data(COSMIC_hg38_demo)}
            ## Remove from input variants those present in COSMIC
            vcf_NoCosmic <- .filter_cosmic_helper(vcf, COSMIC_hg38_demo)
        }
    }

    return(vcf_NoCosmic)
}
