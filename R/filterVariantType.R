#' Filter variants by the specified type
#'
#' This function filters out from the variants of the input vcf file those
#' corresponding to the type specified in the variantType argument: synonymous,
#' nonsynonymous, frameshift, nonsense, not translated or a combination of them.
#'
#' @param vcf a \code{CollapsedVCF} or \code{data.frame} object containing
#' variants
#' @param assembly human genome assembly: "hg19" or "hg38"
#' @param variantType type of variant to remove. Possible values are: synonymous,
#' nonsynonymous, frameshift, nonsense, not translated or a combination of them.
#'
#' @return Returns a \code{CollapsedVCF} or \code{data.frame} object containing
#' those variants which passed the filter.
#'
#' @author Laura Fancello
#'
filterVariantType <- function(vcf, assembly, variantType){

    # Sanity Checks  -----------------------------------------------------------
    ## Check input arguments
    if (!(methods::is(vcf)[1] == "CollapsedVCF")) {
        if (!(methods::is(vcf)[1] == "data.frame")) {
            stop("No valid vcf: please provide a CollapsedVCF object.")
        }
    }
    if (is.null(variantType)) {
        stop("argument 'variantType' is missing, with no default")
    }
    if(assembly == "hg19") {
        BSgenome <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    }
    if(assembly == "hg38") {
        BSgenome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    }
    if ((assembly != "hg19") && (assembly != "hg38")) {
        stop("No valid genome assembly: please specify 'hg19' or 'hg38'")
    }    
    
    # Preprocess input  --------------------------------------------------------
    if (methods::is(vcf)[1] == "CollapsedVCF") {
        
        ## Set UCSC style and allowed chromosome names.
        GenomeInfoDb::seqlevelsStyle(vcf) <- "UCSC"
        GenomeInfoDb::seqlevelsStyle(txdb) <- "UCSC"
        keepchr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                     "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
                     "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                     "chr20", "chr21", "chr22", "chrX", "chrY")
        GenomeInfoDb::seqlevels(vcf) <- keepchr
        GenomeInfoDb::seqlevels(txdb) <- keepchr
    
        ## Expand VCF
        vcf_exp <- VariantAnnotation::expand(vcf)
    
        ## Annotate
        ## NOte: this will annotated SNV mapped in coding regions. The SNV not 
        ## mapped in the exonic reigons will be filter out silently
        GENETIC_CODE <- Biostrings::getGeneticCode(id_or_name2 = "1")
        vcf_annotated <- suppressWarnings(  # This is not good but the warning seem to be generated also 
                                                # in the vignette of the VariantAnnotation package.
            VariantAnnotation::predictCoding(vcf_exp
                                            , txdb
                                            , seqSource=BSgenome
                                            , genetic.code=GENETIC_CODE)
        )
        
       # Remove mutations specified in the variantType argument  --------------
        vcf_filtered <- vcf_annotated[which(!(GenomicRanges::elementMetadata(vcf_annotated)$CONSEQUENCE %in% variantType))]
        #if(getOption("debug")){
        if(debug_env$debug){
            print("Returning variable 'debug_env$variantType' for debugging")
            debug_env$variantType <- vcf_annotated
        }
        
        # CONVERT TO vcf format
        # ------------------
        # We need the output of this function to be a collapsedVCF, an object 
        # from the VariantAnnotation package. This is a trick to convert it
        # 
        # Filter in this direction first to make sure that we do not include frameshifts
        # that fall inside other SNVs. They would be included by mistake
        vcf_out <- IRanges::subsetByOverlaps(vcf, vcf_filtered, type="equal")
        
        return(vcf_out)
        

    } else {
        stop("Error: unexpected vcf file format")
    }
    
    
    return(vcf_out)
}
