#' .remove_noncoding_SNV
#'
#' This helper function removes from a VariantAnnotation VCF object the 
#' SNV mapped in non coding regions
#'
#' @param vcf VariantAnnotation vcf object
#' @param assembly hg19 or hg38
#'
#' @return Granges object
.remove_noncoding_SNV <- function(vcf, assembly){
        
    if(is.null(GenomicRanges::elementMetadata(vcf)$CONSEQUENCE)){
        vcf_exp <- VariantAnnotation::expand(vcf)
        
        # Generate the BSGenome and txdb objects required for the
        # VariantAnnotation::predictCoding function
        if(assembly == "hg19") {
            BSgenome <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
            txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
        }
        if(assembly == "hg38") {
            BSgenome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
            txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
        }
        # Apply the VariantAnnotation::predictCoding function which gives in
        # output coding mutations
        GENETIC_CODE <- Biostrings::getGeneticCode(id_or_name2 = "1")
        vcf_gr <-  suppressWarnings(  # This is not good but the warning seem to be generated also 
            # in the vignette of the VariantAnotatoni package.
            VariantAnnotation::predictCoding(vcf_exp,
                                             txdb,
                                             seqSource=BSgenome, genetic.code=GENETIC_CODE)
        )
        
        # return vcf object filtered
        vcf <- IRanges::subsetByOverlaps(vcf, vcf_gr, type="within")
        
        return(vcf)
        # Return a Granges object
        #return(GenomicRanges::makeGRangesFromDataFrame(as.data.frame(vcf_ok)))
    } 
}