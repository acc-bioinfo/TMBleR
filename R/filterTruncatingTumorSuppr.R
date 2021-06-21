#' Filter truncating variants in tumor suppressors
#'
#' This function removes from the input vcf truncating variants in tumor
#' suppressor genes
#'
#' @param vcf \code{CollapsedVCF} object containing variants
#' @param assembly human genome assembly: hg19 or hg38
#' @param tsList path to file containg list of tumor suppressors. If not
#' provided a list of 1217 tumor suppressors from the TSgene2 database
#' (<https://bioinfo.uth.edu/TSGene/>) is used.
#'
#' @return Returns a \code{data.frame} object with those variants which do not
#' represent truncating mutations in tumor suppressors
#'
#' @author Laura Fancello
#'
filterTruncatingTumorSuppr<- function(vcf, assembly, tsList) {
    
    # Sanity Checks  ----------------------------------------------------------
    # Check input arguments
    if (assembly == "hg19") {
        BSgenome <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    }
    if (assembly == "hg38") {
        BSgenome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    }
    if ((assembly != "hg19") && (assembly != "hg38")) {
        stop("No valid genome assembly: please specify 'hg19' or 'hg38'")
    }

    # Preprocess input  --------------------------------------------------------
    ## Set UCSC style and allowed chromosome names.
    GenomeInfoDb::seqlevelsStyle(txdb) <- "UCSC"
    GenomeInfoDb::seqlevels(txdb, pruning.mode="coarse") <- c("chr1", "chr2", 
                                                              "chr3", "chr4",
                                                "chr5", "chr6", "chr7", "chr8",
                                                "chr9", "chr10", "chr11",
                                                "chr12", "chr13", "chr14",
                                                "chr15", "chr16", "chr17",
                                                "chr18", "chr19", "chr20",
                                                "chr21", "chr22", "chrX",
                                                "chrY")

    ## Annnotate the type of amino acid change for nonsynonymous variants
    vcf_exp <- VariantAnnotation::expand(vcf)
    GENETIC_CODE <- Biostrings::getGeneticCode(id_or_name2 = "1")

    vcf_annotated <- suppressWarnings(
        VariantAnnotation::predictCoding(query =vcf_exp
                                         , subject = txdb
                                         , seqSource = BSgenome
                                         , genetic.code = GENETIC_CODE)
    )
    
    # Remove truncating mutations in tumor suppressors -------------------------
    truncating <- c("nonsense", "frameshift")
    whichTS <- which(GenomicRanges::elementMetadata(vcf_annotated)$GENEID %in% tsList)
    whichTrunc <- which(GenomicRanges::elementMetadata(vcf_annotated)$CONSEQUENCE %in% truncating)
    which_TS_and_Trunc=intersect(whichTS, whichTrunc)
    if(length(which_TS_and_Trunc)>0){
        vcf_filtered <- vcf_annotated[-(which_TS_and_Trunc)]
    }else{
        vcf_filtered <- vcf_annotated
    }
    
    # Retrieve vcf format
    out <- IRanges::subsetByOverlaps(vcf, vcf_filtered)
    return(out)    
}
