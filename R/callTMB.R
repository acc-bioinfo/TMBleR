#' Perform Tumor Mutational Burden (TMB) quantification
#'
#' This function counts the total number of mutations in the object provided in
#' input with variants. It gives in output both the total number of mutations and
#' the number of mutations per Mb (dividing by the sequenced genomic space,
#' extracted from the sequencing design).
#'
##' @param inputForTMB a \code{list} with the following elements: sample name, 
#' filter description, design and a \code{CollapsedVCF} or \code{data.frame} 
#' object containing somatic variants
#' @param genome human genome assembly: hg19 or hg38
#'
#' @return Returns a \code{table} with sample name, filter description, size of
#' the sequenced genomic space and TMB expressed both as total number of
#' mutations in the sequenced space and number of mutations per megabase
#'
#' @author Laura Fancello
#'
callTMB <- function(inputForTMB, genome){

    # Sanity Checks  -----------------------------------------------------------
    ## Check and read input arguments
    if (is.null(inputForTMB)) {
        stop("argument 'inputForTMB' is missing, with no default")
    }
    if (is.null(genome)) {
        stop("argument 'genome' is missing, with no default: please specify 'hg19' or 'hg38'")
    }
    if ((genome != "hg19") && (genome != "hg38")) {
        stop("No valid genome assembly: please indicate 'hg19' or 'hg38'")
    }
    
    # Preprocess input  --------------------------------------------------------
    sample <- inputForTMB$sample
    filter <- inputForTMB$filter
    design <- inputForTMB$design
    design_name <- design@metadata$design_name
    vcf <- inputForTMB$variants

    if (!(methods::is(vcf)[1] == "data.frame") && !(methods::is(vcf)[1] == "CollapsedVCF") && !(methods::is(vcf)[1] == "GRanges")){
        "Input does not contain valid object for variants:
                     please provide valid CollapsedVCF, data.frame or GRanges object"
    }
    
    if (methods::is(vcf)[1] == "CollapsedVCF") {
        if(length(SummarizedExperiment::SummarizedExperiment(vcf))==0){
            vcf <- 0
        }else{
            GenomeInfoDb::seqlevelsStyle(vcf) <- "UCSC"
            vcf <- SummarizedExperiment::rowRanges(vcf)
        }
    }
    if (methods::is(vcf)[1] == "data.frame") {
        if(dim(stats::na.omit(vcf))[1]==0){
            vcf <- 0
        }
            else{
                vcf$CHR <- as.character(as.vector(vcf$CHR))
                vcf$START <- as.character(as.vector(vcf$START))
                vcf$END <- as.character(as.vector(vcf$END))
                vcf$REF <- as.character(as.vector(vcf$REF))
                vcf$ALT <- as.character(as.vector(vcf$ALT))
                vcf <- GenomicRanges::makeGRangesFromDataFrame(vcf, 
                                                           seqnames.field="CHR", 
                                                           start.field="START", 
                                                           end.field="END", 
                                                           keep.extra.columns = FALSE)
                GenomeInfoDb::seqlevelsStyle(vcf) <- "UCSC"
            }
    }
    if (methods::is(vcf)[1] == "GRanges" && nrow(as.data.frame(vcf)) == 0 ) {
        vcf <- 0
    }

    # TMB quantification  -----------------------------------------------------------
    ## Calculate sequencing size using design
    size <- sum(IRanges::width(GenomicRanges::reduce(design))) / 1000000

    ## Calculate TMB
    if(identical(vcf, 0)){
        TMB <- 0
        TMB_perMb <- 0
    }else{
        TMB <- length(vcf)
        TMB_perMb <- round(TMB / size, digits=2)
    }
    
    ## Output a matrix with TMB values
    TMB_results <- c(sample, design_name, filter, size, TMB, TMB_perMb)
    #names(TMB_results) <- c("sample", "design", "filter", "seq_size", "TotNMut", "TMB_perMb")
    #TMBs_matrix <- as.data.frame(matrix(unlist(TMB_results), ncol=5, byrow = TRUE))
    #head <- c("sample", "filter", "seq_size", "Tot_N_Mutations", "TMB_per_Mb")
    #colnames(TMBs_matrix) <- head
    return(TMB_results)
}
