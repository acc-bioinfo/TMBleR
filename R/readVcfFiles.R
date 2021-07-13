#' Read VCF files
#'
#' This function takes in input a (list of) VCF file name(s), reads and process
#' them by setting UCSC style and allowed chromosome names.
#'
#' @param vcfFiles list of one or more \code{character} strings corresponding to
#' the vcf file names. Each element of the list should be named with the sample name. 
#' If no name is provided, the .vcf file name will be used.
#' @param assembly human genome assembly: hg19 or hg38
#'
#' @return a (list of) \code{CollapsedVCF} object(s) containing variants
#' 
#' @examples
#'
#' # Read in input name of vcf files (provide them as a list even if you only have
#' # one file)
#' vcf_files <- list(Horizon5="Horizon5_ExamplePanel.vcf", 
#'                   HorizonFFPEmild="HorizonFFPEmild_ExamplePanel.vcf")
#'                   
#' # For each vcf file, get the absolute path
#' vcf_files <- lapply(vcf_files,
#'                     function(x) system.file("extdata", x, package = "TMBleR", mustWork = TRUE))
#' 
#' # Read in the files
#' vcfs <- readVcfFiles(vcfFiles = vcf_files, assembly = "hg19")
#'
#' @import magrittr
#' @author Laura Fancello
#'
#' @export
readVcfFiles <- function(vcfFiles, assembly){

    # Sanity Checks  -----------------------------------------------------------
    # Check arguments in input
    if(!is.list(vcfFiles)){
        stop("wrong input type: vcfFiles argument needs to be a list")}

    if ((assembly != "hg19") && (assembly != "hg38")) {
        stop("No valid genome specified: please specify 'hg19' or 'hg38'")
    }
    
    for(i in length(vcfFiles)){
        if(!(file.exists(vcfFiles[[i]]))){
            stop("file(s) do not exist: ", vcfFiles[[i]])
        }
    }
    # check if the list is named
    if(is.null(names(vcfFiles)) || length(names(vcfFiles)) != length(vcfFiles) ){
        warning("The list of vcf files was not named. Vcf files names used as 'sample' names by default")
        names(vcfFiles) <- unlist(lapply(vcfFiles, basename))
    }
    
    ## VCF VALIDATOR ----------------------------------------------------------
    ## validate each vcf file
    validator_check <<- lapply(vcfFiles, .validate_vcf, check_genotype=FALSE, AF_softcheck=TRUE)
    
    message("Validate VCF input:")
    print(validator_check)
    
    message("Warnings and Errors exported to 'validator_check' variable")
    
    ## Read vcf files ----------------------------------------------------------
    vcfs <- lapply(vcfFiles, VariantAnnotation::readVcf, assembly)

    
    # Vcf post processing
    for (i in seq_along(vcfs)) {
        
        rows_before <- dim(vcfs[[1]])[1]  # remember how many rows we had before making any change
        chr_before <- GenomeInfoDb::seqlevels(vcfs[[1]])
        
        # Set UCSC style
        GenomeInfoDb::seqlevelsStyle(vcfs[[i]]) <- "UCSC"
        # Filter chromosome names not matching.
        keepchr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                     "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                     "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
                     "chr21", "chr22", "chrX", "chrY")
        GenomeInfoDb::seqlevels(vcfs[[i]], pruning.mode="coarse") <- keepchr
        
        rows_after <- dim(vcfs[[1]])[1]
        
        # Print message in case of missmatch
        if(rows_after != rows_before){
            message(names(vcfs[1]), "\nSome variants in the .vcf file were discarded because did not have a chromosome ID maching the expected format:")
            message("Found: ")
            print(chr_before)
            message("Expected: ") 
            print(keepchr)
            message("Variants discarded: ", rows_before-rows_after, " out of ", rows_before)
        }
    }
    return(vcfs)
}
