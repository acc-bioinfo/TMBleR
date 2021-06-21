#' Provide text and graphical output for Tumor Mutational Burden (TMB)
#' quantification
#'
#' This function reads in input a list (or a list of lists), with sample name,
#' filter description, sn object containing sequencing design and an object
#' containing variants. It quantifies TMB based on these variants and generates
#' in output a text file containing TMB values and a pdf file with a barplot
#' visualization.
#'
#' @param inputForTMB a \code{list} of lists with the following elements: sample name, 
#' filter description, design and a \code{CollapsedVCF} or \code{data.frame} 
#' object containing somatic variants
#' @param assembly human genome assembly: hg19 or hg38
#'
#' @return Returns a data.frame with sample, filter, sequencing size,
#' total number of mutations and number of mutations per megabase.
#'
#' @examples
#'
#' ## Read vcf
#' vcf_files <- list(Horizon5="Horizon5_ExamplePanel.vcf", 
#'                   HorizonFFPEmild="HorizonFFPEmild_ExamplePanel.vcf")
#' vcf_files <- lapply(vcf_files
#'                     , function(x) system.file("extdata", x, package = "TMBleR", mustWork = TRUE))
#' vcfs <- readVcfFiles(vcfFiles = vcf_files , assembly = "hg19")
#' 
#' ## Read design
#' design <- readDesign(system.file("extdata"
#' , "ExamplePanel_GeneIDs.txt"
#' , package = "TMBleR"
#' , mustWork = TRUE)
#' , assembly = "hg19"
#' , ids = "entrezgene_id")
#' 
#' ## Prepare data for TMB quantification, formatting it as required by the 
#' ## applyTMB() function. While this is done automatically within the 
#' ## applyFilters() function, on a non filtered vcf we need to use the function
#' ## here described to get the correct format
#' vcfs_nonfiltered <- applyFilters(vcfs = vcfs, 
#'                                 design = design,
#'                                 assembly = "hg19")
#'                                     
#' ## Perform TMB quantification
#' TMB_res=applyTMB(inputForTMB = vcfs_nonfiltered, assembly = "hg19")
#'
#'
#' @author Laura Fancello
#'
#'
#' @export
applyTMB <- function(inputForTMB, assembly){
    
    # Sanity Checks  ----------------------------------------------------------
    ## Verify input arguments
    if (is.null(inputForTMB)) {
        stop("argument 'inputForTMB' is missing, with no default")
    }
    if (is.null(assembly)) {
        stop("argument 'assembly' is missing, with no default: please specify 'hg19' or 'hg38'")
    }
    if ((assembly != "hg19") && (assembly != "hg38")) {
        stop("No valid genome assembly: please indicate 'hg19' or 'hg38'")
    }

    # TMB quantification  ------------------------------------------------------
    TMBs <- lapply(inputForTMB, callTMB, assembly)

    if (unique(unlist(lapply(TMBs, length))) != 6){
        stop("missing column while computing TMB. Something went wrong. ")
    }

    ## Write results in a table
    TMBs_matrix <- as.data.frame(matrix(unlist(TMBs), ncol=6, byrow = TRUE))
    colnames(TMBs_matrix) <- c("Sample",
                               "Design",
                               "Filter",
                               "Sequencing_Size",
                               "Tot_Number_Mutations",
                               "TMB_per_Mb")
    
    TMBs_matrix$Sequencing_Size=as.numeric(as.vector(TMBs_matrix$Sequencing_Size))
    TMBs_matrix$Tot_Number_Mutations=as.numeric(as.vector(TMBs_matrix$Tot_Number_Mutations))
    TMBs_matrix$TMB_per_Mb=as.numeric(as.vector(TMBs_matrix$TMB_per_Mb))
    
    return(TMBs_matrix)
}


################################################################################
