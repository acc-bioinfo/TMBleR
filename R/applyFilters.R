#' Apply filters on a list of CollapsedVCF objects
#'
#' This function calls the callFilters() function on each element of a list of
#' CollapsedVCF objects and generates in output a list of lists. Each list
#' corresponds to one input element, and consists in four elements: an object
#' containg variants which passed the filter, a character string describing the
#' applied filter (if any), an object containing the sequencing design and a
#' character string with the name of the sample.
#'
#' @param vcfs a \code{list} of one or more \code{CollapsedVCF} object(s)
#' @param assembly human genome assembly: hg19 or hg38
#' @param vaf.cutoff minimum value of variant allele frequency accepted
#' @param remove.cancer logical value 'TRUE' or 'FALSE' indicating whether or
#' not to remove cancer variants (variants described in COSMIC and
#' truncating mutations in tumor suppressors)
#' @param remove.nonexonic logical value 'TRUE' or 'FALSE' indicating whether or
#' not SNV mapped ouside of exons are to be removed. Default is True
#' @param tsList path to file containing list of tumor suppressors. If not
#' provided a list of 1217 tumor suppressors from the TSgene2 database
#' (<http://bioinfo.mc.vanderbilt.edu/TSGene/>) is used by default.
#' @param variantType type of variant to remove. Possible values: synonymous,
#' nonsynonymous, frameshift, nonsense, not translated or a  combination of them
#' specified in a \code{character vector}
#' @param design a \code{GRanges} object containing WES or panel design
#'
#' @return Returns a \code{list} of \code{lists}. Each \code{list} include the
#' following elements: a \code{GRanges}, \code{CollapsedVCF}, \code{data.frame}
#' object containing variants passing the filter, a \code{charcater string}
#' describing applied filter (if any), a \code{GRanges} or \code{character vector}
#' with the sequencing design, a \code{character string} with the name of the
#' sample.
#'
#' @examples
#' 
#' ## Read vcf
#' vcf_files <- list(Horizon5="Horizon5_ExamplePanel.vcf", 
#'                   HorizonFFPEmild="HorizonFFPEmild_ExamplePanel.vcf")
#' vcf_files <- lapply(vcf_files
#'                     , function(x) system.file("extdata", x, package = "TMBleR", mustWork = TRUE))
#' vcfs <- readVcfFiles(vcfFiles = vcf_files, assembly = "hg19")
#' 
#' ## Read design
#' design <- readDesign(system.file("extdata"
#' , "ExamplePanel_GeneIDs.txt"
#' , package = "TMBleR"
#' , mustWork = TRUE)
#' , assembly = "hg19"
#' , ids = "entrezgene_id")
#' 
#' ## Apply filter to remove known cancer variants using the default tumor
#' ## suppressors list
#' vcfs_NoCancer <- applyFilters(vcfs = vcfs
#' , assembly = "hg19"
#' , design = design
#' , remove.cancer = TRUE
#' , tsList = NULL
#' , variantType = NULL)
#'
#' ## Apply filter to remove synonymous mutations
#' vcfs_filtered <- applyFilters(vcfs = vcfs,
#'                              assembly = "hg19",
#'                              design = design,
#'                              remove.cancer = FALSE,
#'                              tsList = NULL,
#'                              variantType = "synonymous")
#'
#'
#' @author Laura Fancello
#'
#' @export
applyFilters <- function(vcfs,
                        assembly,
                        design,
                        vaf.cutoff = 0,
                        remove.nonexonic = TRUE, 
                        remove.cancer = FALSE,
                        tsList = NULL,
                        variantType = c()
                        ){

    # Check input arguments
    if (is.null(vcfs)) {stop("argument \"vcfs\" is missing, with no default")}
    
    for(i in length(vcfs)){
        if (!(methods::is(vcfs[[i]])[1] == "CollapsedVCF")) {
            stop("No valid vcf provided: please provide a CollapsedVCF object")
        }
    }
    
    # QC on VAF
    if (!is.numeric(vaf.cutoff)){stop("vaf.cutoff should be numeric type")}
    if (vaf.cutoff > 1){stop("vaf.cutoff can't be higher than 1")}
    if (vaf.cutoff < 0){stop("vaf.cutoff can't be a negative value")}
        
    if ((assembly != "hg19") && (assembly != "hg38")) {
        stop("No valid genome assembly: please specify 'hg19' or 'hg38'")
    }
    if (is.null(design)) {
        stop("argument 'design' is missing, with no default")
    }
    if (is.null(remove.cancer)) {
        stop("argument \"remove.cancer\" is missing with no default: please specify TRUE or FALSE")
    }

    
    # No Filter called scenario -----------------------------------------------
    if(vaf.cutoff == 0 & is.null(tsList) & is.null(variantType)  & remove.cancer == FALSE & remove.nonexonic == FALSE){
        # Generate input for applyTMB function  
        filter="NoFilter"
        # Apply inputToTMB function on each element of the input list
        vcfs_unfiltered <- lapply(vcfs,
                                  inputToTMB,
                                  design,
                                  filter)
        
        # Append to each list of the output list of lists an element "sample" with
        # the sample name
        for(i in seq_along(vcfs_unfiltered)) {
            vcfs_unfiltered[[i]][["sample"]]=names(vcfs_unfiltered)[i]
        }
        
        return(vcfs_unfiltered)
    }
    
    # Else, Apply filtering ----------------------------------------------------
    # Apply callFilters function on each element of the input list
    vcfs_filtered <- lapply(vcfs,
                            callFilters,
                            assembly,
                            design,
                            vaf.cutoff,
                            remove.cancer,
                            remove.nonexonic,
                            tsList,
                            variantType)

   
    for(i in seq_along(vcfs_filtered)) {
        # Append to each list of the output list of lists an element "sample" with
        # the sample name
        vcfs_filtered[[i]][["sample"]] <- names(vcfs_filtered)[i]
    }
    
    return(vcfs_filtered)
}
