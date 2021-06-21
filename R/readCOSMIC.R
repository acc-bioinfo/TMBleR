#' Read in a COSMIC file and format it for usage in TMBleR
#' 
#' The objective of this script is to convert a COSMIC .vcf file into a R object 
# compatible with TMBleR. For the instruction on what COSMIC file to use and how 
# to acquire it, please refert to the vignette: HOWTO_Import_external_data
#'
#' @param input_file path to CosmicCodingMuts.vcf file downloaded following the vignette instructions
#' @param assembly either hg19 or hg38
#' @param output_file path to where to save the converted .vcf file into a .rda 
#'
#' @importFrom rlang .data
#' @return nothing
#' @export
formatCOSMIC <- function(input_file, assembly, output_file="COSMIC"){
    
    # sanity check ------------------------------------------------------------
    if ((assembly != "hg19") && (assembly != "hg38")) {
        stop("No valid genome specified: please specify 'hg19' or 'hg38'")
    }
    
    if(!file.exists(input_file)){stop("file not found: ", input_file)}
    
    if(tools::file_ext(input_file) != "vcf"){
        stop("file must be a vcf file (have the .vcf file extension): "
             , input_file)
    }
    
    # Check parent dir is writable
    if(file.access(dirname(output_file), mode = 2) == -1){
        stop("File can not be written: ", output_file)
    }
    

    COSMIC <- VariantAnnotation::readVcf(file=input_file, assembly)
    COSMIC <- methods::as(COSMIC, "VRanges")
    
    # Correct formatting ------------------------------------------------------
    GenomeInfoDb::seqlevelsStyle(COSMIC) <- "UCSC"
    GenomeInfoDb::seqlevels(COSMIC, pruning.mode="coarse") <- c(
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
        "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
        "chr13", "chr14", "chr15","chr16", "chr17", "chr18",
        "chr19", "chr20", "chr21","chr22", "chrX", "chrY")
    COSMIC_df <- data.frame(COSMIC) %>%
        dplyr::select(.data$seqnames, .data$start, .data$end, .data$strand, .data$ref, .data$alt) %>%
        dplyr::rename(CHR=.data$seqnames,
                      START=.data$start,
                      END=.data$end,
                      STRAND=.data$strand,
                      REF=.data$ref,
                      ALT=.data$alt)
    
    # COSMIC_hg19 <- data.frame(CHR=COSMIC_df$seqnames,
    #                           START=COSMIC_df$start,
    #                           END=COSMIC_df$end,
    #                           STRAND=COSMIC_df$strand,
    #                           REF=COSMIC_df$ref,
    #                           ALT=COSMIC_df$alt)
    
    # prepara saving
    output_file <- paste0(output_file, ".rda")
    message("Saving: ", output_file)

    if(assembly == "hg19"){
        COSMIC_hg19 <- COSMIC_df                # load in to memory
        save(COSMIC_hg19, file = output_file)   # save to file
        message("COSMIC dataset exported")
    }
    if(assembly == "hg38"){
        COSMIC_hg38 <- COSMIC_df                # load in to memory
        save(COSMIC_hg38, file = output_file)   # save to file
        message("COSMIC dataset exported")
    }
    
    # free up memory
    rm(COSMIC_df)
    rm(COSMIC)
}

# Tests
#formatCOSMIC( input_file = "~/Downloads/CosmicCodingMuts.vcf", "hg19")
#formatCOSMIC( input_file = "~/Downloads/CosmicCodingMuts.vcf", "hg10")
#formatCOSMIC( input_file = "~/Downloads/CosmicCodingMuts", "hg19")
#formatCOSMIC( input_file = "~/Downloads/asdadasd", "hg19")
#formatCOSMIC( input_file = "~/Downloads/CosmicCodingMuts.vcf", "hg19", output_file = "~/Downloads/deleteme")
