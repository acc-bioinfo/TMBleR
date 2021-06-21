#' Filters variants according to input arguments
#'
#' This function calls different filtering functions on the input vcf object
#' according to the input arguments and gives in output a list with the three
#' elements: an object containg variants which passed the filter, a character
#' string describing the applied filter (if any) and an object containing the
#' sequencing design
#'
#' @param vcf \code{CollapsedVCF} object
#' @param assembly human genome assembly: hg19 or hg38
#' @param design a \code{GRanges} object containing WES or panel design
#' @param vaf.cutoff minimum value of variant allele frequency accepted
#' @param remove.cancer logical value 'TRUE' or 'FALSE' indicating whether or
#'  not to remove cancer variants, which is variants described in COSMIC or
#'  truncating mutations in tumor suppressors
#' @param tsList path to file containing list of tumor suppressors. If not
#' provided a list of 1217 tumor suppressors from the TSgene2 database
#' (<https://bioinfo.uth.edu/TSGene/>) is used.
#' @param variantType type of variant to remove: synonymous, frameshift or
#'  nonsense
#' @param remove.nonexonic logical value 'TRUE' or 'FALSE' indicating whether or
#' not to remove SNV mapped in non exonic regions 
#' @return Returns a \code{list} with the following elements: a \code{GRanges},
#' \code{CollapsedVCF}, \code{data.frame} object containing variants passing the
#' filter, a \code{charcater string} describing applied filter (if any), a
#' \code{GRanges} or \code{character vector} with the sequencing design,.
#'
#' @importFrom stats filter
#' @author Laura Fancello
callFilters <- function(vcf,
                        assembly,
                        design,
                        vaf.cutoff,
                        remove.cancer,
                        remove.nonexonic,
                        tsList,
                        variantType
                        ){
    
    # ##########################################################################
    # Apply filtering functions  -----------------------------------------------
    # ##########################################################################
    
    # FOR DEBUGGING
    debug <- debug_env$debug
    vcf_filtered_df <-  cbind(.empty_vcf(), as.data.frame(cbind(filter = as.character())))
  
    .debug_filter <- function(vcf_before_filter, vcf_after_filter, vcf_filtered, filter_label=""){
      tmp_filter <- IRanges::subsetByOverlaps(vcf_before_filter, vcf_after_filter, type="within", invert = TRUE)
      out <- as.data.frame(tmp_filter@rowRanges) %>%
        dplyr::mutate(REF = ".", ALT = ".") %>%
        dplyr::rename(CHR = "seqnames") %>%
        dplyr::rename(START = "start") %>%
        dplyr::rename(END = "end") %>%
        dplyr::select(-.data$paramRangeID, -.data$strand, -.data$width) %>%
        dplyr::mutate(filter = filter_label) %>%
        rbind(.data, vcf_filtered)
      return(out)
    }
    
    
    # REMOVE SNV with WRONG Annotation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Insertions and deletion need to be correctly formatted as described in
    # https://www.internationalgenome.org/wiki/Analysis/Variant%2520Call%2520Format/vcf-variant-call-format-version-40/
    vcf_backup <- vcf  #make a copy
    idx_width_0 = which(GenomicRanges::width(vcf@rowRanges) == 0)
    if(length(idx_width_0)>0){
      vcf <- vcf[-idx_width_0]
      if(debug){
        vcf_filtered_df <- .debug_filter(vcf_before_filter = vcf_backup
                                         , vcf_after_filter = vcf
                                         , vcf_filtered = vcf_filtered_df
                                         , filter_label = "wrong_format")    
      }
    }
    
    
    
    
    # FILTER BY VAF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This is the first filter, We want to have it on top before the remove.nonexonic
    # because we want the VAF filter to be applied to either exonic or not
    if (vaf.cutoff != 0 && nrow(vcf)!=0){
        vcf_backup <- vcf  #make a copy
        vcf <- filterByVAF(vcf, vaf.cutoff)
        if (debug){
          vcf_filtered_df <- .debug_filter(vcf_before_filter = vcf_backup
                                           , vcf_after_filter = vcf
                                           , vcf_filtered = vcf_filtered_df
                                           , filter_label = "vaf.cutoff")
        }
    }

    
    # REMOVE OFF_TARGET 
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Let's remove all the variants mapped outside of the design
    vcf_backup <- vcf  #make a copy
    vcf <- IRanges::subsetByOverlaps(vcf, design, type="within")
    if (debug){
      vcf_filtered_df <- .debug_filter(vcf_before_filter = vcf_backup
                                       , vcf_after_filter = vcf
                                       , vcf_filtered = vcf_filtered_df
                                       , filter_label = "off_target")
    }
    
    # make a copy of the vcf
    vcf_original <- vcf
    vcf_backup <- vcf  
    
    # REMOVE SNV MAPPED in NON CODING regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(methods::is(vcf)[1] == "CollapsedVCF" && nrow(vcf)!=0  ){
      vcf <- .remove_noncoding_SNV(vcf, assembly)  # remove non coding SNV
      if (debug){
        vcf_filtered_df <- .debug_filter(vcf_before_filter = vcf_backup
                                         , vcf_after_filter = vcf
                                         , vcf_filtered = vcf_filtered_df
                                         , filter_label = "non_coding")
      }
    }
    
    # RETRIEVE NON CODING SNV
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Find the variants that we just removed because they were mapped in non coding
    vcf_non_exonic_only_df = NULL
    if(remove.nonexonic == FALSE){
      vcf_non_exonic_only_df = .fetch_filtered_SNV(
          original_vcf = vcf_original
          , filtered_vcf = vcf)  
      
      if (debug){
        # remove the non_coding just filtered out
        vcf_filtered_df <- vcf_filtered_df %>% dplyr::filter(filter != "non_coding")
      }
    }
    
    # FILTER BY CANCER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (remove.cancer == TRUE && nrow(vcf)!=0){
      
      if (is.null(tsList)) {
          # load imported dataset
          utils::data(HumanTumorSuppressors)
          tsList <- HumanTumorSuppressors
      } else{
          tsList <- tsList
      }
      if((!(methods::is(tsList)[1]) == "character")|(!(methods::is(tsList)[2]) == "vector")){
          stop("No valid tsList: please provide a character vector with entrez_gene ids")
      }
      # This also removes all the non-exonic variants
      vcf_backup <- vcf  
      vcf <- removeCancerVariants(vcf, assembly, tsList)
      if (debug){
        vcf_filtered_df <- .debug_filter(vcf_before_filter = vcf_backup
                                         , vcf_after_filter = vcf
                                         , vcf_filtered = vcf_filtered_df
                                         , filter_label = "remove.cancer")
      }
    }

    # FILTER BY VARIANT TYPE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This no longer removes the SNVs mapped in non-coding regions
    # This returns a dataframe
    if (!(is.null(variantType)) && nrow(vcf)!=0) {
      vcf_backup <- vcf  
      vcf <- filterVariantType(vcf, assembly, variantType)
      if (debug){
        vcf_filtered_df <- .debug_filter(vcf_before_filter = vcf_backup
                                         , vcf_after_filter = vcf
                                         , vcf_filtered = vcf_filtered_df
                                         , filter_label = "variant.filter")
      }
      
    }

    # Prepare output VCF FILE - convert it into a dataframe
    # ------------------------------------------------------------------------
    if (nrow(vcf) == 0){
      # if the vcf is empty make sure it's a dataframe of empty shape
      vcf <- .empty_vcf()
    } else {
      # sanity check
      if(class(vcf)[1] == "data.frame"){stop("I was expecting a collapsedVCF object")}
      # convert vcf object into a dataframe
      vcf <- as.data.frame(DelayedArray::rowRanges(vcf)) %>%
        dplyr::rename("CHR" = .data$seqnames ) %>%
        dplyr::rename("START" = .data$start ) %>%
        dplyr::rename("END" = .data$end ) %>%
        dplyr::select(.data$CHR, .data$START, .data$END, .data$REF, .data$ALT)  %>%
        # ALT is returned as a biostring object, let's convert it into a simple chr instead
        dplyr::mutate("ALT" = purrr::map_chr(1:length(.data$ALT), ~ paste(.data$ALT[[.]], collapse="")))
      
    }
    
    # Merge back the non-exonic SNVs
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (remove.nonexonic == FALSE && !is.null(vcf_non_exonic_only_df)){
      vcf <- rbind(vcf_non_exonic_only_df, vcf) 
    }
    
    # Generate character string describing applied filters
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vaf=NULL
    cancer=NULL
    variant=NULL
    if(!(is.null(vaf.cutoff))){
        vaf=paste0("vaf=", vaf.cutoff)
    }else{vaf=""}
    if(remove.cancer==T){
        cancer="NoCancerMuts"
    }else{cancer=""}
    if(!(is.null(variantType))){
        variant=paste(unlist(lapply(variantType, function(x) paste0("No", x))), collapse="")
    }
    filtersDesc=paste0(vaf, cancer, variant)

    # Generate output list
    vcf_filtered <- list(variants=vcf, filter=filtersDesc, design=design)
    
    if (debug){
      print("Returning variable 'debug_env$filters' for debugging")
      debug_env$filters <- vcf_filtered_df
    }
    
    return(vcf_filtered)
}

# Helper function - return empty vcf
.empty_vcf <- function(){
  vcf <- as.data.frame(unique(cbind(as.character(), 
                                     as.numeric(), 
                                     as.numeric(), 
                                     as.character(), 
                                     as.character()))) 
  colnames(vcf) <- c("CHR", "START", "END", "REF", "ALT")
  return(vcf)
}
