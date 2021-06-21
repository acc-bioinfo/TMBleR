#' Plots Tumor Mutational Burden (TMB)
#'
#' This function reads in input the output of applyTMB() function and plots it
#'
#' @param TMB_results a \code{data.frame} with sample, filter, sequencing size,
#' total number of mutations and number of mutations per megabase.
#' @param type a string describing the type of plot to generate: 'barplot', 
#' 'densityplot' or boxplot.
#' @return Returns a pot representing TMB per Mb given in input for of all 
#' samples and with every applied filter separately
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
#' ## Prepare data for TMB quantification, formatting it as required by the 
#' ## applyTMB() function. While this is done automatically within the 
#' ## applyFilters() function, on a non filtered vcf we need to use the function
#' ## here described to get the correct format
#' vcfs_nonfiltered <- applyFilters(vcfs = vcfs, design = design, assembly= "hg19)
#'                                     
#' ## Perform TMB quantification
#' TMB_res=applyTMB(inputForTMB = vcfs_nonfiltered, assembly = "hg19")
#' 
#' ## Plot the TMB values calculated by the applyTMB() function
#' plotTMB(TMB_results = TMB_res, type="barplot")
#'
#'
#' @author Laura Fancello
#'
#' @export
plotTMB <- function(TMB_results, type="barplot"){

    # Sanity Checks  -----------------------------------------------------------
    if((!(type=="barplot"))&(!(type=="densityplot"))&(!(type=="boxplot"))){
	    stop("no valid argument \"type\": please indicate 'barplot' or 'densityplot'")
    }	
    
    # Generate barplot  of TMB values  -----------------------------------------
    if(type=="barplot"){
        p <- ggplot2::ggplot(data=TMB_results, 
            ggplot2::aes(x=as.character(as.vector(.data$Sample)), 
                                        y=as.numeric(as.vector(.data$TMB_per_Mb)), 
                                        fill=as.factor(.data$Filter))) +
            ggplot2::xlab("Sample") +
            ggplot2::ylab("TMB per Mb") +
            ggplot2::geom_bar(stat="identity", color="black", position=ggplot2::position_dodge()) + 
            #ggplot2::theme_bw() + 
            #ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
            #panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
            ggplot2::scale_fill_brewer(palette="Set1", aesthetics="fill", name="Filter") +
            ggplot2::facet_grid(. ~ .data$Design) +
            ggplot2::theme(legend.position="bottom")
    }

    # Generate density plot  of TMB values  ------------------------------------
    if(type=="densityplot"){
        p <- TMB_results %>%
            ggplot2::ggplot(ggplot2::aes(x=.data$TMB_per_Mb, fill=.data$Filter)) + 
            ggplot2::geom_density() +
            ggplot2::facet_grid(.data$Filter ~ .data$Design) +
            ggplot2::theme(
                strip.text.x = ggplot2::element_text(size = 8 ) #, colour="red")
               , strip.text.y = ggplot2::element_text(size = 5 ) # , colour="blue")
               , legend.position="bottom"
               #, legend.title = "Filter"
               #, text = element_text(size=7, colour="green")
               ) +
            ggplot2::scale_fill_brewer(palette="Set1", aesthetics="fill", name="Filter") +
            ggplot2::xlab("TMB per Mb") 
    }
    
    # Generate Boxplot of TMB values  ------------------------------------
    if(type=="boxplot"){
        p <- TMB_results %>% 
            ggplot2::ggplot() + 
            ggplot2::geom_boxplot(ggplot2::aes(x = .data$Filter, y = .data$TMB_per_Mb, fill=.data$Filter)) + 
            ggplot2::facet_grid(.data$Design ~ . , space = "free", scales="free")  +
            ggplot2::scale_fill_brewer(palette="Set1", aesthetics="fill", name="Filter") +
            ggplot2::xlab("Filters") +
            ggplot2::ylab("TMB per Mb") +
            ggplot2::coord_flip() +
            ggplot2::theme(legend.position="bottom")
    }
    
return(p)

}
