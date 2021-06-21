#' Correlate panel-based TMB values with WES-based TMB values
#'
#' This function takes in input TMB values quantified for the same samples from
#' gene panels and from whole exome sequencing datasets and performs correlation
#' analysis. Panel-based TMB values may come from real gene panel sequencing or
#' from a simulated gene panel, generated using the \code{\link{simulatePanel}}
#' function.
#'
#' @param panel.TMB a \code{numeric vector} of length n containing TMB values
#' estimated from a gene panel for n samples
#' @param WES.TMB a \code{numeric vector} of length n containing TMB values
#' estimated from whole exome sequencing for n samples
#' @param corr.coeff  type of correlation coefficient: "pearson", "kendall" or
#' "spearman"
#' @param title.plot a \code{string} describing plot title 
#'
#' @return Returns a pdf file with correlation plot, the Spearman's correlation
#' coefficient with 95% confidence intervals calculated by bootstrapping and the
#' number of samples.
#'
#' @examples
#' # Read vcf files containing somatic mutations identified by WES -------------
#' data(ExampleWESvcfs)
#' 
#' # Read files with WES and gene panel design ---------------------------------
#' data(ExampleWESdesign)
#' data(ExamplePaneldesign)
#' 
#' # Filter for gene panel simulation ------------------------------------------
#' # Remove known cancer mutations
#' vcfs_NoCancer_ForPanel <- applyFilters(vcfs = ExampleWESvcfs, 
#'                                        assembly = "hg19", 
#'                                        design = ExamplePaneldesign, 
#'                                        remove.cancer = TRUE, 
#'                                        tsList = NULL, 
#'                                        variantType = NULL)
#' 
#' # Filter for original WES ---------------------------------------------------
#' # Filter out synonymous mutations
#' vcfs_NoSynonymous_WES <- applyFilters(vcfs = ExampleWESvcfs, 
#'                                       assembly = "hg19", 
#'                                       design = ExampleWESdesign, 
#'                                       remove.cancer = FALSE, 
#'                                       tsList = NULL, 
#'                                       variantType = c("synonymous"))
#' 
#' # Subset the WES dataset so that it will only contain variants in the regions
#' # targeted by the panel you want to simulate
#' SimulatedPanel_NoCancer <- applySimulatePanel(WES = vcfs_NoCancer_ForPanel, 
#'                                               WES.design = ExampleWESdesign,
#'                                               panel.design = ExamplePaneldesign, 
#'                                               assembly = "hg19")
#' 
#' # TMB quantification --------------------------------------------------------
#' # Perform TMB quantification on the simulated panel
#' TMBs_SimulatedPanel <- applyTMB(inputForTMB = SimulatedPanel_NoCancer, assembly = "hg19")
#' # Perform TMB quantification on the original Whole Exome sequencing
#' TMBs_WES <- applyTMB(inputForTMB = vcfs_NoSynonymous_WES, assembly = "hg19")
#' 
#' # Correlate WES-based and simulated panel-based TMB values ------------------
#' cor_plot <- correlateTMBvalues(panel.TMB = TMBs_SimulatedPanel$Tot_Number_Mutations, 
#'                                    WES.TMB = TMBs_WES$Tot_Number_Mutations, 
#'                                    corr.coeff = "spearman",
#'                                    title.plot="Correlation panel-based and WES-based TMB")
#' # Visualize plot of Spearman correlations
#' print(cor_plot)
#'
#' @author Laura Fancello
#'
#' @export
correlateTMBvalues <- function(panel.TMB, WES.TMB, corr.coeff, title.plot){

    # Sanity Checks  ----------------------------------------------------------
    ## Check the input arguments
    if (is.null(panel.TMB)) {
        stop("argument 'panel.TMB' is missing, with no default")
    }
    ## Set the arbitrary minimum number of samples to perform correlation
    ## analysis and get confidence intervals
    if(length(panel.TMB) < 4) {
        stop("Not enough data points for correlation analysis.")
    }
    if(!(methods::is(panel.TMB)[1] == "numeric")) {
        stop("Not valid panel-based TMBs: please provide numeric values.")
    }
    if (is.null(WES.TMB)) {
        stop("argument 'WES.TMB' is missing, with no default")
    }
    ## Set the arbitrary minimum number of samples to perform correlation
    ## analysis and get confidence intervals
    if (length(WES.TMB) < 4) {
        stop("Not enough data points for correlation analysis.")
    }
    if(!(methods::is(WES.TMB)[1] == "numeric")) {
        stop("Not valid WES-based TMBs: please provide numeric values.")
    }

        # Check if the same number of data points for WES and panel-based TMB values
    if(length(panel.TMB) !=  length(WES.TMB)) {
        stop("WES.TMB and panel.TMB have different length")
    }
    
    if (!(corr.coeff == 'spearman') && !(corr.coeff == 'pearson') && !(corr.coeff == 'kendall')) {
        stop("No valid corr.coeff specified: please indicate 'pearson', 'spearman' or 'kendall'")
    }    
    
    # if (is.null(title.plot)) {
    #     stop("No valid title.plot specified: please indicate a titel for the correlation plot")
    # }    
    
    
    # Correlate  ----------------------------------------------------------
    ## Calculate the correlation coefficient
    corr.res <- stats::cor.test(panel.TMB, WES.TMB, method=corr.coeff) 
    
    ## Generate title for correlation plot containing Spearman's correlation
    ## coefficient, 95% confidence intervals and number of data points
    n <- paste0("n=", length(panel.TMB))
    rho <- paste0("corr coeff=",
                  round(as.numeric(as.vector(corr.res$estimate)),
                        digits=2))
    p <- paste0("p-value=",
                round(as.numeric(as.vector(corr.res$p.value)),
                      digits=2))
    
    data <- as.data.frame(cbind(WES.TMB, panel.TMB))
    
    regr <- stats::lm(panel.TMB ~ WES.TMB)
    
    R2 <- summary(regr)$adj.r.squared
    F1 <- summary(regr)$fstatistic
    
        # Plot correlation  ---------------------------------------------------------
    a <- ggplot2::ggplot(data = data, ggplot2::aes(x = WES.TMB, y = panel.TMB)) +
        ggplot2::xlab("WES-based TMB") +
        ggplot2::ylab("Panel-based TMB") +
        ggplot2::labs(title = paste0(title.plot, " ", rho, ", ", p, ", ", n),
                     subtitle = paste0("Linear model: R2=", R2, " F1=", F1)) +
        ggplot2::theme_bw() + 
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
        ggplot2::geom_abline(intercept = regr$coefficients[1], slope = regr$coefficients[2], colour = "blue") +
        ggplot2::geom_point(size = 1)
    
    return(a)
    
}
