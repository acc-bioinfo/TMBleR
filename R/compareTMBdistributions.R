#' Compare TMB distribution in responders and nonresponders
#'
#' This function generates applies the unpaired Mann-Whitney test to compare the
#'  distribution of TMB values between immunotherapy responders and nonresponders.
#'
#' @param dataset a \code{data.frame} object with the following fields:
#' patient identifier, numeric TMB value and response to
#' immunotherapy, expressed as "responders" or "nonresponders"
#' @param TMB name of the \code{data.frame} column containing TMB values

#'
#' @return Returns the statistics from Wilcoxon test
#'
#' @examples
#' ## Compare TMB distribution between immunotherapy responders and nonresponders
#' ## using the unpaired Mann-Whitney test
#' 
#' # Read TMB values and response to immunotherapy
#' data(Hellman_SimulatedFM1Panel_WES)
#' 
#' # Compare TMB distribtions by Wilcoxon test
#' compareTMBdistribution(dataset = Hellman_SimulatedFM1Panel_WES, TMB = "WES.NumMuts")
#' 
#'
#' @author Laura Fancello
#'
#' @export
compareTMBdistribution <- function(dataset,
                             TMB ){
    
    # Sanity Checks  -----------------------------------------------------------
    if (is.null(dataset)) {
        stop("argument 'dataset' is missing, with no default")
    }
    if (is.null(TMB)) {
        stop("argument 'TMB' is missing, with no default")
    }
    
    # Preprocess input  --------------------------------------------------------
    ## Reads input file containing TMB and clinical response
    df <- data.frame(ClinicalResponse=as.factor(dataset$ClinicalResponse),
                     TMB=as.numeric(as.vector(dataset[[TMB]])))
    df <- stats::na.omit(df)

    # Compare distributions  ---------------------------------------------------
    ## Apply unpaired Mann-Whitney test for differences in TMB distribution between
    ## responders and nonresponders
    wilcox_res <- stats::wilcox.test(as.numeric(as.vector(df$TMB))~as.factor(df$ClinicalResponse),
                              data=df,
                              paired=FALSE,
                              alternative="two.sided")
    
    return(wilcox_res)
}