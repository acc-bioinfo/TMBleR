#' Annotate TMB dataframe with clinical outcome
#'
#' Add a new column to the dataframe with the clinical outcome
#'
#' @param TMB_df output from applyTMB()
#' @param ClinicalData Simple dataframe with at least two columns, a
#' Sample column and a ClinicalResponder colum
#'
#' @return annotate dataframe which combines both TMB and clinical responders 
#' @export
#'
#' @examples
#' # Import datasets
#' data(TMB_VanAllen)
#' data("VanAllen_Clinical")
#' # Annotate TMB with clinical reponse
#' TMB_clincal_response <- annotateTMB(TMB_df = TMB_VanAllen, ClinicalData = VanAllen_Clinical )
annotateTMB <- function(TMB_df, ClinicalData){
    
    # sanity check
    if(!methods::is(ClinicalData, "data.frame")) stop("ClinicalData should be a data.frame")
    if(!any(colnames(TMB_df) %in% "Sample")) stop("Sample should be a column in the TMB_df data.frame")
    if(!any(colnames(ClinicalData) %in% "Sample")) stop("Sample should be a column in the ClinicalData data.frame")
    
    # merge clinical data with TMB data into a single dataframe performing an
    # inner join
    out <- merge(ClinicalData, TMB_df, by="Sample", all=F)
    
    return(out)
}