#' Generate and compare two ROC curves of TMB performance in immunotherapy responder
#' classification
#'
#' This function takes in input two sets of TMB values for the same samples and
#' a binary variable representing the clinical response to immunotherapy. It
#' generates two ROC curves showing performance of the two TMB datasets to
#' discriminate between responders and nonresponders, as well as the best TMB
#' cutoff.
#'
#' @param dataset a \code{data.frame} containing numeric vectors of TMB values 
#' and a factor with name 'ClinicalResponse'  and levels 'responder' and 
#' 'nonresponder'
#' @param TMB1 \code{character string} indicating the \code{data.frame} column
#' containing TMB values from the first approach
#' @param TMB2 \code{character string} indicating the \code{data.frame} column
#' containing TMB values from the second approach
#' @param method \code{character string} indicating the method to use: delong,
#' bootstrap or venkatraman. See \code{\link[pROC]{roc.test}} for more details.
#' @param paired a logical indicating whether you want a paired roc.test. If
#' NULL, the paired status will be auto-detected by are.paired. If TRUE but the
#' paired status cannot be assessed by are.paired will produce an error.
#' See \code{\link[pROC]{roc.test}} for more details.
#' @param boot.n for method="bootstrap" and method="venkatraman" only: the
#' number of bootstrap replicates or permutations. Default: 2000.
#' See \code{\link[pROC]{roc.test}} for more details.
#'
#' @return Returns in standard output the result of the test used to compare the
#' two ROC curves and a pdf file containg plots of the ROC curves and their
#' respective AUC, sensitivity and specificity values and the best
#' threshold for sample classification.
#'
#' @examples
#' ## Compare the ROC curves representing the performance of panel-based and
#' ## WES-based TMB to predict clinical response with the default delong method.
#' 
#' # Set the seed to create reproducible train and test sets in generateROC
#' set.seed(949)
#' 
#' 
#' # Read TMB values and response to immunotherapy
#' data(Hellman_SimulatedFM1Panel_WES)
#' 
#' # Compare ROCs
#' compareROC(dataset = Hellman_SimulatedFM1Panel_WES,
#'            TMB1 = "Panel.NumMuts",
#'            TMB2 = "WES.NumMuts",
#'            method = "d",
#'            paired = TRUE)
#'
#' @author Laura Fancello
#'
#' @seealso \code{\link[pROC]{roc.test}}
#'
#' @export
compareROC <- function(dataset,
                       TMB1,
                       TMB2,
                       method,
                       paired,
                       boot.n){
    
    # Sanity Checks  -----------------------------------------------------------
    if (is.null(dataset)) {
        stop("argument 'dataset' is missing, with no default")
    }
    if (is.null(TMB1)) {
        stop("argument 'TMB1' is missing, with no default")
    }
    if (is.null(TMB2)) {
        stop("argument 'TMB2' is missing, with no default")
    }
    if (is.null(method)) {
        stop("argument 'method' is missing, with no default")
    }
    if (is.null(paired)) {
        stop("argument 'paired' is missing, with no default")
    }
    if ((paired != "TRUE") && (paired != "FALSE")) {
        stop("No valid \"paired\" argument: please specify 'TRUE' or 'FALSE'")}
    
    
    # Preprocess input  --------------------------------------------------------
    if(is.null(dataset[[TMB1]])){
        stop("Object \"TMB1\" not found: please specify a vaild column of the input dataset")
    }
    if(is.null(dataset[[TMB2]])){
        stop("Object \"TMB2\" not found: please specify a vaild column of the input dataset")
    }
    
    
    nameTMB1 <- TMB1
    nameTMB2 <- TMB2
    TMB1 <- dataset[[nameTMB1]]
    TMB2 <- dataset[[nameTMB2]]
    
    # Compare ROC curves -------------------------------------------------------
    # Apply test to compare ROCs 
    test.result <- pROC::roc.test(as.factor(dataset$ClinicalResponse),
                                dataset=dataset,
                                TMB1,
                                TMB2,
                                method=method,
                                paired=paired,
                                boot.n=boot.n,
                                alternative="two.sided")
    print(test.result, quote=F)

    # Build ROC curves
    roc1 <- pROC::roc(as.factor(dataset$ClinicalResponse), 
                   TMB1,
                   levels = c("nonresponder", "responder"),
                   ci=TRUE,
                   boot.n=100,
                   ci.alpha=0.95,
                   print.AUC=TRUE,
                   show.thresh=TRUE)
    roc2=pROC::roc(as.factor(dataset$ClinicalResponse),
                   TMB2,
                   levels = c("nonresponder", "responder"),
                   ci=TRUE,
                   boot.n=100,
                   ci.alpha=0.95,
                   print.AUC=TRUE,
                   show.thresh=TRUE
                   )

  
    # Plot ROC curves
    pROC::plot.roc(roc1,
                print.thres="best", # print best threshold according to Youden
                print.thres.pattern = "%.3f (Spec = %.2f, Sens = %.2f)",
                print.thres.cex = .8,
                print.thres.col="orange",
                print.auc=TRUE,
                print.auc.adj=c(-1,12),
                legacy.axes = TRUE,
                col="orange"
    )
    pROC::plot.roc(roc2,
                print.thres="best", # print best threshold according to Youden
                print.thres.pattern = "%.3f (Spec = %.2f, Sens = %.2f)",
                print.thres.cex = .8,
                print.thres.col="blue",
                print.auc=TRUE,
                print.auc.adj=c(-1,14),
                legacy.axes = TRUE,
                add=TRUE,
                col="blue"
    )
    graphics::legend("bottomright", legend=c(nameTMB1, nameTMB2),
           col=c("orange", "blue"), lwd=2)
}
