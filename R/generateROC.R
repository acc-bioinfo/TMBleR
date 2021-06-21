#' Generate ROC curves to evaluate TMB performance in immunotherapy responder
#' classification
#'
#' This function takes in input a data frame with a vector of TMB values and a
#' factor named ClinicalResponse with only these two levels allowed: 'responder'
#' and 'nonresponder'. It gives in output a ROC curve describing TMB performance
#' in responders and nonresponders classification on a training set (random 75%
#' of input data), the best value of TMB cutoff for classification in this 
#' training dataset, a confusion matrix describing the classification performance
#' of this TMB cutoff in a training dataset (25% of input data) and a ROC curve 
#' describing TMB performance in responders and nonresponders classification in 
#' a training set.
#'
#' @param dataset a \code{data.frame} containing a numeric vector of TMB values 
#' and a factor with name 'ClinicalResponse'  and levels 'responder' and 
#' 'nonresponder' and a column with TMB values named "TMB_per_Mb".
#' @param method \code{character string} indicating the method to identify the
#' best TMB cutoff: i.e. Youden, MaxSe, MaxSpSe. For more details see also
#' \code{\link[OptimalCutpoints]{optimal.cutpoints}}
#'
#' @return Returns in output the best TMB cutoff for classification,
#' AUC with 95% C.I. specificity, sensitivity and positive and negative predictive
#' values based on the model built on training data. It also returns the
#' confusion matrix of classification on test data based on the identified TMB
#' cutoff and a the ROC curve on train and test data.
#'
#' @author Laura Fancello
#'
#' @export
generateROC <- function(dataset, method="Youden"){

    # Sanity Checks  -----------------------------------------------------------
    ## Check input arguments
    if (is.null(dataset)) {
        stop("argument 'dataset' is missing, with no default")
    }
    if (!(methods::is(dataset)[1] == "data.frame")) {
        stop("No valid dataset: please provide a data.frame object.")
    }

    # Perform stratified k-fold crossvalidation, governed by criteria such as 
    # ensuring that each fold has the same proportion of observations with a given
    # categorical value which account for classes when 
    # creating folds, to avoid class imbalance
    # folds <- createFolds(factor(data$target), k = 10, list = FALSE)
    # 
    # caret_model <- train(method "glm", trainControl (method = "cv", number = 5, …), data = train, …)
    # 
    # train_control <- trainControl(method="repeatedcv", repeats=5)
    # model <- train(Species~., data=df, trControl=train_control, method="nb")
    # 
    # 
    ## Divide data in training and testing
    inTrain <- caret::createDataPartition(y = as.factor(dataset$ClinicalResponse),
                                          p = .75, # Percent of data in training
                                          list = FALSE)
    training <- dataset[inTrain, ]
    testing  <- dataset[-inTrain, ]

    # Performance analysis on training set  ------------------------------------
    # Identify best TMB cutoff and calculate AUC, specificity, sensitivity,
    # positive and negative predictive values on the training set.
    train_res <- OptimalCutpoints::optimal.cutpoints(X = "TMB_per_Mb",
                                        status = "ClinicalResponse",
                                        tag.healthy = "nonresponder",
                                        methods = method,
                                        data = training,
                                        pop.prev = NULL,
                                        control = OptimalCutpoints::control.cutpoints(),
                                        ci.fit = TRUE,
                                        conf.level = 0.95, 
                                        trace = FALSE)
    
    # Validation on test set  --------------------------------------------------
    ## Use the fitted model to make predictions on new data
    
    # Use best cutoff from model fitted on training to make predictions on test
    cutoff <- train_res[[1]]$Global$optimal.cutoff$cutoff
    testing$pred <- factor( ifelse(testing$TMB_per_Mb >= cutoff,
                                   "responder", "nonresponder") )
    testing$pred = factor(testing$pred, levels=c("responder", "nonresponder"))
    testing$ClinicalResponse = factor(testing$ClinicalResponse, levels=c("responder", "nonresponder"))

    ## Generate confusion matrix using as TMB cutoff the best cutoff determined
    ## on training set
    confMatrix_test <- caret::confusionMatrix(testing$pred, testing$ClinicalResponse)

    ## Plot ROC on test data
    test_res <- OptimalCutpoints::optimal.cutpoints(X = "TMB_per_Mb",
                                                            status = "ClinicalResponse",
                                                            tag.healthy = "nonresponder",
                                                            methods = method,
                                                            data = testing,
                                                            pop.prev = NULL,
                                                            control = OptimalCutpoints::control.cutpoints(),
                                                            ci.fit = TRUE,
                                                            conf.level = 0.95,
                                                            trace = FALSE)
    
    return(list(train_res=train_res, confMatrix_test=confMatrix_test, test_res=test_res))
}
