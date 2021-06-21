context("test_annotateTMB.R")

#' # Import datasets
data("TMB_VanAllen")
data("VanAllen_Clinical")

# Annotate TMB dataset with clinical data
TMB_clinical_response <- annotateTMB(TMB_df = TMB_VanAllen, ClinicalData = VanAllen_Clinical )


# test_that("Positive Test", {
# 
#     p <- generateBoxplot(TMB_clincal_response)
#     
# })

test_that("Negative Test: empy input", {
    
    # This shoul
    tmp <- TMB_clinical_response %>%
        dplyr::filter(Design == "WES") %>%
        dplyr::filter(Filter == "NoCancerMuts") 
    
    expect_error(generateBoxplot(tmp), "TMB_clinical should have at least one row")
    
})



# This should work and generate only two facests
# TMB_clinical_response %>%
#     dplyr::filter(Filter == "NoCancerMuts") %>%
#     generateBoxplot
# 
# # This should also work and only plot one class
# TMB_clinical_response %>%
#     dplyr::filter(ClinicalResponse == "responder") %>%
#     generateBoxplot
# 
# # This only shows one facet
# TMB_clinical_response %>%
#     dplyr::filter(Panel == "WES") %>%
#     generateBoxplot
# 
