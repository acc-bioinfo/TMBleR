context("test_annotateTMB.R")

#' # Import datasets
data("TMB_VanAllen")
data("VanAllen_Clinical")


test_that("Positive Test", {
    # Annotate TMB with clinical reponse
    TMB_clinical_response <- annotateTMB(TMB_df = TMB_VanAllen, ClinicalData = VanAllen_Clinical )
    
    expect_equal(nrow(TMB_clinical_response), 20)
    expect_equal(sum(TMB_clinical_response$Tot_Number_Mutations), 1243)
    expect_equal(table(TMB_clinical_response$ClinicalResponse)[[1]], 10)
})

test_that("Negative Test: wrong input format", {
    # Create errors in the expected columns names 
    TMB_wrong <- TMB_VanAllen %>% dplyr::rename(wrong_name = Sample)
    Clinical_Data_wrong <- VanAllen_Clinical %>% dplyr::rename(wrong_name = Sample)
    
    # Expect column "Sample" in TMB_res
    expect_error(annotateTMB(TMB_df = TMB_wrong
                            , ClinicalData = VanAllen_Clinical )
                 , "Sample should be a column in the TMB_df data.frame")
    
    # Expect column "Sample" in ClinicalData
    expect_error(annotateTMB(TMB_df = TMB_VanAllen
                             , ClinicalData = Clinical_Data_wrong )
                 , "Sample should be a column in the ClinicalData data.frame")
    
})

