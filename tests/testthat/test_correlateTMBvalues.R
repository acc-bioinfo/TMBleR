context("correlateTMBvalues.R")

# Tests
################################################################################
# Read input for the correlateTMBvalues function to test
data(TMBs_SimulatedPanel)
data(TMBs_WES)

# TEST INPUT FORMAT
# ------------------------------------------------------------------------------
test_that("test input: panel.TMB", {
    # missing entry
    expect_error(
        correlateTMBvalues(WES.TMB = TMBs_WES$Tot_Number_Mutations, corr.coeff = "spearman", title.plot="Test plot")
        , "argument \"panel.TMB\" is missing, with no default"
    )
    # Object does not exists
    expect_error(
        correlateTMBvalues(panel.TMB = banana_TMB, WES.TMB = TMBs_WES$Tot_Number_Mutations, corr.coeff = "spearman", title.plot="Test plot")
        , "object 'banana_TMB' not found", fixed=TRUE
    )
})
test_that("test input: WES.TMB", {
    # missing entry
    expect_error(
        correlateTMBvalues(panel.TMB = TMBs_SimulatedPanel$Tot_Number_Mutations, corr.coeff = "spearman", title.plot="Test plot")
        , "argument \"WES.TMB\" is missing, with no default"
    )
    # Object does not exists
    expect_error(
        correlateTMBvalues(panel.TMB = TMBs_SimulatedPanel$Tot_Number_Mutations, WES.TMB = banana_TMB, corr.coeff = "spearman", title.plot="Test plot")
        , "object 'banana_TMB' not found", fixed=TRUE
    )
})
test_that("test input: corr.coeff", {
    # missing entry
    expect_error(
        correlateTMBvalues(panel.TMB = TMBs_SimulatedPanel$Tot_Number_Mutations, WES.TMB = TMBs_WES$Tot_Number_Mutations, title.plot="Test plot")
        , "argument \"corr.coeff\" is missing, with no default"
    )
    # wrong entry
    expect_error(
        correlateTMBvalues(panel.TMB = TMBs_SimulatedPanel$Tot_Number_Mutations, WES.TMB = TMBs_WES$Tot_Number_Mutations, corr.coeff = "banana", title.plot="Test plot")
        , "No valid corr.coeff specified: please indicate 'pearson', 'spearman' or 'kendall'"
    )
})
test_that("test input: title.plot", {
    # missing entry
    expect_error(
        correlateTMBvalues(panel.TMB = TMBs_SimulatedPanel$Tot_Number_Mutations, WES.TMB = TMBs_WES$Tot_Number_Mutations, corr.coeff="spearman")
        , "argument \"title.plot\" is missing, with no default"
    )
})


# TEST OUTPUT OF FUNCTION  
# ----------------------------------------------------------------

# CODE TO GENERATE OUTPUT #
# Calculate Spearman correlation
corr_res <- correlateTMBvalues(TMBs_SimulatedPanel$Tot_Number_Mutations, 
                                   TMBs_WES$Tot_Number_Mutations, 
                                   corr.coeff = "spearman",
                                    title.plot="Test plot")

# TESTS ON OUTPUT #
## Test class of output objects from functions ##
test_that("Class of function output corresponds to the expected one", {
    expect_is(corr_res, "gg")  
})


