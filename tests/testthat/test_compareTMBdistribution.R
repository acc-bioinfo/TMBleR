context("compareTMBdistribution.R")

# Tests
################################################################################
# Read in input TMB values
data(Hellman_SimulatedFM1Panel_WES)

# TEST INPUT FORMAT
# ------------------------------------------------------------------------------
test_that("test input: dataset", {
    # missing entry
    expect_error(
        compareTMBdistribution(TMB = "Panel.NumMuts")
        , "argument \"dataset\" is missing, with no default"
    )
    # Object does not exists
    expect_error(
        compareTMBdistribution(dataset = banana, TMB = "Panel.NumMuts")
        , "object 'banana' not found", fixed=TRUE
    )
})
test_that("test input: TMB", {
    # missing entry
    expect_error(
        compareTMBdistribution(dataset = Hellman_SimulatedFM1Panel_WES)
        , "argument \"TMB\" is missing, with no default"
    )
    # Object does not exists
    expect_error(
        compareTMBdistribution(dataset = banana, TMB = "Panel.NumMuts")
        , "object 'banana' not found", fixed=TRUE
    )
})

# TEST OUTPUT OF FUNCTION  
# ----------------------------------------------------------------

# CODE TO GENERATE OUTPUT #
test_res <- compareTMBdistribution(Hellman_SimulatedFM1Panel_WES, TMB = "Panel.NumMuts")

# TESTS ON OUTPUT #
## Test class of output objects from functions ##
test_that("Class of function output corresponds to the expected one", {
    expect_is(test_res, "htest")
})


    