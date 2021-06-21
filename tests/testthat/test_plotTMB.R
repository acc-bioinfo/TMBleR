context("plotTMB.R")

# Tests
################################################################################

# Read input for the correlateTMBvalues function to test
data(TMB_Horizon)

# TEST INPUT FORMAT
# ------------------------------------------------------------------------------
test_that("test input: TMB_results", {
    # missing entry
    expect_error(
        plotTMB(type="barplot") 
        , "argument \"TMB_results\" is missing, with no default"
    )
    # Object does not exists
    expect_error(
        plotTMB(TMB_results = banana_results, type="barplot")
        , "object 'banana_results' not found", fixed=TRUE
    )
})

test_that("test input: type", {

	# Object does not exists
	expect_error(
		plotTMB(TMB_results = TMB_Horizon, type="bananaplot")
		, "no valid argument \"type\": please indicate 'barplot' or 'densityplot'" 
	)
})
# TEST OUTPUT OF FUNCTION  
# ----------------------------------------------------------------


# Plot results
p <- plotTMB(TMB_results=TMB_Horizon, type="barplot")

## Test class of output objects from functions ##
test_that("Class of function output corresponds to the expected one", {
    expect_is(p, "gg")
})
