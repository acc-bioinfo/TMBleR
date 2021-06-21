context("applyTMB.R")

# Tests
################################################################################
# Read input for the applyTMB function to test
data(vcfs_all)

# TEST INPUT FORMAT
# ------------------------------------------------------------------------------
test_that("test input: inputForTMB", {
    # missing entry
    expect_error(
        applyTMB(assembly = "hg19")
        , "argument \"inputForTMB\" is missing, with no default"
    )
    # Object does not exists
    expect_error(
        applyTMB(inputForTMB = banana_vcf, assembly = "hg19")
        , "object 'banana_vcf' not found", fixed=TRUE
    )
})
test_that("test input: assembly", {
    # missing entry
    expect_error(
        applyTMB(inputForTMB = vcfs_all)
        , "argument \"assembly\" is missing, with no default"
    )
    # wrong entry
    expect_error(
        applyTMB(inputForTMB = vcfs_all, assembly = "banana")
        , "No valid genome assembly: please indicate 'hg19' or 'hg38'"
    )
})

# TEST OUTPUT OF FUNCTION  
# ----------------------------------------------------------------

# CODE TO GENERATE OUTPUT #
# Perform TMB quantification
TMB_vcfs=applyTMB(inputForTMB = vcfs_all, assembly = "hg19")

# TESTS ON OUTPUT #
## Test class of output objects from functions ##
test_that("Class of function output corresponds to the expected one", {
    expect_is(TMB_vcfs, "data.frame")
    expect_is(TMB_vcfs$Sample, "factor")
    expect_is(TMB_vcfs$Filter, "factor")
    expect_is(TMB_vcfs$Sequencing_Size, "numeric")
    expect_is(TMB_vcfs$Tot_Number_Mutations, "numeric")
    expect_is(TMB_vcfs$TMB_per_Mb, "numeric")
})

test_that("Dimensions of function output correspond to the expected ones", {
    expect_equal(dim(TMB_vcfs)[1], 18)
    expect_equal(dim(TMB_vcfs)[2], 6)
})