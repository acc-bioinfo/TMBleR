context("applySimulatePanel.R")

# Tests
################################################################################
# Read input for the applySimulatePanel function to test
data(vcfs_NoCancer_ForPanel)
WES.design <- readDesign(system.file("extdata"
                                     , "ExampleWESdesign_OnlyChr7ForTest.bed"
                                     , package = "TMBleR"
                                     , mustWork = TRUE)
                         , assembly = "hg19")
panel.design <- readDesign(system.file("extdata",
                                        "ExamplePanel_design_OnlyChr5and7ForTest.bed",
                                       package = "TMBleR",
                                       mustWork = TRUE),
                           assembly = "hg19")

# TEST INPUT FORMAT
# ------------------------------------------------------------------------------
test_that("test input: WES", {
    # missing entry
    expect_error(
        applySimulatePanel(WES = banana_WES, WES.design = WES.design, panel.design = panel.design, assembly = "hg19")
        , "object 'banana_WES' not found", fixed=TRUE
    )
    # Object does not exists
    expect_error(
        applySimulatePanel(WES = banana_WES, WES.design = WES.design, panel.design = panel.design, assembly = "hg19")
        , "object 'banana_WES' not found", fixed=TRUE
    )
})
test_that("test input: WES.design", {
    # missing entry
    expect_error(
        applySimulatePanel(WES = vcfs_NoCancer_ForPanel, panel.design = panel.design, assembly = "hg19")
        , "argument \"WES.design\" is missing, with no default"
    )
    # Object does not exists
    expect_error(
        applySimulatePanel(WES = vcfs_NoCancer_ForPanel, WES.design = banana_WESdesign, panel.design = panel.design, assembly = "hg19")
        , "object 'banana_WESdesign' not found", fixed=TRUE
    )
})
test_that("test input: panel.design", {
    # missing entry
    expect_error(
        applySimulatePanel(WES = vcfs_NoCancer_ForPanel, WES.design = WES.design, assembly = "hg19")
        , "argument \"panel.design\" is missing, with no default"
    )
    # Object does not exists
    expect_error(
        applySimulatePanel(WES = vcfs_NoCancer_ForPanel, WES.design = WES.design, panel.design = banana_paneldesign, assembly = "hg19")
        , "object 'banana_paneldesign' not found", fixed=TRUE
    )
})

test_that("test input: assembly", {
    # missing entry
    expect_error(
        applySimulatePanel(WES = vcfs_NoCancer_ForPanel,WES.design = WES.design, panel.design = panel.design)
        , "argument \"assembly\" is missing, with no default"
    )
    # wrong entry
    expect_error(
        applySimulatePanel(WES = vcfs_NoCancer_ForPanel, WES.design = WES.design, panel.design = panel.design, assembly = "banana")
        , "No valid genome specified: please specify 'hg19' or 'hg38'"
    )
})

# TEST OUTPUT OF FUNCTION  
# ----------------------------------------------------------------

# CODE TO GENERATE OUTPUT #
# Read input for the applySimulatePanel function to test
data("vcfs_NoCancer_ForPanel")

# Subset the WES dataset so that it will only contain variants in the regions 
# targeted by the panel you want to simulate
SimulatedPanel_NoCancer <- applySimulatePanel(vcfs_NoCancer_ForPanel, 
                                              WES.design=WES.design,
                                              panel.design=panel.design, 
                                              assembly = "hg19")

# TESTS ON OUTPUT #
## Test class of output objects from functions ##
test_that("Class of function output corresponds to the expected one", {
    expect_is(SimulatedPanel_NoCancer[[1]], "list")
    expect_is(SimulatedPanel_NoCancer[[1]]$variants, "GRanges")
    expect_is(SimulatedPanel_NoCancer[[1]]$filter, "character")
    expect_is(SimulatedPanel_NoCancer[[1]]$design, "GRanges")
    expect_is(SimulatedPanel_NoCancer[[1]]$sample, "character")
})    

## Test dimensions of output objects from functions ##
test_that("Number of variants in output by applySimulatePanel function is not higher than number of variants in input", {
    if(is.null(dim(SimulatedPanel_NoCancer[[1]]$variants)[1])){
        # When no overlap possible between WES and panel design assign 
        # simulated panel size 0
        dim_SimulatedPanel <-0 
    }else{
            dim_SimulatedPanel <- dim(SimulatedPanel_NoCancer[[1]]$variants)[1] 
        }
    expect_true(dim(vcfs_NoCancer_ForPanel[[1]]$variants)[1] >= dim_SimulatedPanel)
    })
test_that("Dimensions of function output correspond to the expected ones", {
    expect_equal(length(SimulatedPanel_NoCancer[[1]]), 4)
    expect_equal(length(SimulatedPanel_NoCancer), 4)
})
    

    