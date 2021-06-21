context("readVcfFiles.R")

# Tests
################################################################################


# Read in input name of vcf files 
# (provide them as a list even if you only have one file)
vcf_files <- list(Horizon5_panel="Horizon5_ExamplePanel_OnlyChr5and7ForTest.vcf", 
                  HorizonFFPEmild="HorizonFFPEmild_ExamplePanel_OnlyChr5and7ForTest.vcf")

# For each vcf file, get the absolute path
vcf_files <- lapply(vcf_files,
                    function(x) system.file("extdata", x, package = "TMBleR", mustWork = TRUE))



# TEST INPUT FORMAT
# ------------------------------------------------------------------------------
test_that("test input: assembly", {
    # missing entry
    expect_error(
        readVcfFiles(vcf_files)
        , "argument \"assembly\" is missing, with no default"
    )
    # wrong entry
    expect_error(
        readVcfFiles(vcf_files, assembly = "banana")
        , "No valid genome specified: please specify 'hg19' or 'hg38'"
    )
})

test_that("test input: vcflist", {
    # missing entry
    expect_error(
        readVcfFiles(assembly="hg19")
        , "argument \"vcfFiles\" is missing, with no default"
    )
    # wrong entry
    expect_error(
        readVcfFiles(vcfFiles = "Horizon5_ExamplePanel_OnlyChr5and7ForTest.vcf", assembly = "hg19")
        , "wrong input type: vcfFiles argument needs to be a list"
    )
    # file does not exists
    expect_error(
        readVcfFiles(vcfFiles = list(apple=system.file("extdata", 
                                                       "banana.vcf",
                                                       package = "TMBleR",
                                                       mustWork = FALSE)
                                     , assembly = "hg38")
                     , "file(s) do not exist", fixed=TRUE)
    )
    
    # Test unnamed list
    unnamed_vcf_files <- list("Horizon5_ExamplePanel_OnlyChr5and7ForTest.vcf"
                              , "HorizonFFPEmild_ExamplePanel_OnlyChr5and7ForTest.vcf")
    unnamed_vcf_files <- lapply(unnamed_vcf_files,
                                function(x) system.file("extdata", x, package = "TMBleR", mustWork = TRUE))
    expect_warning(
        readVcfFiles(vcfFiles = unnamed_vcf_files, assembly = "hg38")
        , "The list of vcf files was not named. Vcf files names used as 'sample' names by default", fixed=TRUE
    )
})

# POSTIVE TEST  ----------------------------------------------------------------
test_that("Positive test", {
    
    vcfs <- readVcfFiles(vcf_files, assembly = "hg19")    
    
    expect_equal(length(vcfs), 2) 
    expect_equal(dim(vcfs$Horizon5), c(46,2))
})

# TEST OUTPUT OF FUNCTION  
# ----------------------------------------------------------------
## Test class of output objects from functions ##
vcfs <- readVcfFiles(vcf_files, assembly = "hg19")    

test_that("Class of function output corresponds to the expected one", {
    # Test readVcfFiles function
    expect_is(vcfs, "list")
    expect_is(vcfs[[1]], "CollapsedVCF")
})
## Test dimensions of output objects from functions ##
test_that("Dimensions of function output correspond to the expected ones", {
    # Number of variants read in input by readVcfFiles
    expect_equal(length(S4Vectors::expand(vcfs[[1]])@rowRanges), 46)
})

