context("readDesign.R")

# Tests
################################################################################

file <- system.file("extdata"
                    , "ExamplePanel_GeneIDs.txt"
                    , package = "TMBleR"
                    , mustWork = TRUE)

# TEST INPUT ARGUMENTS FORMAT
# ------------------------------------------------------------------------------
test_that("test input: filename", {
    # missing filename entry
    expect_error(
        readDesign(assembly = "hg19")
        , "argument \"filename\" is missing, with no default"
    )
    # file does not exists
    expect_error(
        readDesign(filename = "banana.txt", assembly = "hg19")
        , "file(s) do not exist", fixed=TRUE
    )
})

test_that("test inputs", {

    # missing entry
    expect_error(
        readDesign(filename = file)
        , "argument \"assembly\" is missing, with no default"
    )
    # wrong entry
    expect_error(
        readDesign(filename = file, assembly = "banana")
        , "No valid genome specified: please specify 'hg19' or 'hg38'"
    )
})

# TEST BED file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("test input: valid format BED file", {
    # BED file does not contain the minimal information: CHR START END.
    expect_error(
        readDesign(filename = system.file("extdata", "TestBED_LessThan3Columns.bed"
                            , package = "TMBleR"
                            , mustWork = TRUE)
                   , assembly = "hg19")
        , "Invalid BED file: less than 3 columns. Please provide a valid BED file with at least chromosome, start and end columns and check that columns are tab-delimited."
    )
    # Fields are not tab-separated.
    expect_error(
        readDesign(filename = system.file("extdata", "TestBED_NoTabDelim.bed"
                                          , package = "TMBleR"
                                          , mustWork = TRUE)
                  , assembly = "hg19")
        , "Invalid BED file: less than 3 columns. Please provide a valid BED file with at least chromosome, start and end columns and check that columns are tab-delimited."
    )
    # Chromosome field is not described using 'chr' prefix (i.e. chr2)
    expect_error(
        readDesign(filename = system.file("extdata", "TestBED_WrongChr.bed"
                                          , package = "TMBleR"
                                          , mustWork = TRUE)
                   , assembly = "hg19")
        , "Invalid BED file. Please check that all values in column one start with the prefix 'chr'")
    
})


# TEST OUTPUT OF FUNCTION  
# ----------------------------------------------------------------

# Read in input the panel sequencing design
designPanel_txt <- readDesign(system.file("extdata"
                                      , "ExamplePanel_GeneIDs_Only10GenesForTest.txt"
                                      , package = "TMBleR"
                                      , mustWork = TRUE)
                                      , assembly = "hg19"
                                      , ids = "entrezgene_id")

designPanel_bed <- readDesign(system.file("extdata"
                                          , "ExamplePanel_design_OnlyChr5and7ForTest.bed"
                                          , package = "TMBleR"
                                          , mustWork = TRUE)
                                      , assembly = "hg19")

## Test class of output objects from functions ##
test_that("Class of function output corresponds to the expected one", {
    expect_is(designPanel_txt, "GRanges")
    expect_is(designPanel_bed, "GRanges")
})

## Test dimensions of output objects from functions ##
test_that("Dimensions of function output correspond to the expected ones", {
    # Number of ranges of design read in input by readDesign
    expect_equal(length(designPanel_txt), 447)
    expect_equal(length(designPanel_bed), 1008)
})

