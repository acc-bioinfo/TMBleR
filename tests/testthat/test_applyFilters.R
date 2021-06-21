context("applyFilters.R")

# Tests
################################################################################

# Read vcf filenames
vcf_files <- list(Horizon5="Horizon5_ExamplePanel_OnlyChr5and7ForTest.vcf")
# For each vcf file, get the absolute path
vcf_files <- lapply(vcf_files
                    , function(x) system.file("extdata", x, package = "TMBleR", mustWork = TRUE))
# Read in the files
vcfs <- readVcfFiles(vcf_files, assembly = "hg19")


# Read in input the panel sequencing design
designPanel <- readDesign(system.file("extdata"
                                      , "ExamplePanel_GeneIDs_Only10GenesForTest.txt"
                                      , package = "TMBleR"
                                      , mustWork = TRUE)
                                      , assembly = "hg19"
                                      , ids = "entrezgene_id")

# TEST INPUT FORMAT
# ------------------------------------------------------------------------------
test_that("test input: vcfs", {
    # missing entry
    expect_error(
        applyFilters(assembly = "hg19", design = designPanel)
        , "argument \"vcfs\" is missing, with no default"
    )
    # Object does not exists
    expect_error(
        applyFilters(vcfs = banana_vcfs, assembly = "hg19", design = designPanel)
        , "object 'banana_vcfs' not found", fixed=TRUE
    )
})
test_that("test input: assembly", {
    # missing entry
    expect_error(
        applyFilters(vcfs = vcfs, design = designPanel)
        , "argument \"assembly\" is missing, with no default"
    )
    # wrong entry
    expect_error(
        applyFilters(vcfs = vcfs, assembly = "banana", design = designPanel)
        , "No valid genome assembly: please specify 'hg19' or 'hg38'"
    )
})
test_that("test input: design", {
    # missing entry
    expect_error(
        applyFilters(vcfs = vcfs, assembly = "hg19")
        , "argument \"design\" is missing, with no default"
    )
    # Object does not exists
    expect_error(
        applyFilters(vcfs, assembly = "hg19", design = designBanana)
        , "object 'designBanana' not found", fixed=TRUE
    )
})
test_that("test input: vaf.cutoff", {
    # wrong format entry
    expect_error(
        applyFilters(vcfs = vcfs, assembly = "hg19", design = designPanel, vaf.cutoff = "ciao", remove.cancer = FALSE, tsList = NULL, variantType = NULL)
        , "vaf.cutoff should be numeric type"
    )
    # wrong range value entry
    expect_error(
        applyFilters(vcfs = vcfs, assembly = "hg19", design = designPanel, vaf.cutoff = 30)
        , "vaf.cutoff can't be higher than 1"
    )
})
test_that("test input: vaf.cutoff", {
    # wrong format entry
    expect_error(
        applyFilters(vcfs = vcfs, assembly = "hg19", design = designPanel, vaf.cutoff = NULL)
        , "vaf.cutoff should be numeric type"
    )
    
    
})

# test_that("test input: tsList", {
#     # wrong format entry
#     expect_error(
#         applyFilters(vcfs = vcfs
#                      , assembly = "hg19"
#                      , design = designPanel
#                      , remove.cancer = TRUE
#                      , tsList = as.numeric(as.vector(c(1,2,3))), variantType = NULL)
#         , "No valid tsList: please provide a character vector with entrez_gene ids"
#     )
# })



# TEST OUTPUT OF FUNCTION  
# ----------------------------------------------------------------

# CODE TO GENERATE OUTPUT #
# Apply different filters
# Filter 1: remove known cancer mutations (e.g. coding mutations described in 
# COSMIC and truncating mutations in tumor suppressor genes)
vcfs_NoCancer <- suppressWarnings(applyFilters(  vcfs
                              , assembly="hg19"
                              , design=designPanel
                              , remove.cancer=TRUE
                             , tsList=NULL
                             , variantType=NULL))

# Filter 2: remove known cancer mutations (e.g. coding mutations described in 
# COSMIC and truncating mutations in tumor suppressor genes) and synonymous 
# mutations
vcfs_NoCancer_NoSynonymous <- suppressWarnings(applyFilters(vcfs, 
                                                assembly="hg19",
                                                design=designPanel, 
                                                remove.cancer=TRUE, 
                                                tsList=NULL, 
                                                variantType=c("synonymous")))

# Filter 3: remove synonymous mutations
vcfs_NoSynonymous <- applyFilters(vcfs, 
                                       assembly="hg19", 
                                       design=designPanel, 
                                       remove.cancer=FALSE, 
                                       tsList=NULL, 
                                       variantType=c("synonymous"))

# Filter 4: remove synonymous mutations and mutations with variant allele 
# frequency < 0.05
vcfs_NoSynonymous_VAFFilter <- applyFilters(vcfs, 
                                                 assembly="hg19", 
                                                 design=designPanel, 
                                                 vaf.cutoff=0.05, 
                                                 remove.cancer=FALSE, 
                                                 tsList=NULL, 
                                                 variantType=c("synonymous"))

# Filter 5: remove mutations with variant allele frequency < 0.05
vcfs_VAFFilter <- applyFilters(vcfs, 
                                    assembly="hg19", 
                                    design=designPanel, 
                                    vaf.cutoff=0.05,
                                    remove.cancer=FALSE, 
                                    tsList=NULL, 
                                    variantType=NULL)

# Filter 6: reemove known cancer mutations (e.g. coding mutations described in 
# COSMIC and truncating mutations in tumor suppressor genes) and mutations with 
# variant allele frequency < 0.05
vcfs_NoCancer_VAFFilter <- suppressWarnings(applyFilters(vcfs, 
                                             assembly="hg19", 
                                             design=designPanel, 
                                             vaf.cutoff=0.05, 
                                             remove.cancer=TRUE, 
                                             tsList=NULL, 
                                             variantType=NULL))

# Filter 7: remove known cancer mutations (e.g. coding mutations described in 
# COSMIC and truncating mutations in tumor suppressor genes), synonymous 
# mutations and mutations with variant allele frequency < 0.05
vcfs_NoCancer_VAFFilter_NoSynonymous <- suppressWarnings(applyFilters(vcfs, 
                                                          assembly="hg19", 
                                                          design=designPanel, 
                                                          vaf.cutoff=0.05, 
                                                          remove.cancer=TRUE, 
                                                          tsList=NULL, 
                                                          variantType=c("synonymous")))


# TESTS ON OUTPUT #
## Test class of output objects from functions ##
test_that("Class of function output corresponds to the expected one", {
    expect_is(vcfs_NoCancer, "list")
    expect_is(vcfs_NoCancer[[1]], "list")
    expect_is(vcfs_NoCancer[[1]]$variants, "data.frame")
    expect_is(vcfs_NoCancer[[1]]$filter, "character")
    expect_is(vcfs_NoCancer[[1]]$design, "GRanges")
    expect_is(vcfs_NoCancer[[1]]$sample, "character")
    expect_is(vcfs_NoCancer_NoSynonymous, "list")
    expect_is(vcfs_NoCancer_NoSynonymous[[1]], "list")
    expect_is(vcfs_NoCancer_NoSynonymous[[1]]$variants, "data.frame")
    expect_is(vcfs_NoCancer_NoSynonymous[[1]]$filter, "character")
    expect_is(vcfs_NoCancer_NoSynonymous[[1]]$design, "GRanges")
    expect_is(vcfs_NoCancer_NoSynonymous[[1]]$sample, "character")
    expect_is(vcfs_NoSynonymous, "list")
    expect_is(vcfs_NoSynonymous[[1]], "list")
    expect_is(vcfs_NoSynonymous[[1]]$variants, "data.frame")
    expect_is(vcfs_NoSynonymous[[1]]$filter, "character")
    expect_is(vcfs_NoSynonymous[[1]]$design, "GRanges")
    expect_is(vcfs_NoSynonymous[[1]]$sample, "character")
    expect_is(vcfs_NoSynonymous_VAFFilter, "list")
    expect_is(vcfs_NoSynonymous_VAFFilter[[1]], "list")
    expect_is(vcfs_NoSynonymous_VAFFilter[[1]]$variants, "data.frame")
    expect_is(vcfs_NoSynonymous_VAFFilter[[1]]$filter, "character")
    expect_is(vcfs_NoSynonymous_VAFFilter[[1]]$design, "GRanges")
    expect_is(vcfs_NoSynonymous_VAFFilter[[1]]$sample, "character")
    expect_is(vcfs_VAFFilter, "list")
    expect_is(vcfs_VAFFilter[[1]], "list")
    expect_is(vcfs_VAFFilter[[1]]$variants, "data.frame")
    expect_is(vcfs_VAFFilter[[1]]$filter, "character")
    expect_is(vcfs_VAFFilter[[1]]$design, "GRanges")
    expect_is(vcfs_VAFFilter[[1]]$sample, "character")
    expect_is(vcfs_NoCancer_VAFFilter, "list")
    expect_is(vcfs_NoCancer_VAFFilter[[1]], "list")
    expect_is(vcfs_NoCancer_VAFFilter[[1]]$variants, "data.frame")
    expect_is(vcfs_NoCancer_VAFFilter[[1]]$filter, "character")
    expect_is(vcfs_NoCancer_VAFFilter[[1]]$design, "GRanges")
    expect_is(vcfs_NoCancer_VAFFilter[[1]]$sample, "character")
    expect_is(vcfs_NoCancer_VAFFilter_NoSynonymous, "list")
    expect_is(vcfs_NoCancer_VAFFilter_NoSynonymous[[1]], "list")
    expect_is(vcfs_NoCancer_VAFFilter_NoSynonymous[[1]]$variants, "data.frame")
    expect_is(vcfs_NoCancer_VAFFilter_NoSynonymous[[1]]$filter, "character")
    expect_is(vcfs_NoCancer_VAFFilter_NoSynonymous[[1]]$design, "GRanges")
    expect_is(vcfs_NoCancer_VAFFilter_NoSynonymous[[1]]$sample, "character")
})

## Test dimensions of output objects from functions ##
test_that("Number of variants generated in output by applyFilter function is not higher than number of variants in input", {
    # Compare n variants in input expanded VCF object and in output dataframe
    extraction <- vcfs[[1]] %>% 
        S4Vectors::expand() %>%
        SummarizedExperiment::rowRanges() %>%
        length()
    
    expect_true(extraction >= dim(vcfs_NoCancer[[1]]$variants)[1])
    expect_true(extraction >= dim(vcfs_NoCancer_NoSynonymous[[1]]$variants)[1])
    expect_true(extraction >= dim(vcfs_NoSynonymous[[1]]$variants)[1])
    expect_true(extraction >= dim(vcfs_NoSynonymous_VAFFilter[[1]]$variants)[1])
    expect_true(extraction >= dim(vcfs_VAFFilter[[1]]$variants)[1])
    expect_true(extraction >= dim(vcfs_NoCancer_VAFFilter[[1]]$variants)[1])
    expect_true(extraction >= dim(vcfs_NoCancer_VAFFilter_NoSynonymous[[1]]$variants)[1])
})
test_that("Dimensions of function output correspond to the expected ones", {
    # Number of elements in lists generated by applyFilters
    expect_equal(length(vcfs_NoCancer[[1]]), 4)
    expect_equal(length(vcfs_NoCancer_NoSynonymous[[1]]), 4)
    expect_equal(length(vcfs_NoSynonymous[[1]]), 4)
    expect_equal(length(vcfs_NoSynonymous_VAFFilter[[1]]), 4)
    expect_equal(length(vcfs_VAFFilter[[1]]), 4)
    expect_equal(length(vcfs_NoCancer_VAFFilter[[1]]), 4)
    expect_equal(length(vcfs_NoCancer_VAFFilter_NoSynonymous[[1]]), 4)
})


