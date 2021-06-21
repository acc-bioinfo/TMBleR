context("applyFilters_hg38.R")

################################################################################
# Helper function
################################################################################
test_missing <- function(x)
{
    " Check the number of variants (filtered + passed) is equal to those evaluated by the filter (input vcf) "
    if(isFALSE(debug_env$debug)){
        print("Set debug_env$debug = TRUE in '/R/globalVariable.R'")
        debug_env$debug = TRUE
    }
    passed.var <- x[[1]][[1]] %>% nrow
    filtered.var <- debug_env$filters %>% nrow
    input.var <- hg38.vcfs[[1]]@fixed %>% nrow
    list(cat("Filtered:", filtered.var
             ,"\nPassed:", passed.var
             ,"\nInput Variants:", input.var
             ,"\nmissing:", input.var - sum(filtered.var, passed.var)))
    expect_true( input.var == sum(filtered.var, passed.var) )
    
}
# Conditions:
# 1. Set debug_env$debug = TRUE 
# 2. Launch only after the corresponding filter
#test_missing(hg38.vcfs_nonfiltered)
#test_missing(hg38.vcfs_NoCancer)

# Helper fn to check in the vcf file if the FILTER field contains the expected match
# in the original vcf file used as input. If returns the the FILTER field without the
# expected match
find_expected_filter_type <- function(variant_vector=hg38.vcfs$testhg38@fixed$FILTER, expected_filter_type){
    
    idx_filter <- purrr::map_lgl(variant_vector, ~ {
        out <- unlist(strsplit(.x,";", fixed = TRUE))
        out <- unlist(strsplit(out, "="))
        out <- unlist(strsplit(out, "|", fixed=TRUE))
        out <- out == expected_filter_type
        return(any(out))
    })
    
    return(variant_vector[!idx_filter])
}

################################################################################
# PREPARATION 
################################################################################

# Define vcf path and read in vcf files
hg38.vcf_files <- list(testhg38= system.file("extdata", "hg38.test.vcf", package = "TMBleR", mustWork = TRUE))
hg38.vcfs <- readVcfFiles(hg38.vcf_files, assembly = "hg38")

# Read in input the panel sequencing design
designPanel <- readDesign(system.file("extdata"
                                      , "ExamplePanel_design_OnlyChr7ForTest_hg38.bed"
                                      , package = "TMBleR"
                                      , mustWork = TRUE)
                          , assembly = "hg38")

# ------------------------------------------------------------------------------
# Apply different filters
# ------------------------------------------------------------------------------

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter 0: no filters - only filtered due to .bed file & intronic regions
hg38.vcfs_nonfiltered <- applyFilters(vcfs = hg38.vcfs,
                                      assembly = "hg38",
                                      design = designPanel)
#test_missing(hg38.vcfs_nonfiltered)

# This is the number of variants actually left out after filtering
left_over_count = nrow(hg38.vcfs_nonfiltered$testhg38$variants)
# this is the expected number of variants to be filtered
expected_left_over = find_expected_filter_type(expected_filter_type="off_target") %>%
    find_expected_filter_type(expected_filter_type="non_coding") %>% 
    length()

# Finally run the test
test_that("Test expected number of variants left to be filtered: No Filter", {
    expect_equal(left_over_count, expected_left_over)
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter 1: remove known cancer mutations (e.g. coding mutations described in
# COSMIC and truncating mutations in tumor suppressor genes)
hg38.vcfs_NoCancer <-suppressWarnings(applyFilters(vcfs = hg38.vcfs
                                   , assembly="hg38"
                                   , design=designPanel
                                   , remove.cancer=TRUE
                                   , tsList=NULL
                                   , variantType=NULL)
)
#test_missing(hg38.vcfs_NoCancer)
# This is the number of variants actually filtered out
left_over_count = nrow(hg38.vcfs_NoCancer$testhg38$variants)
# this is the expected number of variants to be filtered
expected_left_over = find_expected_filter_type(expected_filter_type="off_target") %>%
    find_expected_filter_type(expected_filter_type="non_coding") %>% 
    find_expected_filter_type(expected_filter_type="cancer") %>% 
    length()

test_that("Test expected number of variants to be filtered: Cancer Filter", {
    expect_equal(left_over_count, expected_left_over)
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter 2: remove known cancer mutations (e.g. coding mutations described in
# COSMIC and truncating mutations in tumor suppressor genes) and synonymous
# mutations
hg38.vcfs_NoCancer_NoSynonymous <- suppressWarnings(applyFilters(vcfs = hg38.vcfs,
                                                                 assembly="hg38",
                                                                 design=designPanel,
                                                                 remove.cancer=TRUE,
                                                                 tsList=NULL,
                                                                 variantType=c("synonymous")))
#test_missing(hg38.vcfs_NoCancer_NoSynonymous)
# This is the number of variants actually filtered out
left_over_count = nrow(hg38.vcfs_NoCancer_NoSynonymous$testhg38$variants)
# this is the expected number of variants to be filtered
expected_filtered_count = find_expected_filter_type(expected_filter_type="off_target") %>%
    find_expected_filter_type(expected_filter_type="non_coding") %>%
    find_expected_filter_type(expected_filter_type="cancer") %>%
    find_expected_filter_type(expected_filter_type="synonimous") %>%
    length()

test_that("Test expected number of variants to be filtered: Cancer Filter + synonimous", {
    expect_equal(left_over_count, expected_filtered_count)
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter 3: remove synonymous mutations
hg38.vcfs_NoSynonymous <- applyFilters(vcfs = hg38.vcfs,
                                       assembly="hg38",
                                       design=designPanel,
                                       remove.cancer=FALSE,
                                       tsList=NULL,
                                       variantType=c("synonymous"))
#test_missing(hg38.vcfs_NoSynonymous)
left_over_count = nrow(hg38.vcfs_NoSynonymous$testhg38$variants)
# this is the expected number of variants to be filtered
expected_filtered_count = find_expected_filter_type(expected_filter_type="off_target") %>%
    find_expected_filter_type(expected_filter_type="non_coding") %>%
    find_expected_filter_type(expected_filter_type="synonimous") %>%
    length()

test_that("Test expected number of variants to be filtered: synonimous", {
    expect_equal(left_over_count, expected_filtered_count)
})


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter 4: remove non-synonymous mutations
hg38.vcfs_NoNonSynonymous <- applyFilters(vcfs = hg38.vcfs,
                                          assembly="hg38",
                                          design=designPanel,
                                          remove.cancer=FALSE,
                                          tsList=NULL,
                                          variantType=c("nonsynonymous"))
#test_missing(hg38.vcfs_NoNonSynonymous)
left_over_count = nrow(hg38.vcfs_NoNonSynonymous$testhg38$variants)
# this is the expected number of variants to be filtered
expected_filtered_count = find_expected_filter_type(expected_filter_type="off_target") %>%
    find_expected_filter_type(expected_filter_type="non_coding") %>%
    find_expected_filter_type(expected_filter_type="nonsynonimous") %>%
    length()

test_that("Test expected number of variants to be filtered: nonsynonimous ", {
    expect_equal(left_over_count, expected_filtered_count)
})


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter 5: remove known cancer mutations (e.g. coding mutations described in
# COSMIC and truncating mutations in tumor suppressor genes) and non-synonymous
# mutations
hg38.vcfs_NoCancerNonSynonymous <-suppressWarnings(applyFilters(vcfs = hg38.vcfs,
                                                assembly="hg38",
                                                design=designPanel,
                                                remove.cancer=TRUE,
                                                tsList=NULL,
                                                variantType=c("nonsynonymous"))
)
#test_missing(hg38.vcfs_NoCancerNonSynonymous)
left_over_count = nrow(hg38.vcfs_NoCancerNonSynonymous$testhg38$variants)
# this is the expected number of variants to be filtered
expected_filtered_count = find_expected_filter_type(expected_filter_type="off_target") %>%
    find_expected_filter_type(expected_filter_type="non_coding") %>%
    find_expected_filter_type(expected_filter_type="cancer") %>%
    find_expected_filter_type(expected_filter_type="nonsynonimous") %>%
    length()

test_that("Test expected number of variants to be filtered: Cancer Filter + nonsynonymous", {
    expect_equal(left_over_count, expected_filtered_count)
})


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter 6: remove known cancer mutations (e.g. coding mutations described in
# COSMIC and truncating mutations in tumor suppressor genes) and synonymous + non-synonymous
# mutations
hg38.vcfs_NoCancerSynonimousNonSynonymous <- suppressWarnings(applyFilters(vcfs = hg38.vcfs,
                                                          assembly="hg38",
                                                          design=designPanel,
                                                          remove.cancer=TRUE,
                                                          tsList=NULL,
                                                          variantType=c("nonsynonymous", "synonymous"))
)
#test_missing(hg38.vcfs_NoCancerSynonimousNonSynonymous)
left_over_count = nrow(hg38.vcfs_NoCancerSynonimousNonSynonymous$testhg38$variants)
# this is the expected number of variants to be filtered
expected_filtered_count = find_expected_filter_type(expected_filter_type="off_target") %>%
    find_expected_filter_type(expected_filter_type="non_coding") %>%
    find_expected_filter_type(expected_filter_type="cancer") %>%
    find_expected_filter_type(expected_filter_type="nonsynonimous") %>%
    find_expected_filter_type(expected_filter_type="synonimous") %>%
    length()

test_that("Test expected number of variants to be filtered: Cancer Filter + nonsynonymous + nonsynonymous", {
    expect_equal(left_over_count, expected_filtered_count)
})


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter 7: remove non-sense mutations
hg38.vcfs_NoNonsense <- applyFilters(vcfs = hg38.vcfs,
                                     assembly="hg38",
                                     design=designPanel,
                                     remove.cancer=FALSE,
                                     tsList=NULL,
                                     variantType=c("nonsense"))
#test_missing(hg38.vcfs_NoNonsense)
left_over_count = nrow(hg38.vcfs_NoNonsense$testhg38$variants)
# this is the expected number of variants to be filtered
expected_filtered_count = find_expected_filter_type(expected_filter_type="off_target") %>%
    find_expected_filter_type(expected_filter_type="non_coding") %>%
    find_expected_filter_type(expected_filter_type="nonsense") %>%
    length()

test_that("Test expected number of variants to be filtered: nonsense", {
    expect_equal(left_over_count, expected_filtered_count)
})


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter 8: remove known cancer mutations (e.g. coding mutations described in
# COSMIC and truncating mutations in tumor suppressor genes) and nonsense
# mutations
hg38.vcfs_NoCancerNonsense <- suppressWarnings(applyFilters(vcfs = hg38.vcfs,
                                           assembly="hg38",
                                           design=designPanel,
                                           remove.cancer=TRUE,
                                           tsList=NULL,
                                           variantType=c("nonsense"))
)
#test_missing(hg38.vcfs_NoCancerNonsense)

left_over_count = nrow(hg38.vcfs_NoCancerNonsense$testhg38$variants)
# this is the expected number of variants to be filtered
expected_filtered_count = find_expected_filter_type(expected_filter_type="off_target") %>%
    find_expected_filter_type(expected_filter_type="non_coding") %>%
    find_expected_filter_type(expected_filter_type="cancer") %>%
    find_expected_filter_type(expected_filter_type="nonsense") %>%
    length()

test_that("Test expected number of variants to be filtered: nonsense + cancer", {
    expect_equal(left_over_count, expected_filtered_count)
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter 9: remove frameshift mutations
hg38.vcfs_NoFrameshift <- applyFilters(vcfs = hg38.vcfs,
                                       assembly="hg38",
                                       design=designPanel,
                                       remove.cancer=FALSE,
                                       tsList=NULL,
                                       variantType=c("frameshift"))
#test_missing(hg38.vcfs_NoFrameshift)
left_over_count = nrow(hg38.vcfs_NoFrameshift$testhg38$variants)
# this is the expected number of variants to be filtered
expected_filtered_count = find_expected_filter_type(expected_filter_type="off_target") %>%
    find_expected_filter_type(expected_filter_type="non_coding") %>%
    find_expected_filter_type(expected_filter_type="frameshift") %>%
    length()

test_that("Test expected number of variants to be filtered: frameshift", {
    expect_equal(left_over_count, expected_filtered_count)
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter 10: remove known cancer mutations (e.g. coding mutations described in
# COSMIC and truncating mutations in tumor suppressor genes) and frameshift
# mutations
hg38.vcfs_NoCancerFrameshift <- suppressWarnings(applyFilters(vcfs = hg38.vcfs,
                                             assembly="hg38",
                                             design=designPanel,
                                             remove.cancer=TRUE,
                                             tsList=NULL,
                                             variantType=c("frameshift"))
)
#test_missing(hg38.vcfs_NoCancerFrameshift)
left_over_count = nrow(hg38.vcfs_NoCancerFrameshift$testhg38$variants)
# this is the expected number of variants to be filtered
expected_filtered_count = find_expected_filter_type(expected_filter_type="off_target") %>%
    find_expected_filter_type(expected_filter_type="non_coding") %>%
    find_expected_filter_type(expected_filter_type="frameshift") %>%
    find_expected_filter_type(expected_filter_type="cancer") %>%
    length()

test_that("Test expected number of variants to be filtered: frameshift", {
    expect_equal(left_over_count, expected_filtered_count)
})















# Filter 11: remove synonymous mutations and mutations with variant allele
# frequency < 0.05
# vcfs_NoSynonymous_VAFFilter <- applyFilters(vcfs,
#                                                  assembly="hg19",
#                                                  design=designPanel,
#                                                  vaf.cutoff=0.05,
#                                                  remove.cancer=FALSE,
#                                                  tsList=NULL,
#                                                  variantType=c("synonymous"))

# # Filter 12: remove mutations with variant allele frequency < 0.05
# vcfs_VAFFilter <- applyFilters(vcfs,
#                                     assembly="hg19",
#                                     design=designPanel,
#                                     vaf.cutoff=0.05,
#                                     remove.cancer=FALSE,
#                                     tsList=NULL,
#                                     variantType=NULL)
#

# Filter 13: remove known cancer mutations (e.g. coding mutations described in
# COSMIC and truncating mutations in tumor suppressor genes) and mutations with
# variant allele frequency < 0.05
# vcfs_NoCancer_VAFFilter <- suppressWarnings(applyFilters(vcfs,
#                                              assembly="hg19",
#                                              design=designPanel,
#                                              vaf.cutoff=0.05,
#                                              remove.cancer=TRUE,
#                                              tsList=NULL,
#                                              variantType=NULL))

# Filter 14: remove known cancer mutations (e.g. coding mutations described in
# COSMIC and truncating mutations in tumor suppressor genes), synonymous
# mutations and mutations with variant allele frequency < 0.05
# vcfs_NoCancer_VAFFilter_NoSynonymous <- suppressWarnings(applyFilters(vcfs,
#                                                           assembly="hg19",
#                                                           design=designPanel,
#                                                           vaf.cutoff=0.05,
#                                                           remove.cancer=TRUE,
#                                                           tsList=NULL,
#                                                           variantType=c("synonymous")))


################################################################################
# TESTS
################################################################################

# TEST WRONG INPUT FORMAT
# ------------------------------------------------------------------------------
test_that("test wrong input: hg38.vcfs", {
    # missing entry
    expect_error(
        applyFilters(assembly = "hg38", design = design_gr)
        , "argument \"vcfs\" is missing, with no default"
    )
    # Object does not exists
    expect_error(
        applyFilters(vcfs = banana_vcfs, assembly = "hg38", design = designPanel)
        , "object 'banana_vcfs' not found", fixed=TRUE
    )
})
test_that("test expect wrong input: assembly", {
    # missing entry
    expect_error(
        applyFilters(vcfs = hg38.vcfs, design = designPanel)
        , "argument \"assembly\" is missing, with no default"
    )
    # wrong entry
    expect_error(
        applyFilters(vcfs = hg38.vcfs, assembly = "banana", design = designPanel)
        , "No valid genome assembly: please specify 'hg19' or 'hg38'"
    )
})
test_that("test wrong input: design", {
    # missing entry
    expect_error(
        applyFilters(vcfs = hg38.vcfs, assembly = "hg38")
        , "argument \"design\" is missing, with no default"
    )
    # Object does not exists
    expect_error(
        applyFilters(vcfs = hg38.vcfs, assembly = "hg38", design = designBanana)
        , "object 'designBanana' not found", fixed=TRUE
    )
})
test_that("test wrong input: vaf.cutoff", {
    # wrong format entry
    expect_error(
        applyFilters(vcfs = hg38.vcfs, assembly = "hg38", design = designPanel, vaf.cutoff = "ciao", remove.cancer = FALSE, tsList = NULL, variantType = NULL)
        , "vaf.cutoff should be numeric type"
    )
    # wrong range value entry
    expect_error(
        applyFilters(vcfs = hg38.vcfs, assembly = "hg38", design = designPanel, vaf.cutoff = 30)
        , "vaf.cutoff can't be higher than 1"
    )
})
test_that("test input: vaf.cutoff", {
    # wrong format entry
    expect_error(
        applyFilters(vcfs = hg38.vcfs, assembly = "hg38", design = designPanel, vaf.cutoff = NULL)
        , "vaf.cutoff should be numeric type"
    )


})


# TEST OUTPUT FORMAT FROM THE FILTER
# ----------------------------------------------------------------
# make sure the output from the filter is in the correct format and consistent
test_that("Class of function output corresponds to the expected one", {
    expect_is(hg38.vcfs_nonfiltered, "list")
    expect_is(hg38.vcfs_NoCancer, "list")
    expect_is(hg38.vcfs_NoCancer[[1]], "list")
    expect_is(hg38.vcfs_NoCancer[[1]]$variants, "data.frame")
    expect_is(hg38.vcfs_NoCancer[[1]]$filter, "character")
    expect_is(hg38.vcfs_NoCancer[[1]]$design, "GRanges")
    expect_is(hg38.vcfs_NoCancer[[1]]$sample, "character")
    expect_is(hg38.vcfs_NoCancer_NoSynonymous, "list")
    expect_is(hg38.vcfs_NoCancer_NoSynonymous[[1]], "list")
    expect_is(hg38.vcfs_NoCancer_NoSynonymous[[1]]$variants, "data.frame")
    expect_is(hg38.vcfs_NoCancer_NoSynonymous[[1]]$filter, "character")
    expect_is(hg38.vcfs_NoCancer_NoSynonymous[[1]]$design, "GRanges")
    expect_is(hg38.vcfs_NoCancer_NoSynonymous[[1]]$sample, "character")
    expect_is(hg38.vcfs_NoSynonymous, "list")
    expect_is(hg38.vcfs_NoSynonymous[[1]], "list")
    expect_is(hg38.vcfs_NoSynonymous[[1]]$variants, "data.frame")
    expect_is(hg38.vcfs_NoSynonymous[[1]]$filter, "character")
    expect_is(hg38.vcfs_NoSynonymous[[1]]$design, "GRanges")
    expect_is(hg38.vcfs_NoSynonymous[[1]]$sample, "character")
    # expect_is(vcfs_NoSynonymous_VAFFilter, "list")
    # expect_is(vcfs_NoSynonymous_VAFFilter[[1]], "list")
    # expect_is(vcfs_NoSynonymous_VAFFilter[[1]]$variants, "data.frame")
    # expect_is(vcfs_NoSynonymous_VAFFilter[[1]]$filter, "character")
    # expect_is(vcfs_NoSynonymous_VAFFilter[[1]]$design, "GRanges")
    # expect_is(vcfs_NoSynonymous_VAFFilter[[1]]$sample, "character")
    # expect_is(vcfs_VAFFilter, "list")
    # expect_is(vcfs_VAFFilter[[1]], "list")
    # expect_is(vcfs_VAFFilter[[1]]$variants, "data.frame")
    # expect_is(vcfs_VAFFilter[[1]]$filter, "character")
    # expect_is(vcfs_VAFFilter[[1]]$design, "GRanges")
    # expect_is(vcfs_VAFFilter[[1]]$sample, "character")
    # expect_is(vcfs_NoCancer_VAFFilter, "list")
    # expect_is(vcfs_NoCancer_VAFFilter[[1]], "list")
    # expect_is(vcfs_NoCancer_VAFFilter[[1]]$variants, "data.frame")
    # expect_is(vcfs_NoCancer_VAFFilter[[1]]$filter, "character")
    # expect_is(vcfs_NoCancer_VAFFilter[[1]]$design, "GRanges")
    # expect_is(vcfs_NoCancer_VAFFilter[[1]]$sample, "character")
    # expect_is(vcfs_NoCancer_VAFFilter_NoSynonymous, "list")
    # expect_is(vcfs_NoCancer_VAFFilter_NoSynonymous[[1]], "list")
    # expect_is(vcfs_NoCancer_VAFFilter_NoSynonymous[[1]]$variants, "data.frame")
    # expect_is(vcfs_NoCancer_VAFFilter_NoSynonymous[[1]]$filter, "character")
    # expect_is(vcfs_NoCancer_VAFFilter_NoSynonymous[[1]]$design, "GRanges")
    # expect_is(vcfs_NoCancer_VAFFilter_NoSynonymous[[1]]$sample, "character")
})

## Test dimensions of output objects from functions ##
test_that("Number of variants generated in output by applyFilter function is not higher than number of variants in input", {
    # Compare n variants in input expanded VCF object and in output dataframe
    extraction <- hg38.vcfs[[1]] %>%
        S4Vectors::expand() %>%
        SummarizedExperiment::rowRanges() %>%
        length()

    expect_true(extraction >= dim(hg38.vcfs_NoCancer[[1]]$variants)[1])
    expect_true(extraction >= dim(hg38.vcfs_NoCancer_NoSynonymous[[1]]$variants)[1])
    expect_true(extraction >= dim(hg38.vcfs_NoSynonymous[[1]]$variants)[1])
    # expect_true(extraction >= dim(vcfs_NoSynonymous_VAFFilter[[1]]$variants)[1])
    # expect_true(extraction >= dim(vcfs_VAFFilter[[1]]$variants)[1])
    # expect_true(extraction >= dim(vcfs_NoCancer_VAFFilter[[1]]$variants)[1])
    # expect_true(extraction >= dim(vcfs_NoCancer_VAFFilter_NoSynonymous[[1]]$variants)[1])
})
test_that("Dimensions of function output correspond to the expected ones", {
    # Number of elements in lists generated by applyFilters
    expect_equal(length(hg38.vcfs_NoCancer[[1]]), 4)
    expect_equal(length(hg38.vcfs_NoCancer_NoSynonymous[[1]]), 4)
    expect_equal(length(hg38.vcfs_NoSynonymous[[1]]), 4)
    # expect_equal(length(vcfs_NoSynonymous_VAFFilter[[1]]), 4)
    # expect_equal(length(vcfs_VAFFilter[[1]]), 4)
    # expect_equal(length(vcfs_NoCancer_VAFFilter[[1]]), 4)
    # expect_equal(length(vcfs_NoCancer_VAFFilter_NoSynonymous[[1]]), 4)
})
test_that("Number of variants filtered correspond to the expected ones", {
    # Number of elements in lists generated by applyFilters
    expect_equal(nrow(hg38.vcfs_nonfiltered[[1]][["variants"]]), 31)
    expect_equal(nrow(hg38.vcfs_NoCancer[[1]][["variants"]]), 14)
    expect_equal(nrow(hg38.vcfs_NoCancer_NoSynonymous[[1]][["variants"]]), 7) 
    expect_equal(nrow(hg38.vcfs_NoSynonymous[[1]][["variants"]]), 24)           
    expect_equal(nrow(hg38.vcfs_NoNonSynonymous[[1]][["variants"]]), 23) 
    expect_equal(nrow(hg38.vcfs_NoCancerNonSynonymous[[1]][["variants"]]), 7)  
    expect_equal(nrow(hg38.vcfs_NoCancerSynonimousNonSynonymous[[1]][["variants"]]), 0)
    expect_equal(nrow(hg38.vcfs_NoNonsense[[1]][["variants"]]), 22)   
    expect_equal(nrow(hg38.vcfs_NoCancerNonsense[[1]][["variants"]]), 14)
    expect_equal(nrow(hg38.vcfs_NoFrameshift[[1]][["variants"]]), 27)
    expect_equal(nrow(hg38.vcfs_NoCancerFrameshift[[1]][["variants"]]), 14)
    # expect_equal(length(vcfs_NoSynonymous_VAFFilter[[1]]), 4)
    # expect_equal(length(vcfs_VAFFilter[[1]]), 4)
    # expect_equal(length(vcfs_NoCancer_VAFFilter[[1]]), 4)
    # expect_equal(length(vcfs_NoCancer_VAFFilter_NoSynonymous[[1]]), 4)
})





context("applyTMB.R")

TMB_vcfs0=applyTMB(inputForTMB = hg38.vcfs_nonfiltered, assembly = "hg38")
TMB_vcfs1=applyTMB(inputForTMB = hg38.vcfs_NoCancer, assembly = "hg38")
TMB_vcfs2=applyTMB(inputForTMB = hg38.vcfs_NoCancer_NoSynonymous, assembly = "hg38")
TMB_vcfs3=applyTMB(inputForTMB = hg38.vcfs_NoSynonymous, assembly = "hg38")
TMB_vcfs4=applyTMB(inputForTMB = hg38.vcfs_NoNonSynonymous, assembly = "hg38")
TMB_vcfs5=applyTMB(inputForTMB = hg38.vcfs_NoCancerNonSynonymous, assembly = "hg38")
TMB_vcfs6=applyTMB(inputForTMB = hg38.vcfs_NoCancerSynonimousNonSynonymous, assembly = "hg38")
TMB_vcfs7=applyTMB(inputForTMB = hg38.vcfs_NoNonsense, assembly = "hg38")
TMB_vcfs8=applyTMB(inputForTMB = hg38.vcfs_NoCancerNonsense, assembly = "hg38")
TMB_vcfs9=applyTMB(inputForTMB = hg38.vcfs_NoFrameshift, assembly = "hg38")
TMB_vcfs10=applyTMB(inputForTMB = hg38.vcfs_NoCancerFrameshift, assembly = "hg38")
# TMB_vcfs11=applyTMB(inputForTMB = vcfs_NoCancer, assembly = "hg19")
# TMB_vcfs12=applyTMB(inputForTMB = vcfs_NoCancer, assembly = "hg19")
# TMB_vcfs13=applyTMB(inputForTMB = vcfs_NoCancer, assembly = "hg19")
# TMB_vcfs14=applyTMB(inputForTMB = vcfs_NoCancer, assembly = "hg19")

# TESTS ON OUTPUT #
## Test class of output objects from functions ##
test_that("Class of function output corresponds to the expected one", {
    expect_is(TMB_vcfs0, "data.frame")
    # expect_is(TMB_vcfs$Sample, "factor")
    # expect_is(TMB_vcfs$Filter, "factor")
    expect_is(TMB_vcfs0$Sequencing_Size, "numeric")
    expect_is(TMB_vcfs0$Tot_Number_Mutations, "numeric")
    expect_is(TMB_vcfs0$TMB_per_Mb, "numeric")
})

test_that("Dimensions of function output correspond to the expected ones", {
    expect_equal(dim(TMB_vcfs0)[1], 1)
    expect_equal(dim(TMB_vcfs0)[2], 6)
})
test_that("TMB results correspond to the expected ones", {
    expect_equal(TMB_vcfs0$TMB_per_Mb,0.19)
    expect_equal(TMB_vcfs1$TMB_per_Mb,0.09)
    expect_equal(TMB_vcfs2$TMB_per_Mb,0.04)
    expect_equal(TMB_vcfs3$TMB_per_Mb,0.15)
    expect_equal(TMB_vcfs4$TMB_per_Mb,0.14)
    expect_equal(TMB_vcfs5$TMB_per_Mb,0.04)
    expect_equal(TMB_vcfs6$TMB_per_Mb,0.00)
    expect_equal(TMB_vcfs7$TMB_per_Mb,0.14)
    expect_equal(TMB_vcfs8$TMB_per_Mb,0.09)
    expect_equal(TMB_vcfs9$TMB_per_Mb,0.17)
    expect_equal(TMB_vcfs10$TMB_per_Mb,0.09)
    # expect_equal(TMB_vcfs11$TMB_per_Mb,0)
    # expect_equal(TMB_vcfs12$TMB_per_Mb,0)
    # expect_equal(TMB_vcfs13$TMB_per_Mb,0)
    # expect_equal(TMB_vcfs14$TMB_per_Mb,0)

})

