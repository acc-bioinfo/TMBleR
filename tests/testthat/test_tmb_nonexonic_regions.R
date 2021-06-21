# Note: This test is not directly connected to any other script.R 
#
# In this script we want to make sure we are testing the following behavious:
# 
# CONDITIONS
# -----------------
# * VCF included in SNV mapped in non-exonic positions;
# * Such positions happened to be included in the target region of the panel;
# 
# EXPECTE BEHAVIOUR
# -----------------
# * The SNV shouldn't not be  removed by any filter and not considered in the TMB estimation

context("test_tmb_nonexonic_regions.R")


# Environment Preparation for TEST using hg19 ================================================================
genome_assembly <- "hg19"

# load a vcf including SNV in non exonic regions 
vcf_files <- list(nonexoni_vcf="Sample1_ExampleWES_chr7.vcf")
# Find it in the file system
vcf_files <- lapply(vcf_files, function(x) 
    system.file("extdata", x, package = "TMBleR", mustWork = TRUE))

# read in the files 
vcfs <- readVcfFiles(vcfFiles = vcf_files, assembly = genome_assembly)

# PREPARE THE DESIGN FILE ------------------------------------------
# Read bed file with the panel sequencing design 
design_df <- data.frame(seqnames = "chr7"
                        , start = "1"
                        , end = "40000")
design_gr <- GenomicRanges::makeGRangesFromDataFrame(design_df)
design_gr@metadata$design_name <- "design_name"

# Check expected intersection
expected_intersection <- IRanges::subsetByOverlaps(
    vcfs$nonexoni_vcf@rowRanges, design_gr) 

# Run filter and Test ----- ----- ----- ----- ----- ----- ----- -----

# This applies no filter - Only SNV not mapped inside the design should be removed
vcfs_noFilter <- applyFilters(  vcfs = vcfs
                                , assembly = genome_assembly
                                , design = design_gr
                                , remove.nonexonic = FALSE)
test_that("test No Filter", {
    expect_equal(length(expected_intersection), nrow(vcfs_noFilter$nonexoni_vcf$variants))
})


vcfs_vafFilter <- applyFilters(  vcfs = vcfs
                                , assembly = genome_assembly
                                , design = design_gr
                                , vaf.cutoff = 0.3
                                , remove.nonexonic = FALSE)
test_that("test Vaf filter", {
    expect_equal(1, nrow(vcfs_vafFilter$nonexoni_vcf$variants))
})

vcfs_vafFilter <- applyFilters(  vcfs = vcfs
                                 , assembly = genome_assembly
                                 , design = design_gr
                                 , vaf.cutoff = 0.5
                                 , remove.nonexonic = FALSE)
test_that("test Vaf filter", {
    expect_equal(0, nrow(vcfs_vafFilter$nonexoni_vcf$variants))
})



                                
# This applies a simple  filter 
vcfs_NoSynonymous <- applyFilters(vcfs = vcfs
                     , assembly = genome_assembly
                     , design = design_gr
                     , variantType = c("synonymous")
                     , remove.nonexonic = FALSE 
                     )
test_that("test No Synonymous", {
    expect_equal(length(expected_intersection), nrow(vcfs_NoSynonymous$nonexoni_vcf$variants))
})

vcfs_NoCancer <- applyFilters(vcfs = vcfs
                                  , assembly = genome_assembly
                                  , design = design_gr
                                  , remove.cancer = TRUE
                                  , remove.nonexonic = FALSE 
                            )
test_that("test No Cancer", {
    expect_equal(length(expected_intersection), nrow(vcfs_NoCancer$nonexoni_vcf$variants))
})





#TMB_res=applyTMB(inputForTMB = vcfs_noFilter, assembly = genome_assembly)                                


# SET UP The TESTS ----------------------------------------------------------



# expected region in the panel design and outside of exonic 
#      chr7 31439
                 
# test_that("test No Synonymous", {
#     n_match = vcfs_NoSynonymous$nonexoni_vcf$variants %>%
#         dplyr::filter(CHR == "chr7") %>%
#         dplyr::filter(START == "31439") %>%
#         nrow() 
#     expect_true(n_match > 0, "The nonexonic SNV was removed while using the No Synonymous filter")
# })
# 
# 
# test_that("test No Filter", {
#     n_match = vcfs_noFilter$nonexoni_vcf$variants %>%
#         dplyr::filter(CHR == "chr7") %>%
#         dplyr::filter(START == "31439") %>%
#         nrow() 
#     expect_true(n_match > 0, "The nonexonic SNV was removed while using the No filter")
# })
#                                 

# 
# 
# library(GenomicRanges)
# library(IRanges)
# IRanges::subsetByOverlaps(
#     makeGRangesFromDataFrame(vcfs_NoSynonymous$nonexoni_vcf$variants),
#     vcfs_NoSynonymous$nonexoni_vcf$design) 


# HG38 TEST ================================================================
genome_assembly <- "hg38"