# Note: This test is not directly connected to any other script.R 
# 
# The objective of this test is to make sure that we are testing the filtering 
# of SNV mapped outside of the design. 
# 

context("test_offtarget_filter.R")


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
                        , start = "299840"
                        , end = "299860")
design_gr <- GenomicRanges::makeGRangesFromDataFrame(design_df)
design_gr@metadata$design_name <- "design_name"

# Check expected intersection
expected_intersection <- IRanges::subsetByOverlaps(
    vcfs$nonexoni_vcf@rowRanges, design_gr) 



# This applies no filter - Only SNV not mapped inside the design should be removed
vcfs_noFilter <- applyFilters(  vcfs = vcfs
                                , assembly = genome_assembly
                                , design = design_gr
                                , remove.nonexonic = FALSE)

test_that("test No Filter", {
    expect_equal(length(expected_intersection), nrow(vcfs_noFilter$nonexoni_vcf$variants))
})

# This applies a simple  filter 
vcfs_NoSynonymous <- applyFilters(vcfs = vcfs
                                  , assembly = genome_assembly
                                  , design = design_gr
                                  , variantType = c("synonymous")
)
test_that("test No Synomymous", {
    # This gets filtered
    expect_equal(nrow(vcfs_NoSynonymous$nonexoni_vcf$variants), 0)
})


vcfs_Noframeshift <- applyFilters(vcfs = vcfs
                                  , assembly = genome_assembly
                                  , design = design_gr
                                  , variantType = c("frameshift")
)
test_that("test No frameshift", {
    # This gets filtered
    expect_equal(nrow(vcfs_Noframeshift$nonexoni_vcf$variants), length(expected_intersection))
})



# 
# vcfs_VAF <- applyFilters(vcfs = vcfs
#                           , assembly = genome_assembly
#                           , design = design_gr
#                           , vaf.cutoff = 0.1
# )
# 


# HG38 TEST ================================================================
genome_assembly <- "hg38"