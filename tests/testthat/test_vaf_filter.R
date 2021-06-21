## THIS TESTING UNIT IS NOT ASSOCIATED WITH ANY R script
## It's objective is to expensively test the VAF filter option
## 
## 
## IN order to make sure the testing of the filters is clean and well organized 
## we will create one vcf file that contains all the SNV that we need to test. 
## The vcf will have a INFO column with a "filter" tag, that will specify what type
##  of filter we are expecting to apply. The list of possible filters is the following:
##  
##      vaf_cutoff = snv filtered out by
##      off_target = regions outside of the panel
##      non_coding = snv mapped outside of coding regions (but inside of the panel design)
##      variant_type = either somatic, frameshift, nonsense, etc..
##      
##  The filter tag will be specified in either the INFO or FORMAT column
##  
##  
##        ##INFO=<ID=filter,Number=1,Type=String,Description="Expected TMBleR Variant filter">									
##        
##        or
##
##        ##FORMAT=<ID=filter,Number=1,Type=String,Description="Expected TMBleR Variant filter">				
##        
##        
##  Note 1: see issue #74, we need to decide whether to use for the AF field 
##  either the AF in the INFO or in the FORMAT field. Ideally, the validator, 
##  if it finds the AF only in the INFO field, and not in the FORMAT, it could
##  copy the AF into the FORMAT field as well (by modifying the header too) and 
##  raise a warning. My preference is the F
## 
## #############################################################################
context("VAF_Filter.R")


# Tests
################################################################################

# PREPARE VCF FILE ----------------------------
# Read vcf filenames
#vcf_files <- list(test_filter="Pat02_vs_NORMAL.vcf")
vcf_files <- list(test_filter="test_filters.vcf")
# For each vcf file, get the absolute path
vcf_files <- lapply(vcf_files
                    , function(x) system.file("extdata", x, package = "TMBleR", mustWork = TRUE))
# Read in the files
vcfs <- readVcfFiles(vcf_files, assembly = "hg19")


# PREPARE DESIGN ----------------------------
#Read in input the panel sequencing design
# designPanel <- readDesign(system.file("extdata"
#                                       , "ExamplePanel_GeneIDs_Only10GenesForTest.txt"
#                                       , package = "TMBleR"
#                                       , mustWork = TRUE)
#                           , assembly = "hg19"
#                           , ids = "entrezgene_id")
# Read bed file with the panel sequencing design 
design_df <- data.frame(seqnames = "chr7"
                        , start = "1"
                        , end = "100000000")
design_gr <- GenomicRanges::makeGRangesFromDataFrame(design_df)
design_gr@metadata$design_name <- "design_name"


###### RUN THE TESTS ###########################################################

# TEST INPUT FORMAT
vcfs_out <- applyFilters(vcfs, 
                         assembly="hg19",
                         design=design_gr, 
                         remove.nonexonic = TRUE,
                         remove.cancer=TRUE, 
                         #tsList=NULL, 
                         #variantType = c("synonymous")
                         
)


# check for insertions
idx_insertion <- which(vcfs$test_filter@assays@data@listData$TAG == "insertion")
