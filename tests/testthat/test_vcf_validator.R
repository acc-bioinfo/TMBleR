# Test suite for the functiong .vcf_validator()
# 
#  VCF validator is private function design to help troubleshoot problems encountered
#  in vcf files in input. It is not designed to validate the entire vcf file but only
#  those parts that are relevant to TMBleR. As a conseguence it will raise
#  warnings and errors based on vcf file format
#  
# ##############################################################################
context("vcf_validator.R")

# Preparation 
################################################################################

# PREPARE VCF FILE ----------------------------
# Read vcf filenames
vcf_files <- list(test_validator="test_vcf_validator.vcf")
# For each vcf file, get the absolute path
vcf_files <- lapply(vcf_files
                    , function(x) system.file("extdata", x, package = "TMBleR", mustWork = TRUE))

# this will apply the validation
validation_out <- .validate_vcf(vcf_files$test_validator)

#vcfs <-  readVcfFiles(vcf_files, assembly = "hg19")


test_that("Check warnings", {
  # CHECK AF DESCRIBED IN THE HEADER BUT NOT IN THE INFO FIELD
  expect_true(any(grepl("line 6: AF is described in the header but it is not defined in the INFO field of this and possibly other variants", validation_out$warnings)))

  # CHECK AF DESCRIBED IN BOTH INFO AND FORMAT FIELDS
  expect_true(any(grepl("line 8: AF is described both in FORMAT and in INFO of this and possibly other variants. FORMAT's AF will be used", validation_out$warnings)))
})

test_that("Check errors", {
  # CHECK CORRECT FIRST LINE
  expect_true(any(grepl("line 1: The vcf file must start with ##fileformat", validation_out$errors, fixed = TRUE)))
  
  # CHECK CORRECT ALT ALLELE
  expect_true(any(grepl("line 6: ALT must be a comma separated list of alternate alleles (A,C,G,T,N,.,*); '' found", validation_out$errors, fixed = TRUE)))
  
  # CHECK CORRECT ALT ALLELE (2)
  expect_true(any(grepl("line 6: the ALT allele is not provided for this and possibly other variants. If this is a deletion, it should be represented according to the 1000 Genome project as shown in the link: https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40]", validation_out$errors, fixed = TRUE)))
  
  # CHECK CORRECT REF ALLELE
  expect_true(any(grepl("line 7: REF must be one of A,C,G,T,N,. (case insensitive); '' found", validation_out$errors, fixed = TRUE)))
  
  # CHECK CORRECT REF ALLELE (2)
  expect_true(any(grepl("line 7: the REF allele is not provided for this and possibly other variants. If this is an insertion, it should be represented according to the 1000 Genome project as shown in the link: https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40]", validation_out$errors, fixed = TRUE)))
})