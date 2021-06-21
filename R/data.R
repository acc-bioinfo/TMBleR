#' Cosmic coding mutations hg19
#'
#' A dataset containing coding mutations from COSMIC for genome assembly hg19. 
#' Since the COSMIC dataset is very large and it requires registration we are
#' bundleing in the package a very small version, just enough to demo the
#' package vignette and test suite. Please refer on the documentation 
#' vignette/HOWTO_Import_external_data.html on how to retrieve and 
#' use COSMIC data with TMBleR

#' data("COSMIC_hg19_demo") will load the demo dataset as a "COSMIC_hg19
#'  object in the global env
#'
#' @format A data frame with 4740646 rows and 6 variables
#' \describe{
#'  \item{CHR}{chromosome name, in UCSC format chrN}
#'  \item{START}{start genomic coordinate, UCSC format}
#'  \item{END}{start genomic coordinate, UCSC format}
#'  \item{STRAND}{strand}
#'  \item{REF}{reference nucleotide}
#'  \item{ALT}{alternative nucleotide}
#' }
#' @source From CosmicCodingMuts.vcf file for hg19 downloaded from
#' https://cancer.sanger.ac.uk/cosmic/download and saved in data-raw as
#' CosmicCodingMuts_hg19.vcf. The dataset was then subset through random 
#' stratified sampling to max 10.000 snv per chromosome.
"COSMIC_hg19_demo"





#' Cosmic coding mutations hg38
#'
#' A dataset containing coding mutations from COSMIC for genome assembly hg38
#' Since the COSMIC dataset is very large and it requires registration we are
#' bundleing in the package a very small version, just enough to demo the
#' package vignette and test suite. Please refer on the documentation 
#' vignette/HOWTO_Import_external_data.html on how to retrieve and 
#' use COSMIC data with TMBleR.
#' 
#' data("COSMIC_hg38_demo") will load the demo dataset as a "COSMIC_hg38
#'  object in the global env
#'
#' @format A data frame with 4668373 rows and 6 variables
#' \describe{
#'  \item{CHR}{chromosome name, in UCSC format chrN}
#'  \item{START}{start genomic coordinate, UCSC format}
#'  \item{END}{start genomic coordinate, UCSC format}
#'  \item{STRAND}{strand}
#'  \item{REF}{reference nucleotide}
#'  \item{ALT}{alternative nucleotide}
#' }
#' @source From CosmicCodingMuts.vcf file for hg38 downloaded from
#' https://cancer.sanger.ac.uk/cosmic/download and saved in data-raw as
#' CosmicCodingMuts_hg38.vcf . The dataset was then subset through random 
#' stratified sampling to max 10.000 snv per chromosome.
"COSMIC_hg38_demo"





#' List of human tumor suppressors
#'
#' List of gene ids known to be human tumor suppressors. It is used as default
#' tsList argument in applyFilters()  if not specified otherwise and the
#' removeCancer argument is set to TRUE. In such a sceneario, truncating mutations in
#' these human tumor suppressor genes are not counted for TMB quantification.
#' This is a corretive strategy usually required for panel-based TMB 
#' quantification to avoid upward skewing bias due to panel enrichment in cancer 
#' genes.
#'
#' @format character vector with 1217 elements
#' @source From the TSGene 2.0 database available at
#' https://bioinfo.uth.edu/TSGene/. See Zhao et al. "TSGene 2.0: an
#' updated literature-based knowledgebase for tumor suppressor genes" NAR 2016.
#' PMID: 26590405
"HumanTumorSuppressors"




#' Dataset containing sequencing design of a WES experiment
#'
#' Sequencing design of the WES experiments used as examples.
#'
#'
#'@format A GRanges object with 199202 ranges and 0 metadata columns
#' \describe{
#'   \item{seqnames}{Chromosome name}
#'   \item{ranges}{Start-End genomic positions}
#'   \item{strand}{'+', '-' or '*' for undefined}
#' }
#' @source ss_v6_merged+exome.bed file containing design of an example WES
"ExampleWESdesign"




#' Dataset containing sequencing design of a panel
#'
#' Sequencing design of the panel used as examples.
#'
#'
#' @format A GRanges object with 8320 ranges and 2 metadata columns
#' \describe{
#'   \item{seqnames}{Chromosome name}
#'   \item{ranges}{Start-End genomic positions}
#'   \item{strand}{'+', '-' or '*' for undefined}
#'   \item{V4}{hgnc gene symbol and start-end cooridnates}
#'   \item{V6}{hgnc gene symbol and Pool (1 or 2)}
#' }
#' @source In house panel design, provided as example.
#' design_file_path = system.file( "data-raw"
#'     , "ExamplePaneldesign.bed"
#'     , package = "TMBleR"
#'     , mustWork = TRUE)
#'     
#' design <- readDesign(
#'     filename = design_file_path
#'     , assembly = "hg19
#'     , ids = "entrezgene_id")
"ExamplePaneldesign"




#' Dataset containing example vcf from WES
#'
#' Mutations in VCF format from WES sequencing of 4 samples
#'
#'
#' @format A list of four elements (corresponding to 4 samples). Each element is a CollapsedVCF object.
#' \describe{
#'   \item{fixed}{dataframe with REF, ALT, QUAL, FILTER}
#'   \item{info}{dataframe with data in INFO field of VCF}
#'   \item{rowRanges}{GRanges object with chromosome, start and end coordinates, strand}
#'   \item{colData}{dataframe with sample name}
#'   \item{assays}{Reference class object of class "ShallowSimpleListAssays", describing fields in INFO of VCF}
#'   \item{elementMetadata}{dataframe}
#'   \item{NAMES}{NULL}
#'   \item{metadata}{dataframe with metadata of VCF}
#' }
#' @source in house examples of WES vcfs.
"ExampleWESvcfs"



#' Dataset containing tumor mutational burden and clinical data from a published
#' study
#'
#' Dataset from Hellman et al., Cancer Cell 2018 containing the original
#' WES-based TMB value and clinical response from the study and a simulated
#' panel-based TMB value calculated in silico by subsetting the WES dataset to
#' only contain genes targeted by a custom panel used as example.
#'
#'
#' @format A data frame with 74 rows and 4 variables
#' \describe{
#'   \item{PatientID}{Patient identifier}
#'   \item{ClinicalResponse}{Clinical response to immunotherapy; responder or
#'   nonresponder}
#'   \item{Panel.NumMuts}{Panel-based tumor mutational burden, in total number of
#'    mutations}
#'   \item{WES.NumMuts}{WES-based tumor mutational burden, in total number of
#'   mutations}
#'   }
#' @source provided as supplementary table by Hellman et al. in "Tumor Mutational
#'  Burden and Efficacy of Nivolumab Monotherapy and in Combination With Ipilimumab
#'  in Small-Cell Lung Cancer". Cancer Cell 2018 (PMID: 29731394) and modified to
#'  add simulated panel TMB
"Hellman_SimulatedFM1Panel_WES"



#' Dataset containing tumor mutational burden and clinical data from a published
#' study 
#'
#' Dataset from Van Allen EM  et al., Science 2015 containing the detailed 
#' clinical and genome characteristics of the indivitual patients enrolled 
#' into the study. The datase includes wherether the patients were classified
#' as either responder or non-responder. 
#'
#' @format A data frame with 110 rows and 2 variables
#' \describe{
#'   \item{Sample}{Patient identifier}
#'   \item{ClinicalResponse}{Clinical response to immunotherapy; responder or
#'   nonresponder}
#'   }
#' @source provided as supplementary file S2 by Van Allen EM et al. in 
#' "Genomic correlates of response to CTLA-4 blockade in metastatic melanoma"
#'  Science 2015 (PMID: 26359337). Only clinical response and patient variable
#'  are saved here
#'  
#'  # load supplementary file from 
#'  VanAllen_Clinical <- read.table("TableS2_Revised.xlsx")
#'  
#'  VanAllen_Clinical %<>%              
#'    dplyr::rename("Sample" = patient) %>%  # rename patient column into Sample
#'    dplyr::mutate(ClinicalResponse = ifelse(VanAllen_Clinical$group == "response"
#'                                            , "responder", "nonresponder")) %>%
#'    dplyr::select(Sample, ClinicalResponse)
"VanAllen_Clinical"



#' Input to the plotTMB function to be used in tests
#'
#' Table containing estimated tumor mutational burden and metadata (sample,
#' filter, size of the sequenced genomic space) calculated using the Horizon 
#' demo dataset
#'
#'
#' @format A data frame with 16 rows and 5 variables
#' \describe{
#'   \item{Sample}{factor, sample identifier}
#'   \item{Filter}{factor, type of filter, if any, applied on input mutations}
#'   \item{Sequencing_Size}{numeric vector, genomic size targeted by the
#'   sequencing experiment (WES or Panel), in Mb}
#'   \item{Tot_Number_Mutations}{numeric vector, tumor mutational burden in the
#'   targeted genomic region}
#'   \item{TMB_per_Mb}{numeric vector, global tumor mutational burden expressed
#'   as number of mutations}
#' }
#' @source generated in house as an example using a custom Panel and toy
#' examples of VCF files from Horizon cell line dataset
"TMB_Horizon"



#' Input to the plotTMB function to be used in tests
#'
#' Table containing estimated tumor mutational burden and metadata (sample,
#' filter, size of the sequenced genomic space) calculated using 4 patients from 
#' the VanAllen et al. dataset.
#'
#'
#' @format A data frame with 16 rows and 5 variables
#' \describe{
#'   \item{Sample}{factor, sample identifier}
#'   \item{Filter}{factor, type of filter, if any, applied on input mutations}
#'   \item{Sequencing_Size}{numeric vector, genomic size targeted by the
#'   sequencing experiment (WES or Panel), in Mb}
#'   \item{Design}{Name of the panel experiment}
#'   \item{Tot_Number_Mutations}{numeric vector, tumor mutational burden in the
#'   targeted genomic region}
#'   \item{TMB_per_Mb}{numeric vector, global tumor mutational burden expressed
#'   as number of mutations}
#' }
#' @source Dataset from Van Allen EM  et al., Science 2015
"TMB_VanAllen"


#' Panel-based TMB to use as input to the correlateTMBvalues and plotCorrelation
#' functions in tests
#'
#' Table containing estimated tumor mutational burden and metadata (smaple,
#' filter, size of the sequenced genomic space) obtained by simulating in silico
#' the TMB that would be estimated using a certain panel design.
#'
#'
#' @format A data frame with 4 rows and 5 variables
#'  \describe{
#'   \item{Sample}{factor, sample identifier}
#'   \item{Filter}{factor, type of filter, if any, applied on input mutations}
#'   \item{Sequencing_Size}{numeric vector, genomic size targeted by the
#'   sequencing experiment (WES or Panel), in Mb}
#'   \item{Tot_Number_Mutations}{numeric vector, tumor mutational burden in the
#'   targeted genomic region}
#'   \item{TMB_per_Mb}{numeric vector, global tumor mutational burden expressed
#'   as number of mutations}
#' }
#' @source generated in house as an example of TMB estimated obtained by
#' simulating in silico the TMB that would be estimated using a certain panel
#' design.
"TMBs_SimulatedPanel"




#' WES-based TMB to use as input to the correlateTMBvalues and plotCorrelation
#' functions in tests
#'
#' Table comntaining estimated tumor mutational burden and metadata (smaple,
#' filter, size of the sequenced genomic space) obtained by WES.
#'
#'
#' @format A data frame with 4 rows and 5 variables
#'  \describe{
#'   \item{Sample}{factor, sample identifier}
#'   \item{Filter}{factor, type of filter, if any, applied on input mutations}
#'   \item{Sequencing_Size}{numeric vector, genomic size targeted by the
#'   sequencing experiment (WES or Panel), in Mb}
#'   \item{Tot_Number_Mutations}{numeric vector, tumor mutational burden in the
#'   targeted genomic region}
#'   \item{TMB_per_Mb}{numeric vector, global tumor mutational burden expressed
#'   as number of mutations}
#' }
#' @source generated in house as an example of TMB estimated obtained by WES.
"TMBs_WES"





#' Data to use as input to the applyTMB functions in tests
#'
#' List of lists, , generated by applyFilters or applyInputToTMB, containing
#' required input to the applyTMB function. Each element of the list corresponds
#' to a sample and is a list with the following elements: sample, filter,
#' variants and design.
#'
#'
#' @format A data frame with 4 rows and 5 variables
#'  \describe{
#'   \item{sample}{character, sample identifier}
#'   \item{Filter}{character, type of filter, if any, applied on input mutations}
#'   \item{design}{GRanges object, regions targeted by the sequencing experiment
#'   (WES or Panel)}
#'   \item{variants}{data.frame, describing mutations by chromosome, start and
#'   end coordinates, REF and ALT alleles}
#' }
#' @source generated in house as an example of TMB estimated obtained by WES.
"vcfs_all"





#' Data to use as input to the applySimulatedPanel function in tests
#'
#' List of lists, generated by applyFilters or applyInputToTMB, containing
#' required input to the applySimulatedPanel function. Each element of the list
#' corresponds to a sample and is a list with the following elements: sample,
#' filter, variants and design.
#'
#'
#' @format A data frame with 4 rows and 5 variables
#'  \describe{
#'   \item{sample}{character, sample identifier}
#'   \item{Filter}{character, type of filter, if any, applied on input mutations}
#'   \item{design}{GRanges object, regions targeted by the sequencing experiment
#'   (WES or Panel)}
#'   \item{variants}{data.frame, describing mutations by chromosome, start and
#'   end coordinates, REF and ALT alleles}
#' }
#' @source generated in house as an example of TMB estimated obtained by WES.
"vcfs_NoCancer_ForPanel"









