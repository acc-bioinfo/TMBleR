#' Read design file
#'
#' Read the file describing the design of the sequencing experiment (i.e. the
#' targeted regions in the sequencing experiment). It can be either a BED file
#' containing genomic coordinates of the targeted genes or a TXT file containing
#' symbols for targeted genes.
#'
#' @param filename name of file describing the sequencing design. It requires
#' ".bed" extension for files containing the genomic coordinates of the targeted
#' regions or ".txt" extension for files containing the gene symbol of targeted
#' genes. BED files must be tab-delimited, with no track description line and it
#' must report at least CHR START and END information.
#' @param assembly human genome assembly: hg19 or hg38
#' @param ids type of gene identifiers, if .txt file with gene list is provided:
#' hgnc_symbol or entrezgene_id.
#' @param name string indicating the name of the design. if not specified, the filename
#' without extension will be used instead. 
#' 
#' @return \code{data.frame} with genomic coordinates or a character vector with
#' symbols of targeted genes
#'
#' @examples
#' #Read in input the panel sequencing design (either a BED file with genomic
#' #coordinates of the targeted genes or TXT file with targeted gene names).
#' design_FromTxt <- readDesign(filename = system.file("extdata",
#'                                  "ExamplePanel_GeneIDs.txt",
#'                                  package = "TMBleR",
#'                                  mustWork = TRUE),
#'                      assembly = "hg19",
#'                      ids = "entrezgene_id")
#'
#' design_FromBed <- readDesign(filename = system.file("extdata",
#'                                 "ExamplePaneldesign.bed",
#'                                  package = "TMBleR",
#'                                  mustWork = TRUE),
#'                               assembly = "hg19")
#'
#'
#' @author Laura Fancello
#' @export
readDesign <- function(filename, assembly, ids, name= NULL){
    
    # Sanity Checks  -----------------------------------------------------------
    if(is.null(filename)){
        stop("argument 'filename' is missing, with no default")
        }
    if(!(file.exists(filename))){
        stop("file(s) do not exist")
    }
    if(is.null(assembly)){
        stop("argument 'assembly' is missing, with no default: please specify 'hg19' or 'hg38'")
    }
    if ((assembly != "hg19") && (assembly != "hg38")) {
        stop("No valid genome specified: please specify 'hg19' or 'hg38'")}
    if(!(grepl("\\.txt$", filename)) && !(grepl("\\.bed$", filename))){
        stop("No valid file name for design: please provide file with .txt 
             or .bed extension")
    }
  
    if(is.null(name)){
      name = basename(filename) %>%  # select only filename
        strsplit(., ".txt") %>% .[[1]] %>% .[1] %>% # remove .txt expension
        strsplit(., ".bed") %>% .[[1]] %>% .[1]     # remove .bed expension
    }
  
  
    # Read and preprocess input  -----------------------------------------------
    # If design is provided in a BED file, read and put it in a GRanges object
    if(grepl("\\.bed$", filename)){
        design <- utils::read.table(filename, header=FALSE, sep="\t")
        if(!(all(grepl("^chr", as.character(as.vector(design[,1])))))){
            stop("Invalid BED file. Please check that all values in column one start with the prefix 'chr'") 
        }
        if(dim(design)[2]<3){
            stop("Invalid BED file: less than 3 columns. Please provide a valid BED file with at least chromosome, start and end columns and check that columns are tab-delimited.")
        }
        design <- design[, c(1:3)]
        colnames(design) <- c("chr", "start", "end")
        design <- GenomicRanges::makeGRangesFromDataFrame(design, seqnames.field = "chr", start.field = "start", end.field = "end", keep.extra.columns = FALSE)
    }
    if(grepl("\\.txt$", filename)){
        # If design is a provided as vector of gene symbols, retrieve their
        # genomic coordinates and put them in a GRanges object

        if(is.null(ids)){
            stop("argument 'ids' is missing, with no default.")
        }
        if ((ids != "entrezgene_id") && (ids != "hgnc_symbol")) {
            stop("No valid argument 'ids': please specify 'entrezgene_id' or 'hgnc_symbol'")
        }
        
        design <- scan(file=filename, what="character")
        if(assembly == "hg19"){
            ensembl <- biomaRt::useEnsembl(biomart="ensembl",
                                           dataset="hsapiens_gene_ensembl",
                                           GRCh=37)
        }
        if(assembly == "hg38"){
            ensembl <- biomaRt::useEnsembl(biomart="ensembl",
                                           dataset="hsapiens_gene_ensembl")
        }
        if ((assembly != "hg19") && (assembly != "hg38")) {
            stop("No valid genome assembly: please specify 'hg19' or 'hg38'")
        }

    # Download exon coordinates for input genes (assuming that in gene panel
    # design only exonic regions and not the whole gene are usually sequenced)
    
    attributes <- c('chromosome_name',
                    'exon_chrom_start',
                    'exon_chrom_end',
                    'strand',
                    'ensembl_exon_id')
    design_coords <- biomaRt::getBM(attributes=attributes,
                                    filters = ids,
                                    values=design,
                                    mart = ensembl)
    design_coords[design_coords$strand == "-1", ]$strand <- "-"
    design_coords[design_coords$strand == "1", ]$strand <- "+"
    design <- GenomicRanges::makeGRangesFromDataFrame(design_coords,
                                                    seqnames.field="chromosome_name",
                                                    start.field="exon_chrom_start",
                                                    end.field="exon_chrom_end",
                                                    strand.field="strand",
                                                    keep.extra.columns=TRUE)
  
    }
    ## Set UCSC style and allowed chromosome names.
    GenomeInfoDb::seqlevelsStyle(design) <- "UCSC"
    keepchr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
                 "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                 "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                 "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
    GenomeInfoDb::seqlevels(design, pruning.mode="coarse") <- keepchr

    # add name of the design as a metadata
    design@metadata$design_name = name
    
    return(design)
}
