filename="./data-raw/ExampleWESdesign.bed"
ExampleWESdesign <- read.table(file=filename, sep="\t", header=FALSE)
ExampleWESdesign <- GenomicRanges::makeGRangesFromDataFrame(ExampleWESdesign,
                                                  seqnames.field="V1",
                                                  start.field="V2",
                                                  end.field="V3",
                                                  strand.field="V5",
                                                  keep.extra.columns=TRUE)
GenomeInfoDb::seqlevelsStyle(ExampleWESdesign) <- "UCSC"
keepchr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
             "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
             "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
             "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
GenomeInfoDb::seqlevels(ExampleWESdesign, pruning.mode="coarse") <- keepchr
usethis::use_data(ExampleWESdesign, internal=FALSE, compress="gzip")


filename="./data-raw/ExamplePaneldesign.bed"
ExamplePaneldesign <- read.table(file=filename, sep="\t", header=FALSE)
ExamplePaneldesign <- GenomicRanges::makeGRangesFromDataFrame(ExamplePaneldesign,
                                                    seqnames.field="V1",
                                                     start.field="V2",
                                                     end.field="V3",
                                                     strand.field="V5",
                                                     keep.extra.columns=TRUE)
GenomeInfoDb::seqlevelsStyle(ExamplePaneldesign) <- "UCSC"
keepchr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
             "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
             "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
             "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
GenomeInfoDb::seqlevels(ExamplePaneldesign, pruning.mode="coarse") <- keepchr
usethis::use_data(ExamplePaneldesign, internal=FALSE, compress="gzip")
