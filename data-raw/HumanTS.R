filename <- "./data-raw/HumanTumorSuppressors_from_TSGene2"
HumanTumorSuppressors <- read.table(file=filename, sep="\t", header=TRUE)
HumanTumorSuppressors <- as.character(as.vector(HumanTumorSuppressors$GeneID))
usethis::use_data(HumanTumorSuppressors, internal=FALSE, compress="gzip")
