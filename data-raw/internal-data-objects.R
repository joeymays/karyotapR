#Internal Data

#h19 cytoband coordinates; cytoBand.hg19.txt from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz on 2022-04-22

cytoband.hg19 <- read.table("cytoBand.hg19.txt", sep= '\t', skip = 1, col.names = c("chromosome", "start.position", "end.position", "cytoband", "stain"))
cytoband.hg19$strand <- c("*")
cytoband.hg19.genomicRanges <- GRanges(seqnames = cytoband.hg19$chromosome,
                                       ranges = IRanges(cytoband.hg19$start.position, cytoband.hg19$end.position),
                                       strand = cytoband.hg19$strand,
                                       cytoband = cytoband.hg19$cytoband)

#usethis::use_data(cytoband.hg19.genomicRanges, overwrite = TRUE, internal = TRUE)

#CO293 panel metadata
#tap.file is arbitrary .h5 tapestri pipeline file using CO293 panel
tap <- createTapestriExperiment(tap.file, panel.id = "CO293", get.cytobands = F, move.non.genome.probes = F, filter.variants = T)
co293.metadata <- as.data.frame(rowData(tap))
co293.metadata <- co293.metadata[,1:4]

usethis::use_data(cytoband.hg19.genomicRanges, co293.metadata, overwrite = TRUE, internal = TRUE)
