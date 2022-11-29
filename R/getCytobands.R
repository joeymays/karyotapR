# Find overlaps between cytobands and amplicon coverage

getCytobands <- function(TapestriExperimentObject){

    message("Adding cytobands from hg19.")

    amplicon.gr <- GRanges(seqnames = Rle(paste0("chr",SingleCellExperiment::rowData(TapestriExperimentObject)$chr)),
                           ranges = IRanges(start = SingleCellExperiment::rowData(TapestriExperimentObject)$start.pos,
                                            end = SingleCellExperiment::rowData(TapestriExperimentObject)$end.pos,
                                            names = SingleCellExperiment::rowData(TapestriExperimentObject)$probe.id),
                           strand = Rle(values = strand("*"),
                                        lengths = nrow(SingleCellExperiment::rowData(TapestriExperimentObject))))


overlap.hits <- findOverlaps(amplicon.gr, cytoband.hg19.genomicRanges)

# Add cytobands as metadata for matches
cytoband.matches <- character(length = queryLength(overlap.hits))
cytoband.matches[] <- NA
cytoband.matches[queryHits(overlap.hits)] <- mcols(cytoband.hg19.genomicRanges)[,"cytoband"][subjectHits(overlap.hits)]


chromosome.arms <- paste0(decode(seqnames(amplicon.gr)), substr(amplicon.gr$cytoband, 1, 1))


mcols(amplicon.gr)$cytoband <- cytoband.matches
mcols(amplicon.gr)$arm <- paste0(decode(seqnames(amplicon.gr)), substr(amplicon.gr$cytoband, 1, 1))


#make sure that amplicons are sorted correctly when applying back to rowData
#ROW DATA STARTED Out OF ORDER CHECK THAT ****

row.data <- SingleCellExperiment::rowData(TapestriExperimentObject)
amplicon.gr.matrix <- as.data.frame(mcols(amplicon.gr)) %>% tibble::rownames_to_column("probe.id")

amplicon.metadata <- merge(row.data, amplicon.gr.matrix, by = "probe.id", sort = F)

amplicon.metadata$chr <- as.factor(amplicon.metadata$chr)



amplicon.metadata$chr <- factor(amplicon.metadata$chr, levels = c(1:22, "X", "Y", "virus_ref", "virus_ref2"))
amplicon.metadata$arm <- factor(amplicon.metadata$arm, levels = gtools::mixedsort(unique(amplicon.metadata$arm)))




}
