#' Retrieve and add chromosome cytobands and chromosome arms to probe metadata of TapestriExperiment object
#'
#' @param tapestri.experiment.object TapestriExperiment object
#' @param genome Character string indicating reference genome to use. Only hg19 is currently supported.
#'
#' @return TapestriExperiment object with rowData updated to include chromosome arms and cytobands
#' @export
#'
#' @examples
#' \dontrun{tapObject <- getCytobands(tapObject, genome = "hg19")}
getCytobands <- function(tapestri.experiment.object, genome = "hg19"){

    genome <- tolower(genome)

    if(genome != "hg19"){
        stop(paste(genome, "not found. Only hg19 is currently supported."))
    }

    if(genome == "hg19"){
        message("Adding cytobands from hg19.")
    }

    #only add chr label for arms to chrs 1-22, X, Y
    chr.vector <- as.character(SingleCellExperiment::rowData(tapestri.experiment.object)$chr)
    chr.vector <- ifelse(chr.vector %in% c(1:22, "X", "Y"), paste0("chr", chr.vector), chr.vector)

    amplicon.gr <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(chr.vector),
                                          ranges = IRanges::IRanges(start = SingleCellExperiment::rowData(tapestri.experiment.object)$start.pos,
                                                                    end = SingleCellExperiment::rowData(tapestri.experiment.object)$end.pos,
                                                                    names = SingleCellExperiment::rowData(tapestri.experiment.object)$probe.id),
                                          strand = S4Vectors::Rle(values = BiocGenerics::strand("*"),
                                                                  lengths = nrow(SingleCellExperiment::rowData(tapestri.experiment.object))))


    overlap.hits <- GenomicRanges::findOverlaps(amplicon.gr, cytoband.hg19.genomicRanges)

    cytoband.matches <- character(length = S4Vectors::queryLength(overlap.hits))
    cytoband.matches[] <- NA
    cytoband.matches[S4Vectors::queryHits(overlap.hits)] <- S4Vectors::mcols(cytoband.hg19.genomicRanges)[,"cytoband"][S4Vectors::subjectHits(overlap.hits)]

    S4Vectors::mcols(amplicon.gr)$cytoband <- cytoband.matches

    chromosome.arms <- ifelse(is.na(amplicon.gr$cytoband), amplicon.gr$cytoband, paste0(S4Vectors::decode(GenomicRanges::seqnames(amplicon.gr)), substr(amplicon.gr$cytoband, 1, 1)))

    S4Vectors::mcols(amplicon.gr)$arm <- chromosome.arms
    S4Vectors::mcols(amplicon.gr)$arm <- factor(chromosome.arms, unique(chromosome.arms))

    row.data <- SingleCellExperiment::rowData(tapestri.experiment.object)
    amplicon.gr.matrix <- as.data.frame(S4Vectors::mcols(amplicon.gr))
    amplicon.gr.matrix$probe.id <- rownames(amplicon.gr.matrix)
    rownames(amplicon.gr.matrix) <- NULL

    amplicon.metadata <- merge(row.data, amplicon.gr.matrix, by = "probe.id", sort = F)

    if(any(amplicon.metadata$probe.id != rowData(tapestri.experiment.object)$probe.id)){
        stop("Something wrong. rowData and new metadata don't line up.")
    }

    SummarizedExperiment::rowData(tapestri.experiment.object)$cytoband <- amplicon.metadata$cytoband
    SummarizedExperiment::rowData(tapestri.experiment.object)$arm <- amplicon.metadata$arm

    return(tapestri.experiment.object)
}
