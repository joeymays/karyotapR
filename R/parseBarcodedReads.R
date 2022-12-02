#' Parse Barcoded Reads
#'
#' @param bam.file Chr string indicating file path of BAM file. `.bai` BAM index file must be in the same location.
#' @param barcode.lookup data.frame where the first column is the barcode identifier/name and the second column is the DNA sequence. Headers are ignored.
#' @param cell.barcode.tag Character string of length 2, indicates cell barcode field in BAM, specified by Tapestri pipeline (currently "RG"). Default "RG".
#' @param contig Chr string of contig or chromosome name for search.
#' @param max.mismatch The maximum and minimum number of mismatching letters allowed. See [Biostrings::matchPattern()].
#' @param with.indels If TRUE then indels are allowed. See [Biostrings::matchPattern()].
#'
#' @return A data.frame of read counts for each specified barcode.
#' @export
#'
#' @rdname parseBarcodedReads
#'
#' @examples
#' \dontrun{counts <- parseBarcodedReadsFromContig(bam.file, barcode.lookup, "virus_ref2")}
parseBarcodedReadsFromContig <- function(bam.file, barcode.lookup, contig, cell.barcode.tag = "RG", max.mismatch = 2, with.indels = T){

    # set bam parameters
    which <- GenomicRanges::GRanges(contig, IRanges::IRanges(1,536870912))
    what <- c("qname", "rname", "isize", "seq")
    param <- Rsamtools::ScanBamParam(which = which, what = what, tag = cell.barcode.tag)

    # check for bam index and throw error if not found
    tryCatch(expr = bam <- Rsamtools::scanBam(file = bam.file, param = param),
             error = function(e){

                error.message.match <- grep(pattern = "failed to load BAM index", x = e$message)
                if(!S4Vectors::isEmpty(error.message.match)){
                    e$message <- paste0(e$message, "\nBAM index .bai file needs to be in the same directory as .bam file.",
                                        "\nRun `Rsamtools::indexBam(bam.filename)` to generate BAM index file if it cannot be found.")
                }

                 stop(e$message)
             })

    # filter bam by cell barcode tag
    bam.filter <- data.frame(qname = bam[[1]]$qname, rname = bam[[1]]$rname, tlen = bam[[1]]$isize, seq  = bam[[1]]$seq, tag = bam[[1]]$tag[[cell.barcode.tag]])

    # parse barcode lookup table by column index
    barcode.lut.vector <- barcode.lookup[,2]
    names(barcode.lut.vector) <- barcode.lookup[,1]
    barcode.set <- Biostrings::DNAStringSet(x = barcode.lut.vector, use.names = T)

    # match input barcodes to reads and convert to logical
    sequence.matches <- lapply(barcode.set, function(x){
        as.logical(Biostrings::vcountPattern(pattern = x, subject = bam.filter$seq, max.mismatch = max.mismatch, with.indels = with.indels))
    })

    # match barcoded reads to cell barcodes, count, combine
    n <- names(sequence.matches)
    names(n) <- names(sequence.matches)

    sequence.match.counts <- lapply(n, function(x){
        counts <- as.data.frame(table(bam.filter$tag[sequence.matches[[x]]]))
        colnames(counts) <- c("cell.barcode", x)
        return(counts)
    })

    sequence.match.counts <- Reduce(function(x, y) merge(x, y, by="cell.barcode", all = T), sequence.match.counts)
    rownames(sequence.match.counts) <- sequence.match.counts$cell.barcode
    sequence.match.counts[is.na(sequence.match.counts)] <- 0

    return(sequence.match.counts)
}

#' Parse Barcoded Reads
#'
#' `parseBarcodedReads()` and `parseBarcodedReadsFromContig()` match exogenous DNA barcode sequences from plasmid transduction to their associated cell barcodes and outputs them as a table of counts.
#' `parseBarcodedReads()` is a shortcut for `parseBarcodedReadsFromContig()`, allowing the user to specify 'gRNA' or 'sample.barcode'.
#'
#' @param bam.file Chr string indicating file path of BAM file. `.bai` BAM index file must be in the same location.
#' @param barcode.lookup data.frame where the first column is the barcode identifier/name and the second column is the DNA sequence. Headers are ignored.
#' @param probe Chr string, either "gRNA" or "sample.barcode" to parse counts from grnaCounts or sampleBarcodeCounts alternative experiments, respectively.
#' @param ... Arguments to pass on to `parseBarcodedReadsFromContig()`.
#' @param tapestri.experiment.object TapestriExperiment
#'
#' @return A data.frame of read counts for each specified barcode
#' @export
#'
#' @examples
#' \dontrun{counts <- parseBarcodedReads(tapestri.experiment.object,
#' bam.file, barcode.lookup, "gRNA")}
parseBarcodedReads <- function(tapestri.experiment.object, bam.file, barcode.lookup, probe, ...){

    probe <- tolower(probe)

    if(probe == "grna"){
        contig <- as.character(rowData(altExp(tapestri.experiment.object, "grnaCounts"))[grnaProbe(tapestri.experiment.object),"chr"])
    } else if(probe == "sample.barcode"){
        contig <- as.character(rowData(altExp(tapestri.experiment.object, "sampleBarcodeCounts"))[barcodeProbe(tapestri.experiment.object),"chr"])
    } else {
        stop("Probe not recognized. Try probe = 'grna' or 'sample.barcode'.")
    }

    parseBarcodedReadsFromContig(bam.file, barcode.lookup, contig = contig, ...)

}
