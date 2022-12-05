#' Parse Barcoded Reads
#'
#' @param bam.file Chr string indicating file path of BAM file. `.bai` BAM index file must be in the same location.
#' @param barcode.lookup data.frame where the first column is the barcode identifier/name and the second column is the DNA sequence. Headers are ignored.
#' @param cell.barcode.tag Character string of length 2, indicates cell barcode field in BAM, specified by Tapestri pipeline (currently "RG"). Default "RG".
#' @param contig Chr string of contig or chromosome name to search for barcodes in. Can be a vector of more than one contig to expand search space.
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

    # match barcoded reads to cell barcodes and count. Returns NULL if no matches.
    sequence.match.ids <- names(sequence.matches)
    names(sequence.match.ids) <- names(sequence.matches)

    sequence.match.counts <- lapply(sequence.match.ids, function(x){

        bam.matches <- bam.filter$tag[sequence.matches[[x]]]

        if(S4Vectors::isEmpty(bam.matches)){
            message(paste0("No read matches found for ID '", x,"'."))
            return(NULL)
        }

        counts <- as.data.frame(table(bam.matches))
        counts$bam.matches <- as.character(counts$bam.matches)
        colnames(counts) <- c("cell.barcode", x)
        return(counts)
    })

    # stop if no matches
    if(all(sapply(sequence.match.counts, is.null))){
        stop("No matches found for any barcode IDs.")
    }

    # remove NULL items from sequence.match.counts, removes barcode IDs with no matches
    sequence.match.counts <- sequence.match.counts[!sapply(sequence.match.counts, is.null)]

    # combine count tables
    sequence.match.counts <- Reduce(function(x, y) merge(x, y, by="cell.barcode", all = T), sequence.match.counts)
    rownames(sequence.match.counts) <- sequence.match.counts$cell.barcode
    sequence.match.counts[is.na(sequence.match.counts)] <- 0

    return(sequence.match.counts)
}

#' Parse Barcoded Reads
#'
#' `parseBarcodedReads()` and `parseBarcodedReadsFromContig()` match exogenous DNA barcode sequences to their associated cell barcodes and saves them to the colData (cell barcode metadata) of the input TapestriExperiment object.
#' `parseBarcodedReads()` is a shortcut for `parseBarcodedReadsFromContig()`, allowing the user to specify 'gRNA' or 'sample.barcode', and assign the result to the input object..
#'
#' @param bam.file Chr string indicating file path of BAM file. `.bai` BAM index file must be in the same location.
#' @param barcode.lookup data.frame where the first column is the barcode identifier/name and the second column is the DNA sequence. Headers are ignored.
#' @param probe.tag Chr string, either "gRNA" or "sample.barcode" to parse counts from grnaCounts or sampleBarcodeCounts alternative experiments, respectively.
#' @param ... Arguments to pass on to `parseBarcodedReadsFromContig()`.
#' @param TapestriExperiment TapestriExperiment object
#' @param return.table If TRUE, returns table of read counts per barcode. If FALSE, returns TapestriExperiment. Default FALSE.
#'
#' @return An updated TapestriExperiment object with read counts added to the colData slot. If `return.table == TRUE`, a data.frame of read counts for each specified barcode.
#' @export
#'
#' @examples
#' \dontrun{counts <- parseBarcodedReads(TapestriExperiment,
#' bam.file, barcode.lookup, "gRNA")}
parseBarcodedReads <- function(TapestriExperiment, bam.file, barcode.lookup, probe.tag, return.table = F, ...){

    probe.tag <- tolower(probe.tag)

    if(probe.tag == "grna"){
        contig <- as.character(rowData(altExp(TapestriExperiment, "grnaCounts"))[grnaProbe(TapestriExperiment),"chr"])
    } else if(probe.tag == "sample.barcode"){
        contig <- as.character(rowData(altExp(TapestriExperiment, "sampleBarcodeCounts"))[barcodeProbe(TapestriExperiment),"chr"])
    } else {
        stop(paste0("Probe tag '", probe.tag, "' not recognized. Try probe = 'grna' or 'sample.barcode'."))
    }

    result <- parseBarcodedReadsFromContig(bam.file, barcode.lookup, contig = contig, ...)

    if(return.table){
        return(result)
    } else {
        # get and merge colData
        cell.data <- as.data.frame(SingleCellExperiment::colData(TapestriExperiment))
        updated.cell.data <- merge(cell.data, result, all.x = T, sort = F)

        # set NAs to 0
        ids <- setdiff(colnames(result), "cell.barcode")
        updated.cell.data[,ids][is.na(updated.cell.data[,ids])] <- 0
        rownames(updated.cell.data) <- updated.cell.data$cell.barcode

        # update TapestriExperiment
        SummarizedExperiment::colData(TapestriExperiment) <- S4Vectors::DataFrame(updated.cell.data)

        return(TapestriExperiment)
    }
}

#' Call sample labels based on metadata counts.
#'
#' `callSampleLables()` determines sample labels by comparing auxiliary count data, most likely generated from barcoded reads (see[parseBarcodedReads()]).
#' Labels are dictated by whichever `coldata.labels` element has the highest number of counts.
#' By default, ties are broken by choosing whichever label has the lowest index position (`ties.method = "first"`).
#' Samples with 0 counts for all labels are labeled "none".
#'
#' @param TapestriExperiment A TapestriExperiment object.
#' @param coldata.labels A chr vector of column labels corresponding to colData.
#' @param method A chr string indicating call method. Only "max" currently supported, calls based on whichever label has the most counts.
#' @param ties.method A chr string passed to `max.col()` indicating how to break ties. Default "first".
#'
#' @return A chr vectpr of sample labels.
#' @export
#'
#' @examples
#' \dontrun{sample.calls <- determineSampleLables(TapestriExperiment, c("g7", "gNC"))}
callSampleLables <- function(TapestriExperiment, coldata.labels, method = "max", ties.method = "first"){

    if(method != "max"){
        stop("Method not recognized. Only 'max' currently supported.")
    } else {

        # check for missing colData labels
        if(any(!coldata.labels %in% colnames(SingleCellExperiment::colData(TapestriExperiment)))){
            stop(paste0(coldata.labels[!coldata.labels %in% colnames(SingleCellExperiment::colData(TapestriExperiment))], " not found in colData."))
        }

        # subset colData
        coldata.subset <- as.data.frame(SingleCellExperiment::colData(TapestriExperiment)[,coldata.labels])

        # check if numeric
        if(any(!apply(coldata.subset, 2, is.numeric))){
            stop("Selected coldata.labels are not numeric.")
        }

        # make calls
        sample.calls <- coldata.labels[max.col(coldata.subset, ties.method = ties.method)]
        names(sample.calls) <- rownames(coldata.subset)
        sample.calls[rowSums(coldata.subset) == 0] <- "none"

        return(sample.calls)
    }
}
