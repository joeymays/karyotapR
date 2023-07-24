#' Parse Barcoded Reads
#'
#' @param bam.file File path of BAM file. `.bai` BAM index file must be in the same location (can be generated using [`Rsamtools::indexBam`]).
#' @param barcode.lookup `data.frame` where the first column is the barcode identifier/name and the second column is the DNA sequence. Headers are ignored.
#' @param cell.barcode.tag Character of length 2, indicates cell barcode field in BAM, specified by Tapestri pipeline (currently "RG"). Default "RG".
#' @param contig Character, contig or chromosome name to search for barcodes in. Can be a vector of more than one contig to expand search space.
#' @param max.mismatch Numeric, the maximum and minimum number of mismatching letters allowed. Default 2.
#' @param with.indels If `TRUE`, then indels are allowed. Default `FALSE`.
#'
#' @return A data.frame of read counts for each specified barcode.
#' @export
#'
#' @rdname countBarcodedReads
#' @order 2
#'
#' @seealso [Biostrings::matchPattern()]
#'
#' @examples
#' \dontrun{counts <- countBarcodedReadsFromContig(bam.file, barcode.lookup, "virus_ref2")}
countBarcodedReadsFromContig <- function(bam.file, barcode.lookup, contig, cell.barcode.tag = "RG", max.mismatch = 2, with.indels = FALSE){

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
                 cli::cli_abort(e$message)
             })

    # filter bam by cell barcode tag
    bam.filter <- data.frame(qname = bam[[1]]$qname, rname = bam[[1]]$rname, tlen = bam[[1]]$isize, seq  = bam[[1]]$seq, tag = bam[[1]]$tag[[cell.barcode.tag]])

    # parse barcode lookup table by column index
    barcode.lut.vector <- barcode.lookup[,2]
    names(barcode.lut.vector) <- barcode.lookup[,1]
    barcode.set <- Biostrings::DNAStringSet(x = barcode.lut.vector, use.names = TRUE)

    # match input barcodes to reads and convert to logical
    # progress bar
    sequence.matches <- purrr::map(barcode.set, function(x){
        as.logical(Biostrings::vcountPattern(pattern = x, subject = bam.filter$seq, max.mismatch = max.mismatch, with.indels = with.indels))
    }, .progress = "Matching barcodes...")

    # match barcoded reads to cell barcodes and count. Returns NULL if no matches.
    sequence.match.ids <- names(sequence.matches)
    names(sequence.match.ids) <- names(sequence.matches)

    sequence.match.counts <- lapply(sequence.match.ids, function(x){

        bam.matches <- bam.filter$tag[sequence.matches[[x]]]

        if(S4Vectors::isEmpty(bam.matches)){
            cli::cli_alert_info("No read matches found for ID {.q {x}}.")
            return(NULL)
        }

        counts <- as.data.frame(table(bam.matches))
        counts$bam.matches <- as.character(counts$bam.matches)
        colnames(counts) <- c("cell.barcode", x)
        cli::cli_alert_info("{sum(counts[,2])} read matches found for ID {.q {x}}.")
        return(counts)
    })

    # stop if no matches
    if(all(vapply(sequence.match.counts, is.null, FUN.VALUE = logical(1)))){
        cli::cli_abort("No matches found for any barcode IDs.")
    }

    # remove NULL items from sequence.match.counts, removes barcode IDs with no matches
    sequence.match.counts <- sequence.match.counts[!vapply(sequence.match.counts, is.null, logical(1))]

    # combine count tables
    sequence.match.counts <- Reduce(function(x, y) merge(x, y, by="cell.barcode", all = TRUE), sequence.match.counts)
    rownames(sequence.match.counts) <- sequence.match.counts$cell.barcode
    sequence.match.counts[is.na(sequence.match.counts)] <- 0

    return(sequence.match.counts)
}

#' @name countBarcodedReads
#'
#' @title Get read counts from barcoded reads
#'
#' @description `countBarcodedReads()` and `countBarcodedReadsFromContig()` match exogenous DNA barcode sequences to their associated
#' cell barcodes and saves them to the `colData` (cell barcode metadata) of `TapestriExperiment`.
#' `countBarcodedReads()` is a shortcut for `countBarcodedReadsFromContig()`, allowing the user to specify 'gRNA' or 'barcode'
#' to use the `grnaCounts` or `barcodeCounts` `altExp` slots.
#' The entries in the `barcode.lookup` table do not have to be present in the sample,
#' allowing users to keep one master table/file of available barcode sequences for use in all experiments.
#'
#' @param bam.file File path of BAM file. `.bai` BAM index file must be in the same location (can be generated using [`Rsamtools::indexBam`]).
#' @param barcode.lookup `data.frame`, first column is the barcode identifier/name and the second column is the DNA sequence. Headers are ignored.
#' @param probe Character, either "gRNA" or "barcode" to parse counts from `grnaCounts` or `barcodeCounts` `altExp` slots, respectively.
#' @param ... Arguments to pass on to `countBarcodedReadsFromContig()`.
#' @param TapestriExperiment `TapestriExperiment` object
#' @param return.table Logical, if `TRUE`, returns table of read counts per barcode. If `FALSE`, returns `TapestriExperiment.` Default `FALSE`.
#' @param max.mismatch Numeric, the maximum and minimum number of mismatching letters allowed. Default 2.
#' @param with.indels If `TRUE`, then indels are allowed. Default `FALSE`.
#'
#' @return `TapestriExperiment` with barcoded read counts added to `colData`.
#' @export
#'
#' @concept barcoded reads
#'
#' @rdname countBarcodedReads
#' @order 1
#'
#' @seealso [Rsamtools::indexBam()]
#'
#' @examples
#' \dontrun{counts <- countBarcodedReads(TapestriExperiment,
#' bam.file, barcode.lookup, "gRNA")}
countBarcodedReads <- function(TapestriExperiment, bam.file, barcode.lookup, probe, return.table = FALSE, max.mismatch = 2, with.indels = FALSE, ...){

    probe <- tolower(probe)

    # get contig
    if(probe == "grna"){
        contig <- as.character(rowData(altExp(TapestriExperiment, "grnaCounts"))[grnaProbe(TapestriExperiment),"chr"])
    } else if(probe == "barcode"){
        contig <- as.character(rowData(altExp(TapestriExperiment, "barcodeCounts"))[barcodeProbe(TapestriExperiment),"chr"])
    } else {
        cli::cli_abort("{.var probe} {.q {probe}} not recognized. Try {.var probe} = {.q grna} or {.q barcode.}")
    }

    result <- countBarcodedReadsFromContig(bam.file, barcode.lookup, contig = contig, max.mismatch = max.mismatch, with.indels = with.indels, ...)

    if(return.table){
        return(result)
    } else {
        # get existing colData
        existing.cell.data <- as.data.frame(SingleCellExperiment::colData(TapestriExperiment))

        # drop columns if already exist to allow overwriting
        existing.cell.data <- existing.cell.data[,!colnames(existing.cell.data) %in% setdiff(colnames(result), "cell.barcode")]

        # merge result and existing colData
        updated.cell.data <- merge(existing.cell.data, result, by = "cell.barcode", all.x = TRUE, sort = FALSE)

        # set NAs to 0
        ids <- setdiff(colnames(result), "cell.barcode")
        updated.cell.data[,ids][is.na(updated.cell.data[,ids])] <- 0

        # reorder to match colData
        rownames(updated.cell.data) <- updated.cell.data$cell.barcode
        updated.cell.data <- updated.cell.data[rownames(SingleCellExperiment::colData(TapestriExperiment)),]

        # update TapestriExperiment
        SummarizedExperiment::colData(TapestriExperiment) <- S4Vectors::DataFrame(updated.cell.data)

        return(TapestriExperiment)
    }
}

#' Call sample labels based on feature counts
#'
#' `callSampleLables()` assigns labels (stored as `colData` column) to cells using feature count data in `colData`.
#' This is most useful for assigning barcode labels based on barcoded reads (see [countBarcodedReads]).
#' For `method = max`, labels are dictated by whichever `input.features` column has the highest number of counts.
#' By default, ties are broken by choosing whichever label has the lowest index position (`ties.method = "first"`).
#' Samples with 0 counts for all `input.features` columns are labeled according to `neg.label`.
#' If only one feature column is used, labels are assigned to cells with counts > `min.count.threshold`, and `neg.label` otherwise.
#'
#' @param TapestriExperiment A `TapestriExperiment` object.
#' @param input.features Character vector, column names in `colData` to evaluate.
#' @param method Character, call method. Only "max" currently supported, calls based on whichever `input.features` column has the most counts.
#' @param ties.method Character, passed to `max.col()` indicating how to break ties. Default "first".
#' @param neg.label Character, label for samples with no counts. Default `NA`.
#' @param output.feature Character, column name to use for the call output. Default "sample.call".
#' @param return.table Logical, if `TRUE`, returns a data.frame of the sample.calls. If `FALSE` (default), returns updated `TapestriExperiment` object.
#' @param min.count.threshold Numeric, minimum number of counts per cell to use for call. Default 1.
#'
#' @return A `TapestriExperiment` object with sample calls added to `colData` column `sample.name`. If `return.table == TRUE`, a `data.frame` of sample calls.
#' @export
#'
#' @concept barcoded reads
#'
#' @examples
#' \dontrun{
#' TapestriExperiment <- callSampleLables(TapestriExperiment,
#'   input.features = c("g7", "gNC"),
#'   output.feature = "sample.grna"
#' )
#' }
callSampleLables <- function(TapestriExperiment, input.features, output.feature = "sample.call", return.table = FALSE, neg.label = NA, method = "max", ties.method = "first", min.count.threshold = 1) {
  if (method != "max") {
    cli::cli_abort("Method not recognized. Only 'max' currently supported.")
  } else {
    # check for missing colData labels
    if (any(!input.features %in% colnames(SingleCellExperiment::colData(TapestriExperiment)))) {
      cli::cli_abort("{.q {input.features[!input.features %in% colnames(SingleCellExperiment::colData(TapestriExperiment))]}} not found in {.var colData}.")
    }

    # subset colData
    existing.cell.data <- as.data.frame(SingleCellExperiment::colData(TapestriExperiment))
    coldata.subset <- existing.cell.data[, input.features, drop = FALSE]

    # check if numeric
    if (any(!apply(coldata.subset, 2, is.numeric))) {
      cli::cli_abort("Selected input.features are not numeric.")
    }

    # apply min threshold, set counts to 0
    if(min.count.threshold < 1){
        cli::cli_abort("min.count.threshold must be positive.")
    } else {
        coldata.subset[coldata.subset < min.count.threshold] <- 0
    }

    # make calls
    sample.calls <- input.features[max.col(coldata.subset, ties.method = ties.method)]
    sample.calls <- data.frame(cell.barcode = rownames(coldata.subset), sample.call = sample.calls)
    sample.calls[rowSums(coldata.subset) == 0, "sample.call"] <- neg.label # set label if no call is made
    sample.calls$sample.call <- as.factor(sample.calls$sample.call)
    rownames(sample.calls) <- sample.calls$cell.barcode
    colnames(sample.calls)[2] <- output.feature

    if (return.table) {
      return(sample.calls)
    } else {
      # drop from existing data if already exist to allow overwriting
      existing.cell.data <- existing.cell.data[, !colnames(existing.cell.data) %in% setdiff(colnames(sample.calls), "cell.barcode")]

      # merge result and existing colData
      updated.cell.data <- merge(existing.cell.data, sample.calls, by = "cell.barcode", all.x = TRUE, sort = FALSE)

      # reorder to match colData
      rownames(updated.cell.data) <- updated.cell.data$cell.barcode
      updated.cell.data <- updated.cell.data[rownames(SingleCellExperiment::colData(TapestriExperiment)), ]

      # update TapestriExperiment
      SummarizedExperiment::colData(TapestriExperiment) <- S4Vectors::DataFrame(updated.cell.data)

      return(TapestriExperiment)
    }
  }
}

