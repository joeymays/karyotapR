#' Normalize raw counts
#'
#' Normalizes raw counts from `counts` slot in `TapestriExperiment` and returns the object with normalized counts in the `normcounts` slot.
#' "mb" is the only method supported currently.
#' This method performs the same normalization scheme as in Mission Bio's mosaic package for python.
#' Counts for each barcode are normalized relative to their barcode's mean and probe counts are normalized relative to their probe's median.
#' Also calculates the standard deviation for each probe using normalized counts and adds it to `rowData`.
#'
#' @param TapestriExperiment `TapestriExperiment` object.
#' @param method Character, normalization method. Default "mb".
#'
#' @return `TapestriExperiment` object with normalized counts added to `normcounts` slot.
#' @export
#'
#' @importFrom stats median sd
#'
#' @examples
#' \dontrun{TapestriExperiment <- normalizeCounts(TapestriExperiment)}

normalizeCounts <- function(TapestriExperiment, method = "mb"){

    method <- tolower(method)

    if(method != "mb"){
        message("Only method 'mb' is currently supported. Running as `method = 'mb'.")
    }

    raw.count.matrix <- SummarizedExperiment::assay(TapestriExperiment, "counts")

    # get "good barcodes", barcodes that have at least 10% the counts of the 11th barcode ranked for highest number of counts
    barcode.sums <- apply(raw.count.matrix, MARGIN = 2, sum)
    barcode.sorted <- sort(barcode.sums, decreasing = T)
    good.barcodes <- barcode.sums > (barcode.sorted[11] / 10)

    # normalize barcodes relative to barcode means
    barcode.means <- apply(raw.count.matrix, 2, mean) + 1
    read.counts.normal <- sweep(raw.count.matrix, 2, barcode.means, "/")

    # normalize probes relative to probe median. medians calculated using "good barcodes"
    probe.medians <- apply(read.counts.normal[,good.barcodes], 1, median) + 0.05
    read.counts.normal <- sweep(read.counts.normal, 1, probe.medians, "/")

    # scale all counts by 2X for diploid baseline
    read.counts.normal <- read.counts.normal * 2

    # add to normalized counts slot
    normcounts(TapestriExperiment) <- read.counts.normal

    # calculate norm count SD for each amplicon and add to rowData
    current.probe.data <- SummarizedExperiment::rowData(TapestriExperiment)
    current.probe.order <- rownames(current.probe.data)
    new.probe.data <- apply(read.counts.normal, 1, sd)
    names(new.probe.data) <- rownames(read.counts.normal)
    new.probe.data <- new.probe.data[current.probe.order]
    SummarizedExperiment::rowData(TapestriExperiment)$norm.count.sd <-  new.probe.data

    return(TapestriExperiment)
}
