#' @param ploidy.all Numeric, sets all entries of `ploidy` column in output. Default 2.
#' @param sample.label.all Character, sets all entries of `sample.label` column in output. Default "cluster".
#'
#' @return `data.frame` with 3 columns named `arm`, `ploidy`, and `sample.label`
#' @export
#'
#' @describeIn getPloidy generates a `data.frame` template for `control.ploidy` in `getPloidy()`.
#' @order 2
#'
generateControlPloidyTemplate <- function(ploidy.all = 2, sample.label.all = "cluster") {
  ploidy.template <- data.frame(
    arm = c(
      "chr1p", "chr1q", "chr2p",
      "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q",
      "chr6p", "chr6q", "chr7p", "chr7q", "chr8p", "chr8q",
      "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q",
      "chr12p", "chr12q", "chr13q", "chr14q", "chr15q",
      "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q",
      "chr19p", "chr19q", "chr20p", "chr20q", "chr21q", "chr22q",
      "chrXp", "chrXq"
    ),
    ploidy = ploidy.all,
    sample.label = sample.label.all
  )
  rownames(ploidy.template) <- ploidy.template$arm

  return(ploidy.template)
}

#' @name getPloidy
#'
#' @title Calculate ploidy for each cell-chromosome pair using baseline control
#'
#' @description `getPloidy()` transforms the normalized count matrix of a `TapestriExperiment` object
#' into ploidy values based on a reference subset of cell barcodes and given ploidy value (e.g. 2 for diploid).
#' This is practically used to set the median ploidy of a usually diploid reference/control
#' cell population to a known ploidy value, e.g. 2, and then calculate the ploidy for all the
#' cells relative to that control population. This occurs individually for each probe,
#' such that the result is one ploidy value per cell barcode per probe.
#' `control.ploidy` is a `data.frame` lookup table used to indicate the ploidy value and cell barcodes
#' to use as the reference. A template for `control.ploidy` can be generated using [`generateControlPloidyTemplate()`].
#'
#' The `control.ploidy` data.frame should include 3 columns named `arm`, `ploidy`, and `sample.label`.
#' `arm` is chromosome arms from chr1p through chrXq, `ploidy` is the reference ploidy value, and `sample.label` is the
#' value corresponding to the `colData` column given in `sample.category` to indicate the cell barcode subset to use to set the ploidy.
#' This is best used in a workflow where the cells are clustered first, and then one cluster is used as the reference population
#' the other clusters. This also allows for the baseline ploidy to be set for each chromosome individually in the case where the
#' reference population is not completely diploid.
#'
#' @param TapestriExperiment `TapestriExperiment` object.
#' @param control.ploidy A `data.frame` with columns `arm`, `ploidy`, and `sample.label`. See details.
#' @param sample.category Character, `colData` column to use for subsetting cell.barcodes. Default "cluster".
#'
#' @return `TapestriExperiment` object with ploidy values in `ploidy` assay slot.
#' @export
#'
#' @rdname getPloidy
#' @order 1
#'
#' @examples
#' \dontrun{
#' control.ploidy <- generateControlPloidyTemplate()
#' TapestriExperiment <- getPloidy(TapestriExperiment,
#'   control.ploidy,
#'   sample.category = "cluster"
#' )
#' }
getPloidy <- function(TapestriExperiment, control.ploidy, sample.category = "cluster") {
  sample.category <- tolower(sample.category)

  # error checks
  if (!sample.category %in% colnames(SummarizedExperiment::colData(TapestriExperiment))) {
    stop(paste0("sample.category '", sample.category, "' not found in colData"))
  }

  if (any(!unique(control.ploidy$sample.label) %in% unique(SummarizedExperiment::colData(TapestriExperiment)[, sample.category]))) {
    stop(paste0("control.ploidy sample.label elements not found in colData. Check control.ploidy."))
  }

  counts.mat <- SummarizedExperiment::assay(TapestriExperiment, "normcounts")

  # get median normalized counts for each amplicon based on control.ploidy
  probe.table <- as.data.frame(SummarizedExperiment::rowData(TapestriExperiment))[, c("probe.id", "arm")]
  probe.table <- merge(probe.table, control.ploidy, by = "arm", all.x = T, sort = F)
  rownames(probe.table) <- probe.table$probe.id

  sample.category.lookup <- SummarizedExperiment::colData(TapestriExperiment)[, sample.category, drop = F]

  # define function for calculating median from cell subset
  getProbeMedian <- function(idx) {
    probe.info <- probe.table[idx, ]
    probe.median <- median(counts.mat[probe.info$probe.id, rownames(sample.category.lookup)[sample.category.lookup[, 1] == probe.info$sample.label]])
    return(probe.median)
  }

  probe.medians <- lapply(seq_len(nrow(probe.table)), getProbeMedian)
  probe.medians <- unlist(probe.medians)
  names(probe.medians) <- probe.table$probe.id

  if (any(probe.medians == 0)) {
    stop(paste0(names(probe.medians[probe.medians == 0]), " control cell median equal to 0. Filter out prior to proceeding.\n"))
  }

  probe.medians <- probe.medians[rownames(SummarizedExperiment::rowData(TapestriExperiment))] # reorder based on rowData
  counts.ploidy <- sweep(x = counts.mat, 1, probe.medians, "/") # normalize relative to medians
  probe.table <- probe.table[rownames(SummarizedExperiment::rowData(TapestriExperiment)), ] # reorder based on rowData
  counts.ploidy <- sweep(x = counts.ploidy, 1, probe.table$ploidy, "*") # scale to control ploidy

  SummarizedExperiment::assay(TapestriExperiment, "ploidy") <- counts.ploidy

  return(TapestriExperiment)
}

#' Smooth ploidy values by chromosome and chromosome arm
#'
#' `smoothPloidy()` takes `ploidy` slot values for probes on a chromosome and smoothes them by median (default) for each chromosome
#' and chromosome arm, resulting in one ploidy value per chromosome and chromosome arm for each cell barcode.
#' Values are then discretized into integers by conventional rounding.
#' Smoothed ploidy and discretized smoothed ploidy values are stored as `smoothPloidy` and `DiscretePloidy` assays,
#' in `altExp` slots `smoothedPloidyByChrom` for chromosome-level smoothing, and `smoothedPloidyByArm` for chromosome arm-level smoothing.
#'
#' @param TapestriExperiment `TapestriExperiment` object.
#' @param method Character, smoothing method. Supports median (default) or mean.
#'
#' @importFrom rlang .data
#'
#' @return `TapestriExperiment` with `smoothPloidy` and `DiscretePloidy` assays in `altExp` slots `smoothedPloidyByChrom` and `smoothedPloidyByArm`.
#' @export
#'
#' @examples
#' \dontrun{
#' TapestriExperiment <- smoothPloidy(TapestriExperiment)
#' }
smoothPloidy <- function(TapestriExperiment, method = "median") {
  method <- tolower(method)

  if (method == "median") {
    smooth.func <- stats::median
  } else if (method == "mean") {
    smooth.func <- mean
  } else {
    stop(paste0("method '", method, "' not recognized. Please use mean or median."))
  }

  ploidy.counts <- SummarizedExperiment::assay(TapestriExperiment, "ploidy")

  ploidy.tidy <- ploidy.counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("probe.id") %>%
    tidyr::pivot_longer(cols = !tidyr::matches("probe.id"), names_to = "cell.barcode", values_to = "ploidy") %>%
    dplyr::left_join(as.data.frame(SummarizedExperiment::rowData(TapestriExperiment)[, c("probe.id", "chr", "arm")]), by = "probe.id")

  smoothed.ploidy.chr <- ploidy.tidy %>%
    dplyr::group_by(.data$cell.barcode, .data$chr) %>%
    dplyr::summarize(smooth.ploidy = smooth.func(.data$ploidy), .groups = "drop") %>%
    tidyr::pivot_wider(id_cols = dplyr::all_of("chr"), values_from = dplyr::all_of("smooth.ploidy"), names_from = dplyr::all_of("cell.barcode")) %>%
    tibble::column_to_rownames("chr")

  smoothed.ploidy.arm <- ploidy.tidy %>%
    dplyr::group_by(.data$cell.barcode, .data$arm) %>%
    dplyr::summarize(smooth.ploidy = smooth.func(.data$ploidy), .groups = "drop") %>%
    tidyr::pivot_wider(id_cols = dplyr::all_of("arm"), values_from = dplyr::all_of("smooth.ploidy"), names_from = dplyr::all_of("cell.barcode")) %>%
    tibble::column_to_rownames("arm")

  discrete.ploidy.chr <- round(smoothed.ploidy.chr, 0)
  discrete.ploidy.arm <- round(smoothed.ploidy.arm, 0)


  smoothed.ploidy.chr <- SingleCellExperiment::SingleCellExperiment(list(
    smoothedPloidy = smoothed.ploidy.chr,
    discretePloidy = discrete.ploidy.chr
  ))

  smoothed.ploidy.arm <- SingleCellExperiment::SingleCellExperiment(list(
    smoothedPloidy = smoothed.ploidy.arm,
    discretePloidy = discrete.ploidy.arm
  ))

  smoothed.ploidy.chr <- .TapestriExperiment(smoothed.ploidy.chr)
  smoothed.ploidy.arm <- .TapestriExperiment(smoothed.ploidy.arm)

  SingleCellExperiment::altExp(TapestriExperiment, "smoothedPloidyByChrom", withDimnames = T) <- smoothed.ploidy.chr
  SingleCellExperiment::altExp(TapestriExperiment, "smoothedPloidyByArm", withDimnames = T) <- smoothed.ploidy.arm

  return(TapestriExperiment)
}
