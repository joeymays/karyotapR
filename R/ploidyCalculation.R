#' Generate Control Ploidy Template
#'
#' `generateControlPloidyTemplate()` generates a data.frame template for use in [`getPloidy()`].
#'
#' @param ploidy.all Numeric, sets all entries of `ploidy` in the output. Default 2.
#' @param sample.label.all Chr string, sets all entries of `sample.label` in the output. Default "cluster".
#'
#' @return data.frame with 3 columns named `arm`, `ploidy`, and `sample.label`
#' @export
#'
#' @describeIn getPloidy
#'
#' @examples
#' \dontrun{control.ploidy <- generateControlPloidyTemplate()}
generateControlPloidyTemplate <- function(ploidy.all = 2, sample.label.all = "cluster"){

    ploidy.template <- data.frame(arm = c("chr1p","chr1q","chr2p",
                                          "chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                          "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q",
                                          "chr9p","chr9q","chr10p","chr10q","chr11p","chr11q",
                                          "chr12p","chr12q","chr13q","chr14q","chr15q",
                                          "chr16p","chr16q","chr17p","chr17q","chr18p","chr18q",
                                          "chr19p","chr19q","chr20p","chr20q","chr21q","chr22q",
                                          "chrXp","chrXq"),
                                  ploidy = ploidy.all,
                                  sample.label = sample.label.all)
    rownames(ploidy.template) <- ploidy.template$arm

    return(ploidy.template)
}

#' Calculate ploidy from control cell population
#'
#' `getPloidy()` transforms the normalized count matrix of a TapestriExperiment object into ploidy values based on a reference subset of cell barcodes and given ploidy value (e.g. 2 for diploid).
#' This is practically used to set the median ploidy of a control, usually diploid, cell population to a known ploidy value, e.g. 2, and then calculate the ploidy for all the cells relative to that control population.
#' This occurs individually for each probe.
#' `control.ploidy` is a lookup table (data.frame) to indicate the ploidy value and cell barcodes to use as the reference.
#' A template for `control.ploidy` can be generated using [`generateControlPloidyTemplate()`].
#' The `control.ploidy` data.frame should include 3 columns named `arm`, `ploidy`, and `sample.label`.
#' `arm` is chromosome arms from chr1p through chrXq, `ploidy` is the reference ploidy value, and `sample.label` is the value corresponding to the colData column given in `coldata.set` to indicate the cell barcode subset to use to set the ploidy.
#'
#' @param TapestriExperiment TapestriExperiment object.
#' @param control.ploidy A data.frame with columns `arm`, `ploidy`, and `sample.label`. See description.
#' @param coldata.set A chr string indicating colData column to use for subsetting cell.barcodes. Default "cluster".
#'
#' @return TapestriExperiment object with ploidy values in `ploidy` assay slot.
#' @export
#'
#' @examples
#' \dontrun{TapestriExperiment <- getPloidy(TapestriExperiment,
#' control.ploidy, coldata.set = "cluster")}
getPloidy <- function(TapestriExperiment, control.ploidy, coldata.set = "cluster"){

    coldata.set <- tolower(coldata.set)

    # error checks
    if(!coldata.set %in% colnames(SummarizedExperiment::colData(TapestriExperiment))){
        stop(paste0("coldata.set '", coldata.set, "' not found in colData"))
    }

    if(any(!unique(control.ploidy$sample.label) %in% unique(SummarizedExperiment::colData(TapestriExperiment)[,coldata.set]))){
        stop(paste0("control.ploidy sample.label elements not found in colData. Check control.ploidy."))
    }

    counts.mat <- SummarizedExperiment::assay(TapestriExperiment, "normcounts")

    # get median normalized counts for each amplicon based on control.ploidy
    probe.table <- as.data.frame(SummarizedExperiment::rowData(TapestriExperiment))[, c("probe.id", "arm")]
    probe.table <- merge(probe.table, control.ploidy, by = "arm", all.x = T, sort = F)
    rownames(probe.table) <- probe.table$probe.id

    coldata.set.lookup <- SummarizedExperiment::colData(TapestriExperiment)[,coldata.set, drop = F]

    # define function for calculating median from cell subset
    getProbeMedian <- function(idx){
        probe.info <- probe.table[idx,]
        probe.median <- median(counts.mat[probe.info$probe.id, rownames(coldata.set.lookup)[coldata.set.lookup[,1] == probe.info$sample.label]])
        return(probe.median)
    }

    probe.medians <- lapply(seq_len(nrow(probe.table)), getProbeMedian)
    probe.medians <- unlist(probe.medians)
    names(probe.medians) <- probe.table$probe.id

    if(any(probe.medians == 0)){
        stop(paste0(names(probe.medians[probe.medians == 0]), " control cell median equal to 0. Filter out prior to proceeding.\n"))
    }

    probe.medians <- probe.medians[rownames(SummarizedExperiment::rowData(TapestriExperiment))] # reorder based on rowData
    counts.ploidy <- sweep(x = counts.mat, 1, probe.medians, "/") # normalize relative to medians
    probe.table <- probe.table[rownames(SummarizedExperiment::rowData(TapestriExperiment)),] # reorder based on rowData
    counts.ploidy <- sweep(x = counts.ploidy, 1, probe.table$ploidy, "*") # scale to control ploidy

    SummarizedExperiment::assay(TapestriExperiment, "ploidy") <- counts.ploidy

    return(TapestriExperiment)
}

#' Smooth Ploidy Values by Chromosome and Chromosome Arm
#'
#' `smoothPloidy()` takes ploidy values across probes and smooths them by median (default) for each chromosome and chromosome arm, resulting in one ploidy value per chromosome/arm per cell barcode.
#' Values are then discretized into integers by conventional rounding.
#' Results are stored as alternate experiments in the TapestriObject, with smoothPloidy and DiscretePloidy assays.
#'
#' @param TapestriExperiment TapestriExperiment object.
#' @param method Chr string indicating smoothing method. Supports median (default) or mean.
#'
#' @importFrom rlang .data
#'
#' @return TapestriExperiment object with smoothed ploidy values in altExp slot
#' @export
#'
#' @examples
#' \dontrun{TapestriExperiment <- smoothPloidy(TapestriExperiment)}
smoothPloidy <- function(TapestriExperiment, method = "median"){

    method <- tolower(method)

    if(method == "median"){
        smooth.func <- stats::median
    } else if(method == "mean"){
        smooth.func <- mean
    } else {
        stop(paste0("method '", method, "' not recognized. Please use mean or median."))
    }

    ploidy.counts <- SummarizedExperiment::assay(TapestriExperiment, "ploidy")

    ploidy.tidy <- ploidy.counts %>% as.data.frame() %>% tibble::rownames_to_column("probe.id")  %>%
        tidyr::pivot_longer(cols = !tidyr::matches("probe.id"), names_to = "cell.barcode", values_to = "ploidy") %>%
        dplyr::left_join(as.data.frame(SummarizedExperiment::rowData(TapestriExperiment)[,c("probe.id", "chr", "arm")]), by = "probe.id")

    smoothed.ploidy.chr <- ploidy.tidy %>% dplyr::group_by(.data$cell.barcode, .data$chr) %>% dplyr::summarize(smooth.ploidy = smooth.func(.data$ploidy), .groups = "drop") %>%
        tidyr::pivot_wider(id_cols = dplyr::all_of("chr"), values_from = dplyr::all_of("smooth.ploidy"), names_from = dplyr::all_of("cell.barcode")) %>% tibble::column_to_rownames("chr")

    smoothed.ploidy.arm <- ploidy.tidy %>% dplyr::group_by(.data$cell.barcode, .data$arm) %>% dplyr::summarize(smooth.ploidy = smooth.func(.data$ploidy), .groups = "drop") %>%
        tidyr::pivot_wider(id_cols = dplyr::all_of("arm"), values_from = dplyr::all_of("smooth.ploidy"), names_from = dplyr::all_of("cell.barcode")) %>% tibble::column_to_rownames("arm")

    discrete.ploidy.chr <- round(smoothed.ploidy.chr, 0)
    discrete.ploidy.arm <- round(smoothed.ploidy.arm, 0)


    smoothed.ploidy.chr <- SingleCellExperiment::SingleCellExperiment(list(smoothedPloidy = smoothed.ploidy.chr,
                                                                           discretePloidy = discrete.ploidy.chr))

    smoothed.ploidy.arm <- SingleCellExperiment::SingleCellExperiment(list(smoothedPloidy = smoothed.ploidy.arm,
                                                                           discretePloidy = discrete.ploidy.arm))

    smoothed.ploidy.chr <- .TapestriExperiment(smoothed.ploidy.chr)
    smoothed.ploidy.arm <- .TapestriExperiment(smoothed.ploidy.arm)

    SingleCellExperiment::altExp(TapestriExperiment, "smoothedPloidyByChrom", withDimnames = T) <- smoothed.ploidy.chr
    SingleCellExperiment::altExp(TapestriExperiment, "smoothedPloidyByArm", withDimnames = T) <- smoothed.ploidy.arm

    return(TapestriExperiment)
}

























