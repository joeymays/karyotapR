#' Print the top-left corner of a matrix
#'
#' Outputs up to 5 rows and columns of the input matrix object (with `rownames` and `colnames`) to get a quick look without filling the console.
#'
#' @param input.mat A matrix-like object.
#'
#' @return A matrix-like object matching input class, subset to a maximum of 5 rows and columns.
#' @export
#'
#' @examples
#' \dontrun{
#' corner(assay(TapestriObject))
#' }
corner <- function(input.mat) {
  if (nrow(input.mat) > 4) {
    row.out <- 5
  } else {
    row.out <- nrow(input.mat)
  }

  if (ncol(input.mat) > 4) {
    col.out <- 5
  } else {
    col.out <- ncol(input.mat)
  }

  input.mat[1:row.out, 1:col.out]
}

#' Get tidy-style data from `TapestriExperiment` objects
#'
#' `getTidyData()` pulls data from the indicated `assay` and/or `altExp` slot(s), and rearranges it into tidy format.
#' `colData` (cell metadata) from the top-level/main experiment is included.
#' `rowData` (probe metadata) from the indicated `assay` and/or `altExp` slot(s) is included.
#' Attempts are made to sort by "chr" and "start.pos" columns if they are present to simplify plotting and other downstream operations.
#'
#' @param TapestriExperiment `TapestriExperiment` object.
#' @param alt.exp Character, `altExp` slot to use. `NULL` (default) uses top-level/main experiment.
#' @param assay Character, `assay` slot to use. `NULL` (default) uses first-indexed assay (often "counts").
#' @param feature.id.as.factor Logical, if `TRUE` (default), the feature.id column is returned as a factor.
#'
#' @return A `tibble` of tidy data with corresponding metadata from `colData` and `rowData`.
#' @export
#'
#' @examples \dontrun{getTidyData(TapestriObject, alt.exp = "alleleFrequency")}
getTidyData <- function(TapestriExperiment, alt.exp = NULL, assay = NULL, feature.id.as.factor = TRUE){

    if(is.null(alt.exp)){
        target.exp <- TapestriExperiment
    } else {
        target.exp <- altExp(TapestriExperiment, alt.exp)
    }

    if(is.null(assay)){
        assay <- SummarizedExperiment::assayNames(target.exp)[1]
    } else if(!assay %in% SummarizedExperiment::assayNames(target.exp)){
        stop(paste0("assay '", assay, "' not found. Available assays are:\n",
                    paste0(SummarizedExperiment::assayNames(target.exp), collapse = ", ")))
    }

    tidy.data <- SummarizedExperiment::assay(target.exp, assay)

    tidy.data <- as.data.frame(tidy.data) %>% tibble::rownames_to_column("feature.id") %>% dplyr::as_tibble() %>%
        tidyr::pivot_longer(cols = !tidyr::matches("feature.id"), names_to = "cell.barcode", values_to = assay)

    tidy.data <- tidy.data %>% dplyr::left_join(as.data.frame(SummarizedExperiment::colData(TapestriExperiment)), by = "cell.barcode")

    # get rowdata and add rownames as a column to target for join with feature.id
    rowdata.to.join <- as.data.frame(SummarizedExperiment::rowData(target.exp)) %>% tibble::rownames_to_column("feature.id")

    tidy.data <- tidy.data %>% dplyr::left_join(rowdata.to.join,
                                                by = "feature.id",
                                                suffix = c(".bc", ".feature"))

    # attempt to sort by chr and start.pos if present
    if(all(c("chr", "start.pos") %in% colnames(tidy.data))){
        tidy.data <- tidy.data %>% dplyr::arrange(.data$chr, .data$start.pos)
    } else if("chr" %in% colnames(tidy.data)){
        tidy.data <- tidy.data %>% dplyr::arrange(.data$chr)
    }

    # make feature.id a factor to enable sorting in visualization
    if(feature.id.as.factor){
        tidy.data$feature.id <- factor(tidy.data$feature.id, levels = unique(tidy.data$feature.id))
    }

    return(tidy.data)
}

#' Check for available `altExp` and `assay` slots
#'
#' Checks if assay is available in given `altExp.` If not, lists available assays in error.
#' If assay is `NULL`, sets assay to first assay in `mainExp` or `altExp.`
#'
#' @param TapestriExperiment `TapestriExperiment` object
#' @param alt.exp Character, `altExp` slot to use. `NULL` (default) uses top-level/main experiment.
#' @param assay Character, `assay` slot to use. `NULL` (default) uses first-indexed assay (often "counts").
#'
#' @return Character, assay selection
#' @noRd
#'
.SelectAssay <- function(TapestriExperiment, alt.exp = NULL, assay = NULL){

    if(is.null(alt.exp)){
        target.exp <- TapestriExperiment
    } else {
        target.exp <- altExp(TapestriExperiment, alt.exp)
    }

    if(is.null(assay)){
        assay <- SummarizedExperiment::assayNames(target.exp)[1]
    } else if(!assay %in% SummarizedExperiment::assayNames(target.exp)){
        stop(paste0("assay '", assay, "' not found. Available assays are:\n",
                    paste0(SummarizedExperiment::assayNames(target.exp), collapse = ", ")))
    }

    return(assay)
}
