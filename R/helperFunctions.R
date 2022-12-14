#' Glimpse the top corner of a matrix
#'
#' Outputs up to 5 rows and columns of the input matrix object to get a quick look without filling the console.
#'
#' @param input.mat a matrix-like object
#'
#' @return matrix-like object matching input class
#' @export
#'
#' @examples
#' \dontrun{corner(assay(TapestriObject))}
corner <- function(input.mat){
    if(nrow(input.mat) > 4){
        row.out = 5
    } else { row.out = nrow(input.mat)}

    if(ncol(input.mat) > 7){
        col.out = 8
    } else { col.out = ncol(input.mat)}

    print(input.mat[1:row.out, 1:col.out])
}

#' Get tidy-style data from TapestriExperiment objects
#'
#' `getTidyData()` pulls the matrix from the indicated assay and/or altExp slot(s), and morphs it into tidy format.
#' ColData from the top-level/main experiment is merged.
#' RowData from the indicated assay and/or altExp slot(s) is merged.
#' Attempts are made to sort by "chr" and "start.pos" columns if they are present in order to
#'
#' @param TapestriExperiment A TapestriExperiment object
#' @param alt.exp Chr string indicating altExp slot to pull from. `NULL` (default) pulls from top-level/main experiment.
#' @param assay Chr string indicating assay slot to pull from. `NULL` (default) pulls from first-indexed assay (often "counts").
#' @param feature.id.as.factor Logical indicating whether feature.id column should be coerced into a factor. In practice, used to preserve order for downstream plotting. Default T.
#'
#' @return A tibble of tidy data with corresponding metadata from colData and rowData.
#' @export
#'
#' @examples \dontrun{getTidyData(TapestriObject, alt.exp = "alleleFrequency")}
getTidyData <- function(TapestriExperiment, alt.exp = NULL, assay = NULL, feature.id.as.factor = T){

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
