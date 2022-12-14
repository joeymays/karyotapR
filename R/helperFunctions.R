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

getTidyPloidy <- function(TapestriExperiment){

    ploidy.tidy <- SummarizedExperiment::assay(TapestriExperiment, "ploidy")

    ploidy.tidy <- as.data.frame(ploidy.tidy) %>% tibble::rownames_to_column("probe.id") %>% dplyr::as_tibble() %>%
        tidyr::pivot_longer(cols = !tidyr::matches("probe.id"), names_to = "cell.barcode", values_to = "ploidy")

    ploidy.tidy <- ploidy.tidy %>% dplyr::left_join(as.data.frame(SummarizedExperiment::colData(TapestriExperiment)), by = "cell.barcode")

    ploidy.tidy <- ploidy.tidy %>% dplyr::left_join(as.data.frame(SummarizedExperiment::rowData(TapestriExperiment)), by = "probe.id")

    #order factor by probe.id order for plotting
    ploidy.tidy$probe.id <- factor(ploidy.tidy$probe.id, levels = SummarizedExperiment::rowData(TapestriExperiment)$probe.id)

    return(ploidy.tidy)
}
