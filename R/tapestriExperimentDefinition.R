#' TapestriExperiment Class Definition
#'
#' @slot barcodeProbe character.
#' @slot grnaProbe character.
#'
#' @return TapestriExperiment
#' @export
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @examples
#' x <- new("TapestriExperiment")
.TapestriExperiment <- setClass(
    "TapestriExperiment",
    contains = "SingleCellExperiment",
    slots = c(barcodeProbe = "character", grnaProbe = "character")
)
