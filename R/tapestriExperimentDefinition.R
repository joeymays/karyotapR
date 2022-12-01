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
    slots = c(barcodeProbe = "character",
              grnaProbe = "character")
)


#' @describeIn .TapestriExperiment Show method for Tapestri Experiment
#' @param TapestriExperiment A TapestriExperiment object
#' @param object An R object
#'
#' @export
#' @importMethodsFrom SingleCellExperiment show
setMethod("show", "TapestriExperiment", function(object) {
    callNextMethod()
    cat(
        "barcodeProbe: ", object@barcodeProbe, "\n",
        "grnaProbe: ", object@grnaProbe, "\n",
        sep = ""
    )
})

setGeneric("barcodeProbe", function(x) standardGeneric("barcodeProbe"))
setMethod("barcodeProbe", "TapestriExperiment", function(x) x@barcodeProbe)
setGeneric("barcodeProbe<-", function(x, value) standardGeneric("barcodeProbe<-"))
setMethod("barcodeProbe<-", "TapestriExperiment", function(x, value) {
    x@barcodeProbe <- value
    validObject(x)
    x
})

setGeneric("grnaProbe", function(x) standardGeneric("grnaProbe"))
setMethod("grnaProbe", "TapestriExperiment", function(x) x@grnaProbe)
setGeneric("grnaProbe<-", function(x, value) standardGeneric("grnaProbe<-"))
setMethod("grnaProbe<-", "TapestriExperiment", function(x, value) {
    x@grnaProbe <- value
    validObject(x)
    x
})


S4Vectors::setValidity2("TapestriExperiment", function(object) {
    msg <- NULL

    if (length(grnaProbe(object)) > 1) {
        msg <- c(msg, "'grnaProbe' should have length equal to 1")
    }
    if (length(barcodeProbe(object)) > 1) {
        msg <- c(
            msg, "'barcodeProbe' should have length equal to 1"
        )
    }

    if (length(msg)) {
        msg
    } else TRUE
})
