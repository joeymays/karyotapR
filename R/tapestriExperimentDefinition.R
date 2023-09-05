#' TapestriExperiment Class Definition
#'
#' @slot barcodeProbe character.
#' @slot grnaProbe character.
#' @slot gmmParams list.
#'
#' @return TapestriExperiment object
#' @export
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @examples
#' tapExpObject <- new("TapestriExperiment")
.TapestriExperiment <- setClass(
  "TapestriExperiment",
  contains = "SingleCellExperiment",
  slots = c(
    barcodeProbe = "character",
    grnaProbe = "character",
    gmmParams = "list"
  )
)


#' @describeIn .TapestriExperiment Show method for `TapestriExperiment`
#' @param TapestriExperiment A `TapestriExperiment` object
#' @param object An R object
#'
#' @concept build experiment
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
    S4Vectors::coolcat("gmmParams(%d): %s\n", names(object@gmmParams))
})


#' @name Custom Slot Getters and Setters
#'
#' @title Getter and Setter functions for `TapestriExperiment` slots
#'
#' @description Get and set custom slots in `TapestriExperiment`. Slots include
#' `barcodeProbe` for a sample barcode probe ID
#' and `grnaProbe` for a gRNA-associated probe ID. These are used as shortcuts for
#' [moveNonGenomeProbes()] and [countBarcodedReads()].
#' `gmmParams` holds parameters and metadata for GMM copy number calling models.
#'
#' @param x A `TapestriExperiment` object
#'
#' @concept build experiment
#'
#' @export
#' @rdname slotGettersSetters
#' @aliases barcodeProbe
#' 
#' @return For the getter methods `barcodeProbe`, `grnaProbe`, and `gmmParams`, the value
#' of the given slot is returned. For the setter methods `barcodeProbe` and `grnaProbe`,
#' a `TapestriExperiment` object is returned with modifications made to the given slot.
#' 
#' @examples
#' tap.object <- newTapestriExperimentExample() # example TapestriExperiment object
#' barcodeProbe(tap.object) <- "Probe01"
#' barcodeProbe(tap.object)
#' 
#' grnaProbe(tap.object) <- "Probe02"
#' grnaProbe(tap.object)
#' 
#' gmmParams(tap.object)
setGeneric("barcodeProbe", function(x) standardGeneric("barcodeProbe"))

#' @param TapestriExperiment A `TapestriExperiment` object
#'
#' @export
#' @describeIn slotGettersSetters barcodeProbe getter
setMethod("barcodeProbe", "TapestriExperiment", function(x) x@barcodeProbe)

#' @param x A `TapestriExperiment` object
#'
#' @param value Character, probe ID to assign to slot
#'
#' @export
#' @rdname slotGettersSetters
setGeneric("barcodeProbe<-", function(x, value) standardGeneric("barcodeProbe<-"))

#' @param TapestriExperiment A `TapestriExperiment` object
#'
#' @export
#' @describeIn slotGettersSetters barcodeProbe setter
setMethod("barcodeProbe<-", "TapestriExperiment", function(x, value) {
  x@barcodeProbe <- value
  validObject(x)
  x
})


#' @param x A `TapestriExperiment` object
#'
#' @export
#' @rdname slotGettersSetters
setGeneric("grnaProbe", function(x) standardGeneric("grnaProbe"))

#' @param TapestriExperiment A `TapestriExperiment` object
#'
#' @export
#' @describeIn slotGettersSetters grnaProbe getter
setMethod("grnaProbe", "TapestriExperiment", function(x) x@grnaProbe)

#' @param x A `TapestriExperiment` object
#'
#' @param value Character, probe ID to assign to slot
#'
#' @export
#' @rdname slotGettersSetters
setGeneric("grnaProbe<-", function(x, value) standardGeneric("grnaProbe<-"))

#' @param TapestriExperiment A `TapestriExperiment` object
#'
#' @export
#' @describeIn slotGettersSetters grnaProbe setter
setMethod("grnaProbe<-", "TapestriExperiment", function(x, value) {
  x@grnaProbe <- value
  validObject(x)
  x
})


#' @param x A `TapestriExperiment` object
#'
#' @export
#' @rdname slotGettersSetters
setGeneric("gmmParams", function(x) standardGeneric("gmmParams"))

#' @param TapestriExperiment A `TapestriExperiment` object
#'
#' @export
#' @describeIn slotGettersSetters gmmParams getter
setMethod("gmmParams", "TapestriExperiment", function(x) x@gmmParams)




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
  } else {
    TRUE
  }
})
