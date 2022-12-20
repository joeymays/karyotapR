#' Move non-genome probes counts and metadata to altExp slots
#'
#' `moveNonGenomeProbes()` takes the probe IDs corresponding to `grnaProbe` and `barcodeProbe` slots of the `TapestriExperiment` object,
#' as well as probes on chrY, and moves them to their own `altExp` slots in the object.
#' This allows those counts and associated metadata to be manipulated separately without interfering with the other probes in the panel.
#' [SingleCellExperiment::splitAltExps()] can be used for manual specification of probes to move to `altExp` slots if the shortcut slots are not used.
#'
#' `moveNonGenomeProbes()` moves probes corresponding to the specified tags to `altExp` (alternative experiment) slots in the `TapestriExperiment` object.
#' These probes should be those which do not correspond to a chromosome and therefore would not be used to call copy number variants.
#' The exception is probes on chromosome Y; chrY does not generally experience CNVs, so we move it to an `altExp` for separate analysis.
#' Probes corresponding to the `barcodeProbe` and `grnaProbe` slots, which are specified by the `panel.id` shortcut or manually (see [Custom Slot Getters and Setters]),
#' are automatically moved to `altExp` by this operation as well.
#' If such probes are not present, the function will generate a message but not throw an error, so it is always safe to run by default.
#' This is run automatically by default and with default behavior as part of [createTapestriExperiment()].
#'
#' @param TapestriExperiment `TapestriExperiment` object.
#' @param move.non.genome.probes Character vector, a combination of "grna", "sample.barcode", and/or "Y" (default), or `FALSE`.
#'
#' @return A `TapestriExperiment` with `altExp` slots filled with counts and metadata for the specified probes.
#' @export
#'
#' @seealso [SingleCellExperiment::splitAltExps()] for manual specification of probes to move to `altExp` slots.
#'
#' @examples
#' \dontrun{
#' tapObject <- moveNonGenomeProbes(tapObject)
#' }
moveNonGenomeProbes <- function(TapestriExperiment, move.non.genome.probes = c("grna", "sample.barcode", "Y")) {
  if (identical(move.non.genome.probes, FALSE)) {
    stop("No non-genomic probes specified.")
  }

  move.non.genome.probes <- tolower(move.non.genome.probes)

  if (any(!move.non.genome.probes %in% c("grna", "sample.barcode", "y"))) {
    stop(paste(move.non.genome.probes, "<- Non-genomic probes not recognized. Use a combination of 'grna', 'sample.barcode', and/or 'Y'."))
  }

  feature.type <- rep("CNV", nrow(TapestriExperiment))

  barcodeProbe <- TapestriExperiment@barcodeProbe
  grnaProbe <- TapestriExperiment@grnaProbe

  if ("grna" %in% move.non.genome.probes && grnaProbe != "not specified") {
    probe.index <- which(rownames(TapestriExperiment) == grnaProbe)

    if (S4Vectors::isEmpty(probe.index)) {
      message("gRNA probe ID not found in TapestriExperiment object.")
    } else {
      feature.type[probe.index] <- "grnaCounts"
      message(paste("Moving gRNA probe", rownames(TapestriExperiment)[probe.index], "to altExp slot 'grnaCounts'."))
    }
  }

  if ("sample.barcode" %in% move.non.genome.probes && barcodeProbe != "not specified") {
    probe.index <- which(rownames(TapestriExperiment) == barcodeProbe)

    if (S4Vectors::isEmpty(probe.index)) {
      message("Sample Barcode probe ID not found in TapestriExperiment object.")
    } else {
      feature.type[probe.index] <- "sampleBarcodeCounts"
      message(paste("Moving sample barcode probe", rownames(TapestriExperiment)[probe.index], "to altExp slot 'sampleBarcodeCounts'."))
    }
  }

  if ("y" %in% move.non.genome.probes) {
    probe.index <- which(rowData(TapestriExperiment)$chr == "Y")

    if (S4Vectors::isEmpty(probe.index)) {
      message("ChrY probe ID(s) not found in TapestriExperiment object.")
    } else {
      feature.type[probe.index] <- "chrYCounts"
      message(paste("Moving chrY probe(s)", paste(rownames(TapestriExperiment)[probe.index], collapse = ", "), "to altExp slot 'chrYCounts'."))
    }
  }

  if (all(feature.type == "CNV")) {
    message("No non-genomic probe IDs found.")
    return(TapestriExperiment)
  }

  TapestriExperiment <- SingleCellExperiment::splitAltExps(TapestriExperiment, feature.type, ref = "CNV")

  return(TapestriExperiment)
}
