#' Move specialty probe counts to alternative experiment slots
#'
#' `splitSpecialtyProbes()` takes the probe IDs corresponding to `grnaProbe` and `barcodeProbe` slots of the `TapestriExperiment` object, as well as probes on chrY, and moves them to their own `altExp` (alternative experiment) slots in the object.
#' This allows those counts and associated metadata to be manipulated separately without interfering with the other probes in the panel.
#'
#' @param tapestri.experiment.object TapestriExperiment object
#' @param separate.specialty.probes Chr vector including a combination of "grna", "sample.barcode", and or "Y" (default), or FALSE.
#'
#' @return A `TapestriExperiment` with `altExp` slots filled with counts and metadata for the specified probes.
#' @export
#'
#' @seealso [SingleCellExperiment::splitAltExps()] for manual specification of probes to move to `altExp` slots.
#'
#' @examples
#' \dontrun{tapObject <- splitSpecialtyProbes(tapObject)}
splitSpecialtyProbes <- function(tapestri.experiment.object, separate.specialty.probes = c("grna", "sample.barcode", "Y")){

    if(identical(separate.specialty.probes, F)){
        stop("No specialty probes specified.")
    }

    separate.specialty.probes <- tolower(separate.specialty.probes)

    if(any(!separate.specialty.probes %in% c("grna", "sample.barcode", "y"))){
        stop(paste(separate.specialty.probes, "<- Specialty probes not recognized. Use a combination of 'grna', 'sample.barcode', and/or 'Y'."))
    }

    feature.type <- rep("normal", nrow(tapestri.experiment.object))

    barcodeProbe <- tapestri.experiment.object@barcodeProbe
    grnaProbe <- tapestri.experiment.object@grnaProbe

    if("grna" %in% separate.specialty.probes && grnaProbe != "not specified"){

        probe.index <- which(rownames(tapestri.experiment.object) == grnaProbe)

        if(S4Vectors::isEmpty(probe.index)){
            message("gRNA probe ID not found in TapestriExperiment object.")
        } else {
            feature.type[probe.index] <- "grnaCounts"
            message(paste("Moving gRNA probe", rownames(tapestri.experiment.object)[probe.index], "to altExp slot 'grnaCounts'."))
        }
    }

    if("sample.barcode" %in% separate.specialty.probes && barcodeProbe != "not specified"){

        probe.index <- which(rownames(tapestri.experiment.object) == barcodeProbe)

        if(S4Vectors::isEmpty(probe.index)){
            message("Sample Barcode probe ID not found in TapestriExperiment object.")
        } else {

            feature.type[probe.index] <- "sampleBarcodeCounts"
            message(paste("Moving sample barcode probe", rownames(tapestri.experiment.object)[probe.index], "to altExp slot 'sampleBarcodeCounts'."))
        }
    }

    if("y" %in% separate.specialty.probes){

        probe.index <- which(rowData(tapestri.experiment.object)$chr == "Y")

        if(S4Vectors::isEmpty(probe.index)){
            message("ChrY probe ID(s) not found in TapestriExperiment object.")
        } else {

            feature.type[probe.index] <- "chrYCounts"
            message(paste("Moving chrY probe(s)", paste(rownames(tapestri.experiment.object)[probe.index], collapse = ", "), "to altExp slot 'chrYCounts'."))
        }
    }

    if(all(feature.type == "normal")){
        stop("No specialty probe IDs found. Aborting.")
    }

    tapestri.experiment.object <- SingleCellExperiment::splitAltExps(tapestri.experiment.object, feature.type, ref = "normal")

    return(tapestri.experiment.object)
}
