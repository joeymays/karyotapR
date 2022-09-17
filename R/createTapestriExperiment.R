#' Create TapestriExperiment object from Tapestri Pipeline output
#'
#' @param h5.filename file path for .h5 file form Tapestri Pipeline output
#' @param panel.id Tapestri panel name, CO261 and CO293 supported only.
#'
#' @return TapestriExperiment
#' @export
#'
#' @import rhdf5
#' @import S4Vectors
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @examples
#' \dontrun{x <- CreateTapestriExperiment("myh5file.h5", "CO293")}
CreateTapestriExperiment <- function(h5.filename, panel.id = ""){

    if(!panel.id %in% c("CO261", "CO293")){
        stop(paste("panel.id", panel.id, "is not recognized. Please use CO261 or CO293."))
    }

    tapestri.h5 <- rhdf5::H5Fopen(file.path(h5.filename))

    read.counts.raw <- t(matrix(data = tapestri.h5$"/assays/dna_read_counts/layers/read_counts",
                                ncol=length(tapestri.h5$'/assays/dna_read_counts/ca/id'),
                                byrow = T))

    sce <- SingleCellExperiment::SingleCellExperiment(list(read.counts.raw = read.counts.raw),
                                                      colData = S4Vectors::DataFrame(cell.barcode = tapestri.h5$'/assays/dna_read_counts/ra/barcode',
                                                                                     row.names = tapestri.h5$'/assays/dna_read_counts/ra/barcode'),
                                                      rowData = S4Vectors::DataFrame(probe.id = tapestri.h5$'/assays/dna_read_counts/ca/id',
                                                                                     row.names = tapestri.h5$'/assays/dna_read_counts/ca/id'))

    tapestri.object <- .TapestriExperiment(sce)

    #barcode and gRNA provided in function call are currently orphaned
    if(panel.id == "CO293"){
        tapestri.object@barcodeProbe = "AMPL205334"
        tapestri.object@grnaProbe = "AMPL205666"
    } else if(panel.id == "CO261") {
        tapestri.object@barcodeProbe = "CO261_AMP1"
        tapestri.object@grnaProbe = NULL
    }

    rhdf5::H5Fclose(tapestri.h5)

    return(tapestri.object)

}
